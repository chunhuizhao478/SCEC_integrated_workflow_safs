"""introspect.py -- recover a project descriptor from a SHIPPED deck.

This is what makes the reproduction exercise's E2 rung tractable, and it is independently
useful to anyone adopting the workflow with decks they already have.

It reads what the deck actually contains -- nc axes, easi `!ConstantMap` values, the
literals inside the `rs_muw` and `Tnuc_s` Lua, and `parameters.par` -- and emits a DRAFT
descriptor plus the build parameters it inferred.

Everything it cannot infer is marked `UNKNOWN`, never guessed.  A silently-wrong default
here would poison the whole exercise: the descriptor is what the rebuild is driven from,
so a guessed value would be "reproduced" perfectly and mean nothing.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

__all__ = ["introspect_deck", "DeckFacts", "UNKNOWN"]

UNKNOWN = "UNKNOWN"


class DeckFacts(dict):
    """What was recovered, plus what could not be."""

    @property
    def unknown(self) -> list[str]:
        return sorted(k for k, v in self.items() if v is UNKNOWN or v == UNKNOWN)

    def report(self) -> str:
        w = max((len(k) for k in self), default=0)
        lines = []
        for k in sorted(self):
            v = self[k]
            mark = "  ?" if (v is UNKNOWN or v == UNKNOWN) else "   "
            lines.append(f"{mark} {k:<{w}} : {v}")
        if self.unknown:
            lines.append(f"\n  {len(self.unknown)} field(s) could NOT be inferred: "
                         f"{self.unknown}")
        return "\n".join(lines)


def _axes(nc):
    from netCDF4 import Dataset
    with Dataset(str(nc)) as ds:
        return tuple(np.asarray(ds.variables[a][:], float) for a in ("x", "y", "z"))


def _box(nc) -> dict:
    x, y, z = _axes(nc)
    return {"xmin": float(x[0]), "xmax": float(x[-1]),
            "ymin": float(y[0]), "ymax": float(y[-1]),
            "zmin": float(z[0]), "zmax": float(z[-1]),
            "dx": float(np.diff(x).mean()), "dz": float(np.diff(z).mean())}


def introspect_deck(deck_dir) -> DeckFacts:
    """Read a shipped deck and report everything recoverable."""
    d = Path(deck_dir)
    if not d.is_dir():
        raise FileNotFoundError(f"deck not found: {d}")
    f = DeckFacts()
    f["deck"] = d.name

    # ---- grids, straight off the nc axes ----------------------------------------
    for kind, pat in (("stress", "*stress*.nc"), ("material", "*material*.nc"),
                      ("friction", "*friction*.nc"), ("plasticity", "*plasticity*.nc"),
                      ("thermal", "*thermal*.nc")):
        hit = next(d.glob(pat), None)
        f[f"{kind}_nc"] = hit.name if hit else UNKNOWN
        if hit:
            b = _box(hit)
            f[f"{kind}_grid"] = (f"x[{b['xmin']:.0f},{b['xmax']:.0f}] "
                                 f"y[{b['ymin']:.0f},{b['ymax']:.0f}] "
                                 f"z[{b['zmin']:.0f},{b['zmax']:.0f}] "
                                 f"dx={b['dx']:.0f} dz={b['dz']:.0f}")
            if kind == "stress":
                f["stress_box"] = b

    mesh = next(d.glob("*.puml.h5"), None)
    f["mesh"] = mesh.name if mesh else UNKNOWN

    # ---- the stress design, from the FILENAME and then verified against the field --
    s = f.get("stress_nc")
    if s and s != UNKNOWN:
        m = re.search(r"_k(\d+(?:\.\d+)?)", s)
        f["k_from_filename"] = float(m.group(1)) if m else UNKNOWN
        # The Phase 2 naming rule: the freeze tag appears only when the freeze is ON.
        mz = re.search(r"_freeze(\d+(?:\.\d+)?)m", s)
        f["freeze_from_filename"] = float(mz.group(1)) if mz else 0.0
        f["freeze_note"] = ("no _freeze tag in the name -> the shallow freeze is OFF for "
                            "this deck; VERIFY against the field, do not trust the name")

    # ---- easi yamls: constants and the Lua literals -------------------------------
    for y in sorted(d.glob("*.yaml")):
        t = y.read_text()
        for key in ("rs_b", "rs_sl0", "rs_a", "rs_srW"):
            m = re.search(rf"^\s*{key}:\s*([-\d.eE+]+)\s*$", t, re.M)
            if m and f"{key}_constant" not in f:
                f[f"{key}_constant"] = float(m.group(1))
        if "rs_muw" in t and "LuaMap" in t:
            f.update(_parse_rs_muw(t))
        if "Tnuc_s" in t and "LuaMap" in t:
            f.update(_parse_tnuc(t))

    # ---- parameters.par -----------------------------------------------------------
    par = d / "parameters.par"
    if par.is_file():
        t = par.read_text()
        active = "\n".join(ln for ln in t.splitlines()
                           if not ln.strip().startswith("!"))
        for key, cast in (("Plasticity", int), ("Tv", float), ("FreqCentral", float),
                          ("FreqRatio", float), ("EndTime", float), ("FL", int)):
            m = re.search(rf"^\s*{key}\s*=\s*([-\d.eE+]+)", active, re.M)
            f[key] = cast(m.group(1)) if m else UNKNOWN
    else:
        f["parameters_par"] = UNKNOWN

    # ---- plasticity angles, from the FIELD (the yaml prose may lie) ---------------
    pnc = next(d.glob("*plasticity*.nc"), None)
    if pnc:
        from deckbuild.asagi import read_asagi
        _, _, _, pf, _ = read_asagi(pnc, fields=["bulkFriction"])
        vals = np.unique(pf["bulkFriction"])
        angles = sorted(float(np.degrees(np.arctan(v))) for v in vals)
        f["phi_deg_from_field"] = [round(a, 4) for a in angles[:8]]
        f["phi_note"] = ("read from bulkFriction = tan(phi), NOT from the yaml comment; "
                         "a deck's prose has been wrong before")
    return f


def _parse_rs_muw(t: str) -> dict:
    """Recover the strike frame and the f_w design from the emitted Lua literals."""
    out = {}
    m = re.search(r'x\["x"\]\s*-\s*([-\d.eE+]+)\)\s*\*\s*([-\d.eE+]+)\s*\+\s*'
                  r'\(x\["y"\]\s*-\s*([-\d.eE+]+)\)\s*\*\s*([-\d.eE+]+)', t)
    if m:
        ox, cx, oy, cy = (float(g) for g in m.groups())
        out["strike_origin_xy"] = (ox, oy)
        # su = (sin az, cos az) -> az = atan2(cx, cy), degrees east of north
        out["strike_azimuth_deg"] = round(float(np.degrees(np.arctan2(cx, cy))) % 360.0, 6)
    base = re.search(r"return \{ rs_muw = ([-\d.eE+]+)", t)
    if base:
        out["fw_base"] = float(base.group(1))
    tr = re.findall(r"f_w ([-\d.eE+]+) -> ([-\d.eE+]+) across s = ([-\d.eE+]+) - "
                    r"([-\d.eE+]+) km", t)
    if tr:
        out["fw_values"] = [float(tr[0][0])] + [float(a[1]) for a in tr]
        out["fw_boundaries_s_km"] = [(float(a[2]), float(a[3])) for a in tr]
        out["fw_n_regions"] = len(out["fw_values"])
    return out


def _parse_tnuc(t: str) -> dict:
    out = {}
    xs = re.search(r'dx\s*=\s*x\["x"\]\s*-\s*([-\d.eE+]+)', t)
    ys = re.search(r'dy\s*=\s*x\["y"\]\s*-\s*([-\d.eE+]+)', t)
    zs = re.search(r'dz\s*=\s*x\["z"\]\s*([-+]\s*[\d.eE+]+)', t)
    if xs and ys:
        z = 0.0
        if zs:
            raw = zs.group(1).replace(" ", "")
            z = -float(raw[1:]) if raw[0] == "+" else float(raw)
        out["hypocenter_xyz"] = (float(xs.group(1)), float(ys.group(1)), z)
    r = re.search(r"R2\s*=\s*([-\d.eE+]+)\s*\*\s*([-\d.eE+]+)", t)
    if r:
        out["nucleation_radius_m"] = float(r.group(1))
    a = re.search(r"return \{ Tnuc_s = ([-\d.eE+]+)", t)
    if a:
        out["nucleation_amplitude"] = float(a.group(1))
    return out
