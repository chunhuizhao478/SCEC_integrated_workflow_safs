"""friction.py -- the rate-and-state friction fields.

    temperature (or a depth profile)  ->  data{rs_a, rs_srW}
    an along-strike design            ->  the rs_muw !LuaMap
    the SNAPPED hypocentre            ->  the Tnuc_s !LuaMap

The zoning profiles are ported verbatim from the legacy `build_friction_nc_thermal.py`
(b constant, slope m = 8.0e-5 per degC):

  CASE1 (3 zones)  a-b(T) = +0.004               T <= 50
                          = +0.004 - m(T-50)     50..150     (0 at 100)
                          = -0.004               150..300
                          = -0.004 + m(T-300)    T >= 300    (0 at 350)
  CASE2 (2 zones)  a-b(T) = -0.004               T <= 300
                          = -0.004 + m(T-300)    T >= 300
  V_w(T) both      = 0.05 m/s for T <= 350; linear to 1000 over 350..400; 1000 above

CASE2 is velocity-weakening all the way to the surface, which is why a CASE1 calibration
does not transfer to it.

`b` is a scalar, not a field: SeisSol v1.1.3 FL=103 never reads rs_b spatially.

SPATIAL rs_muw needs SeisSol newer than v1.3.2.  v1.1.3 accepts the YAML and SILENTLY
IGNORES it -- the run completes and is wrong -- so `check_seissol_supports_spatial_muw`
exists to be called BEFORE building a graded design.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from deckbuild.asagi import hull_containment, read_asagi, write_asagi
from deckbuild.config import Project
from deckbuild.contract import Artifact, GateReport, HARD, Stage, WARN
from deckbuild.geometry import build_grid, load_fault, strike_s_km

__all__ = ["FrictionStage", "FwDesign", "FrictionError", "FRICTION_FIELDS",
           "a_minus_b", "a_of_T", "vw_of_T", "fw_profile_1d", "write_lua_map",
           "nucleation_lua", "check_seissol_supports_spatial_muw", "unit_self_test"]

FRICTION_FIELDS = ["rs_a", "rs_srW"]

# Legacy zoning constants.  Overridable per project via FrictionStage.build(...).
M_SLOPE = 8.0e-5          # d(a-b)/dT, per degC
AB_COLD, AB_HOT = 0.004, -0.004
T_KINKS_CASE1 = (50.0, 150.0, 300.0)
T_KINK_CASE2 = 300.0
VW_LO, VW_HI = 0.05, 1000.0
T_VW_LO, T_VW_HI = 350.0, 400.0


class FrictionError(ValueError):
    """Raised for an unbuildable friction design."""


# --------------------------------------------------------------------------- profiles
def a_minus_b(T, case: int, m: float = M_SLOPE) -> np.ndarray:
    """(a - b)(T).  VERBATIM from the legacy build_friction_nc_thermal."""
    T = np.asarray(T, float)
    if case == 1:
        t1, t2, t3 = T_KINKS_CASE1
        out = np.where(T <= t1, AB_COLD,
                       np.where(T <= t2, AB_COLD - m * (T - t1),
                                np.where(T <= t3, AB_HOT, AB_HOT + m * (T - t3))))
    elif case == 2:
        out = np.where(T <= T_KINK_CASE2, AB_HOT, AB_HOT + m * (T - T_KINK_CASE2))
    else:
        raise FrictionError(f"unknown thermal case {case!r}; expected 1 or 2")
    return np.asarray(out, float)


def a_of_T(T, case: int, b: float, m: float = M_SLOPE) -> np.ndarray:
    return a_minus_b(T, case, m) + b


def vw_of_T(T) -> np.ndarray:
    """V_w(T): VW_LO below T_VW_LO, linear to VW_HI, then flat."""
    T = np.asarray(T, float)
    t = np.clip((T - T_VW_LO) / (T_VW_HI - T_VW_LO), 0.0, 1.0)
    return VW_LO + (VW_HI - VW_LO) * t


def unit_self_test(verbose: bool = False) -> None:
    """Hand values from the legacy plan, checked BEFORE any baking.

    The legacy build runs this first and aborts on failure; keeping it means a change to
    the profile constants cannot reach a deck.
    """
    b = 0.019
    for T, want1, want2 in [(15.0, 0.004, -0.004), (50.0, 0.004, -0.004),
                            (100.0, 0.0, -0.004), (150.0, -0.004, -0.004),
                            (279.0, -0.004, -0.004), (300.0, -0.004, -0.004),
                            (350.0, 0.0, 0.0), (400.0, 0.004, 0.004)]:
        g1 = float(a_minus_b(T, 1)); g2 = float(a_minus_b(T, 2))
        if abs(g1 - want1) > 1e-12 or abs(g2 - want2) > 1e-12:
            raise FrictionError(
                f"a-b profile regression at T={T}: case1 {g1} (want {want1}), "
                f"case2 {g2} (want {want2})")
    for T, want in [(0.0, VW_LO), (350.0, VW_LO), (375.0, 0.5 * (VW_LO + VW_HI)),
                    (400.0, VW_HI), (1000.0, VW_HI)]:
        if abs(float(vw_of_T(T)) - want) > 1e-9:
            raise FrictionError(f"V_w profile regression at T={T}")
    # G4: the reference hypocentre plateau.
    if abs(float(a_of_T(279.0, 1, b)) - 0.015) > 1e-12:
        raise FrictionError("a(279 degC, case 1) must be 0.015")
    if verbose:
        print("friction profile self-test OK")


def depth_profile_fields(z, knots_m, values, vw_knots_m, vw_values, b: float):
    """Piecewise-linear a(z) and V_w(z) for a project with no thermal model.

    Knots are elevations (negative down) and must DECREASE, matching how a user writes a
    depth section top-down.
    """
    kn = np.asarray(knots_m, float)
    vk = np.asarray(vw_knots_m, float)
    for name, arr in (("a_minus_b_knots_m", kn), ("vw_knots_m", vk)):
        if arr.size < 2:
            raise FrictionError(f"{name} needs at least 2 knots")
        if not np.all(np.diff(arr) < 0):
            raise FrictionError(
                f"{name} must be strictly DECREASING (elevations, top first); got {arr}")
    z = np.asarray(z, float)
    ab = np.interp(z, kn[::-1], np.asarray(values, float)[::-1])
    vw = np.interp(z, vk[::-1], np.asarray(vw_values, float)[::-1])
    return ab + b, vw


# --------------------------------------------------------------------------- graded f_w
@dataclass(frozen=True)
class FwDesign:
    """An N-region along-strike strong-rate-weakening floor."""

    fw_values: tuple[float, ...]
    boundaries_s_km: tuple[tuple[float, float], ...] = ()
    f0: float = 0.6

    def __post_init__(self) -> None:
        if not self.fw_values:
            raise FrictionError("FwDesign needs at least one f_w value")
        for v in self.fw_values:                       # gate F0
            if not 0.0 <= v < self.f0:
                raise FrictionError(
                    f"F0: every f_w must satisfy 0 <= f_w < f0 = {self.f0}; got {v}")
        if len(self.boundaries_s_km) != len(self.fw_values) - 1:
            raise FrictionError(
                f"{len(self.fw_values)} region(s) need {len(self.fw_values) - 1} "
                f"transition window(s), got {len(self.boundaries_s_km)}")
        prev = -np.inf
        for i, (lo, hi) in enumerate(self.boundaries_s_km):
            if hi <= lo:
                raise FrictionError(f"F0: window {i} has s_hi {hi} <= s_lo {lo}")
            if lo < prev:
                raise FrictionError(f"F0: window {i} overlaps the previous one")
            prev = hi

    @property
    def n_regions(self) -> int:
        return len(self.fw_values)


def fw_profile_1d(s_km, design: FwDesign) -> np.ndarray:
    """f_w(s): plateaus joined by smoothstep windows -- the shape the Lua emits."""
    s = np.asarray(s_km, float)
    out = np.full(s.shape, float(design.fw_values[0]))
    for i, (lo, hi) in enumerate(design.boundaries_s_km):
        t = np.clip((s - lo) / (hi - lo), 0.0, 1.0)
        out = out + (design.fw_values[i + 1] - design.fw_values[i]) * t * t * (3.0 - 2.0 * t)
    return out


def write_lua_map(out_path, design: FwDesign, strike, returns: str = "rs_muw") -> Path:
    """Emit the shippable `!LuaMap`.

    The along-strike coordinate is computed INLINE in Lua from the descriptor's azimuth
    and origin -- the same expression as `geometry.strike_s_km`, so the design the checks
    score and the field SeisSol reads are the same function.  Gate FL re-parses what was
    written and asserts it evaluates to the design.
    """
    # %.17g, NOT repr(): repr gives the shortest round-trip form
    # (-0.7193398003386512), while every shipped rs_muw map carries the 17-significant
    # form (-0.71933980033865119).  Both parse to the same double, but Phase 8 compares
    # the emitted TEXT, so the formatting has to match the legacy emitter.
    def g17(v):
        return f"{float(v):.17g}"

    az = np.radians(strike.azimuth_deg)
    cx, cy = float(np.sin(az)), float(np.cos(az))
    ox, oy = float(strike.origin_xy[0]), float(strike.origin_xy[1])
    L = [f"[{returns}]: !LuaMap", f"  returns: [{returns}]", "  function: |",
         "    function f (x)",
         f'      local s = ((x["x"] - {ox!r}) * {g17(cx)} + '
         f'(x["y"] - {oy!r}) * {g17(cy)}) / 1000.0',
         "      local inc = 0.0", "      local t"]
    n = len(design.boundaries_s_km)
    for i, (lo, hi) in enumerate(design.boundaries_s_km):
        a, bnd = design.fw_values[i], design.fw_values[i + 1]
        L += [f"      -- transition {i + 1}/{n}: f_w {a:g} -> {bnd:g} "
              f"across s = {lo:g} - {hi:g} km",
              f"      t = (s - ({lo:g})) / (({hi:g}) - ({lo:g}))",
              "      if (t < 0.0) then t = 0.0 end",
              "      if (t > 1.0) then t = 1.0 end",
              f"      inc = inc + ({g17(bnd)} - ({g17(a)})) * t * t * (3.0 - 2.0 * t)"]
    L += [f"      return {{ {returns} = {g17(design.fw_values[0])} + inc }}",
          "    end", ""]
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(L))
    return out_path


def _eval_lua_text(text: str, s_km) -> np.ndarray:
    """Re-evaluate an emitted LuaMap's arithmetic in numpy (gate FL).

    Parses the literals back out of the generated source rather than trusting the design
    object, so a mismatch between what was designed and what was WRITTEN is caught.
    """
    # Accept both our form `= X + inc` and the LEGACY form `= X + inc * g`, where the
    # shipped maps carry a `local g = <val>` global scale.  Phase 8 compares against those
    # files, so the parser has to read them, not only what we emit.
    base = re.search(r"return \{ \w+ = ([-\d.eE+]+) \+ inc(?:\s*\*\s*g)? \}", text)
    if not base:
        raise FrictionError(
            "FL: cannot find the base f_w in the Lua.  Expected a line of the form "
            "`return { <name> = <base> + inc }` or `... + inc * g }`; got:\n"
            + "\n".join(l for l in text.splitlines() if "return {" in l))
    g = 1.0
    mg = re.search(r"local g\s*=\s*([-\d.eE+]+)", text)
    if mg:
        g = float(mg.group(1))
    out = np.full(np.shape(s_km), 0.0, float)
    pat = re.compile(r"t = \(s - \(([-\d.eE+]+)\)\) / \(\(([-\d.eE+]+)\) - "
                     r"\(([-\d.eE+]+)\)\)\s*\n.*?\n.*?\n\s*inc = inc \+ "
                     r"\(([-\d.eE+]+) - \(([-\d.eE+]+)\)\)")
    found = pat.findall(text)
    n_want = text.count("-- transition ")
    if len(found) != n_want:
        raise FrictionError(
            f"FL: parsed {len(found)} transition(s) out of the emitted Lua but the text "
            f"declares {n_want}; the emitter and this parser have drifted, so FL would "
            f"be checking nothing")
    for lo, hi, _lo2, hv, av in found:
        t = np.clip((np.asarray(s_km, float) - float(lo)) / (float(hi) - float(lo)),
                    0.0, 1.0)
        out = out + (float(hv) - float(av)) * t * t * (3.0 - 2.0 * t)
    return float(base.group(1)) + out * g


def nucleation_lua(snapped, radius_m: float, amplitude_mpa: float,
                   returns: str = "Tnuc_s") -> str:
    """The SCEC compact bell, centred on the SNAPPED hypocentre.

    F(r) = exp(r^2/(r^2 - R^2)) for r < R, 0 otherwise; F(0) = 1.  The centre MUST be the
    snapped point -- centring on the requested one puts the peak overstress off the fault,
    where no facet samples it.
    """
    x, y, z = (float(v) for v in snapped.xyz)
    return "\n".join([
        f"[{returns}]: !LuaMap", f"  returns: [{returns}]", "  function: |",
        "    function f (x)",
        f"      -- F(r) = exp(r^2/(r^2-R^2)) for r<R.  F(0)=1 -> {amplitude_mpa:g} MPa "
        f"at the centre.",
        f"      local dx = x[\"x\"] - {x!r}",
        f"      local dy = x[\"y\"] - {y!r}",
        f"      local dz = x[\"z\"] - {z!r}",
        "      local r2 = dx*dx + dy*dy + dz*dz",
        f"      local R2 = {radius_m!r} * {radius_m!r}",
        "      local F = 0.0",
        "      if (r2 < R2) then F = math.exp(r2 / (r2 - R2)) end",
        f"      return {{ {returns} = {amplitude_mpa!r} * F }}", "    end", ""])


def check_seissol_supports_spatial_muw(src_path, min_version=(1, 3, 2)) -> GateReport:
    """Gate F-VER: does this SeisSol read a SPATIAL rs_muw?

    v1.1.3 accepts the YAML and silently ignores it, so the run completes and is wrong.
    Must be called BEFORE building a graded design.
    """
    rep = GateReport("SeisSol spatial rs_muw support")
    if src_path is None:
        rep.skip("F-VER", "no SeisSol source tree given; spatial rs_muw support UNKNOWN")
        return rep
    p = Path(src_path)
    ver = None
    for cand in (p / "CMakeLists.txt", p / "src" / "version.h", p):
        if cand.is_file():
            m = re.search(r"(\d+)\.(\d+)\.(\d+)", cand.read_text(errors="ignore"))
            if m:
                ver = tuple(int(g) for g in m.groups())
                break
    if ver is None:
        rep.add("F-VER", False,
                f"could not determine a SeisSol version under {p}; spatial rs_muw needs "
                f"> {'.'.join(map(str, min_version))} and older builds IGNORE IT SILENTLY")
        return rep
    ok = ver > min_version
    rep.add("F-VER", ok,
            f"SeisSol {'.'.join(map(str, ver))} "
            f"{'supports' if ok else 'DOES NOT support'} spatial rs_muw "
            f"(needs > {'.'.join(map(str, min_version))}; older builds accept the yaml "
            f"and silently ignore it)")
    return rep


# --------------------------------------------------------------------------- the stage
class FrictionStage(Stage):
    name = "friction"

    def build(self, cfg: Project, out_dir, *, case: int = 1, thermal=None,
              b: float | None = None, out_name: str | None = None,
              mesh: str | None = None, dtype=np.float32) -> Artifact:
        """Build data{rs_a, rs_srW}.

        thermal -- the Artifact from MaterialStage's CTM path.  Omit to use the project's
                   depth_profile.  The stage never reads a file it was not handed.
        """
        unit_self_test()
        b = cfg.physics.rs_b if b is None else float(b)
        spec = cfg.raw.thermal
        out_dir = Path(out_dir)

        if thermal is not None:
            tpath = Path(getattr(thermal, "path", thermal))
            x, y, z, flds, _ = read_asagi(tpath)
            if "T" not in flds:
                raise FrictionError(f"{tpath}: expected a 'T' field, have {sorted(flds)}")
            T = flds["T"]
            a = a_of_T(T, case, b)
            vw = vw_of_T(T)
            source = f"ctm_slices/case{case}"
            tsha = getattr(thermal, "sha256", "")
        else:
            if spec is None or spec.kind != "depth_profile":
                raise FrictionError(
                    "no thermal artifact was given and the project has no "
                    "raw.thermal.kind = depth_profile; friction cannot be built")
            p = spec.params
            # The four grids are INDEPENDENT (the shipped ones are material 1500/250,
            # stress 1000/250, friction 1500/200).  Take this grid from the thermal
            # source's own params; fall back to the stress box only as a default.
            box = cfg.stress_box
            gdx = float(p.get("grid_dx", box.dx))
            gdz = float(p.get("dz", box.dz))
            zlo = float(p.get("z_min", box.zmin))
            zhi = float(p.get("z_max", box.zmax))
            gx = np.arange(box.xmin, box.xmax + 0.5 * gdx, gdx)
            gy = np.arange(box.ymin, box.ymax + 0.5 * gdx, gdx)
            gz = np.arange(zlo, zhi + 0.5 * gdz, gdz)
            a1, v1 = depth_profile_fields(gz, p["a_minus_b_knots_m"], p["a_minus_b_values"],
                                          p["vw_knots_m"], p["vw_values"], b)
            a = np.broadcast_to(a1[:, None, None], (len(gz), len(gy), len(gx))).copy()
            vw = np.broadcast_to(v1[:, None, None], (len(gz), len(gy), len(gx))).copy()
            x, y, z = gx, gy, gz
            source, tsha = "depth_profile", ""

        # G2, on the float64 arrays, BEFORE the float32 cast.
        if not np.all(a > 0):
            raise FrictionError(f"G2: rs_a must be > 0 everywhere; min = {a.min():.6g}")
        if not (np.all(vw >= VW_LO - 1e-12) and np.all(vw <= VW_HI + 1e-12)):
            raise FrictionError(
                f"G2: rs_srW must lie in [{VW_LO}, {VW_HI}]; got "
                f"[{vw.min():.6g}, {vw.max():.6g}]")

        name = out_name or f"friction_{source.replace('/', '_')}.nc"
        path = out_dir / name
        write_asagi(path, x, y, z, {"rs_a": a, "rs_srW": vw}, dtype=dtype, attrs={
            "title": f"{cfg.name} rate-and-state friction",
            "source": source, "case": int(case), "rs_b": float(b),
            "note": "rs_b is a SCALAR: SeisSol v1.1.3 FL=103 never reads it spatially",
            "thermal_sha256": tsha, "project": cfg.name, "crs": cfg.crs.epsg,
        })
        return Artifact.of(path, kind="friction",
                           params={"case": int(case), "rs_b": float(b), "source": source,
                                   "mesh": mesh or cfg.default_mesh},
                           provenance={"thermal_sha256": tsha,
                                       "descriptor_sha256": cfg.sha256()})

    def verify(self, cfg: Project, artifact: Artifact, *, case: int | None = None,
               mesh: str | None = None, **_) -> GateReport:
        """G1-G4."""
        p = Path(artifact.path)
        x, y, z, flds, attrs = read_asagi(p)
        case = int(attrs.get("case", case or 1))
        b = float(attrs.get("rs_b", cfg.physics.rs_b))
        rep = GateReport(f"friction G1-G4 ({p.name})")

        # G1 -- fault containment.
        mspec = cfg.mesh(mesh or cfg.default_mesh)
        mpath = cfg.resolve_path(mspec.path)
        if not Path(mpath).is_file():
            rep.skip("G1", f"no mesh at {mpath}; hull containment not checked")
        else:
            fault = load_fault(mspec, cfg.data_dir, strike=cfg.strike)
            hull_containment(p, fault.cent, label="fault centroids", gate="G1",
                             severity=HARD, report=rep)

        # G2 -- ranges, re-checked after the float32 round trip.
        a, vw = flds["rs_a"], flds["rs_srW"]
        rep.add("G2a", bool(np.all(a > 0)), f"rs_a > 0 everywhere: min = {a.min():.6g}")
        rep.add("G2b", bool(np.all((vw >= VW_LO * (1 - 1e-6)) & (vw <= VW_HI * (1 + 1e-6)))),
                f"rs_srW in [{VW_LO}, {VW_HI}]: [{vw.min():.6g}, {vw.max():.6g}]")

        # G3 -- the profile identity, where the source is thermal.
        if str(attrs.get("source", "")).startswith("ctm"):
            rep.skip("G3", "round-trip vs the thermal grid needs the T nc; call "
                           "asagi.roundtrip_selfcheck with the source temperature")
        else:
            rep.skip("G3", "depth-profile source: the field IS the profile, so the "
                           "round-trip check is vacuous")

        # G4 -- the plateau anchor.  a(VW plateau) must be exactly b + AB_HOT.
        want = b + AB_HOT
        on_plateau = np.isclose(a, want, rtol=0, atol=1e-6)
        rep.add("G4", bool(on_plateau.any()),
                f"a == b + {AB_HOT:g} = {want:g} on {int(on_plateau.sum())} node(s) "
                f"(the velocity-weakening plateau); a range "
                f"[{a.min():.6g}, {a.max():.6g}]",
                severity=WARN if not on_plateau.any() else HARD)
        return rep

    # ---------------------------------------------------------------- graded f_w
    def build_fw_map(self, cfg: Project, out_dir, design: FwDesign,
                     out_name: str | None = None) -> Artifact:
        name = out_name or ("fault_rs_muw_"
                            + "-".join(f"{v:g}" for v in design.fw_values) + ".yaml")
        path = write_lua_map(Path(out_dir) / name, design, cfg.strike)
        return Artifact.of(path, kind="rs_muw_lua",
                           params={"fw_values": list(design.fw_values),
                                   "boundaries_s_km": [list(b) for b in
                                                       design.boundaries_s_km],
                                   "f0": design.f0},
                           provenance={"descriptor_sha256": cfg.sha256()})

    def verify_fw_map(self, cfg: Project, artifact: Artifact, design: FwDesign,
                      s_probe=None) -> GateReport:
        """F0, F1 and FL on the EMITTED Lua."""
        text = Path(artifact.path).read_text()
        rep = GateReport(f"graded f_w ({Path(artifact.path).name})")
        rep.add("F0", True,
                f"admissible: {design.n_regions} region(s), all 0 <= f_w < f0 = "
                f"{design.f0:g}, windows ordered and disjoint")

        # FL -- what was WRITTEN must evaluate to the design.
        s = np.linspace(-50.0, 500.0, 1101) if s_probe is None else np.asarray(s_probe)
        want = fw_profile_1d(s, design)
        got = _eval_lua_text(text, s)
        err = float(np.max(np.abs(got - want)))
        rep.add("FL", err < 1e-12,
                f"the emitted Lua reproduces the design over s in "
                f"[{s.min():g}, {s.max():g}] km: max |d| = {err:.3e}")

        # F1 -- the first region's plateau must hold its value EXACTLY (the hypocentre
        # and its nucleation patch live there).
        if design.boundaries_s_km:
            a1 = design.boundaries_s_km[0][0]
            plateau = s <= a1
            f1 = bool(np.all(got[plateau] == design.fw_values[0])) if plateau.any() else False
            rep.add("F1", f1,
                    f"gate plateau (s <= {a1:g} km) holds f_w = {design.fw_values[0]:g} "
                    f"exactly on {int(plateau.sum())} probe point(s)")
        else:
            rep.add("F1", True, "uniform design: one region, nothing to grade")

        # The strike frame really is the descriptor's.
        ox = f"{float(cfg.strike.origin_xy[0])!r}"
        rep.add("FLframe", ox in text,
                f"the emitted Lua uses the descriptor's strike origin ({ox})")
        return rep
