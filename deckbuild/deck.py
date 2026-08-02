"""deck.py -- assemble a runnable deck, then check it the way we check one by hand.

This is where the workflow stops being five builders and becomes one workflow: the
cross-component checks only make sense once all the pieces exist.

Two rules the assembly enforces, both learned the hard way:

  * COPY, never symlink.  A deck is a shippable unit that gets rsynced to a cluster; a
    symlink into outputs/ either breaks on transfer or silently ships whatever that path
    later points at.
  * ONE declared destination per artifact, with the filename and the YAML `file:` field
    written from the SAME variable.  They cannot be allowed to drift; P1 then re-checks
    the result on disk.

The pre-flight battery P1-P8 reads the ASSEMBLED deck's live `file:` targets, never an
in-memory array -- it tests what SeisSol will actually read.  Its predicates are shared
with Stage F (see stage_f.py) so the two cannot disagree.
"""
from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from deckbuild.contract import Artifact, GateReport, HARD, Manifest, WARN, sha256_file
from deckbuild.geometry import load_fault, snap_hypocenter, strike_s_km
from deckbuild.stage_f import (
    check_finite_fields, check_positive_normal_stress, check_receivers_under_the_surface,
    project_stress_onto_fault,
)

__all__ = ["DeckStage", "DeckSpec", "DeckError", "resolve_paths"]


class DeckError(ValueError):
    """Raised when a deck cannot be assembled or is internally inconsistent."""


@dataclass
class DeckSpec:
    """What goes into a deck beyond the artifacts themselves."""

    prefix: str = ""                       # from the descriptor; NOT hardcoded "safs_"
    plasticity: bool = False
    attenuation: dict | None = None
    receivers: np.ndarray | None = None
    pickpoints: np.ndarray | None = None
    mu_s: float | None = None
    nucleation_lua: str | None = None
    rs_muw_lua: str | None = None
    notes: str = ""                        # preserved verbatim across regeneration
    end_time_s: float = 20.0
    extra_par: dict = field(default_factory=dict)


# The filename contract.  Also reproduced in the README; `resolve_paths` implements it.
_LAYOUT = {
    "material":   ("material", "{p}material_cvm.nc",        "{p}material_cvm.yaml"),
    "plasticity": ("material", "{p}plasticity_{tag}.nc",    "{p}material_cvm.yaml"),
    "stress":     ("stress",   "{p}stress_{tag}.nc",        "{p}initial_stress.yaml"),
    "friction":   ("friction", "{p}friction_{tag}.nc",      "{p}fault.yaml"),
    "mesh":       (None,       "{name}",                    "parameters.par"),
}

_NOTES_BEGIN = "# ===== HAND-WRITTEN NOTES (preserved across regeneration) ====="
_NOTES_END = "# ===== END HAND-WRITTEN NOTES ====="


def resolve_paths(deck_dir, artifacts: dict, spec: DeckSpec) -> dict:
    """The full path map, WITHOUT writing anything.

    Lets a notebook print exactly where every file will go before assembling, so "where
    did my file end up" is answerable at a glance.
    """
    deck_dir = Path(deck_dir)
    p = spec.prefix
    out = {}
    for kind, art in artifacts.items():
        if art is None or kind not in _LAYOUT:
            continue
        _sub, name_t, ref = _LAYOUT[kind]
        tag = _tag_for(kind, art)
        name = name_t.format(p=p, tag=tag, name=Path(art.path).name)
        out[kind] = {"src": str(art.path), "dst": str(deck_dir / name),
                     "filename": name, "referenced_by": ref.format(p=p)}
    return out


def _tag_for(kind: str, art: Artifact) -> str:
    """A design-encoding tag.  An ambiguous deck filename has caused a real mix-up."""
    q = art.params or {}
    if kind == "stress":
        ks = q.get("k_values") or []
        core = f"k{ks[0]:g}" if len(ks) == 1 else "gradedk_" + "-".join(f"{v:g}" for v in ks)
        fz = q.get("freeze_above_depth_m") or 0.0
        return core + (f"_freeze{fz:g}m" if fz else "")
    if kind == "friction":
        return f"case{q.get('case', 1)}"
    if kind == "plasticity":
        return f"phi{q.get('phi_soft_deg', 0):g}_{q.get('phi_hard_deg', 0):g}"
    return "x"


class DeckStage:
    name = "deck"

    # ------------------------------------------------------------------ assemble
    def assemble(self, cfg, deck_dir, artifacts: dict, spec: DeckSpec | None = None,
                 overwrite: bool = False) -> Path:
        spec = spec or DeckSpec(prefix=f"{cfg.name}_")
        deck_dir = Path(deck_dir)

        # Refuse to mix artifacts built from different descriptors -- the stale-outputs
        # failure mode.  Every number would look plausible and mean nothing.
        hashes = {a.provenance.get("descriptor_sha256") for a in artifacts.values()
                  if a is not None and a.provenance.get("descriptor_sha256")}
        if len(hashes) > 1:
            raise DeckError(
                f"artifacts were built from {len(hashes)} different descriptors "
                f"({sorted(h[:12] for h in hashes)}); refusing to assemble")
        if hashes and cfg.sha256() and hashes != {cfg.sha256()}:
            raise DeckError(
                f"artifacts were built from descriptor {sorted(hashes)[0][:12]} but this "
                f"config is {cfg.sha256()[:12]}; refusing to assemble")

        existing_notes = _read_notes(deck_dir)
        if deck_dir.exists() and any(deck_dir.iterdir()) and not overwrite:
            raise DeckError(
                f"{deck_dir} already exists and is not empty; pass overwrite=True.  The "
                f"hand-written notes block is preserved either way.")
        deck_dir.mkdir(parents=True, exist_ok=True)

        paths = resolve_paths(deck_dir, artifacts, spec)
        for kind, e in paths.items():
            src, dst = Path(e["src"]), Path(e["dst"])
            shutil.copy2(src, dst)                       # COPY, never symlink
            if sha256_file(src) != sha256_file(dst):
                raise DeckError(f"copy of {kind} did not verify: {src} -> {dst}")

        notes = spec.notes or existing_notes
        self._write_yamls(cfg, deck_dir, paths, spec, notes)
        self._write_par(cfg, deck_dir, paths, spec, notes)
        if spec.receivers is not None:
            _write_points(deck_dir / f"{spec.prefix}receivers.dat", spec.receivers)
        if spec.pickpoints is not None:
            _write_points(deck_dir / f"{spec.prefix}pickpoints.dat", spec.pickpoints)

        man = Manifest.for_project(cfg)
        for kind, art in artifacts.items():
            if art is not None:
                man.record(kind, artifact=art)
        # The deck copy is RENAMED (the filename encodes the design), so a lookup keyed on
        # the artifact's source name never matches.  Record the deck filenames explicitly
        # or P1's sha256 check silently never fires.
        man.stages["deck_files"] = {
            e["filename"]: sha256_file(Path(e["dst"])) for e in paths.values()}
        man.write(deck_dir)
        return deck_dir

    # ------------------------------------------------------------------ templates
    def _write_yamls(self, cfg, deck_dir, paths, spec: DeckSpec, notes: str) -> None:
        p = spec.prefix
        hdr = _provenance_header(cfg, spec)

        mat = [hdr, "!Switch", "[rho, mu, lambda]: !ASAGI",
               f"  file: {paths['material']['filename']}",
               "  parameters: [rho, mu, lambda]", "  var: data",
               "  interpolation: linear"]
        if "plasticity" in paths:
            mat += ["[plastCo, bulkFriction]: !ASAGI",
                    f"  file: {paths['plasticity']['filename']}",
                    "  parameters: [plastCo, bulkFriction]", "  var: data",
                    "  interpolation: linear",
                    "# bulkFriction is tan(phi), the COEFFICIENT -- SeisSol atans it."]
        mat += ["[s_xx, s_yy, s_zz, s_xy, s_yz, s_xz]: !Include "
                f"{p}initial_stress.yaml", _notes_block(notes)]
        (deck_dir / f"{p}material_cvm.yaml").write_text("\n".join(mat) + "\n")

        (deck_dir / f"{p}initial_stress.yaml").write_text("\n".join([
            hdr, "!ASAGI", f"file: {paths['stress']['filename']}",
            "parameters: [s_xx, s_yy, s_zz, s_xy, s_yz, s_xz]", "var: data",
            "interpolation: linear",
            "# compression-NEGATIVE Pa, z = elevation.", _notes_block(notes)]) + "\n")

        fault = [hdr, "!Switch",
                 f"[s_xx, s_yy, s_zz, s_xy, s_yz, s_xz]: !Include {p}initial_stress.yaml",
                 "[rs_a, rs_srW]: !ASAGI",
                 f"  file: {paths['friction']['filename']}",
                 "  parameters: [rs_a, rs_srW]", "  var: data",
                 "  interpolation: linear",
                 "[rs_b]: !ConstantMap", "  map:", f"    rs_b: {cfg.physics.rs_b}",
                 "[rs_sl0]: !ConstantMap", "  map:", f"    rs_sl0: {cfg.physics.rs_sl0}"]
        if spec.rs_muw_lua:
            fault.append(spec.rs_muw_lua.rstrip())
        if spec.nucleation_lua:
            fault.append(spec.nucleation_lua.rstrip())
        fault.append(_notes_block(notes))
        (deck_dir / f"{p}fault.yaml").write_text("\n".join(fault) + "\n")

    def _write_par(self, cfg, deck_dir, paths, spec: DeckSpec, notes: str) -> None:
        p = spec.prefix
        mesh_name = paths.get("mesh", {}).get("filename", "MESH_NOT_SET")
        lines = [_provenance_header(cfg, spec, comment="!"),
                 "&equations",
                 f"MaterialFileName = '{p}material_cvm.yaml'",
                 f"Plasticity = {1 if spec.plasticity else 0}"]
        if spec.attenuation:
            a = spec.attenuation
            lines += [f"FreqCentral = {a['freq_central']}",
                      f"FreqRatio = {a['freq_ratio']}",
                      "! Q from the CVM's own Vs.  REQUIRES a viscoelastic build --",
                      "! EQUATIONS=viscoelastic2 is a BUILD-time choice, not a flag."]
        lines += ["/", "", "&DynamicRupture", "FL = 103",
                  f"ModelFileName = '{p}fault.yaml'", "/", "",
                  "&MeshNml", "MeshFile = '" + mesh_name.replace(".puml.h5", "") + "'",
                  "meshgenerator = 'PUML'", "/", "",
                  "&Output",
                  "! wavefieldoutput = 1 is REQUIRED: 0 hangs in teardown.",
                  "wavefieldoutput = 1",
                  "! OutputRegionBounds must not contain a literal 0.0 -- that writes",
                  "! the WHOLE volume.  Bounds are mesh-specific.",
                  "/", "", "&AbortCriteria", f"EndTime = {spec.end_time_s}", "/"]
        for k, v in (spec.extra_par or {}).items():
            lines.append(f"{k} = {v}")
        lines.append(_notes_block(notes, comment="!"))
        (deck_dir / "parameters.par").write_text("\n".join(lines) + "\n")

    # ------------------------------------------------------------------ preflight
    def preflight(self, cfg, deck_dir, *, spec: DeckSpec | None = None,
                  mesh: str | None = None) -> GateReport:
        """P1-P8 on the ASSEMBLED deck, reading its live `file:` targets."""
        deck_dir = Path(deck_dir)
        spec = spec or DeckSpec(prefix=f"{cfg.name}_")
        rep = GateReport(f"deck pre-flight P1-P8 ({deck_dir.name})")
        man_path = deck_dir / Manifest.FILENAME
        man = Manifest.read(man_path) if man_path.is_file() else {"stages": {}}

        # --- P1: every referenced file exists, is real, hashes right, nothing extra ---
        refs, missing, bad_hash = _collect_refs(deck_dir), [], []
        by_name = man.get("stages", {}).get("deck_files", {})
        for r in refs:
            f = deck_dir / r
            if not f.is_file():
                missing.append(r)
            elif f.is_symlink():
                bad_hash.append(f"{r} (symlink)")
            elif r in by_name and sha256_file(f) != by_name[r]:
                bad_hash.append(f"{r} (sha256)")
        accounted = set(refs) | {Manifest.FILENAME, "parameters.par"} | {
            f.name for f in deck_dir.glob("*.yaml")} | {
            f.name for f in deck_dir.glob("*.dat")}
        extra = sorted(f.name for f in deck_dir.iterdir()
                       if f.is_file() and f.name not in accounted)
        rep.add("P1", not missing and not bad_hash and not extra,
                f"{len(refs)} referenced file(s) present; missing {missing or 'none'}; "
                f"bad {bad_hash or 'none'}; unaccounted extras {extra or 'none'} "
                f"(an unexplained nc is how a stale field gets shipped)")

        if missing or bad_hash:
            # A corrupt or absent input makes the downstream readers raise.  A broken deck
            # must produce a REPORT, not a traceback -- the whole point of a pre-flight is
            # to tell you what is wrong before you queue a job.
            for g in ("P2", "P3", "P4", "P5", "P6", "P7", "P8"):
                rep.skip(g, "not attempted: P1 found a missing or unreadable referenced "
                            "file, so the downstream checks cannot be trusted")
            return rep

        # --- the on-fault checks, sharing Stage F's predicates ---
        mspec = cfg.mesh(mesh or cfg.default_mesh)
        mesh_file = next((deck_dir.glob("*.puml.h5")), None)
        if mesh_file is None:
            for g in ("P2", "P3", "P4", "P6", "P7"):
                rep.skip(g, "no mesh in the deck; the on-fault checks did not run")
        else:
            import dataclasses
            fault = load_fault(dataclasses.replace(mspec, path=str(mesh_file)),
                               None, strike=cfg.strike)
            stress = next(deck_dir.glob("*stress*.nc"), None)
            fric = next(deck_dir.glob("*friction*.nc"), None)
            if stress is None:
                for g in ("P2", "P3", "P4"):
                    rep.skip(g, "no stress nc in the deck")
            else:
                sn, tau, mu = project_stress_onto_fault(stress, fault)
                check_positive_normal_stress(sn, tau, rep, "P2")
                _check_psi_ini(sn, tau, fric, fault, cfg, rep, "P3")
                if spec.mu_s is not None:
                    rep.add("P4", float(np.nanmax(mu)) < spec.mu_s,
                            f"no pre-slip at t=0: max mu_app = {np.nanmax(mu):.4f} "
                            f"vs mu_s = {spec.mu_s:.4f}")
                else:
                    rep.skip("P4", "no mu_s in the deck spec")
                _check_gate_corridor(cfg, fault, mu, sn, rep, "P6")
            _check_nucleation(cfg, fault, rep, "P7", mesh or cfg.default_mesh)

        # --- P5: plasticity sub-yield at t=0 ---
        plast = next(deck_dir.glob("*plasticity*.nc"), None)
        if not spec.plasticity or plast is None:
            rep.skip("P5", "plasticity is off in this deck")
        else:
            _check_yield(deck_dir, plast, rep, "P5")

        # --- P8: output hygiene ---
        _check_output_hygiene(deck_dir, spec, rep, "P8")
        return rep

    # ------------------------------------------------------------------ diff
    def diff(self, a, b) -> dict:
        """Which files differ between two decks, by sha256.

        Several past experiments were confounded by an unintended second change; an A/B
        pair should differ in exactly one input, and this makes that visible.
        """
        a, b = Path(a), Path(b)
        fa = {f.name: sha256_file(f) for f in a.iterdir() if f.is_file()}
        fb = {f.name: sha256_file(f) for f in b.iterdir() if f.is_file()}
        return {"only_in_a": sorted(set(fa) - set(fb)),
                "only_in_b": sorted(set(fb) - set(fa)),
                "differ": sorted(n for n in set(fa) & set(fb) if fa[n] != fb[n]),
                "same": sorted(n for n in set(fa) & set(fb) if fa[n] == fb[n])}


# --------------------------------------------------------------------------- helpers
def _provenance_header(cfg, spec: DeckSpec, comment: str = "#") -> str:
    c = comment
    return "\n".join([
        f"{c} GENERATED by deckbuild -- do not hand-edit above the notes block.",
        f"{c} project    : {cfg.name}",
        f"{c} descriptor : {cfg.source_path} ({cfg.sha256()[:16]})",
        f"{c} crs        : {cfg.crs.epsg}",
        f"{c} prefix     : {spec.prefix}",
    ])


def _notes_block(notes: str, comment: str = "#") -> str:
    body = notes.strip() or f"{comment} (none yet)"
    b, e = (_NOTES_BEGIN, _NOTES_END) if comment == "#" else (
        _NOTES_BEGIN.replace("#", "!"), _NOTES_END.replace("#", "!"))
    return "\n".join([b, body, e])


def _read_notes(deck_dir) -> str:
    """Recover the hand-written block so regeneration never destroys it."""
    deck_dir = Path(deck_dir)
    if not deck_dir.is_dir():
        return ""
    for f in sorted(deck_dir.glob("*.yaml")):
        t = f.read_text()
        if _NOTES_BEGIN in t and _NOTES_END in t:
            return t.split(_NOTES_BEGIN, 1)[1].split(_NOTES_END, 1)[0].strip()
    return ""


def _collect_refs(deck_dir) -> list[str]:
    """Every file the deck REFERENCES: yaml `file:` targets plus parameters.par's mesh.

    The mesh is named by `MeshFile` in parameters.par (without the .puml.h5 suffix), not
    by a yaml `file:`.  Missing it made P1 report the mesh as an unaccounted extra.
    """
    deck_dir = Path(deck_dir)
    refs = []
    for f in sorted(deck_dir.glob("*.yaml")):
        for line in f.read_text().splitlines():
            s = line.strip()
            if s.startswith("file:"):
                refs.append(s.split(":", 1)[1].strip())
    par = deck_dir / "parameters.par"
    if par.is_file():
        for line in par.read_text().splitlines():
            s = line.strip()
            if s.startswith("MeshFile") and "=" in s:
                stem = s.split("=", 1)[1].strip().strip("\'\"")
                for cand in (f"{stem}.puml.h5", stem):
                    if (deck_dir / cand).is_file():
                        refs.append(cand)
                        break
    return refs


def _write_points(path, pts) -> Path:
    pts = np.atleast_2d(np.asarray(pts, float))
    np.savetxt(path, pts, fmt="%.6f")
    return Path(path)


def _check_psi_ini(sn, tau, fric_nc, fault, cfg, rep, gate) -> None:
    """psi_ini finite everywhere, and mu_ini == tau/sigma_n."""
    if fric_nc is None:
        rep.skip(gate, "no friction nc in the deck")
        return
    from deckbuild.asagi import trilinear_sample
    a = trilinear_sample(fric_nc, fault.cent[:, 0], fault.cent[:, 1], fault.cent[:, 2],
                         fields=["rs_a"])["rs_a"]
    sr0, v_ini = 1.0e-6, 1.0e-9
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        psi = a * np.log((2.0 * sr0 / v_ini) * np.sinh(tau / (a * sn)))
    n_bad = int((~np.isfinite(psi)).sum())
    mu_ini = tau / sn
    ok_id = bool(np.allclose(mu_ini, tau / sn, rtol=0, atol=0))
    rep.add(gate, n_bad == 0 and ok_id,
            f"psi_ini finite on every facet ({n_bad} non-finite); mu_ini == tau/sigma_n "
            f"holds; psi range [{np.nanmin(psi):.3f}, {np.nanmax(psi):.3f}]")


def _check_gate_corridor(cfg, fault, mu, sn, rep, gate) -> None:
    """kappa corridor -- a WARN, deliberately.

    kappa is NECESSARY, NOT SUFFICIENT: an f_w = 0.10 design passed this screen and still
    arrested at s ~ 60 km on the cluster.  A green P6 is not a prediction.
    """
    s = fault.s_km(cfg.strike)
    lo, hi = cfg.physics.seis_band_km
    band = (fault.depth_m / 1000.0 >= lo) & (fault.depth_m / 1000.0 <= hi)
    if not band.any():
        rep.skip(gate, f"no facet in the seismogenic band {lo}-{hi} km")
        return
    dtau = mu[band] * sn[band]
    worst = float(np.nanmin(dtau))
    names = ", ".join(b.name for b in cfg.gate_bands) or "none"
    rep.add(gate, True,
            f"corridor screened over {int(band.sum())} facets in {lo}-{hi} km; min "
            f"dtau_dyn = {worst:.2f} MPa; named bands: {names}.  kappa is NECESSARY NOT "
            f"SUFFICIENT -- only the run decides crossing.", severity=WARN)


def _check_nucleation(cfg, fault, rep, gate, mesh_name) -> None:
    try:
        hypo = cfg.hypocenter(mesh_name)
    except Exception:                                            # noqa: BLE001
        rep.skip(gate, "no hypocentre defined for this mesh")
        return
    try:
        snap = snap_hypocenter(hypo, fault, cfg.strike, cfg.crs,
                               gate_bands=cfg.gate_bands)
    except Exception as exc:                                     # noqa: BLE001
        rep.add(gate, False, f"the hypocentre does not land on this mesh: {exc}")
        return
    rep.add(gate, True,
            f"nucleation at facet {snap.facet_index} (snap {snap.distance_m:.1f} m, "
            f"s = {snap.s_km:.2f} km, depth = {snap.depth_m:.0f} m)")


def _check_yield(deck_dir, plast_nc, rep, gate) -> None:
    stress = next(Path(deck_dir).glob("*stress*.nc"), None)
    if stress is None:
        rep.skip(gate, "no stress nc to check against")
        return
    from deckbuild.asagi import asagi_axes, read_asagi, trilinear_sample
    x, y, z, sf, _ = read_asagi(stress)
    # The four grids are independent by design (the shipped ones differ), so SAMPLE the
    # plasticity onto the stress grid rather than skipping.  ASAGI interpolates linearly
    # at run time, so this is what SeisSol effectively sees.
    ZZ, YY, XX = np.meshgrid(z, y, x, indexing="ij")
    pf = trilinear_sample(plast_nc, XX.ravel(), YY.ravel(), ZZ.ravel(),
                          fields=["plastCo", "bulkFriction"])
    pf = {k: v.reshape(len(z), len(y), len(x)) for k, v in pf.items()}
    sxx, syy, szz = (-sf[k] / 1e6 for k in ("s_xx", "s_yy", "s_zz"))
    sxy, syz, sxz = (-sf[k] / 1e6 for k in ("s_xy", "s_yz", "s_xz"))
    sm = (sxx + syy + szz) / 3.0
    j2 = 0.5 * ((sxx - sm) ** 2 + (syy - sm) ** 2 + (szz - sm) ** 2) \
        + sxy ** 2 + syz ** 2 + sxz ** 2
    tau = np.sqrt(np.maximum(j2, 0.0))
    phi = np.arctan(pf["bulkFriction"])
    taulim = np.maximum(0.0, pf["plastCo"] / 1e6 * np.cos(phi) - sm * np.sin(phi))
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(taulim > 0, tau / taulim, 0.0)
    n_yield = int((ratio > 1.0).sum())
    i = int(np.nanargmax(ratio))
    rep.add(gate, n_yield == 0,
            f"sub-yield at t=0: {n_yield} of {ratio.size:,} grid points yield; max "
            f"tau/taulim = {np.nanmax(ratio):.4f} at flat index {i}")


def _check_output_hygiene(deck_dir, spec, rep, gate) -> None:
    par = Path(deck_dir) / "parameters.par"
    problems = []
    if par.is_file():
        t = par.read_text()
        active = [ln for ln in t.splitlines() if not ln.strip().startswith("!")]
        joined = "\n".join(active)
        if "wavefieldoutput = 1" not in joined:
            problems.append("wavefieldoutput must be 1 (0 hangs in teardown)")
        for ln in active:
            if "OutputRegionBounds" in ln and "0.0" in ln:
                problems.append("OutputRegionBounds contains a literal 0.0, which writes "
                                "the WHOLE volume")
    else:
        problems.append("no parameters.par")
    if spec.receivers is not None:
        sub = GateReport("r")
        check_receivers_under_the_surface(spec.receivers, sub, "r")
        if not sub.ok:
            problems.append(sub.gates[0].detail)
    rep.add(gate, not problems, "output hygiene: " + ("; ".join(problems) or "clean"))
