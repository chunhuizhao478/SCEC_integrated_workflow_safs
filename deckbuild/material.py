"""material.py -- the rock-property model and everything derived from it.

    raw velocity slices (or a 1-D model)  ->  data{rho, mu, lambda}
                                          ->  plasticity data{plastCo, bulkFriction}
                                          ->  the Sv(z) profile the stress closure needs
    raw temperature slices                ->  data{T}

Two ordering rules are load-bearing and look like inefficiencies, so they are called out
where they happen:

  1. Moduli are formed AT THE SOURCE NODES (mu = rho Vs^2, lambda = rho(Vp^2 - 2Vs^2)) and
     only THEN resampled onto the uniform z axis.  Converting after resampling would store
     something that is not the piecewise-linear interpolant of the node moduli.
  2. `plastCo = cohesion_factor * mu` is LINEAR in mu, so ASAGI's linear interpolation of
     plastCo equals the factor times the interpolated mu SeisSol uses for the elastic
     material -- the cohesion an element sees is exactly consistent with its shear modulus.
     `bulkFriction` is the one thresholded field.

SeisSol wants the friction COEFFICIENT tan(phi), not degrees: Model/Plasticity.h computes
atan(bulkFriction).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from deckbuild.asagi import read_asagi, write_asagi
from deckbuild.config import Project
from deckbuild.contract import Artifact, GateReport, HARD, Stage, WARN
from deckbuild.geometry import load_fault

__all__ = ["MaterialStage", "MaterialArtifacts", "PlasticitySpec", "AttenuationSpec",
           "MaterialError", "derive_plasticity", "sv_profile_from_material"]

G_ACCEL = 9.81
RHO_W = 1000.0


class MaterialError(ValueError):
    """Raised for an unbuildable material model."""


@dataclass(frozen=True)
class PlasticitySpec:
    """Roten et al. (2014) Drucker-Prager.

    The PAPER's angles are 35/45.  The SAFS production decks use 30/40.  Neither may be a
    silent default, so the paper values are the default and a deck that wants otherwise
    must say so.
    """

    phi_soft_deg: float = 35.0
    phi_hard_deg: float = 45.0
    vs_threshold: float = 2500.0
    cohesion_factor: float = 1.0e-4


@dataclass(frozen=True)
class AttenuationSpec:
    """Q from the CVM's own Vs.  This stage only records the numbers and emits the YAML
    block; SeisSol computes the Q fields from mu/rho at run time."""

    qs_over_vs: float = 0.05
    qp_over_qs: float = 2.0
    freq_central: float = 0.5
    freq_ratio: float = 100.0


@dataclass
class MaterialArtifacts:
    material: Artifact | None = None
    plasticity: Artifact | None = None
    thermal: Artifact | None = None
    sv_profile: Artifact | None = None
    attenuation: dict | None = None


# --------------------------------------------------------------------------- derivations
def derive_plasticity(rho, mu, spec: PlasticitySpec):
    """(plastCo, bulkFriction) pointwise on the material grid."""
    rho = np.asarray(rho, float)
    mu = np.asarray(mu, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        vs = np.sqrt(np.where(rho > 0, mu / rho, np.nan))
    tan_soft = float(np.tan(np.radians(spec.phi_soft_deg)))
    tan_hard = float(np.tan(np.radians(spec.phi_hard_deg)))
    bulk = np.where(vs < spec.vs_threshold, tan_soft, tan_hard)
    return spec.cohesion_factor * mu, bulk, vs


def sv_profile_from_material(material_nc, g: float = G_ACCEL, rho_w: float = RHO_W,
                             pore: str = "hydrostatic"):
    """Lateral-mean Sv_eff(depth) from a material nc.

    Sv_total is the trapezoidal integral of the LATERAL-MEAN rho(z) downward from the top
    of the grid; Sv_eff subtracts the pore pressure.  This is a property of the material,
    so it is derived here and handed to the stress stage rather than shipped as a file.
    """
    if pore != "hydrostatic":
        raise MaterialError(f"unknown pore-pressure model {pore!r} (only 'hydrostatic')")
    _, _, z, flds, _ = read_asagi(material_nc, fields=["rho"])
    rho_mean = flds["rho"].reshape(len(z), -1).mean(axis=1)     # per z level
    order = np.argsort(-np.asarray(z))                          # top-down
    ztop = np.asarray(z)[order]
    rho_td = rho_mean[order]

    # SEA-LEVEL referenced: depth = -z, so depth 0 means z = 0.  That is exactly what the
    # stress closure asks for (it evaluates Sv at depth = -z_level).  Referencing to the
    # GRID TOP instead would offset the entire stress field by z_max -- +3250 m for the
    # shipped SAFS CVM, since extend_z_top exists precisely to cover topography.  The
    # eigenvalue gate V2 cannot catch that: it checks the RATIO k, which is invariant
    # under an overall scale of Sv.
    depth = -ztop                                               # 0 at SEA LEVEL
    # Integrate the overburden downward from the actual grid top, then re-reference.
    span = ztop[0] - ztop                                       # 0 at the grid top
    col = np.concatenate([[0.0], np.cumsum(
        0.5 * (rho_td[1:] + rho_td[:-1]) * g * np.diff(span))]) / 1.0e6
    pp = rho_w * g * np.maximum(depth, 0.0) / 1.0e6
    return depth, col - pp


# --------------------------------------------------------------------------- readers
def _layered_1d(cfg: Project, p):
    """A laterally uniform 1-D model.  No community model needed."""
    d = np.asarray(p["depth_m"], float)
    if not np.all(np.diff(d) > 0):
        raise MaterialError(f"layered_1d depth_m must be strictly increasing, got {d}")
    box = cfg.stress_box
    dx = float(p.get("grid_dx", box.dx))
    dz = float(p.get("dz", box.dz))
    z_min = float(p.get("z_min", box.zmin))
    z_max = float(p.get("z_max", box.zmax))
    gx = np.arange(box.xmin, box.xmax + 0.5 * dx, dx)
    gy = np.arange(box.ymin, box.ymax + 0.5 * dx, dx)
    gz = np.arange(z_min, z_max + 0.5 * dz, dz)
    depth = np.clip(-gz, d[0], d[-1])
    vp = np.interp(depth, d, np.asarray(p["vp"], float))
    vs = np.interp(depth, d, np.asarray(p["vs"], float))
    rho = np.interp(depth, d, np.asarray(p["rho"], float))
    # Rule 1: moduli at the (interpolated) nodes, from vp/vs/rho -- never from resampled
    # moduli.  For a 1-D model the two coincide, but the ORDER is kept the same as the
    # slice reader so the two paths cannot drift.
    mu = rho * vs ** 2
    lam = rho * (vp ** 2 - 2.0 * vs ** 2)
    if not np.all(lam > 0):
        bad = int((lam <= 0).sum())
        raise MaterialError(
            f"M1: lambda = rho(Vp^2 - 2Vs^2) <= 0 at {bad} depth(s); this is a bad "
            f"velocity model, not a code bug (Vp must exceed sqrt(2)*Vs)")
    shape = (len(gz), len(gy), len(gx))

    def bc(a):
        return np.broadcast_to(a[:, None, None], shape).copy()
    return gx, gy, gz, {"rho": bc(rho), "mu": bc(mu), "lambda": bc(lam)}


def _cvm_slices(cfg: Project, p, path_glob=None):
    """Raw CVM horizontal slices -> {rho, mu, lambda} on a uniform grid.

    RULE 1 in action: the moduli are formed AT THE SOURCE NODES and only then
    z-resampled.  interp_slices is therefore asked for vp/vs/rho, and the conversion
    happens on the SLICE stack before the vertical resample -- not after.
    """
    from deckbuild.rawslices import interp_slices
    paths = sorted(Path(cfg.data_dir).glob(path_glob)) if path_glob else []
    if not paths:
        raise MaterialError(
            f"no CVM slices matched {path_glob!r} under {cfg.data_dir}")
    cols = list(p.get("columns", ["lon", "lat", "vp", "vs", "density"]))[2:]
    gx, gy, gz, f3, info = interp_slices(
        paths, cfg.crs, cols,
        grid_dx=float(p.get("grid_dx", 1500.0)),
        z_min=float(p.get("z_min", -45000.0)), z_max=float(p.get("z_max", 100.0)),
        dz=float(p.get("dz", 250.0)),
        extend_z_top=float(p.get("extend_z_top", 100.0)),
        expect_spacing_m=p.get("expect_spacing_m"))
    vp, vs, rho = (f3[c] for c in cols)
    mu = rho * vs ** 2
    lam = rho * (vp ** 2 - 2.0 * vs ** 2)
    if not np.all(lam > 0):
        bad = int((lam <= 0).sum())
        i = int(np.argmin(lam))
        raise MaterialError(
            f"M1: lambda = rho(Vp^2 - 2Vs^2) <= 0 at {bad} node(s); worst at flat index "
            f"{i} (Vp={vp.ravel()[i]:.1f}, Vs={vs.ravel()[i]:.1f}, "
            f"rho={rho.ravel()[i]:.1f}).  This is a bad velocity model, not a code bug.")
    return gx, gy, gz, {"rho": rho, "mu": mu, "lambda": lam}


def _ctm_slices(cfg: Project, p, path_glob=None):
    """Raw CTM temperature slices -> {T} on a uniform grid."""
    from deckbuild.rawslices import interp_slices
    paths = sorted(Path(cfg.data_dir).glob(path_glob)) if path_glob else []
    if not paths:
        raise MaterialError(f"no CTM slices matched {path_glob!r} under {cfg.data_dir}")
    cols = list(p.get("columns", ["Lon", "Lat", "Temperature"]))[2:]
    gx, gy, gz, f1, info = interp_slices(
        paths, cfg.crs, cols,
        grid_dx=float(p.get("grid_dx", 1500.0)),
        z_min=float(p.get("z_min", -21000.0)), z_max=float(p.get("z_max", 0.0)),
        dz=float(p.get("dz", 200.0)),
        extend_z_top=float(p.get("extend_z_top", 200.0)),
        expect_spacing_m=p.get("expect_spacing_m", 200.0))
    return gx, gy, gz, {"T": f1[cols[0]]}


VELOCITY_READERS = {"layered_1d": _layered_1d, "cvm_slices": _cvm_slices}
THERMAL_READERS = {"ctm_slices": _ctm_slices}


# --------------------------------------------------------------------------- the stage
class MaterialStage(Stage):
    name = "material"

    def build(self, cfg: Project, out_dir, *, plasticity: PlasticitySpec | None = None,
              attenuation: AttenuationSpec | None = None, thermal: bool = True,
              dtype=np.float32, prefix: str = "") -> MaterialArtifacts:
        out_dir = Path(out_dir)
        spec = cfg.raw.velocity
        if spec is None:
            raise MaterialError("this project has no raw.velocity source")
        if spec.kind not in VELOCITY_READERS:
            raise MaterialError(
                f"velocity kind {spec.kind!r} is not implemented yet; have "
                f"{sorted(VELOCITY_READERS)}")
        rd = VELOCITY_READERS[spec.kind]
        gx, gy, gz, flds = (rd(cfg, dict(spec.params), spec.path)
                            if spec.kind == "cvm_slices" else rd(cfg, dict(spec.params)))

        out = MaterialArtifacts()
        mpath = out_dir / f"{prefix}material.nc"
        write_asagi(mpath, gx, gy, gz, flds, dtype=dtype, attrs={
            "title": f"{cfg.name} elastic material", "source": spec.kind,
            "units": "rho kg/m^3; mu, lambda Pa", "project": cfg.name,
            "crs": cfg.crs.epsg})
        out.material = Artifact.of(mpath, kind="material",
                                   params={"source": spec.kind},
                                   provenance={"descriptor_sha256": cfg.sha256()})

        if plasticity is not None:
            pc, bf, _vs = derive_plasticity(flds["rho"], flds["mu"], plasticity)
            ppath = out_dir / (f"{prefix}plasticity_phi"
                               f"{plasticity.phi_soft_deg:g}_{plasticity.phi_hard_deg:g}.nc")
            write_asagi(ppath, gx, gy, gz, {"plastCo": pc, "bulkFriction": bf},
                        dtype=dtype, attrs={
                            "title": "Roten et al. (2014) Drucker-Prager",
                            "phi_soft_deg": plasticity.phi_soft_deg,
                            "phi_hard_deg": plasticity.phi_hard_deg,
                            "vs_threshold_m_s": plasticity.vs_threshold,
                            "cohesion_factor": plasticity.cohesion_factor,
                            "note": "bulkFriction is tan(phi), the COEFFICIENT -- "
                                    "SeisSol computes atan() of it"})
            out.plasticity = Artifact.of(ppath, kind="plasticity",
                                         params={"phi_soft_deg": plasticity.phi_soft_deg,
                                                 "phi_hard_deg": plasticity.phi_hard_deg,
                                                 "vs_threshold": plasticity.vs_threshold,
                                                 "cohesion_factor":
                                                     plasticity.cohesion_factor})

        # The Sv profile: cached on the material's own hash, so a rebuilt velocity model
        # can never be paired with a stale profile.
        depth, sv_eff = sv_profile_from_material(mpath)
        spath = out_dir / f"{prefix}sv_profile_{out.material.sha256[:12]}.npz"
        np.savez(spath, depth_m=depth, sv_eff_mpa=sv_eff, reference="sea_level")
        out.sv_profile = Artifact.of(spath, kind="sv_profile",
                                     params={"pore": "hydrostatic", "g": G_ACCEL},
                                     provenance={"material_sha256": out.material.sha256})

        if thermal and cfg.raw.thermal is not None and cfg.raw.thermal.kind == "ctm_slices":
            ts = cfg.raw.thermal
            tx, ty, tz, tf = _ctm_slices(cfg, dict(ts.params), ts.path)
            tpath = out_dir / f"{prefix}thermal_T.nc"
            write_asagi(tpath, tx, ty, tz, tf, dtype=dtype, attrs={
                "title": f"{cfg.name} temperature", "source": "ctm_slices",
                "units": "degC", "project": cfg.name, "crs": cfg.crs.epsg})
            out.thermal = Artifact.of(tpath, kind="thermal",
                                      params={"source": "ctm_slices"},
                                      provenance={"descriptor_sha256": cfg.sha256()})

        if attenuation is not None:
            out.attenuation = {"qs_over_vs": attenuation.qs_over_vs,
                               "qp_over_qs": attenuation.qp_over_qs,
                               "freq_central": attenuation.freq_central,
                               "freq_ratio": attenuation.freq_ratio}
        return out

    def verify(self, cfg: Project, artifact, *, mesh: str | None = None,
               plasticity_artifact=None, **_) -> GateReport:
        """M1-M5."""
        art = artifact.material if isinstance(artifact, MaterialArtifacts) else artifact
        p = Path(art.path)
        x, y, z, flds, _ = read_asagi(p)
        rho, mu, lam = flds["rho"], flds["mu"], flds["lambda"]
        rep = GateReport(f"material M1-M5 ({p.name})")

        rep.add("M1", bool(np.all(lam > 0)),
                f"lambda > 0 everywhere: min = {lam.min():.6g} Pa")
        finite = np.isfinite(rho) & np.isfinite(mu) & np.isfinite(lam)
        rep.add("M2", bool(np.all(finite)),
                f"no NaN/Inf after interpolation: {int((~finite).sum())} bad node(s)")

        mspec = cfg.mesh(mesh or cfg.default_mesh)
        mpath = cfg.resolve_path(mspec.path)
        if not Path(mpath).is_file():
            rep.skip("M3", f"no mesh at {mpath}; hull containment not checked")
        else:
            from deckbuild.asagi import hull_containment
            fault = load_fault(mspec, cfg.data_dir, strike=cfg.strike)
            # HARD over the fault (where the physics is), WARN in the far field (ASAGI
            # edge-clamps it, and the shipped PREFERRED CVM legitimately falls short there).
            hull_containment(p, fault.cent, label="fault bbox", gate="M3",
                             severity=HARD, report=rep)

        spec = cfg.raw.velocity
        if spec is not None and spec.kind == "layered_1d":
            from deckbuild.asagi import roundtrip_selfcheck
            pp = dict(spec.params)
            dk = np.asarray(pp["depth_m"], float)
            rho_k = np.asarray(pp["rho"], float)
            vs_k = np.asarray(pp["vs"], float)

            def ev(pts):                      # mu = rho Vs^2 from the 1-D model
                dep = np.clip(-pts[:, 2], dk[0], dk[-1])
                return np.interp(dep, dk, rho_k) * np.interp(dep, dk, vs_k) ** 2
            # The 1-D model is piecewise linear in DEPTH but mu is a PRODUCT of two such
            # interpolants, so it is piecewise QUADRATIC -- trilinear sampling of it is
            # not exact.  Scale the tolerance to the field (mu ~ 1e10 Pa).
            roundtrip_selfcheck(p, ev, field="mu", n=400,
                                median_tol=2e8, max_tol=2e9, gate="M4", report=rep,
                                kink_straddling=lambda q: np.ones(len(q), bool))
        else:
            rep.skip("M4", f"no cheap evaluator for reader {getattr(spec, 'kind', None)!r}; "
                           f"run asagi.roundtrip_selfcheck against the raw stack")

        with np.errstate(divide="ignore", invalid="ignore"):
            vs = np.sqrt(np.where(rho > 0, mu / rho, np.nan))
        lo, hi = float(np.nanmin(vs)), float(np.nanmax(vs))
        rep.add("M5", 100.0 < lo and hi < 8000.0,
                f"Vs in a plausible range: [{lo:.1f}, {hi:.1f}] m/s",
                severity=WARN)

        if plasticity_artifact is not None:
            _, _, _, pf, pattrs = read_asagi(plasticity_artifact.path)
            bf = pf["bulkFriction"]
            want = {float(np.tan(np.radians(pattrs["phi_soft_deg"]))),
                    float(np.tan(np.radians(pattrs["phi_hard_deg"])))}
            got = set(np.round(np.unique(bf), 6))
            ok = got <= {round(w, 6) for w in want}
            rep.add("M6", ok,
                    f"bulkFriction takes only tan(phi_soft), tan(phi_hard) = "
                    f"{sorted(round(w, 6) for w in want)}; found {sorted(got)}")
            rep.add("M7", bool(np.all(pf["plastCo"] >= 0)),
                    f"plastCo >= 0: min = {pf['plastCo'].min():.6g} Pa")
        return rep
