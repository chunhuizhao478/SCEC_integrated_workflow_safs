"""stress.py -- the initial stress field.

    orientation x Sv(z) x closure k  ->  data{s_xx..s_xz}, compression-NEGATIVE Pa

The closure is Andersonian C1 (legacy `magnitudes_C1` + `build_tensor_andersonian`),
copied expression for expression:

    sig2 = Sv_eff                          (Anderson strike-slip: sigma2 is VERTICAL)
    sig3 = Sv_eff / ((1 - R) k + R)
    sig1 = k sig3
    sigma = sig1 eH eH^T + sig3 eh eh^T + sig2 ev ev^T,  eH = (sin az, cos az, 0)

That vertical-sigma2 assumption is a STRIKE-SLIP assumption.  A normal or thrust fault
needs a different recipe (the plan adds an "analytic" stress kind in Phase 9).

`k` is a COLUMN property, exactly like R and the azimuth, so a per-column `k_map`
broadcasts through `magnitudes_C1` the same way a scalar does -- that is what makes
graded-k a two-line change rather than a new code path.

Two shallow treatments, both bug fixes, both preserved verbatim:
  * DAYLIGHT PATCH -- fill z >= 0 with the shallowest sub-surface slice, so a daylighting
    facet does not get sigma_n = tau_0 = 0 and hence psi_ini = -inf ("Inf/NaN in
    energies" at init).
  * SHALLOW FREEZE -- every node shallower than `freeze_above_depth_m` carries its OWN
    column's recipe evaluated AT that depth.  sigma_n and tau freeze together (the C1
    tensor scales with Sv_eff at fixed k), so mu_app and the lateral structure are
    unchanged and BELOW the depth the field is bit-identical.  Without it the trace band
    sits at sigma_n < 1 MPa and slides for the whole run (446 m of slip in one real case).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from deckbuild.asagi import read_asagi, write_asagi
from deckbuild.config import ConfigError, Project
from deckbuild.contract import Artifact, GateReport, HARD, Stage, WARN
from deckbuild.geometry import UP, build_grid, load_fault, strike_s_km
from deckbuild.orientation import read_orientation

__all__ = ["StressStage", "KDesign", "StressError", "STRESS_FIELDS",
           "magnitudes_C1", "build_tensor_andersonian", "k_profile_1d", "sv_eff_at"]

STRESS_FIELDS = ["s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz"]
COMP_IJ = {"s_xx": (0, 0), "s_yy": (1, 1), "s_zz": (2, 2),
           "s_xy": (0, 1), "s_yz": (1, 2), "s_xz": (0, 2)}


class StressError(ValueError):
    """Raised for an unbuildable stress design."""


# --------------------------------------------------------------------------- physics
def magnitudes_C1(Sv_eff, R, k):
    """Closure C1.  VERBATIM from the legacy `magnitudes_C1`.

    Ordered sig1 >= sig2 >= sig3 > 0 for k >= 1, R in [0, 1].  `k` may be a scalar or a
    per-column array -- the algebra is pure elementwise, so both broadcast identically.
    """
    sig2 = Sv_eff
    sig3 = Sv_eff / ((1.0 - R) * k + R)
    sig1 = k * sig3
    return sig1, sig2, sig3


def build_tensor_andersonian(az_deg, sig1, sig2, sig3):
    """One principal axis exactly vertical.  VERBATIM from the legacy function.

    e_H = (sin az, cos az, 0) carries sig1; e_h = (cos az, -sin az, 0) carries sig3;
    e_v = (0, 0, 1) carries sig2.  az is SHmax, degrees east of north.
    """
    az = np.radians(az_deg)
    z0 = np.zeros_like(az)
    eH = np.stack([np.sin(az), np.cos(az), z0], axis=1)
    eh = np.stack([np.cos(az), -np.sin(az), z0], axis=1)
    ev = np.broadcast_to(UP, eH.shape)
    return (sig1[:, None, None] * np.einsum("ni,nj->nij", eH, eH)
            + sig3[:, None, None] * np.einsum("ni,nj->nij", eh, eh)
            + sig2[:, None, None] * np.einsum("ni,nj->nij", ev, ev))


def sv_eff_at(depth_m, depth_grid, sv_grid):
    """Sv_eff (MPa) at positive-down depths, clamped.  legacy: sv_total_at."""
    return np.interp(np.clip(depth_m, depth_grid[0], depth_grid[-1]), depth_grid, sv_grid)


# --------------------------------------------------------------------------- graded k
@dataclass(frozen=True)
class KDesign:
    """An N-region along-strike closure-ratio design.

    k_values      one k per region, ordered along increasing s.
    boundaries    exactly len(k_values) - 1 smoothstep windows (s_lo, s_hi) in km.
    """

    k_values: tuple[float, ...]
    boundaries_s_km: tuple[tuple[float, float], ...] = ()

    def __post_init__(self) -> None:
        if not self.k_values:
            raise StressError("KDesign needs at least one k value")
        if any(k < 1.0 for k in self.k_values):
            raise StressError(
                f"every k must be >= 1 (the C1 closure assumes sigma1 >= sigma3); "
                f"got {self.k_values}")
        if len(self.boundaries_s_km) != len(self.k_values) - 1:
            raise StressError(
                f"{len(self.k_values)} region(s) need "
                f"{len(self.k_values) - 1} transition window(s), got "
                f"{len(self.boundaries_s_km)}")
        prev_hi = -np.inf
        for i, (lo, hi) in enumerate(self.boundaries_s_km):
            if hi <= lo:
                raise StressError(
                    f"transition window {i} has s_hi {hi} <= s_lo {lo}")
            if lo < prev_hi:
                raise StressError(
                    f"transition window {i} (s {lo}-{hi}) overlaps the previous one; "
                    f"windows must be ordered and disjoint")
            prev_hi = hi

    @property
    def n_regions(self) -> int:
        return len(self.k_values)

    @property
    def is_uniform(self) -> bool:
        return len(self.k_values) == 1


def k_profile_1d(s_km, design: KDesign) -> np.ndarray:
    """k(s) as region plateaus joined by smoothstep windows.

    smoothstep t*t*(3 - 2t) is the same shape the shipped rs_muw Lua uses, so a k design
    and an f_w design grade identically.
    """
    s = np.asarray(s_km, float)
    out = np.full(s.shape, float(design.k_values[0]))
    for i, (lo, hi) in enumerate(design.boundaries_s_km):
        t = np.clip((s - lo) / (hi - lo), 0.0, 1.0)
        out = out + (design.k_values[i + 1] - design.k_values[i]) * t * t * (3.0 - 2.0 * t)
    return out


# --------------------------------------------------------------------------- the stage
class StressStage(Stage):
    name = "stress"

    def build(self, cfg: Project, out_dir, *, sv_profile, k=None, design=None,
              freeze_above_depth_m: float = 500.0, daylight_patch=None,
              daylight_min_depth_m=None, mesh: str | None = None,
              out_name: str | None = None, dtype=np.float32) -> Artifact:
        """Build the stress nc.

        sv_profile -- REQUIRED.  The Artifact (or (depth_m, sv_eff_mpa) tuple) produced by
                      MaterialStage.  Making it an argument rather than a default path
                      means a rebuilt velocity model cannot silently fail to reach the
                      stress field, and the artifact records which material it came from.
        """
        if (k is None) == (design is None):
            raise StressError("give exactly one of k= (scalar) or design= (KDesign)")
        if k is not None:
            design = KDesign(k_values=(float(k),))
        if design.is_uniform and design.k_values[0] < 1.0:
            raise StressError(f"k must be >= 1, got {design.k_values[0]}")

        mesh_name = mesh or cfg.default_mesh
        mspec = cfg.mesh(mesh_name)
        if daylight_patch is None:
            daylight_patch = mspec.daylights
        if daylight_min_depth_m is None:
            daylight_min_depth_m = mspec.daylight_min_depth_m
        if daylight_patch and cfg.stress_box.zmax > 0.0 and daylight_min_depth_m <= 0.0:
            # The legacy patch reads z_topo(x, y) off the mesh's free surface and fills
            # each column with ITS OWN overburden.  This implementation writes one
            # constant slab, which is correct only for a flat-top domain.  Refuse rather
            # than silently write the approximation on a topographic one.
            raise StressError(
                "daylight_patch on a domain with topography (stress_box.zmax = "
                f"{cfg.stress_box.zmax:g} > 0) needs a per-column topographic overburden, "
                "which this build does not yet compute.  Set "
                "meshes.<name>.daylight_min_depth_m > 0 to use the constant-slab "
                "approximation deliberately, or keep zmax <= 0 for a flat top.")

        depth_grid, sv_grid, sv_sha = _sv_arrays(sv_profile)
        gx, gy, gz = build_grid(cfg.stress_box)

        # Per-column orientation and k, on the (x-major) column list -- the same order
        # the legacy step2 used: meshgrid(gx, gy, indexing='ij').ravel().
        XX, YY = np.meshgrid(gx, gy, indexing="ij")
        colx, coly = XX.ravel(), YY.ravel()
        ofield = read_orientation(cfg.raw.orientation, cfg.crs, cfg.data_dir)
        az_col, R_col, n_fb = ofield.at(colx, coly)
        s_col = strike_s_km(colx, coly, cfg.strike)
        k_col = k_profile_1d(s_col, design)

        comps = _assemble(gz, depth_grid, sv_grid, az_col, R_col, k_col,
                          len(gx), len(gy), freeze_above_depth_m)

        if daylight_patch:
            _apply_daylight_patch(comps, gz, depth_grid, sv_grid, az_col, R_col, k_col,
                                  len(gx), len(gy), daylight_min_depth_m)

        out_dir = Path(out_dir)
        name = out_name or _default_name(design, freeze_above_depth_m)
        path = out_dir / name
        attrs = {
            "title": f"{cfg.name} initial stress, Andersonian C1",
            "convention": "compression-NEGATIVE Pa, z = elevation",
            "closure": "C1: sig2 = Sv_eff (vertical); sig3 = Sv_eff/((1-R)k+R); sig1 = k*sig3",
            "k_values": ", ".join(f"{v:g}" for v in design.k_values),
            "boundaries_s_km": str([list(b) for b in design.boundaries_s_km]),
            "freeze_above_depth_m": float(freeze_above_depth_m),
            "daylight_patch": int(bool(daylight_patch)),
            "daylight_min_depth_m": float(daylight_min_depth_m),
            "orientation_kind": ofield.kind,
            "orientation_nearest_fallback_columns": int(n_fb),
            "sv_profile_sha256": sv_sha,
            "project": cfg.name,
            "crs": cfg.crs.epsg,
        }
        write_asagi(path, gx, gy, gz, comps, attrs=attrs, dtype=dtype)
        return Artifact.of(path, kind="stress",
                           params={"k_values": list(design.k_values),
                                   "boundaries_s_km": [list(b) for b in
                                                       design.boundaries_s_km],
                                   "freeze_above_depth_m": float(freeze_above_depth_m),
                                   "daylight_patch": bool(daylight_patch),
                                   "daylight_min_depth_m": float(daylight_min_depth_m),
                                   "mesh": mesh_name},
                           provenance={"sv_profile_sha256": sv_sha,
                                       "orientation_kind": ofield.kind,
                                       "descriptor_sha256": cfg.sha256()})

    # ---------------------------------------------------------------- verify
    def verify(self, cfg: Project, artifact: Artifact, *, design=None,
               freeze_above_depth_m=None, **_) -> GateReport:
        """V1-V4, names preserved from the legacy graded_k.verify.

        V1 (regional bit-identity against a per-region constant-k build) needs the whole
        build pipeline and is exercised in the stage's own tests, where the reference
        builds are scratch; here it is reported as skipped unless a design is supplied.
        """
        p = Path(artifact.path)
        x, y, z, flds, attrs = read_asagi(p)
        rep = GateReport(f"stress V1-V4 ({p.name})")

        design = design or _design_from_attrs(attrs)
        if freeze_above_depth_m is None:
            freeze_above_depth_m = float(attrs.get("freeze_above_depth_m", 0.0))

        rep.skip("V1", "not run HERE: regional bit-identity needs the Sv profile and a "
                       "scratch dir.  Call StressStage.verify_regions(...) and merge its "
                       "report with rep.extend() for the complete V1-V4 battery.")

        # V2 -- the WRITTEN tensor's horizontal eigenvalue ratio must equal k(s).
        sig = _tensor_stack(flds)                       # (nz, ny, nx, 3, 3), comp-POSITIVE
        XX, YY = np.meshgrid(x, y, indexing="ij")
        s_col = strike_s_km(XX.ravel(), YY.ravel(), cfg.strike).reshape(len(x), len(y))
        k_want = k_profile_1d(s_col, design)            # (nx, ny)

        # Use a mid-depth slice: shallow levels may be frozen or zero.
        kz = _probe_level(z, sig, freeze_above_depth_m)
        if kz is None:
            rep.add("V2", False, "no level with a non-zero tensor to probe")
        else:
            h = sig[kz][:, :, :2, :2]                   # (ny, nx, 2, 2)
            w = np.linalg.eigvalsh(h)                   # ascending
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = w[:, :, 1] / w[:, :, 0]
            got = ratio.T                               # -> (nx, ny)
            good = np.isfinite(got)
            err = np.abs(got[good] - k_want[good])
            worst = float(err.max()) if err.size else np.inf
            rep.add("V2", worst < 2e-3,
                    f"eig sigma1/sigma3 vs k(s) at z = {z[kz]:.0f} m: max |d| = "
                    f"{worst:.2e} over {int(good.sum())} columns "
                    f"(k range {k_want.min():.3f}-{k_want.max():.3f})")

        # V3 -- the freeze.
        if freeze_above_depth_m > 0:
            zf = -abs(freeze_above_depth_m)
            above = np.flatnonzero(z > zf)
            if above.size == 0:
                rep.skip("V3", f"no grid level above the freeze depth {zf:.0f} m")
            else:
                kref = int(np.argmin(np.abs(z - zf)))
                same = all(np.array_equal(flds[f][i], flds[f][kref])
                           for f in STRESS_FIELDS for i in above)
                rep.add("V3", same,
                        f"freeze: {above.size} level(s) above {zf:.0f} m "
                        f"{'all equal' if same else 'DIFFER from'} the z = {z[kref]:.0f} m "
                        f"slice")
        else:
            rep.skip("V3", "freeze_above_depth_m = 0 (the shallow freeze is off)")

        # V4 -- the on-fault screen.  Requires the mesh; skipped when it is absent.
        mesh_spec = cfg.mesh(str(attrs.get("mesh") or cfg.default_mesh)
                             if attrs.get("mesh") else None)
        mesh_path = cfg.resolve_path(mesh_spec.path)
        if not Path(mesh_path).is_file():
            # ABSENT is a legitimate skip; PRESENT-BUT-UNUSABLE is a failure.  Collapsing
            # the two would let a wrong fault_bc or a corrupt file pass as "ok".
            rep.skip("V4", f"no mesh at {mesh_path}; the on-fault screen did not run")
            return rep
        try:
            fault = load_fault(mesh_spec, cfg.data_dir, strike=cfg.strike)
        except Exception as exc:                        # noqa: BLE001
            rep.add("V4", False,
                    f"mesh {mesh_path} exists but could not be read: "
                    f"{type(exc).__name__}: {exc}")
            return rep

        from deckbuild.asagi import hull_containment
        hull_containment(p, fault.cent, label="fault facets", gate="V4hull",
                         severity=HARD, report=rep)
        sn, tau, mu = _project_onto_fault(p, fault)
        finite = np.isfinite(mu)
        rep.add("V4a", bool(np.all(sn[finite] > 0)),
                f"sigma_n > 0 on every facet: min = {np.nanmin(sn):.3f} MPa "
                f"(tau_0 min = {np.nanmin(tau):.4f} MPa)")
        rep.add("V4b", bool(finite.any()),
                f"mu_app finite on {int(finite.sum())}/{len(mu)} facets; "
                f"max mu_app = {np.nanmax(mu):.4f}", severity=WARN)
        return rep

    def verify_regions(self, cfg: Project, artifact: Artifact, *, sv_profile,
                       design: KDesign, tmp_dir, **kw) -> GateReport:
        """V1: each region must be bit-identical to a constant-k build of that region.

        The reference builds are SCRATCH -- they are written into `tmp_dir` and are not
        deliverables.  This is the check that proves grading did not perturb a region.
        """
        p = Path(artifact.path)
        _, _, _, graded, _ = read_asagi(p)
        x, y, z, _, attrs = read_asagi(p)
        XX, YY = np.meshgrid(x, y, indexing="ij")
        s_col = strike_s_km(XX.ravel(), YY.ravel(), cfg.strike).reshape(len(x), len(y))

        rep = GateReport(f"stress V1 ({p.name})")
        edges = _region_edges(design)
        for i, kv in enumerate(design.k_values):
            lo, hi = edges[i]
            plateau = (s_col >= lo) & (s_col <= hi)      # (nx, ny)
            if not plateau.any():
                rep.add(f"V1r{i + 1}", True,
                        f"region {i + 1}/{design.n_regions} (k={kv:g}) has no column at "
                        f"its plateau -- a window this narrow is a taper, not a region",
                        severity=WARN)
                continue
            ref = self.build(cfg, tmp_dir, sv_profile=sv_profile, k=kv,
                             freeze_above_depth_m=float(
                                 attrs.get("freeze_above_depth_m", 0.0)),
                             daylight_patch=bool(attrs.get("daylight_patch", 0)),
                             daylight_min_depth_m=float(
                                 attrs.get("daylight_min_depth_m", 0.0)),
                             out_name=f"_v1_ref_k{kv:g}.nc", **kw)
            _, _, _, cst, _ = read_asagi(ref.path)
            mask = plateau.T                              # -> (ny, nx)
            same = all(np.array_equal(graded[f][:, mask], cst[f][:, mask])
                       for f in STRESS_FIELDS)
            rep.add(f"V1r{i + 1}", same,
                    f"region {i + 1}/{design.n_regions} (k={kv:g}, s {lo:.1f}-{hi:.1f} km, "
                    f"{int(plateau.sum())} columns) "
                    f"{'bit-identical to' if same else 'DIFFERS from'} a constant-k build")
        return rep


# --------------------------------------------------------------------------- helpers
def _sv_arrays(sv_profile):
    """Accept an Artifact (npz), a path, or a (depth_m, sv_eff_mpa) tuple."""
    if isinstance(sv_profile, tuple) and len(sv_profile) == 2:
        d, s = (np.asarray(v, float) for v in sv_profile)
        return d, s, ""
    path = Path(getattr(sv_profile, "path", sv_profile))
    if not path.is_file():
        raise StressError(f"Sv profile not found: {path}")
    z = np.load(path)
    for dk, sk in (("depth_m", "sv_eff_mpa"), ("depth_grid", "sv_total"),
                   ("depth", "sv")):
        if dk in z and sk in z:
            sha = getattr(sv_profile, "sha256", "")
            return np.asarray(z[dk], float), np.asarray(z[sk], float), sha
    raise StressError(
        f"{path}: expected arrays (depth_m, sv_eff_mpa); have {sorted(z.files)}")


def _assemble(gz, depth_grid, sv_grid, az_col, R_col, k_col, nx, ny,
              freeze_above_depth_m):
    """Build every z level.  Returns {field: (nz, ny, nx)} compression-NEGATIVE Pa."""
    nz = len(gz)
    comps = {f: np.zeros((nz, ny, nx)) for f in STRESS_FIELDS}
    zf = -abs(freeze_above_depth_m) if freeze_above_depth_m else None
    for iz, zz in enumerate(gz):
        # The shallow freeze: a node above the freeze depth uses its OWN column's recipe
        # evaluated AT the freeze depth.  sigma_n and tau scale together, so mu_app and
        # the lateral structure are untouched.
        z_eval = zz if (zf is None or zz <= zf) else zf
        depth = -z_eval
        if depth <= 0:
            continue                       # above sea level with no freeze -> zero
        sv = float(sv_eff_at(np.array([depth]), depth_grid, sv_grid)[0])
        if sv <= 0:
            continue
        s1, s2, s3 = magnitudes_C1(np.full(nx * ny, sv), R_col, k_col)
        sig = build_tensor_andersonian(az_col, s1, s2, s3)     # comp-POSITIVE, MPa
        for f, (a, b) in COMP_IJ.items():
            comps[f][iz] = (-sig[:, a, b] * 1.0e6).reshape(nx, ny).T
    return comps


def _apply_daylight_patch(comps, gz, depth_grid, sv_grid, az_col, R_col, k_col,
                          nx, ny, min_depth_m):
    """Fill z >= 0 with the shallowest sub-surface slice.

    Without this a daylighting facet gets sigma_n = tau_0 = 0, and the rate-and-state
    state init psi = a*ln[(2 sr0/V_ini) sinh(tau_0/(a sigma_n))] is -inf -- SeisSol aborts
    with "Inf/NaN in energies" before the first step.
    """
    above = np.flatnonzero(gz >= 0.0)
    if above.size == 0:
        return
    depth = max(float(min_depth_m), float(-gz[gz < 0].max()) if (gz < 0).any() else 0.0)
    if depth <= 0:
        return
    sv = float(sv_eff_at(np.array([depth]), depth_grid, sv_grid)[0])
    s1, s2, s3 = magnitudes_C1(np.full(nx * ny, sv), R_col, k_col)
    sig = build_tensor_andersonian(az_col, s1, s2, s3)
    for f, (a, b) in COMP_IJ.items():
        slab = (-sig[:, a, b] * 1.0e6).reshape(nx, ny).T
        for iz in above:
            comps[f][iz] = slab


def _default_name(design: KDesign, freeze_m: float) -> str:
    if design.is_uniform:
        core = f"k{design.k_values[0]:g}"
    else:
        core = "gradedk_" + "-".join(f"{v:g}" for v in design.k_values)
    tag = f"_freeze{freeze_m:g}m" if freeze_m else ""
    return f"stress_andersonian_{core}{tag}.nc"


def _design_from_attrs(attrs) -> KDesign:
    ks = tuple(float(v) for v in str(attrs.get("k_values", "")).split(",") if v.strip())
    if not ks:
        raise StressError("cannot recover the k design from the nc attributes")
    import ast
    b = attrs.get("boundaries_s_km", "[]")
    bounds = tuple(tuple(float(v) for v in w) for w in ast.literal_eval(str(b)))
    return KDesign(k_values=ks, boundaries_s_km=bounds)


def _region_edges(design: KDesign):
    """(plateau_lo, plateau_hi) in s for each region -- the flat part, excluding tapers."""
    out = []
    for i in range(design.n_regions):
        lo = -np.inf if i == 0 else design.boundaries_s_km[i - 1][1]
        hi = np.inf if i == design.n_regions - 1 else design.boundaries_s_km[i][0]
        out.append((lo, hi))
    return out


def _tensor_stack(flds):
    """{field: (nz,ny,nx)} comp-NEGATIVE Pa -> (nz,ny,nx,3,3) comp-POSITIVE MPa."""
    shape = next(iter(flds.values())).shape
    sig = np.zeros(shape + (3, 3))
    for f, (a, b) in COMP_IJ.items():
        sig[..., a, b] = -flds[f] / 1.0e6
        sig[..., b, a] = sig[..., a, b]
    return sig


def _probe_level(z, sig, freeze_m):
    """A z index safely below the freeze depth with a non-zero tensor."""
    zf = -abs(freeze_m) if freeze_m else 0.0
    for kz in range(len(z)):
        if z[kz] <= zf - 1.0 and np.abs(sig[kz]).max() > 0:
            return kz
    for kz in range(len(z)):
        if np.abs(sig[kz]).max() > 0:
            return kz
    return None


def _project_onto_fault(nc_path, fault):
    """Sample the stress nc at facet centroids -> (sigma_n, tau, mu), comp-POSITIVE MPa."""
    from deckbuild.asagi import asagi_axes, trilinear_sample
    from deckbuild.geometry import resolve_tractions
    _, _, zg = asagi_axes(nc_path)
    smp = trilinear_sample(nc_path, fault.cent[:, 0], fault.cent[:, 1],
                           np.clip(fault.cent[:, 2], zg.min(), zg.max()))
    sig = np.zeros((len(fault), 3, 3))
    for f, (a, b) in COMP_IJ.items():
        sig[:, a, b] = -smp[f] / 1.0e6
        sig[:, b, a] = sig[:, a, b]
    return resolve_tractions(sig, fault.strikes, fault.dips, fault.normals)
