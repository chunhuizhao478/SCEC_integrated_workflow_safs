"""geometry.py -- the fault frame, shared by every stage.

Ported from the legacy `generate_stress_nc_from_raw.py` / `project_csm_stress_to_vtu.py`.
The expressions are copied character-for-character; only the signatures changed, so that
the strike frame, the BC map and the harmonisation hint come from the project descriptor
instead of module constants.

Conventions (unchanged from the legacy code, and load-bearing):
  * `z` is elevation, positive up.
  * Facet normals are harmonised into one half-space using the strike hint, so `sigma_n`
    has a consistent sign across the whole fault.
  * Tandem basis: `s = up x n` (strike), `d = s x n` (down-dip).
  * `resolve_tractions` returns compression-POSITIVE sigma_n (the tensor passed in must
    already be the effective, compression-positive one).
  * `s_km` increases toward the strike azimuth, measured from the descriptor's origin.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from deckbuild.config import CRS, ConfigError, GridBox, Hypocenter, MeshSpec, StrikeFrame
from deckbuild.contract import GateReport, HARD, WARN

__all__ = [
    "Fault", "SnappedHypocenter", "GeometryError",
    "strike_s_km", "load_fault", "fault_trace", "resolve_tractions", "build_grid",
    "snap_hypocenter", "project_point", "EPS", "UP",
]

EPS = 1.0e-9                      # legacy: generate_stress_nc_from_raw.EPS
UP = np.array([0.0, 0.0, 1.0])    # legacy: generate_stress_nc_from_raw.UP
# SeisSol tet face -> local vertex indices.  legacy: FACE_VERTS
FACE_VERTS = ((1, 0, 2), (0, 1, 3), (1, 2, 3), (2, 0, 3))


class GeometryError(ValueError):
    """Raised for a malformed mesh or an un-snappable hypocentre."""


# --------------------------------------------------------------------------- the frame
def strike_s_km(x, y, strike: StrikeFrame) -> np.ndarray:
    """Along-strike distance in km from the descriptor's origin.

    VERBATIM from legacy `strike_distance_km`: `su = (sin az, cos az)`, dot with the
    offset, divide by 1000.  A sign flip here silently mirrors every along-strike design,
    so this expression must not be "simplified".  With the SAFS frame it reproduces the
    coefficients that appear literally in every shipped `rs_muw` LuaMap.
    """
    su = np.array([np.sin(np.radians(strike.azimuth_deg)),
                   np.cos(np.radians(strike.azimuth_deg))])
    q = np.stack([np.asarray(x).ravel() - strike.origin_xy[0],
                  np.asarray(y).ravel() - strike.origin_xy[1]], axis=1)
    return (q @ su) / 1000.0


def build_grid(box: GridBox):
    """(gx, gy, gz) for a GridBox.  legacy: _common.build_grid."""
    gx = np.arange(box.xmin, box.xmax + 0.5 * box.dx, box.dx)
    gy = np.arange(box.ymin, box.ymax + 0.5 * box.dx, box.dx)
    gz = np.arange(box.zmin, box.zmax + 0.5 * box.dz, box.dz)
    return gx, gy, gz


# --------------------------------------------------------------------------- the mesh
def _load_puml(path: Path):
    """geometry / connect / per-face BC codes from a PUML HDF5 mesh.

    VERBATIM from legacy `load_puml`, except that an unsupported boundary-format raises
    instead of calling sys.exit (a library must not kill the interpreter).
    """
    import h5py
    with h5py.File(str(path), "r") as f:
        geom = f["geometry"][:]
        conn = f["connect"][:].astype(np.int64)
        fmt = f.attrs.get("boundary-format", "i32")
        if isinstance(fmt, bytes):
            fmt = fmt.decode()
        bnd = f["boundary"][:]
    bits = {"i32": 8, "i64": 16}.get(fmt)
    if np.asarray(bnd).ndim == 2:
        bc = np.asarray(bnd).astype(np.int64)
    elif bits is None:
        raise GeometryError(
            f"{path}: unsupported boundary-format {fmt!r}; expected 'i32' or 'i64'")
    else:
        b = bnd.astype(np.int64)
        mask = (1 << bits) - 1
        bc = np.stack([(b >> (bits * k)) & mask for k in range(4)], axis=1)
    return geom, conn, bc


def _extract_fault(geom, conn, bc, fault_bc: int, path):
    """Deduplicated fault triangles + compacted point array.  legacy: extract_fault."""
    tris = np.concatenate(
        [conn[bc[:, k] == fault_bc][:, FACE_VERTS[k]] for k in range(4)],
        axis=0)
    if len(tris) == 0:
        present = sorted(int(v) for v in np.unique(bc))
        raise GeometryError(
            f"{path}: mesh has no BC-{fault_bc} (dynamic rupture) faces; BC codes present "
            f"are {present}.  Check MeshSpec.fault_bc and MeshSpec.tag_to_bc.")
    _, keep = np.unique(np.sort(tris, axis=1), axis=0, return_index=True)
    tris = tris[np.sort(keep)]
    used, inv = np.unique(tris.ravel(), return_inverse=True)
    return geom[used], inv.reshape(tris.shape)


def _triangle_geometry(points, tris):
    """legacy: triangle_geometry.  Degenerate facets get NaN normals."""
    p = points[tris]
    centroids = p.mean(axis=1)
    cross = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    twice_area = np.linalg.norm(cross, axis=1)
    areas = 0.5 * twice_area
    degen = twice_area < EPS
    normals = cross / np.where(degen, 1.0, twice_area)[:, None]
    normals[degen] = np.nan
    return centroids, normals, areas


def _harmonise_normals(normals, strike_hint_az_deg):
    """Flip facet normals into one half-space.  legacy: harmonise_normals."""
    az = np.radians(strike_hint_az_deg)
    n_global = np.array([-np.cos(az), np.sin(az), 0.0])
    out = normals.copy()
    finite = ~np.isnan(normals).any(axis=1)
    flip = finite & (out @ n_global < 0)
    out[flip] = -out[flip]
    return out, n_global


def _tandem_basis(normals):
    """s = up x n (strike), d = s x n (down-dip).  legacy: tandem_basis."""
    s_raw = np.cross(np.broadcast_to(UP, normals.shape), normals)
    s_norm = np.linalg.norm(s_raw, axis=1)
    degen = ~(s_norm > np.sqrt(2e-6))
    strikes = np.full_like(normals, np.nan)
    dips = np.full_like(normals, np.nan)
    good = ~degen
    strikes[good] = s_raw[good] / s_norm[good, None]
    dips[good] = np.cross(strikes[good], normals[good])
    return strikes, dips, degen


@dataclass(frozen=True, eq=False)
class Fault:
    """Per-facet fault geometry.  The currency between geometry and every stage.

    eq=False: the generated __eq__ would compare numpy arrays field by field and then take
    their truth value, which raises ValueError for any two distinct faults.  Identity
    semantics are the honest default for a bulk array container; use `same_geometry_as`
    for a real comparison.
    """

    cent: np.ndarray        # (N, 3) facet centroids
    normals: np.ndarray     # (N, 3) harmonised unit normals
    strikes: np.ndarray     # (N, 3) unit strike vectors
    dips: np.ndarray        # (N, 3) unit down-dip vectors
    areas: np.ndarray       # (N,)   facet areas, m^2
    mesh_name: str = ""
    n_degenerate: int = 0   # facets with a NaN basis (zero area or ~horizontal)

    def __len__(self) -> int:
        return len(self.cent)

    @property
    def depth_m(self) -> np.ndarray:
        """Positive-down depth of each facet centroid."""
        return -self.cent[:, 2]

    def s_km(self, strike: StrikeFrame) -> np.ndarray:
        return strike_s_km(self.cent[:, 0], self.cent[:, 1], strike)

    def bbox(self) -> tuple[np.ndarray, np.ndarray]:
        return self.cent.min(axis=0), self.cent.max(axis=0)

    def same_geometry_as(self, other: "Fault", tol: float = 0.0) -> bool:
        """True when both faults have identical facet centroids and normals.

        tol=0 is an exact comparison -- what a mesh-identity check wants.  Phase 5's
        "fault unchanged vs the deployed mesh" gate and Phase 8's reproduction diff both
        need this, and both need it not to raise.
        """
        if len(self) != len(other):
            return False
        if tol == 0.0:
            return (np.array_equal(self.cent, other.cent)
                    and np.array_equal(self.normals, other.normals))
        return (np.allclose(self.cent, other.cent, rtol=0, atol=tol)
                and np.allclose(self.normals, other.normals, rtol=0, atol=tol))


def load_fault(mesh: MeshSpec, data_dir: str | Path | None = None,
               strike: StrikeFrame | None = None, mesh_name: str = "") -> Fault:
    """Load the fault facets of a PUML mesh in the Tandem (strike, dip) basis.

    strike -- supplies the harmonisation hint.  Required: without it the normal signs are
              arbitrary and sigma_n would flip sign across the fault.
    """
    if strike is None:
        raise GeometryError(
            "load_fault needs a StrikeFrame for normal harmonisation; pass cfg.strike")
    path = Path(mesh.path)
    if not path.is_absolute() and data_dir is not None:
        path = Path(data_dir) / path
    if not path.is_file():
        raise GeometryError(f"mesh not found: {path.resolve()}")

    geom, conn, bc = _load_puml(path)
    pts, tris = _extract_fault(geom, conn, bc, mesh.fault_bc, path)
    cent, normals_raw, areas = _triangle_geometry(pts, tris)
    normals, _n_global = _harmonise_normals(normals_raw, strike.hint_deg)
    strikes, dips, degen = _tandem_basis(normals)
    return Fault(cent=cent, normals=normals, strikes=strikes, dips=dips, areas=areas,
                 mesh_name=mesh_name or path.stem, n_degenerate=int(degen.sum()))


def fault_trace(fault: Fault, depth_km: float | None = None, tol_km: float = 0.5):
    """(x, y) of the fault facets, for map overlays.

    The fault is a 3-D surface, so its map trace depends on depth.  Pass `depth_km` to get
    a clean curve at that depth; None returns every facet.  legacy: _common.fault_trace.
    """
    fx, fy, fz = fault.cent[:, 0], fault.cent[:, 1], fault.cent[:, 2]
    if depth_km is None:
        return fx, fy
    m = np.abs(-fz / 1000.0 - depth_km) <= tol_km
    return fx[m], fy[m]


def resolve_tractions(sigma_eff, strikes, dips, normals):
    """Project per-facet effective tensors onto each facet.  legacy: resolve_tractions.

    Returns (sigma_n_eff, tau_mag, mu).  sigma_n_eff is compression-POSITIVE; mu is NaN
    where sigma_n_eff <= 0 (a facet in tension has no meaningful apparent friction).
    """
    t = np.einsum("kij,kj->ki", sigma_eff, normals)
    sigma_n_eff = np.einsum("ki,ki->k", normals, t)
    tau_strike = np.einsum("ki,ki->k", strikes, t)
    tau_dip = np.einsum("ki,ki->k", dips, t)
    tau_mag = np.hypot(tau_strike, tau_dip)
    mu = np.full(len(t), np.nan)
    ok = sigma_n_eff > 0
    mu[ok] = tau_mag[ok] / sigma_n_eff[ok]
    return sigma_n_eff, tau_mag, mu


# --------------------------------------------------------------------------- hypocentre
@dataclass(frozen=True)
class SnappedHypocenter:
    """The user's requested point, moved onto the real fault."""

    xyz: tuple[float, float, float]
    requested_xyz: tuple[float, float, float]
    distance_m: float
    facet_index: int
    s_km: float
    depth_m: float
    report: GateReport

    def __repr__(self) -> str:
        return (f"SnappedHypocenter(xyz={self.xyz}, moved {self.distance_m:.1f} m, "
                f"s={self.s_km:.2f} km, depth={self.depth_m:.1f} m)")


def project_point(hypo: Hypocenter, crs: CRS) -> tuple[float, float, float]:
    """The requested point in projected coordinates.

    Geographic input is converted with pyproj (always_xy).  `depth_m` is positive-down and
    is referenced to SEA LEVEL, i.e. z = -depth_m.

    DEVIATION from the plan's wording ("converted to z with the local topo"): the local
    topographic elevation lives on the mesh's free surface, which this function does not
    have.  Sea-level referencing is what the SAFS descriptors already use (their z values
    are elevations), and it is unambiguous.  A topo-relative option would need the free
    surface passed in and is not required by any current project.
    """
    if hypo.is_projected:
        return float(hypo.x), float(hypo.y), float(hypo.z)
    if not hypo.is_geographic:
        raise ConfigError(
            "hypocentre has neither a complete projected (x, y, z) nor a complete "
            "geographic (lon, lat, depth_m) point")
    from pyproj import Transformer
    tf = Transformer.from_crs(crs.geographic_epsg, crs.epsg, always_xy=True)
    x, y = tf.transform(float(hypo.lon), float(hypo.lat))
    return float(x), float(y), -float(hypo.depth_m)


def snap_hypocenter(hypo: Hypocenter, fault: Fault, strike: StrikeFrame, crs: CRS,
                    gate_bands=(), gate: str = "H1") -> SnappedHypocenter:
    """Move the requested hypocentre onto the nearest fault facet centroid.

    Snapping is a named, reported, GATED operation, not a silent nearest-neighbour lookup:
    every downstream number (the Tnuc_s Lua centre, S_E, the friction region the
    nucleation sits in) depends on this point.

    Snaps to a facet CENTROID, not a vertex: the stress and friction fields are evaluated
    per facet, so a facet is the addressable unit.

    HARD-fails past `hypo.snap_tol_m`.  A large snap means the point belongs to a
    different fault strand, or the CRS is wrong -- both are user errors that must not be
    silently absorbed.
    """
    from scipy.spatial import cKDTree

    req = project_point(hypo, crs)
    if len(fault) == 0:
        raise GeometryError("cannot snap a hypocentre onto a fault with no facets")

    tree = cKDTree(fault.cent)
    dist, idx = tree.query(np.asarray(req, float)[None, :], k=1, workers=-1)
    dist = float(np.atleast_1d(dist)[0])
    idx = int(np.atleast_1d(idx)[0])
    snapped = tuple(float(v) for v in fault.cent[idx])

    s_km = float(strike_s_km(snapped[0], snapped[1], strike)[0])
    depth_m = -snapped[2]

    rep = GateReport(f"hypocentre snap ({hypo.label or 'unlabelled'})")
    ok = dist <= hypo.snap_tol_m
    lo, hi = fault.bbox()
    rep.add(gate, ok,
            f"requested ({req[0]:.1f}, {req[1]:.1f}, {req[2]:.1f}) -> facet {idx} at "
            f"({snapped[0]:.1f}, {snapped[1]:.1f}, {snapped[2]:.1f}); moved "
            f"{dist:.1f} m (tol {hypo.snap_tol_m:.0f} m); s = {s_km:.2f} km, "
            f"depth = {depth_m:.1f} m")
    if not ok:
        rep.add(f"{gate}ctx", False,
                f"fault bbox x[{lo[0]:.0f}, {hi[0]:.0f}] y[{lo[1]:.0f}, {hi[1]:.0f}] "
                f"z[{lo[2]:.0f}, {hi[2]:.0f}] -- a snap this large usually means the "
                f"point is on a different strand, or the CRS ({crs.epsg}) is wrong",
                severity=HARD)

    # Did the snap move the point across a NAMED band?  That silently changes which
    # friction region the nucleation sits in, which no downstream gate would catch.
    s_req = float(strike_s_km(req[0], req[1], strike)[0])
    for band in gate_bands:
        in_req = band.s_start_km <= s_req <= band.s_end_km
        in_snap = band.s_start_km <= s_km <= band.s_end_km
        if in_req != in_snap:
            rep.add(f"{gate}band", False,
                    f"the snap moved the hypocentre {'out of' if in_req else 'into'} "
                    f"band {band.name!r} (s {s_req:.2f} -> {s_km:.2f} km); this changes "
                    f"which friction region the nucleation sits in",
                    severity=WARN)

    out = SnappedHypocenter(xyz=snapped, requested_xyz=req, distance_m=dist,
                            facet_index=idx, s_km=s_km, depth_m=depth_m, report=rep)
    rep.raise_if_failed()
    return out
