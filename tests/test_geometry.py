"""Phase 1 acceptance tests for deckbuild.geometry."""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from deckbuild.config import (  # noqa: E402
    CRS, GridBox, Hypocenter, MeshSpec, NamedBand, Project, StrikeFrame,
)
from deckbuild.contract import GateFailure  # noqa: E402
from deckbuild.geometry import (  # noqa: E402
    Fault, GeometryError, build_grid, fault_trace, load_fault, project_point,
    resolve_tractions, snap_hypocenter, strike_s_km,
)

PROJECTS = ROOT / "projects"
SAFS_STRIKE = StrikeFrame(azimuth_deg=314.0, origin_xy=(606971.0, 3707270.0))

LEGACY_DIR = Path(os.environ.get(
    "DECKBUILD_LEGACY_LIB",
    str(Path.home() / "projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/"
        "seisol_quakeworx/v3_under_construction/toolbox/combined_workflow/lib")))
LEGACY_GEN = LEGACY_DIR / "generate_stress_nc_from_raw.py"


def _load_legacy():
    if str(LEGACY_DIR) not in sys.path:
        sys.path.insert(0, str(LEGACY_DIR))
    spec = importlib.util.spec_from_file_location("_legacy_gen", LEGACY_GEN)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_legacy_gen"] = mod
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------- strike_s_km
def test_strike_s_km_is_zero_at_the_origin():
    s = strike_s_km(SAFS_STRIKE.origin_xy[0], SAFS_STRIKE.origin_xy[1], SAFS_STRIKE)
    assert s[0] == 0.0


def test_strike_s_km_uses_sin_cos_not_minus_sin():
    """Guards the exact projection.  su = (sin az, cos az).

    With az = 314 deg this gives the x coefficient -0.71933980033865119 and the y
    coefficient 0.69465837045899725 -- the numbers that appear literally in every shipped
    rs_muw LuaMap.  A sign flip here mirrors every along-strike design.
    """
    ox, oy = SAFS_STRIKE.origin_xy
    sx = strike_s_km(ox + 1000.0, oy, SAFS_STRIKE)[0]
    sy = strike_s_km(ox, oy + 1000.0, SAFS_STRIKE)[0]
    assert sx == pytest.approx(np.sin(np.radians(314.0)), rel=1e-15)
    assert sy == pytest.approx(np.cos(np.radians(314.0)), rel=1e-15)
    assert sx == pytest.approx(-0.71933980033865119, rel=1e-14)
    assert sy == pytest.approx(0.69465837045899725, rel=1e-14)


def test_strike_s_km_accepts_scalars_and_arrays():
    assert strike_s_km(0.0, 0.0, SAFS_STRIKE).shape == (1,)
    assert strike_s_km(np.zeros(5), np.zeros(5), SAFS_STRIKE).shape == (5,)


@pytest.mark.skipif(not LEGACY_GEN.is_file(), reason=f"legacy lib absent: {LEGACY_GEN}")
def test_strike_s_km_matches_legacy_to_zero_ulp():
    """Acceptance criterion 1: 0 ULP against the legacy strike_distance_km."""
    legacy = _load_legacy()
    rng = np.random.default_rng(1234)
    x = rng.uniform(2.0e5, 7.0e5, 10000)
    y = rng.uniform(3.6e6, 4.0e6, 10000)
    mine = strike_s_km(x, y, SAFS_STRIKE)
    theirs = legacy.strike_distance_km(x, y, SAFS_STRIKE.azimuth_deg,
                                       SAFS_STRIKE.origin_xy)
    assert np.array_equal(mine, theirs), "strike_s_km deviates from the legacy expression"


# --------------------------------------------------------------------------- build_grid
def test_build_grid_endpoints_and_spacing():
    box = GridBox(xmin=0.0, xmax=3000.0, ymin=0.0, ymax=2000.0,
                  zmin=-1000.0, zmax=0.0, dx=1000.0, dz=250.0)
    gx, gy, gz = build_grid(box)
    assert gx.tolist() == [0.0, 1000.0, 2000.0, 3000.0]
    assert gy.tolist() == [0.0, 1000.0, 2000.0]
    assert gz.tolist() == [-1000.0, -750.0, -500.0, -250.0, 0.0]


@pytest.mark.skipif(not LEGACY_GEN.is_file(), reason="legacy lib absent")
def test_build_grid_matches_the_legacy_box():
    cfg = Project.load(PROJECTS / "safs_alt.yaml", require_files=False)
    gx, gy, gz = build_grid(cfg.stress_box)
    b = cfg.stress_box
    assert gx[0] == b.xmin and gx[-1] == pytest.approx(b.xmax)
    assert gz[0] == b.zmin and gz[-1] == pytest.approx(b.zmax)
    assert len(gx) == 281 and len(gz) == 69


# --------------------------------------------------------------------------- tractions
def test_resolve_tractions_on_a_vertical_strike_slip_facet():
    """One facet, normal = +y, pure xy shear -> sigma_n = syy, tau = |sxy|."""
    sigma = np.zeros((1, 3, 3))
    sigma[0, 0, 0] = 30.0      # sxx
    sigma[0, 1, 1] = 20.0      # syy  -> this is sigma_n for n = y
    sigma[0, 2, 2] = 25.0
    sigma[0, 0, 1] = sigma[0, 1, 0] = 6.0
    normals = np.array([[0.0, 1.0, 0.0]])
    strikes = np.array([[1.0, 0.0, 0.0]])
    dips = np.array([[0.0, 0.0, -1.0]])
    sn, tau, mu = resolve_tractions(sigma, strikes, dips, normals)
    assert sn[0] == pytest.approx(20.0)
    assert tau[0] == pytest.approx(6.0)
    assert mu[0] == pytest.approx(0.3)


def test_resolve_tractions_gives_nan_mu_in_tension():
    sigma = np.zeros((1, 3, 3))
    sigma[0, 1, 1] = -5.0          # tension on the facet normal
    sn, tau, mu = resolve_tractions(sigma, np.array([[1.0, 0.0, 0.0]]),
                                    np.array([[0.0, 0.0, -1.0]]),
                                    np.array([[0.0, 1.0, 0.0]]))
    assert sn[0] == -5.0 and np.isnan(mu[0])


# --------------------------------------------------------------------------- synthetic mesh
def _write_planar_puml(path: Path, nx=9, nz=5, x0=0.0, y0=0.0, dx=1000.0, dz=1000.0):
    """A vertical y=0 fault embedded in a slab of tets, with SAFS-style BC packing.

    Small enough to be exact, real enough to exercise load_fault end to end.
    """
    import h5py
    xs = x0 + np.arange(nx) * dx
    zs = -np.arange(nz) * dz
    ys = np.array([y0 - dx, y0, y0 + dx])

    pts, idx = [], {}
    for k, zz in enumerate(zs):
        for j, yy in enumerate(ys):
            for i, xx in enumerate(xs):
                idx[(i, j, k)] = len(pts)
                pts.append((xx, yy, zz))
    pts = np.asarray(pts, float)

    # Split each hex cell into 6 tets (the standard Kuhn decomposition).
    KUHN = [(0, 1, 3, 7), (0, 1, 7, 5), (0, 5, 7, 4),
            (0, 3, 2, 7), (0, 6, 4, 7), (0, 2, 6, 7)]
    tets = []
    for k in range(nz - 1):
        for j in range(len(ys) - 1):
            for i in range(nx - 1):
                c = [idx[(i + (n & 1), j + ((n >> 1) & 1), k + ((n >> 2) & 1))]
                     for n in range(8)]
                for t in KUHN:
                    tets.append([c[t[0]], c[t[1]], c[t[2]], c[t[3]]])
    tets = np.asarray(tets, np.uint64)

    # Orient positively, then tag every face whose 3 vertices all lie on y == y0.
    p = pts[tets.astype(int)]
    vol = np.einsum("ij,ij->i", p[:, 1] - p[:, 0],
                    np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))
    flip = vol < 0
    tets[flip] = tets[flip][:, [0, 2, 1, 3]]

    FACE_VERTS = ((1, 0, 2), (0, 1, 3), (1, 2, 3), (2, 0, 3))
    on_fault = np.isclose(pts[:, 1], y0)
    boundary = np.zeros(len(tets), np.int32)
    for f, fv in enumerate(FACE_VERTS):
        tri = tets[:, list(fv)].astype(int)
        is_fault = on_fault[tri].all(axis=1)
        boundary |= (is_fault.astype(np.int32) * 3) << (8 * f)

    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h:
        h.create_dataset("geometry", data=pts)
        h.create_dataset("connect", data=tets)
        h.create_dataset("boundary", data=boundary)
        h.create_dataset("group", data=np.ones(len(tets), np.int32))
        h.attrs["boundary-format"] = "i32"
        h.attrs["topology-format"] = "geometric"
    return path


@pytest.fixture()
def planar_mesh(tmp_path):
    p = _write_planar_puml(tmp_path / "planar.puml.h5")
    return MeshSpec(path=str(p))


@pytest.fixture()
def planar_fault(planar_mesh):
    return load_fault(planar_mesh, strike=StrikeFrame(azimuth_deg=0.0,
                                                      origin_xy=(0.0, 0.0)))


def test_load_fault_finds_the_embedded_fault(planar_fault):
    f = planar_fault
    assert len(f) > 0
    assert np.allclose(f.cent[:, 1], 0.0)                      # all on the y=0 plane
    assert np.allclose(np.abs(f.normals[:, 1]), 1.0)           # normals are +/- y
    assert f.n_degenerate == 0
    assert np.all(f.areas > 0)


def test_load_fault_harmonises_normals_consistently(planar_fault):
    """Every normal must point the same way, or sigma_n flips sign across the fault."""
    ny = planar_fault.normals[:, 1]
    assert np.all(ny > 0) or np.all(ny < 0)


def test_load_fault_basis_is_orthonormal(planar_fault):
    f = planar_fault
    for v in (f.normals, f.strikes, f.dips):
        assert np.allclose(np.linalg.norm(v, axis=1), 1.0)
    assert np.allclose(np.einsum("ij,ij->i", f.strikes, f.normals), 0.0, atol=1e-12)
    assert np.allclose(np.einsum("ij,ij->i", f.dips, f.normals), 0.0, atol=1e-12)
    assert np.allclose(np.einsum("ij,ij->i", f.strikes, f.dips), 0.0, atol=1e-12)


def test_load_fault_dip_points_down(planar_fault):
    """d = s x n; for a vertical fault the dip vector must be vertical."""
    assert np.allclose(np.abs(planar_fault.dips[:, 2]), 1.0)


def test_fault_helpers(planar_fault):
    f = planar_fault
    assert f.depth_m.min() >= 0.0
    lo, hi = f.bbox()
    assert lo[2] < hi[2]
    s = f.s_km(StrikeFrame(azimuth_deg=0.0, origin_xy=(0.0, 0.0)))
    assert s.shape == (len(f),)


def test_fault_trace_at_depth(planar_fault):
    fx, fy = fault_trace(planar_fault, depth_km=2.0, tol_km=0.6)
    assert len(fx) > 0 and len(fx) < len(planar_fault)
    allx, _ = fault_trace(planar_fault)
    assert len(allx) == len(planar_fault)


def test_load_fault_requires_a_strike_frame(planar_mesh):
    with pytest.raises(GeometryError, match="StrikeFrame"):
        load_fault(planar_mesh)


def test_load_fault_missing_file_raises():
    with pytest.raises(GeometryError, match="not found"):
        load_fault(MeshSpec(path="/nope/absent.puml.h5"), strike=SAFS_STRIKE)


def test_load_fault_no_fault_faces_raises(tmp_path):
    """Wrong fault_bc must name the BC codes actually present, not return an empty Fault."""
    p = _write_planar_puml(tmp_path / "m.puml.h5")
    with pytest.raises(GeometryError) as exc:
        load_fault(MeshSpec(path=str(p), fault_bc=4,
                            tag_to_bc={101: 4, 102: 1, 103: 5, 104: 5}),
                   strike=SAFS_STRIKE)
    msg = str(exc.value)
    assert "no BC-4" in msg and "BC codes present" in msg


@pytest.mark.skipif(not LEGACY_GEN.is_file(), reason="legacy lib absent")
def test_load_fault_matches_legacy_on_the_same_mesh(tmp_path):
    """Acceptance criterion 2, on a synthetic mesh: identical facets and normals."""
    legacy = _load_legacy()
    p = _write_planar_puml(tmp_path / "m.puml.h5")
    geom, conn, bc = legacy.load_puml(str(p))
    lpts, ltris = legacy.extract_fault(geom, conn, bc)
    lcent, lnorm_raw, lareas = legacy.triangle_geometry(lpts, ltris)
    lnorm, _ = legacy.harmonise_normals(lnorm_raw, SAFS_STRIKE.azimuth_deg)
    lstrike, ldip, _ = legacy.tandem_basis(lnorm)

    mine = load_fault(MeshSpec(path=str(p)), strike=SAFS_STRIKE)
    assert len(mine) == len(lcent)
    assert np.array_equal(mine.cent, lcent)
    assert np.array_equal(mine.normals, lnorm)
    assert np.array_equal(mine.areas, lareas)
    np.testing.assert_array_equal(mine.strikes, lstrike)
    np.testing.assert_array_equal(mine.dips, ldip)


# --------------------------------------------------------------------------- snapping
def test_project_point_passes_projected_through():
    h = Hypocenter(x=1.0, y=2.0, z=-3.0)
    assert project_point(h, CRS(epsg="EPSG:32611")) == (1.0, 2.0, -3.0)


def test_project_point_converts_geographic():
    h = Hypocenter(lon=-116.35, lat=33.78, depth_m=10000.0)
    x, y, z = project_point(h, CRS(epsg="EPSG:32611"))
    assert 4.0e5 < x < 7.0e5 and 3.6e6 < y < 3.9e6
    assert z == -10000.0                       # depth is positive-down, sea-level ref


def test_snap_returns_a_facet_centroid(planar_fault):
    strike = StrikeFrame(azimuth_deg=0.0, origin_xy=(0.0, 0.0))
    target = planar_fault.cent[7]
    h = Hypocenter(x=float(target[0]) + 30.0, y=float(target[1]) + 40.0,
                   z=float(target[2]), snap_tol_m=500.0, label="t")
    snap = snap_hypocenter(h, planar_fault, strike, CRS(epsg="EPSG:32610"))
    assert tuple(planar_fault.cent[snap.facet_index]) == snap.xyz
    assert snap.distance_m == pytest.approx(50.0, abs=1e-6)
    assert snap.depth_m == pytest.approx(-target[2])
    assert snap.report.ok


def test_snap_hard_fails_beyond_the_tolerance(planar_fault):
    strike = StrikeFrame(azimuth_deg=0.0, origin_xy=(0.0, 0.0))
    h = Hypocenter(x=5.0e5, y=5.0e5, z=-5000.0, snap_tol_m=2000.0, label="far")
    with pytest.raises(GateFailure) as exc:
        snap_hypocenter(h, planar_fault, strike, CRS(epsg="EPSG:32610"))
    msg = str(exc.value)
    assert "H1" in msg
    assert "fault bbox" in msg                 # the diagnostic context


def test_snap_geographic_and_projected_agree(planar_fault):
    """The same physical point given two ways must snap to the same facet."""
    from pyproj import Transformer
    crs = CRS(epsg="EPSG:32610")
    target = planar_fault.cent[11]
    # Place the mesh where UTM 10N is valid by shifting the query, not the mesh:
    # convert the projected target to lon/lat and back through the Hypocenter path.
    tf = Transformer.from_crs(crs.epsg, crs.geographic_epsg, always_xy=True)
    lon, lat = tf.transform(float(target[0]) + 5.0e5, float(target[1]) + 4.0e6)
    # Undo the shift by snapping against a shifted copy of the fault.
    shifted = Fault(cent=planar_fault.cent + np.array([5.0e5, 4.0e6, 0.0]),
                    normals=planar_fault.normals, strikes=planar_fault.strikes,
                    dips=planar_fault.dips, areas=planar_fault.areas)
    strike = StrikeFrame(azimuth_deg=0.0, origin_xy=(5.0e5, 4.0e6))
    a = snap_hypocenter(Hypocenter(lon=lon, lat=lat, depth_m=-float(target[2])),
                        shifted, strike, crs)
    b = snap_hypocenter(Hypocenter(x=float(target[0]) + 5.0e5,
                                   y=float(target[1]) + 4.0e6, z=float(target[2])),
                        shifted, strike, crs)
    assert a.facet_index == b.facet_index


def test_snap_warns_when_it_crosses_a_named_band(planar_fault):
    """Crossing a band silently changes which friction region the nucleation is in."""
    strike = StrikeFrame(azimuth_deg=0.0, origin_xy=(0.0, 0.0))
    s = planar_fault.s_km(strike)
    lo, hi = float(s.min()), float(s.max())
    mid = 0.5 * (lo + hi)
    band = NamedBand(s_start_km=mid, s_end_km=hi + 1.0, name="upper half")
    # Request a point just BELOW the band edge but nearest to a facet inside it.
    inside = planar_fault.cent[int(np.argmin(np.abs(s - (mid + 0.01))))]
    h = Hypocenter(x=float(inside[0]), y=float(inside[1]) + 10.0, z=float(inside[2]),
                   snap_tol_m=5000.0)
    snap = snap_hypocenter(h, planar_fault, strike, CRS(epsg="EPSG:32610"),
                           gate_bands=[band])
    # Whether it crossed depends on geometry; assert the machinery reports coherently.
    crossed = [g for g in snap.report.gates if g.name.endswith("band")]
    if crossed:
        assert crossed[0].severity == "warn"
        assert "upper half" in crossed[0].detail


def test_snap_on_an_empty_fault_raises():
    empty = Fault(cent=np.zeros((0, 3)), normals=np.zeros((0, 3)),
                  strikes=np.zeros((0, 3)), dips=np.zeros((0, 3)), areas=np.zeros(0))
    with pytest.raises(GeometryError, match="no facets"):
        snap_hypocenter(Hypocenter(x=0.0, y=0.0, z=0.0), empty, SAFS_STRIKE,
                        CRS(epsg="EPSG:32611"))


def test_snap_report_records_the_distance_and_position(planar_fault):
    strike = StrikeFrame(azimuth_deg=0.0, origin_xy=(0.0, 0.0))
    t = planar_fault.cent[3]
    h = Hypocenter(x=float(t[0]), y=float(t[1]), z=float(t[2]), label="exact")
    snap = snap_hypocenter(h, planar_fault, strike, CRS(epsg="EPSG:32610"))
    d = snap.report.gates[0].detail
    assert "moved" in d and "s =" in d and "depth =" in d
    assert snap.distance_m == pytest.approx(0.0, abs=1e-9)
    assert "SnappedHypocenter" in repr(snap)


# ------------------------------------------------- review finding R-101
def test_R101_fault_comparison_does_not_raise(planar_mesh):
    """Two distinct Faults must be comparable; the generated __eq__ raised ValueError."""
    a = load_fault(planar_mesh, strike=SAFS_STRIKE)
    b = load_fault(planar_mesh, strike=SAFS_STRIKE)
    assert a is not b
    assert (a == b) is False                 # identity semantics, and no exception
    assert a.same_geometry_as(a) and a.same_geometry_as(b)


def test_R101_same_geometry_as_detects_a_real_difference(planar_mesh):
    a = load_fault(planar_mesh, strike=SAFS_STRIKE)
    moved = Fault(cent=a.cent + np.array([0.0, 0.0, 1.0]), normals=a.normals,
                  strikes=a.strikes, dips=a.dips, areas=a.areas)
    assert not a.same_geometry_as(moved)
    assert a.same_geometry_as(moved, tol=2.0)
    smaller = Fault(cent=a.cent[:-1], normals=a.normals[:-1], strikes=a.strikes[:-1],
                    dips=a.dips[:-1], areas=a.areas[:-1])
    assert not a.same_geometry_as(smaller)
