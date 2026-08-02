"""Phase 1 acceptance tests for deckbuild.asagi and the comparison harness."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from deckbuild.asagi import (  # noqa: E402
    AsagiError, asagi_axes, asagi_fields, axis_spacing, hull_containment, read_asagi,
    roundtrip_selfcheck, trilinear_sample, write_asagi,
)
from nc_compare import compare_nc  # noqa: E402


# --------------------------------------------------------------------------- fixtures
def make_grid(nx=7, ny=5, nz=4, dx=1000.0, dz=250.0, x0=0.0, y0=0.0, z0=-3000.0):
    x = np.arange(nx, dtype=np.float64) * dx + x0
    y = np.arange(ny, dtype=np.float64) * dx + y0
    z = np.arange(nz, dtype=np.float64) * dz + z0
    return x, y, z


def linear_field(x, y, z, a=2.0, b=-3.0, c=0.5, d=7.0):
    """f = a*X + b*Y + c*Z + d on the STORED (nz, ny, nx) grid.

    Linear, so its trilinear interpolant is exact everywhere -- the round-trip check can
    then use a tight tolerance without any kink caveat.
    """
    Z, Y, X = np.meshgrid(z, y, x, indexing="ij")
    return a * X + b * Y + c * Z + d


@pytest.fixture()
def nc(tmp_path):
    x, y, z = make_grid()
    f = linear_field(x, y, z)
    p = write_asagi(tmp_path / "g.nc", x, y, z, {"f": f, "g": 2.0 * f},
                    attrs={"note": "test"})
    return p, x, y, z, f


# --------------------------------------------------------------------------- axes
def test_axis_spacing_accepts_a_uniform_axis():
    a = np.arange(10, dtype=np.float64) * 250.0
    assert axis_spacing(a, "z") == pytest.approx(250.0)


def test_axis_spacing_rejects_float32():
    with pytest.raises(AsagiError, match="float64"):
        axis_spacing(np.arange(5, dtype=np.float32), "x")


def test_axis_spacing_rejects_non_increasing():
    a = np.array([0.0, 1000.0, 500.0, 2000.0])
    with pytest.raises(AsagiError, match="strictly increasing"):
        axis_spacing(a, "x")


def test_axis_spacing_rejects_non_equidistant():
    """ASAGI assumes uniform spacing and does not check -- so we must."""
    a = np.array([0.0, 1000.0, 2000.0, 3500.0])
    with pytest.raises(AsagiError) as exc:
        axis_spacing(a, "y")
    msg = str(exc.value)
    assert "not equidistant" in msg and "'y'" in msg


def test_axis_spacing_rejects_too_short():
    with pytest.raises(AsagiError, match=">= 2"):
        axis_spacing(np.array([1.0]), "z")


def test_arange_drift_is_tolerated():
    """np.arange accumulates ~1e-16 relative drift; that must not trip the check."""
    a = np.arange(-45000.0, 3251.0, 250.0)
    assert axis_spacing(a, "z") == pytest.approx(250.0)


# --------------------------------------------------------------------------- write/read
def test_write_read_roundtrip_is_exact(nc):
    p, x, y, z, f = nc
    rx, ry, rz, flds, attrs = read_asagi(p)
    assert np.array_equal(rx, x) and np.array_equal(ry, y) and np.array_equal(rz, z)
    assert sorted(flds) == ["f", "g"]
    # float32 storage: compare against the float32 cast, which must be exact.
    assert np.array_equal(flds["f"], f.astype(np.float32).astype(float))
    assert attrs["note"] == "test"


def test_written_dims_are_z_y_x(nc):
    from netCDF4 import Dataset
    p = nc[0]
    with Dataset(str(p)) as ds:
        assert ds.variables["data"].dimensions == ("z", "y", "x")
        assert ds.variables["x"].dtype == np.float64


def test_asagi_fields_lists_members(nc):
    assert asagi_fields(nc[0]) == ["f", "g"]


def test_asagi_axes(nc):
    p, x, y, z, _ = nc
    ax, ay, az = asagi_axes(p)
    assert np.array_equal(ax, x) and np.array_equal(az, z)


def test_write_rejects_transposed_field(tmp_path):
    """(nx,ny,nz) instead of (nz,ny,nx) is easy and catastrophic; it must not be written."""
    x, y, z = make_grid()
    bad = np.zeros((len(x), len(y), len(z)))
    with pytest.raises(AsagiError) as exc:
        write_asagi(tmp_path / "b.nc", x, y, z, {"f": bad})
    msg = str(exc.value)
    assert "(7, 5, 4)" in msg and "(4, 5, 7)" in msg
    assert "np.transpose" in msg          # the actionable hint


def test_write_rejects_wrong_shape(tmp_path):
    x, y, z = make_grid()
    with pytest.raises(AsagiError, match="expected"):
        write_asagi(tmp_path / "b.nc", x, y, z, {"f": np.zeros((3, 3, 3))})


def test_write_rejects_no_fields(tmp_path):
    x, y, z = make_grid()
    with pytest.raises(AsagiError, match="at least one field"):
        write_asagi(tmp_path / "b.nc", x, y, z, {})


def test_write_rejects_non_finite(tmp_path):
    x, y, z = make_grid()
    f = linear_field(x, y, z)
    f[0, 0, 0] = np.nan
    with pytest.raises(AsagiError, match="non-finite"):
        write_asagi(tmp_path / "b.nc", x, y, z, {"f": f})


def test_write_rejects_non_equidistant_axis(tmp_path):
    x, y, z = make_grid()
    x = x.copy()
    x[-1] += 500.0
    with pytest.raises(AsagiError, match="not equidistant"):
        write_asagi(tmp_path / "b.nc", x, y, z,
                    {"f": np.zeros((len(z), len(y), len(x)))})


# --------------------------------------------------------------------------- sampling
def test_trilinear_is_exact_on_a_linear_field(nc):
    p, x, y, z, _ = nc
    rng = np.random.default_rng(0)
    n = 500
    qx = rng.uniform(x[0], x[-1], n)
    qy = rng.uniform(y[0], y[-1], n)
    qz = rng.uniform(z[0], z[-1], n)
    got = trilinear_sample(p, qx, qy, qz, fields=["f"])["f"]
    want = 2.0 * qx - 3.0 * qy + 0.5 * qz + 7.0
    # float32 storage of values ~1e4 gives ~1e-3 absolute resolution.
    assert np.allclose(got, want, rtol=0, atol=5e-2)


def test_trilinear_hits_node_values_exactly(nc):
    p, x, y, z, f = nc
    got = trilinear_sample(p, [x[2]], [y[1]], [z[3]], fields=["f"])["f"]
    assert got[0] == pytest.approx(float(f[3, 1, 2].astype(np.float32)), rel=0, abs=1e-3)


def test_trilinear_edge_clamps_and_reports(nc):
    p, x, y, z, _ = nc
    qx = np.array([x[0] - 5000.0, x[2]])
    qy = np.array([y[0], y[1]])
    qz = np.array([z[0], z[1]])
    out, clamped = trilinear_sample(p, qx, qy, qz, fields=["f"], return_clamped=True)
    assert clamped.tolist() == [True, False]
    assert np.isfinite(out["f"]).all()


def test_trilinear_rejects_ragged_queries(nc):
    p = nc[0]
    with pytest.raises(AsagiError, match="same length"):
        trilinear_sample(p, [0.0, 1.0], [0.0], [0.0])


def test_read_rejects_unknown_field(nc):
    with pytest.raises(AsagiError, match="no such field"):
        read_asagi(nc[0], fields=["nope"])


# --------------------------------------------------------------------------- gates
def test_hull_containment_passes_inside(nc):
    p, x, y, z, _ = nc
    pts = np.array([[x[1], y[1], z[1]], [x[-2], y[-2], z[-2]]])
    rep = hull_containment(p, pts, label="fault")
    assert rep.ok
    assert "margin" in rep.gates[0].detail


def test_hull_containment_fails_outside_and_counts(nc):
    p, x, y, z, _ = nc
    pts = np.array([[x[0] - 1.0, y[1], z[1]], [x[1], y[1], z[1]]])
    rep = hull_containment(p, pts)
    assert not rep.ok
    assert "1 point(s) OUTSIDE" in rep.gates[0].detail


def test_hull_containment_rejects_bad_shape(nc):
    with pytest.raises(AsagiError, match=r"\(N, 3\)"):
        hull_containment(nc[0], np.zeros((4, 2)))


def test_roundtrip_selfcheck_passes_on_a_linear_field(nc):
    p = nc[0]

    def ev(pts):
        return 2.0 * pts[:, 0] - 3.0 * pts[:, 1] + 0.5 * pts[:, 2] + 7.0

    rep = roundtrip_selfcheck(p, ev, field="f", n=200, median_tol=1e-1, max_tol=1e-1)
    assert rep.ok


def test_roundtrip_selfcheck_fails_on_a_wrong_evaluator(nc):
    p = nc[0]
    rep = roundtrip_selfcheck(p, lambda pts: np.zeros(len(pts)), field="f", n=100)
    assert not rep.ok


def test_roundtrip_selfcheck_rejects_bad_evaluator_shape(nc):
    with pytest.raises(AsagiError, match="evaluator returned shape"):
        roundtrip_selfcheck(nc[0], lambda pts: np.zeros((len(pts), 2)), field="f", n=10)


# --------------------------------------------------------------------------- compare_nc
def test_compare_identical_files(tmp_path):
    x, y, z = make_grid()
    f = linear_field(x, y, z)
    a = write_asagi(tmp_path / "a.nc", x, y, z, {"f": f})
    b = write_asagi(tmp_path / "b.nc", x, y, z, {"f": f})
    d = compare_nc(a, b)
    assert d.data_identical
    assert d.first_diff is None


def test_compare_detects_a_single_changed_value(tmp_path):
    x, y, z = make_grid()
    f = linear_field(x, y, z)
    g = f.copy()
    g[2, 1, 3] += 1.0
    a = write_asagi(tmp_path / "a.nc", x, y, z, {"f": f})
    b = write_asagi(tmp_path / "b.nc", x, y, z, {"f": g})
    d = compare_nc(a, b)
    assert not d.data_identical
    assert d.field_diffs["f"][2] == 1                 # exactly one differing value
    assert d.first_diff[0] == "f"
    assert d.field_diffs["f"][0] == pytest.approx(1.0, abs=1e-3)


def test_compare_detects_a_shifted_axis(tmp_path):
    x, y, z = make_grid()
    f = linear_field(x, y, z)
    a = write_asagi(tmp_path / "a.nc", x, y, z, {"f": f})
    b = write_asagi(tmp_path / "b.nc", x + 1.0, y, z, {"f": f})
    d = compare_nc(a, b)
    assert not d.data_identical and d.axis_diffs


def test_compare_detects_a_missing_field(tmp_path):
    x, y, z = make_grid()
    f = linear_field(x, y, z)
    a = write_asagi(tmp_path / "a.nc", x, y, z, {"f": f, "g": f})
    b = write_asagi(tmp_path / "b.nc", x, y, z, {"f": f})
    d = compare_nc(a, b)
    assert not d.data_identical and d.missing_in_b == ["g"]


def test_compare_attrs_can_be_ignored(tmp_path):
    x, y, z = make_grid()
    f = linear_field(x, y, z)
    a = write_asagi(tmp_path / "a.nc", x, y, z, {"f": f}, attrs={"built": "monday"})
    b = write_asagi(tmp_path / "b.nc", x, y, z, {"f": f}, attrs={"built": "tuesday"})
    assert not compare_nc(a, b).data_identical
    d = compare_nc(a, b, ignore_attrs=("built",))
    assert d.data_identical
    assert "built" in d.attr_diffs        # still reported


def test_compare_reports_a_missing_file(tmp_path):
    x, y, z = make_grid()
    a = write_asagi(tmp_path / "a.nc", x, y, z, {"f": linear_field(x, y, z)})
    d = compare_nc(a, tmp_path / "nope.nc")
    assert d.error and not d.data_identical


def test_summary_is_informative(tmp_path):
    x, y, z = make_grid()
    f = linear_field(x, y, z)
    g = f.copy(); g[0, 0, 0] += 5.0
    a = write_asagi(tmp_path / "a.nc", x, y, z, {"f": f})
    b = write_asagi(tmp_path / "b.nc", x, y, z, {"f": g})
    s = compare_nc(a, b).summary()
    assert "data_identical=False" in s and "first difference" in s


# ------------------------------------------------- review findings R-102 / R-103
def test_R102_large_smooth_excursion_is_a_hard_fail(nc):
    """The legacy G3 max gate hard-fails an excursion with no kink excuse.

    Before the fix the max gate was unconditionally WARN, so it could never fail and a
    genuine interpolation bug in a smooth region would have passed.
    """
    rep = roundtrip_selfcheck(nc[0], lambda pts: np.zeros(len(pts)), field="f",
                              n=50, median_tol=1e9, max_tol=1e-6)
    assert not rep.ok
    bad = [g for g in rep.gates if g.name == "G3max"]
    assert bad and bad[0].severity == "hard"
    assert "no kink_straddling attestation" in bad[0].detail


def test_R102_kink_attestation_downgrades_to_warn(nc):
    rep = roundtrip_selfcheck(nc[0], lambda pts: np.zeros(len(pts)), field="f",
                              n=50, median_tol=1e9, max_tol=1e-6,
                              kink_straddling=lambda pts: np.ones(len(pts), bool))
    assert rep.ok
    g = [g for g in rep.gates if g.name == "G3max"][0]
    assert g.severity == "warn" and "straddling a kink" in g.detail


def test_R102_partial_attestation_still_fails(nc):
    """If even one excursion is NOT excused, the gate must fail."""
    rep = roundtrip_selfcheck(
        nc[0], lambda pts: np.zeros(len(pts)), field="f", n=50,
        median_tol=1e9, max_tol=1e-6,
        kink_straddling=lambda pts: np.arange(len(pts)) > 0)   # first one unexcused
    assert not rep.ok


def test_R105_nonfinite_can_be_opted_into(tmp_path):
    x, y, z = make_grid()
    f = linear_field(x, y, z)
    f[0, 0, 0] = np.nan
    with pytest.raises(AsagiError, match="allow_nonfinite=True"):
        write_asagi(tmp_path / "a.nc", x, y, z, {"f": f})
    p = write_asagi(tmp_path / "b.nc", x, y, z, {"f": f}, allow_nonfinite=True)
    assert np.isnan(read_asagi(p)[3]["f"][0, 0, 0])


# --- acceptance criterion 5: 0 ULP vs legacy on a SHIPPED nc -----------------
import os  # noqa: E402

_DECKS = Path(os.environ.get("DECKBUILD_DECKS",
                             str(Path.home() / "Downloads" / "seisol_quakeworx")))
_SHIPPED = next(_DECKS.glob("*/safs_stress_andersonian_k1.7.nc"), None) \
    if _DECKS.is_dir() else None
_LEGACY_GEN = Path(os.environ.get(
    "DECKBUILD_LEGACY_LIB",
    str(Path.home() / "projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/"
        "seisol_quakeworx/v3_under_construction/toolbox/combined_workflow/lib"))
) / "generate_stress_nc_from_raw.py"


@pytest.mark.skipif(_SHIPPED is None or not _LEGACY_GEN.is_file(),
                    reason="shipped deck or legacy lib not available")
def test_R103_trilinear_matches_legacy_on_the_shipped_nc():
    """Acceptance criterion 5, on the real file: 0 ULP at 1000 fixed-seed points."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("_lg", _LEGACY_GEN)
    lg = importlib.util.module_from_spec(spec)
    sys.modules["_lg"] = lg
    spec.loader.exec_module(lg)

    x, y, z = asagi_axes(_SHIPPED)
    rng = np.random.default_rng(12345)
    qx = rng.uniform(x[0], x[-1], 1000)
    qy = rng.uniform(y[0], y[-1], 1000)
    qz = rng.uniform(z[0], z[-1], 1000)
    mine = trilinear_sample(_SHIPPED, qx, qy, qz)
    theirs = lg.trilinear_sample(str(_SHIPPED), qx, qy, qz)
    assert set(mine) == set(theirs)
    for f in theirs:
        assert np.array_equal(mine[f], theirs[f]), f"{f} deviates from the legacy sampler"
