"""Phase 2 acceptance tests for deckbuild.stress and deckbuild.orientation."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from deckbuild.asagi import read_asagi  # noqa: E402
from deckbuild.config import Project, SourceSpec, CRS  # noqa: E402
from deckbuild.geometry import strike_s_km  # noqa: E402
from deckbuild.orientation import (  # noqa: E402
    OrientationError, csm_tensors_tension, read_orientation,
)
from deckbuild.stress import (  # noqa: E402
    KDesign, STRESS_FIELDS, StressError, StressStage, build_tensor_andersonian,
    k_profile_1d, magnitudes_C1, sv_eff_at,
)

PROJECTS = ROOT / "projects"


@pytest.fixture()
def demo(tmp_path):
    """The demo project, whose orientation is a constant and needs no downloads."""
    return Project.load(PROJECTS / "demo_planar.yaml", require_files=False,
                        data_dir=tmp_path / "data")


@pytest.fixture()
def sv():
    """A simple lithostat minus hydrostat: (depth_m, Sv_eff MPa)."""
    d = np.linspace(0.0, 30000.0, 301)
    return d, (2700.0 - 1000.0) * 9.81 * d / 1.0e6


# --------------------------------------------------------------------------- closure
def test_magnitudes_C1_ordering_and_ratio():
    sv = np.array([100.0])
    for k in (1.0, 1.4, 1.7, 2.5):
        for R in (0.0, 0.3, 0.5, 1.0):
            s1, s2, s3 = magnitudes_C1(sv, np.array([R]), k)
            assert s1[0] >= s2[0] >= s3[0] > 0, (k, R)
            assert s1[0] / s3[0] == pytest.approx(k, rel=1e-12)


def test_magnitudes_C1_k_may_be_a_per_column_array():
    """k is a COLUMN property; an array must broadcast exactly like a scalar."""
    sv = np.full(4, 50.0)
    R = np.full(4, 0.4)
    ka = np.array([1.2, 1.5, 1.8, 2.1])
    s1a, s2a, s3a = magnitudes_C1(sv, R, ka)
    for i, kv in enumerate(ka):
        s1, s2, s3 = magnitudes_C1(sv[i:i + 1], R[i:i + 1], kv)
        assert (s1a[i], s2a[i], s3a[i]) == pytest.approx((s1[0], s2[0], s3[0]), rel=0)


def test_k_equal_one_is_isotropic_horizontal():
    s1, s2, s3 = magnitudes_C1(np.array([80.0]), np.array([0.5]), 1.0)
    assert s1[0] == pytest.approx(s3[0])


def test_build_tensor_is_symmetric_with_vertical_sigma2():
    az = np.array([30.0, 314.0])
    s1, s2, s3 = np.array([3.0, 3.0]), np.array([2.0, 2.0]), np.array([1.0, 1.0])
    t = build_tensor_andersonian(az, s1, s2, s3)
    assert np.allclose(t, np.transpose(t, (0, 2, 1)))
    # sigma2 is exactly vertical: the (2,2) entry is s2 and the off-diagonals vanish.
    assert np.allclose(t[:, 2, 2], s2)
    assert np.allclose(t[:, 0, 2], 0.0) and np.allclose(t[:, 1, 2], 0.0)
    w = np.linalg.eigvalsh(t)
    assert np.allclose(np.sort(w, axis=1), np.array([[1.0, 2.0, 3.0]] * 2))


def test_build_tensor_puts_sigma1_along_shmax():
    az = 30.0
    t = build_tensor_andersonian(np.array([az]), np.array([3.0]), np.array([2.0]),
                                 np.array([1.0]))
    eH = np.array([np.sin(np.radians(az)), np.cos(np.radians(az)), 0.0])
    assert (eH @ t[0] @ eH) == pytest.approx(3.0)


def test_sv_eff_at_clamps():
    d = np.array([0.0, 1000.0, 2000.0])
    s = np.array([0.0, 17.0, 34.0])
    assert sv_eff_at(np.array([-500.0]), d, s)[0] == 0.0
    assert sv_eff_at(np.array([9999.0]), d, s)[0] == 34.0
    assert sv_eff_at(np.array([500.0]), d, s)[0] == pytest.approx(8.5)


# --------------------------------------------------------------------------- k design
def test_kdesign_rejects_bad_input():
    with pytest.raises(StressError, match="at least one"):
        KDesign(k_values=())
    with pytest.raises(StressError, match=">= 1"):
        KDesign(k_values=(0.9,))
    with pytest.raises(StressError, match="transition window"):
        KDesign(k_values=(1.4, 1.8))                       # missing the window
    with pytest.raises(StressError, match="s_hi"):
        KDesign(k_values=(1.4, 1.8), boundaries_s_km=((50.0, 20.0),))
    with pytest.raises(StressError, match="overlaps"):
        KDesign(k_values=(1.4, 1.8, 2.0),
                boundaries_s_km=((20.0, 60.0), (40.0, 80.0)))


def test_k_profile_plateaus_and_smoothsteps():
    d = KDesign(k_values=(1.4, 1.8), boundaries_s_km=((100.0, 140.0),))
    assert k_profile_1d(np.array([0.0]), d)[0] == pytest.approx(1.4)
    assert k_profile_1d(np.array([300.0]), d)[0] == pytest.approx(1.8)
    assert k_profile_1d(np.array([120.0]), d)[0] == pytest.approx(1.6)   # smoothstep(0.5)
    s = np.linspace(0.0, 300.0, 601)
    k = k_profile_1d(s, d)
    assert np.all(np.diff(k) >= -1e-12)                     # monotone
    assert k.min() == pytest.approx(1.4) and k.max() == pytest.approx(1.8)


def test_k_profile_three_regions():
    d = KDesign(k_values=(1.4, 1.9, 1.6),
                boundaries_s_km=((50.0, 70.0), (150.0, 170.0)))
    assert k_profile_1d(np.array([0.0, 100.0, 300.0]), d) == pytest.approx(
        [1.4, 1.9, 1.6])


def test_uniform_design_is_flat():
    d = KDesign(k_values=(1.7,))
    assert np.allclose(k_profile_1d(np.linspace(-500, 500, 50), d), 1.7)


# --------------------------------------------------------------------------- orientation
def test_constant_orientation(demo):
    f = read_orientation(demo.raw.orientation, demo.crs, demo.data_dir)
    assert f.is_uniform and f.kind == "constant"
    az, R, n = f.at([1.0, 2.0], [3.0, 4.0])
    assert az.tolist() == [45.0, 45.0] and R.tolist() == [0.5, 0.5] and n == 0


def test_constant_orientation_validates():
    with pytest.raises(OrientationError, match="azimuth_deg"):
        read_orientation(SourceSpec(kind="constant", params={}), CRS(epsg="EPSG:32610"))
    with pytest.raises(OrientationError, match=r"\[0, 1\]"):
        read_orientation(SourceSpec(kind="constant",
                                    params={"azimuth_deg": 10.0, "shape_ratio": 2.0}),
                         CRS(epsg="EPSG:32610"))


def test_unknown_orientation_kind_raises():
    with pytest.raises(OrientationError, match="unknown orientation kind"):
        read_orientation(SourceSpec(kind="tea_leaves"), CRS(epsg="EPSG:32610"))


def test_csm_tensor_assembly_is_symmetric():
    S = np.arange(6.0)[None, :]
    T = csm_tensors_tension(S)
    assert np.allclose(T, T.transpose(0, 2, 1))
    assert T[0, 0, 0] == 0.0 and T[0, 2, 2] == 5.0


def test_azimuth_interpolation_wraps_at_180():
    """SHmax is an AXIS: averaging 179 and 1 must give 0, not 90."""
    from deckbuild.orientation import OrientationField
    f = OrientationField(cx=np.array([0.0, 1000.0, 0.0, 1000.0]),
                         cy=np.array([0.0, 0.0, 1000.0, 1000.0]),
                         az_deg=np.array([179.0, 1.0, 179.0, 1.0]),
                         R=np.full(4, 0.5), kind="test")
    az, R, _ = f.at([500.0], [500.0])
    assert min(az[0], 180.0 - az[0]) < 1.0, f"got {az[0]}, expected ~0/180"


# --------------------------------------------------------------------------- build
def _build(demo, sv, tmp_path, **kw):
    return StressStage().build(demo, tmp_path, sv_profile=sv, **kw)


def test_build_writes_all_six_components(demo, sv, tmp_path):
    art = _build(demo, sv, tmp_path, k=1.7, freeze_above_depth_m=0.0)
    x, y, z, flds, attrs = read_asagi(art.path)
    assert sorted(flds) == sorted(STRESS_FIELDS)
    assert attrs["k_values"] == "1.7"
    assert art.kind == "stress" and len(art.sha256) == 64


def test_build_is_compression_negative(demo, sv, tmp_path):
    """SeisSol/easi wants compression NEGATIVE inside the nc."""
    art = _build(demo, sv, tmp_path, k=1.7, freeze_above_depth_m=0.0)
    _, _, z, flds, _ = read_asagi(art.path)
    deep = np.flatnonzero(z < -5000.0)
    assert np.all(flds["s_zz"][deep] < 0.0)
    assert np.all(flds["s_xx"][deep] < 0.0)


def test_built_field_reproduces_the_closure_ratio(demo, sv, tmp_path):
    """The written tensor's horizontal eigenvalue ratio must be exactly k."""
    art = _build(demo, sv, tmp_path, k=1.7, freeze_above_depth_m=0.0)
    _, _, z, flds, _ = read_asagi(art.path)
    kz = int(np.argmin(np.abs(z - (-8000.0))))
    h = np.zeros(flds["s_xx"][kz].shape + (2, 2))
    h[..., 0, 0] = -flds["s_xx"][kz] / 1e6
    h[..., 1, 1] = -flds["s_yy"][kz] / 1e6
    h[..., 0, 1] = h[..., 1, 0] = -flds["s_xy"][kz] / 1e6
    w = np.linalg.eigvalsh(h)
    assert np.allclose(w[..., 1] / w[..., 0], 1.7, rtol=2e-3)


def test_build_requires_exactly_one_of_k_or_design(demo, sv, tmp_path):
    with pytest.raises(StressError, match="exactly one"):
        _build(demo, sv, tmp_path)
    with pytest.raises(StressError, match="exactly one"):
        _build(demo, sv, tmp_path, k=1.5, design=KDesign(k_values=(1.5,)))


def test_build_rejects_k_below_one(demo, sv, tmp_path):
    with pytest.raises(StressError, match=">= 1"):
        _build(demo, sv, tmp_path, k=0.5)


def test_sv_profile_is_required(demo, tmp_path):
    with pytest.raises(TypeError):
        StressStage().build(demo, tmp_path, k=1.7)          # no sv_profile


def test_sv_profile_missing_file_raises(demo, tmp_path):
    with pytest.raises(StressError, match="Sv profile not found"):
        _build(demo, tmp_path / "nope.npz", tmp_path, k=1.7)


def test_sv_profile_accepts_an_npz(demo, sv, tmp_path):
    p = tmp_path / "sv.npz"
    np.savez(p, depth_m=sv[0], sv_eff_mpa=sv[1])
    art = _build(demo, p, tmp_path, k=1.7, freeze_above_depth_m=0.0)
    assert Path(art.path).is_file()


def test_sv_profile_npz_without_the_right_keys_raises(demo, tmp_path):
    p = tmp_path / "bad.npz"
    np.savez(p, wrong=np.zeros(3))
    with pytest.raises(StressError, match="expected arrays"):
        _build(demo, p, tmp_path, k=1.7)


def test_default_name_encodes_the_design(demo, sv, tmp_path):
    a = _build(demo, sv, tmp_path, k=1.7, freeze_above_depth_m=0.0)
    b = _build(demo, sv, tmp_path, k=1.7, freeze_above_depth_m=500.0)
    assert Path(a.path).name == "stress_andersonian_k1.7.nc"
    assert Path(b.path).name == "stress_andersonian_k1.7_freeze500m.nc"
    g = _build(demo, sv, tmp_path,
               design=KDesign(k_values=(1.4, 1.8), boundaries_s_km=((10.0, 30.0),)),
               freeze_above_depth_m=0.0)
    assert "gradedk_1.4-1.8" in Path(g.path).name


# --------------------------------------------------------------------------- freeze
def test_freeze_makes_shallow_levels_identical(demo, sv, tmp_path):
    art = _build(demo, sv, tmp_path, k=1.7, freeze_above_depth_m=500.0)
    _, _, z, flds, _ = read_asagi(art.path)
    above = np.flatnonzero(z > -500.0)
    kref = int(np.argmin(np.abs(z - (-500.0))))
    assert above.size > 0
    for f in STRESS_FIELDS:
        for i in above:
            assert np.array_equal(flds[f][i], flds[f][kref])


def test_below_the_freeze_depth_the_field_is_bit_identical(demo, sv, tmp_path):
    """The freeze must touch ONLY the shallow band -- that is its whole contract."""
    a = _build(demo, sv, tmp_path, k=1.7, freeze_above_depth_m=0.0,
               out_name="unfrozen.nc")
    b = _build(demo, sv, tmp_path, k=1.7, freeze_above_depth_m=500.0,
               out_name="frozen.nc")
    _, _, z, fa, _ = read_asagi(a.path)
    _, _, _, fb, _ = read_asagi(b.path)
    deep = np.flatnonzero(z <= -500.0)
    for f in STRESS_FIELDS:
        assert np.array_equal(fa[f][deep], fb[f][deep]), f


def test_freeze_preserves_mu_app(demo, sv, tmp_path):
    """sigma_n and tau freeze TOGETHER, so the apparent friction is unchanged."""
    art = _build(demo, sv, tmp_path, k=1.7, freeze_above_depth_m=500.0)
    _, _, z, flds, _ = read_asagi(art.path)
    kz = int(np.argmin(np.abs(z - (-250.0))))       # inside the frozen band
    ratio = flds["s_xy"][kz] / flds["s_yy"][kz]
    kdeep = int(np.argmin(np.abs(z - (-4000.0))))
    ratio_deep = flds["s_xy"][kdeep] / flds["s_yy"][kdeep]
    assert np.allclose(ratio, ratio_deep, rtol=1e-5)


# --------------------------------------------------------------------------- gates
def test_verify_v2_passes_on_a_uniform_build(demo, sv, tmp_path):
    art = _build(demo, sv, tmp_path, k=1.7, freeze_above_depth_m=0.0)
    rep = StressStage().verify(demo, art)
    v2 = [g for g in rep.gates if g.name == "V2"][0]
    assert v2.passed, v2.detail


def test_verify_v3_reports_skip_when_the_freeze_is_off(demo, sv, tmp_path):
    art = _build(demo, sv, tmp_path, k=1.7, freeze_above_depth_m=0.0)
    rep = StressStage().verify(demo, art)
    v3 = [g for g in rep.gates if g.name == "V3"][0]
    assert v3.severity == "skip" and not v3.passed      # a skip is never a pass


def test_verify_v3_passes_with_the_freeze_on(demo, sv, tmp_path):
    art = _build(demo, sv, tmp_path, k=1.7, freeze_above_depth_m=500.0)
    rep = StressStage().verify(demo, art)
    v3 = [g for g in rep.gates if g.name == "V3"][0]
    assert v3.passed, v3.detail


def test_verify_v2_detects_a_graded_design(demo, sv, tmp_path):
    d = KDesign(k_values=(1.4, 2.0), boundaries_s_km=((10.0, 30.0),))
    art = _build(demo, sv, tmp_path, design=d, freeze_above_depth_m=0.0)
    rep = StressStage().verify(demo, art)
    v2 = [g for g in rep.gates if g.name == "V2"][0]
    assert v2.passed, v2.detail
    assert "1.400-2.000" in v2.detail or "1.4" in v2.detail


def test_verify_v1_regions_are_bit_identical(demo, sv, tmp_path):
    """V1: grading must not perturb a region's plateau."""
    d = KDesign(k_values=(1.4, 2.0), boundaries_s_km=((10.0, 30.0),))
    art = _build(demo, sv, tmp_path, design=d, freeze_above_depth_m=0.0,
                 out_name="graded.nc")
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    rep = StressStage().verify_regions(demo, art, sv_profile=sv, design=d,
                                       tmp_dir=scratch)
    assert rep.ok, "\n".join(rep.lines())
    assert len([g for g in rep.gates if g.name.startswith("V1r")]) == 2


def test_verify_v1_catches_a_perturbed_region(demo, sv, tmp_path):
    """A graded build whose plateau k does NOT match must fail V1."""
    built = KDesign(k_values=(1.4, 2.0), boundaries_s_km=((10.0, 30.0),))
    art = _build(demo, sv, tmp_path, design=built, freeze_above_depth_m=0.0,
                 out_name="graded.nc")
    claimed = KDesign(k_values=(1.4, 1.5), boundaries_s_km=((10.0, 30.0),))
    scratch = tmp_path / "scratch2"
    scratch.mkdir()
    rep = StressStage().verify_regions(demo, art, sv_profile=sv, design=claimed,
                                       tmp_dir=scratch)
    assert not rep.ok


def test_verify_v4_skips_without_a_mesh(demo, sv, tmp_path):
    art = _build(demo, sv, tmp_path, k=1.7, freeze_above_depth_m=0.0)
    rep = StressStage().verify(demo, art)
    v4 = [g for g in rep.gates if g.name == "V4"]
    assert v4 and v4[0].severity == "skip"


def test_check_outputs_flags_the_v1_scratch_builds(demo, sv, tmp_path):
    art = _build(demo, sv, tmp_path, k=1.7, freeze_above_depth_m=0.0)
    (tmp_path / "_v1_ref_k1.4.nc").write_text("scratch")
    rep = StressStage().check_outputs(tmp_path, keep=[art.path])
    assert rep["ok"] and any("_v1_ref" in e for e in rep["extra"])


# ------------------------------------------------- review findings R-201..R-203
def test_R201_daylight_patch_refuses_the_flat_top_approximation_on_topography(sv, tmp_path):
    """The legacy patch is PER-COLUMN topographic; this build writes a constant slab.

    On a domain with topography those differ, so the build must refuse rather than
    silently write an approximation that Phase 8 could never reproduce.
    """
    import yaml
    raw = yaml.safe_load((PROJECTS / "demo_planar.yaml").read_text())
    raw["stress_box"]["zmax"] = 3000.0                 # give it topography
    raw["meshes"]["planar"]["daylights"] = True
    raw["meshes"]["planar"]["daylight_min_depth_m"] = 0.0
    p = tmp_path / "topo.yaml"
    p.write_text(yaml.safe_dump(raw))
    cfg = Project.load(p, require_files=False, data_dir=tmp_path)
    with pytest.raises(StressError, match="topograph"):
        StressStage().build(cfg, tmp_path, sv_profile=sv, k=1.7)


def test_R201_constant_slab_can_be_chosen_deliberately(sv, tmp_path):
    import yaml
    raw = yaml.safe_load((PROJECTS / "demo_planar.yaml").read_text())
    raw["stress_box"]["zmax"] = 3000.0
    raw["meshes"]["planar"]["daylights"] = True
    raw["meshes"]["planar"]["daylight_min_depth_m"] = 250.0   # explicit opt-in
    p = tmp_path / "topo2.yaml"
    p.write_text(yaml.safe_dump(raw))
    cfg = Project.load(p, require_files=False, data_dir=tmp_path)
    art = StressStage().build(cfg, tmp_path, sv_profile=sv, k=1.7,
                              freeze_above_depth_m=0.0)
    _, _, z, f, _ = read_asagi(art.path)
    above = np.flatnonzero(z >= 0.0)
    assert above.size > 0
    assert np.all(f["s_zz"][above] < 0.0), "daylighting levels must carry real stress"


def test_R202_corrupt_mesh_fails_v4_rather_than_skipping(demo, sv, tmp_path):
    """An ABSENT mesh is a legitimate skip; a PRESENT but unreadable one is a failure."""
    bad = Path(demo.data_dir) / "demo_planar.puml.h5"
    bad.parent.mkdir(parents=True, exist_ok=True)
    bad.write_bytes(b"not an hdf5 file at all")
    art = StressStage().build(demo, tmp_path, sv_profile=sv, k=1.7,
                              freeze_above_depth_m=0.0)
    rep = StressStage().verify(demo, art)
    v4 = [g for g in rep.gates if g.name == "V4"][0]
    assert v4.severity == "hard" and not v4.passed
    assert "could not be read" in v4.detail


def test_R203_graded_k_varies_the_field_analytically(demo, sv, tmp_path):
    """R=0.5 with az=45 makes s_xx k-INVARIANT -- s_xy is the component that moves.

    sigma_xx = 0.5*sig3*(1+k) and sig3 = 2Sv/(k+1), so sigma_xx == Sv for EVERY k.
    A lateral-variation probe on s_xx would pass even if grading were dropped entirely.
    """
    d = KDesign(k_values=(1.2, 2.4), boundaries_s_km=((20.0, 30.0),))
    art = StressStage().build(demo, tmp_path, sv_profile=sv, design=d,
                              freeze_above_depth_m=0.0, out_name="graded_xy.nc")
    _, _, z, f, _ = read_asagi(art.path)
    kz = int(np.argmin(np.abs(z - (-10000.0))))
    Sv = 1700.0 * 9.81 * 10000.0 / 1e6

    # the degeneracy itself, asserted so nobody "fixes" the test back
    assert np.ptp(f["s_xx"][kz]) == pytest.approx(0.0, abs=1.0)

    row = f["s_xy"][kz].mean(axis=1) / 1e6
    assert row[0] == pytest.approx(-Sv * (1.2 - 1) / (1.2 + 1), abs=0.05)
    assert row[-1] == pytest.approx(-Sv * (2.4 - 1) / (2.4 + 1), abs=0.05)
    assert abs(np.ptp(row)) > 40.0
