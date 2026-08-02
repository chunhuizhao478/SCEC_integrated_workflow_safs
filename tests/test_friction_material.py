"""Phase 3 + Phase 4 acceptance tests: friction and material."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from deckbuild.asagi import read_asagi  # noqa: E402
from deckbuild.config import Project, StrikeFrame  # noqa: E402
from deckbuild.friction import (  # noqa: E402
    AB_HOT, FRICTION_FIELDS, FrictionError, FrictionStage, FwDesign, M_SLOPE, VW_HI,
    VW_LO, a_minus_b, a_of_T, check_seissol_supports_spatial_muw, depth_profile_fields,
    fw_profile_1d, nucleation_lua, unit_self_test, vw_of_T, write_lua_map, _eval_lua_text,
)
from deckbuild.material import (  # noqa: E402
    AttenuationSpec, MaterialError, MaterialStage, PlasticitySpec, derive_plasticity,
    sv_profile_from_material,
)

PROJECTS = ROOT / "projects"


@pytest.fixture()
def demo(tmp_path):
    return Project.load(PROJECTS / "demo_planar.yaml", require_files=False,
                        data_dir=tmp_path / "data")


# ============================================================ Phase 3: friction profiles
def test_profile_self_test_passes():
    unit_self_test()


def test_a_minus_b_case1_kinks():
    """0 at 100 degC and again at 350; the VW plateau is -0.004."""
    assert a_minus_b(15.0, 1) == pytest.approx(0.004)
    assert a_minus_b(100.0, 1) == pytest.approx(0.0, abs=1e-15)
    assert a_minus_b(200.0, 1) == pytest.approx(-0.004)
    assert a_minus_b(350.0, 1) == pytest.approx(0.0, abs=1e-15)
    assert a_minus_b(400.0, 1) == pytest.approx(0.004)


def test_a_minus_b_case2_is_vw_to_the_surface():
    """CASE2's trace is fully velocity-WEAKENING -- why CASE1 calibration cannot transfer."""
    assert a_minus_b(0.0, 2) == pytest.approx(-0.004)
    assert a_minus_b(15.0, 2) == pytest.approx(-0.004)
    assert a_minus_b(1.0, 1) > 0 > a_minus_b(1.0, 2)


def test_a_minus_b_slope_is_the_legacy_one():
    assert (a_minus_b(120.0, 1) - a_minus_b(110.0, 1)) / 10.0 == pytest.approx(-M_SLOPE)


def test_unknown_case_raises():
    with pytest.raises(FrictionError, match="unknown thermal case"):
        a_minus_b(100.0, 3)


def test_vw_profile():
    assert vw_of_T(0.0) == pytest.approx(VW_LO)
    assert vw_of_T(350.0) == pytest.approx(VW_LO)
    assert vw_of_T(375.0) == pytest.approx(0.5 * (VW_LO + VW_HI))
    assert vw_of_T(400.0) == pytest.approx(VW_HI)
    assert vw_of_T(9999.0) == pytest.approx(VW_HI)


def test_g4_anchor():
    assert a_of_T(279.0, 1, 0.019) == pytest.approx(0.015)
    assert vw_of_T(279.0) == pytest.approx(0.05)


def test_depth_profile_requires_decreasing_knots():
    with pytest.raises(FrictionError, match="DECREASING"):
        depth_profile_fields(np.array([-100.0]), [0.0, -1000.0][::-1], [0.0, 1.0],
                             [0.0, -1000.0], [0.05, 0.05], 0.019)


def test_depth_profile_interpolates():
    z = np.array([0.0, -1500.0, -3000.0])
    a, vw = depth_profile_fields(z, [0.0, -3000.0], [0.004, -0.004],
                                 [0.0, -3000.0], [0.05, 1000.0], 0.019)
    assert a[0] == pytest.approx(0.023) and a[-1] == pytest.approx(0.015)
    assert a[1] == pytest.approx(0.019)
    assert vw[1] == pytest.approx(0.5 * (0.05 + 1000.0))


# ============================================================ Phase 3: graded f_w
def test_fwdesign_f0_admissibility():
    with pytest.raises(FrictionError, match="F0"):
        FwDesign(fw_values=(0.7,), f0=0.6)
    with pytest.raises(FrictionError, match="F0"):
        FwDesign(fw_values=(-0.1,))
    with pytest.raises(FrictionError, match="transition window"):
        FwDesign(fw_values=(0.0, 0.05))
    with pytest.raises(FrictionError, match="s_hi"):
        FwDesign(fw_values=(0.0, 0.05), boundaries_s_km=((30.0, 10.0),))


def test_fw_profile_plateaus():
    d = FwDesign(fw_values=(0.0, 0.045), boundaries_s_km=((150.0, 160.0),))
    assert fw_profile_1d(np.array([0.0]), d)[0] == 0.0
    assert fw_profile_1d(np.array([300.0]), d)[0] == pytest.approx(0.045)
    assert fw_profile_1d(np.array([155.0]), d)[0] == pytest.approx(0.0225)


def test_lua_map_roundtrips_exactly(tmp_path):
    """FL: the EMITTED text must evaluate to the design, not just the design object."""
    strike = StrikeFrame(azimuth_deg=314.0, origin_xy=(606971.0, 3707270.0))
    d = FwDesign(fw_values=(0.0, 0.045, 0.03), boundaries_s_km=((150.0, 160.0),
                                                                (176.0, 182.0)))
    p = write_lua_map(tmp_path / "m.yaml", d, strike)
    s = np.linspace(-50.0, 400.0, 901)
    assert np.max(np.abs(_eval_lua_text(p.read_text(), s) - fw_profile_1d(s, d))) < 1e-12


def test_lua_map_embeds_the_descriptor_strike_frame(tmp_path):
    """A wrong frame here silently mirrors or shifts every along-strike design."""
    strike = StrikeFrame(azimuth_deg=314.0, origin_xy=(606971.0, 3707270.0))
    p = write_lua_map(tmp_path / "m.yaml",
                      FwDesign(fw_values=(0.0, 0.05), boundaries_s_km=((10.0, 20.0),)),
                      strike)
    t = p.read_text()
    assert "606971.0" in t and "3707270.0" in t
    # the coefficients that appear in every shipped rs_muw map
    assert "-0.71933980033865119" in t and "0.69465837045899725" in t


def test_fw_stage_gates(demo, tmp_path):
    st = FrictionStage()
    d = FwDesign(fw_values=(0.0, 0.05), boundaries_s_km=((10.0, 20.0),))
    art = st.build_fw_map(demo, tmp_path, d)
    rep = st.verify_fw_map(demo, art, d)
    assert rep.ok, "\n".join(rep.lines())
    assert {g.name for g in rep.gates} >= {"F0", "F1", "FL", "FLframe"}


def test_fw_verify_catches_a_mismatched_design(demo, tmp_path):
    st = FrictionStage()
    built = FwDesign(fw_values=(0.0, 0.05), boundaries_s_km=((10.0, 20.0),))
    art = st.build_fw_map(demo, tmp_path, built)
    claimed = FwDesign(fw_values=(0.0, 0.09), boundaries_s_km=((10.0, 20.0),))
    assert not st.verify_fw_map(demo, art, claimed).ok


def test_seissol_version_gate():
    """v1.1.3 accepts a spatial rs_muw and SILENTLY IGNORES it."""
    assert check_seissol_supports_spatial_muw(None).skipped
    import tempfile
    for ver, want in (("1.1.3", False), ("1.3.2", False), ("1.4.0", True)):
        d = Path(tempfile.mkdtemp())
        (d / "CMakeLists.txt").write_text(f"project(SeisSol VERSION {ver})")
        rep = check_seissol_supports_spatial_muw(d)
        assert rep.ok is want, (ver, rep.gates[0].detail)


def test_nucleation_lua_uses_the_snapped_point():
    class S:
        xyz = (604446.944, 3704576.3853, -10067.9819)
    t = nucleation_lua(S(), radius_m=2000.0, amplitude_mpa=75.0)
    assert "604446.944" in t and "-10067.9819" in t
    assert "math.exp(r2 / (r2 - R2))" in t and "75.0" in t


# ============================================================ Phase 3: the nc
def test_friction_build_from_depth_profile(demo, tmp_path):
    art = FrictionStage().build(demo, tmp_path)
    x, y, z, f, attrs = read_asagi(art.path)
    assert sorted(f) == sorted(FRICTION_FIELDS)
    assert attrs["source"] == "depth_profile"
    assert np.all(f["rs_a"] > 0)
    assert f["rs_srW"].min() >= VW_LO * (1 - 1e-6)


def test_friction_refuses_without_a_source(demo, tmp_path, monkeypatch):
    from deckbuild.config import RawSources
    object.__setattr__(demo, "raw", RawSources(orientation=demo.raw.orientation,
                                               velocity=demo.raw.velocity, thermal=None))
    with pytest.raises(FrictionError, match="no thermal artifact"):
        FrictionStage().build(demo, tmp_path)


def test_friction_gates(demo, tmp_path):
    st = FrictionStage()
    art = st.build(demo, tmp_path)
    rep = st.verify(demo, art)
    assert rep.ok, "\n".join(rep.lines())
    assert {g.name for g in rep.gates} >= {"G2a", "G2b", "G4"}


def test_friction_g2_rejects_nonpositive_a(demo, tmp_path):
    """rs_a <= 0 is unphysical and must never reach a deck."""
    import yaml
    raw = yaml.safe_load((PROJECTS / "demo_planar.yaml").read_text())
    raw["raw"]["thermal"]["params"]["a_minus_b_values"] = [-0.5, -0.5, -0.5, -0.5]
    p = tmp_path / "bad.yaml"
    p.write_text(yaml.safe_dump(raw))
    cfg = Project.load(p, require_files=False, data_dir=tmp_path)
    with pytest.raises(FrictionError, match="G2"):
        FrictionStage().build(cfg, tmp_path)


def test_friction_from_a_thermal_nc(demo, tmp_path):
    """The CTM path: build a T grid, then friction from it."""
    from deckbuild.asagi import write_asagi
    from deckbuild.geometry import build_grid
    gx, gy, gz = build_grid(demo.stress_box)
    T = np.broadcast_to((-gz / 1000.0 * 25.0)[:, None, None],
                        (len(gz), len(gy), len(gx))).copy()
    tnc = write_asagi(tmp_path / "T.nc", gx, gy, gz, {"T": T})

    class A:
        path = tnc
        sha256 = "deadbeef"
    art = FrictionStage().build(demo, tmp_path, case=1, thermal=A(), out_name="fc1.nc")
    _, _, z, f, attrs = read_asagi(art.path)
    assert attrs["source"] == "ctm_slices/case1"
    # At 20 km depth T = 500 degC -> V_w is on the hot plateau.
    deep = int(np.argmin(np.abs(z + 20000.0)))
    assert f["rs_srW"][deep].max() == pytest.approx(VW_HI, rel=1e-5)


# ============================================================ Phase 4: material
def test_material_build_and_gates(demo, tmp_path):
    st = MaterialStage()
    out = st.build(demo, tmp_path, plasticity=PlasticitySpec(),
                   attenuation=AttenuationSpec())
    assert out.material and out.plasticity and out.sv_profile
    rep = st.verify(demo, out, plasticity_artifact=out.plasticity)
    assert rep.ok, "\n".join(rep.lines())
    assert out.attenuation["freq_central"] == 0.5


def test_material_moduli_are_consistent(demo, tmp_path):
    out = MaterialStage().build(demo, tmp_path)
    _, _, z, f, _ = read_asagi(out.material.path)
    vs = np.sqrt(f["mu"] / f["rho"])
    assert 1500.0 < vs.min() and vs.max() < 4200.0
    assert np.all(f["lambda"] > 0)


def test_material_rejects_a_bad_velocity_model(tmp_path):
    """Vp <= sqrt(2) Vs gives lambda <= 0 -- a bad model, and the message must say so."""
    import yaml
    raw = yaml.safe_load((PROJECTS / "demo_planar.yaml").read_text())
    raw["raw"]["velocity"]["params"]["vp"] = [2000.0, 2000.0, 2000.0, 2000.0]
    raw["raw"]["velocity"]["params"]["vs"] = [1800.0, 1800.0, 1800.0, 1800.0]
    p = tmp_path / "badv.yaml"
    p.write_text(yaml.safe_dump(raw))
    cfg = Project.load(p, require_files=False, data_dir=tmp_path)
    with pytest.raises(MaterialError, match="M1"):
        MaterialStage().build(cfg, tmp_path)


def test_plasticity_is_tan_phi_not_degrees():
    """SeisSol computes atan(bulkFriction), so we must store the COEFFICIENT."""
    rho = np.array([2500.0, 2500.0])
    mu = np.array([2500.0 * 2000.0 ** 2, 2500.0 * 3000.0 ** 2])   # Vs 2000 / 3000
    pc, bf, vs = derive_plasticity(rho, mu, PlasticitySpec(30.0, 40.0, 2500.0, 1e-4))
    assert vs == pytest.approx([2000.0, 3000.0])
    assert bf[0] == pytest.approx(np.tan(np.radians(30.0)))
    assert bf[1] == pytest.approx(np.tan(np.radians(40.0)))
    assert bf[0] == pytest.approx(0.5773502, abs=1e-6)
    assert pc == pytest.approx(1e-4 * mu)


def test_plasticity_default_is_the_paper_not_the_deck():
    """35/45 is Roten's; the SAFS decks use 30/40.  Neither may be silent."""
    s = PlasticitySpec()
    assert (s.phi_soft_deg, s.phi_hard_deg) == (35.0, 45.0)


def test_cohesion_is_linear_in_mu():
    """Linearity is WHY ASAGI-interpolated plastCo stays consistent with mu."""
    rho = np.full(3, 2500.0)
    mu = np.array([1e9, 2e9, 3e9])
    pc, _, _ = derive_plasticity(rho, mu, PlasticitySpec())
    assert pc[1] == pytest.approx(0.5 * (pc[0] + pc[2]))


def test_sv_profile_is_monotone_and_physical(demo, tmp_path):
    out = MaterialStage().build(demo, tmp_path)
    d = np.load(out.sv_profile.path)
    depth, sv = d["depth_m"], d["sv_eff_mpa"]
    assert depth[0] == 0.0 and sv[0] == pytest.approx(0.0)
    assert np.all(np.diff(sv) > 0)
    # ~(rho - rho_w) g z: at 10 km with rho ~ 2700 that is ~165 MPa.
    i = int(np.argmin(np.abs(depth - 10000.0)))
    assert 120.0 < sv[i] < 200.0


def test_sv_profile_cache_key_tracks_the_material(demo, tmp_path):
    """A rebuilt material must invalidate the profile -- the name carries its hash."""
    a = MaterialStage().build(demo, tmp_path / "a")
    import yaml
    raw = yaml.safe_load((PROJECTS / "demo_planar.yaml").read_text())
    raw["raw"]["velocity"]["params"]["grid_dx"] = 4000.0
    p = tmp_path / "coarse.yaml"
    p.write_text(yaml.safe_dump(raw))
    cfg2 = Project.load(p, require_files=False, data_dir=tmp_path)
    b = MaterialStage().build(cfg2, tmp_path / "b")
    assert a.material.sha256 != b.material.sha256
    assert Path(a.sv_profile.path).name != Path(b.sv_profile.path).name


def test_cvm_reader_with_no_matching_slices_says_so(demo, tmp_path):
    """The cvm_slices reader IS implemented (Phase 8 R-801); an empty glob must be loud."""
    import yaml
    raw = yaml.safe_load((PROJECTS / "demo_planar.yaml").read_text())
    raw["raw"]["velocity"]["kind"] = "cvm_slices"
    raw["raw"]["velocity"]["path"] = "raw/cvm/CVM_*.csv"
    p = tmp_path / "cvm.yaml"
    p.write_text(yaml.safe_dump(raw))
    cfg = Project.load(p, require_files=False, data_dir=tmp_path)
    with pytest.raises(MaterialError, match="no CVM slices matched"):
        MaterialStage().build(cfg, tmp_path)


# ------------------------------------------------- review findings R-301..R-304
import yaml  # noqa: E402


def _variant(tmp_path, **mut):
    raw = yaml.safe_load((PROJECTS / "demo_planar.yaml").read_text())
    for path, val in mut.items():
        node = raw
        keys = path.split(".")
        for kk in keys[:-1]:
            node = node[kk]
        node[keys[-1]] = val
    p = tmp_path / "v.yaml"
    p.write_text(yaml.safe_dump(raw))
    return Project.load(p, require_files=False, data_dir=tmp_path / "data")


def test_R301_sv_profile_is_sea_level_referenced(tmp_path):
    """depth 0 must mean z = 0, whatever the material grid's z_max is.

    Grid-top referencing offset the WHOLE stress field by z_max (+3250 m for the shipped
    SAFS CVM).  No gate caught it: V2 checks the eigenvalue RATIO k, which is invariant
    under an overall scale of Sv.
    """
    cfg = _variant(tmp_path, **{"raw.velocity.params.z_max": 500.0})
    out = MaterialStage().build(cfg, tmp_path)
    d = np.load(out.sv_profile.path)
    depth, sv = d["depth_m"], d["sv_eff_mpa"]
    assert str(d["reference"]) == "sea_level"

    i = int(np.argmin(np.abs(depth - 0.0)))
    assert depth[i] == pytest.approx(0.0)
    # At sea level there is 500 m of rock ABOVE, so Sv_eff is already positive there.
    assert sv[i] > 0.0

    j = int(np.argmin(np.abs(depth - 10000.0)))
    # ~ (rho - rho_w) g z with rho ~ 2700 plus the 500 m of extra overburden.
    assert 150.0 < sv[j] < 185.0


def test_R301_sea_level_reference_is_independent_of_the_grid_top(tmp_path):
    """The SAME physical depth must give the same Sv whatever z_max is."""
    a = MaterialStage().build(_variant(tmp_path, **{"raw.velocity.params.z_max": 0.0}),
                              tmp_path / "a")
    b = MaterialStage().build(_variant(tmp_path, **{"raw.velocity.params.z_max": 500.0}),
                              tmp_path / "b")
    da, db = np.load(a.sv_profile.path), np.load(b.sv_profile.path)
    for probe in (2000.0, 10000.0):
        va = np.interp(probe, da["depth_m"], da["sv_eff_mpa"])
        vb = np.interp(probe, db["depth_m"], db["sv_eff_mpa"])
        # b carries 500 m more overburden, so it is HIGHER -- but by a bounded, physical
        # amount, not by a wholesale shift of the depth axis.
        assert vb > va
        assert (vb - va) < 15.0, f"at {probe} m: {va:.2f} vs {vb:.2f} MPa"


def test_R302_friction_grid_is_settable_independently(tmp_path):
    """The four grids are independent; friction must not silently inherit the stress one."""
    cfg = _variant(tmp_path, **{"raw.thermal.params.dz": 200.0,
                                "raw.thermal.params.grid_dx": 5000.0})
    art = FrictionStage().build(cfg, tmp_path)
    _, _, z, _, _ = read_asagi(art.path)
    assert float(np.diff(z)[0]) == pytest.approx(200.0)
    assert float(np.diff(z)[0]) != cfg.stress_box.dz


def test_R303_FL_parser_fails_loudly_if_it_drifts_from_the_emitter(tmp_path):
    """A regex that silently matches nothing would make FL check nothing."""
    strike = StrikeFrame(azimuth_deg=0.0, origin_xy=(0.0, 0.0))
    d = FwDesign(fw_values=(0.0, 0.05), boundaries_s_km=((10.0, 20.0),))
    p = write_lua_map(tmp_path / "m.yaml", d, strike)
    text = p.read_text().replace("inc = inc +", "inc = inc  +")   # break the parser
    with pytest.raises(FrictionError, match="drifted"):
        _eval_lua_text(text, np.array([0.0, 15.0, 30.0]))


def test_R304_M4_runs_for_the_layered_reader(demo, tmp_path):
    """M4 is the material's round-trip guard; it must actually run, not always skip."""
    st = MaterialStage()
    out = st.build(demo, tmp_path)
    rep = st.verify(demo, out)
    m4 = [g for g in rep.gates if g.name == "M4"]
    assert m4, "M4 missing entirely"
    assert m4[0].severity != "skip", "M4 must run for layered_1d"
    assert m4[0].passed, m4[0].detail


def test_eval_lua_reads_the_LEGACY_form_with_a_global_scale(tmp_path):
    """The shipped maps use `return { rs_muw = X + inc * g }` with a `local g`.

    Our emitter writes `X + inc`.  Phase 8 compares against the shipped files, so the
    parser has to read BOTH -- it previously raised "cannot find the base f_w" on every
    real deck.
    """
    strike = StrikeFrame(azimuth_deg=314.0, origin_xy=(606971.0, 3707270.0))
    d = FwDesign(fw_values=(0.0, 0.05), boundaries_s_km=((10.0, 20.0),))
    ours = write_lua_map(tmp_path / "m.yaml", d, strike).read_text()
    legacy = ours.replace("+ inc }", "+ inc * g }").replace(
        "      local inc = 0.0", "      local g = 1.0\n      local inc = 0.0")
    s = np.linspace(-20.0, 60.0, 201)
    assert np.allclose(_eval_lua_text(legacy, s), _eval_lua_text(ours, s), atol=0)


def test_eval_lua_honours_the_global_scale(tmp_path):
    strike = StrikeFrame(azimuth_deg=0.0, origin_xy=(0.0, 0.0))
    d = FwDesign(fw_values=(0.0, 0.06), boundaries_s_km=((10.0, 20.0),))
    ours = write_lua_map(tmp_path / "m.yaml", d, strike).read_text()
    halved = ours.replace("+ inc }", "+ inc * g }").replace(
        "      local inc = 0.0", "      local g = 0.5\n      local inc = 0.0")
    s = np.array([100.0])                       # well past the transition
    assert _eval_lua_text(halved, s)[0] == pytest.approx(0.03)


def test_eval_lua_unparseable_names_the_offending_line(tmp_path):
    with pytest.raises(FrictionError, match="return \\{"):
        _eval_lua_text("function f(x)\n  return { rs_muw = wat }\nend\n",
                       np.array([0.0]))
