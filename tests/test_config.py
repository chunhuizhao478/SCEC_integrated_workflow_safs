"""Phase 0 acceptance tests for deckbuild.config."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from deckbuild.config import (  # noqa: E402
    ConfigError, GridBox, Hypocenter, MeshSpec, PhysicsDefaults, Project, StrikeFrame,
    as_dict,
)

PROJECTS = ROOT / "projects"


# --------------------------------------------------------------------------- helpers
def load_raw(stem: str) -> dict:
    with open(PROJECTS / f"{stem}.yaml") as fh:
        return yaml.safe_load(fh)


def write_tmp(tmp_path: Path, raw: dict, name: str = "p.yaml") -> Path:
    p = tmp_path / name
    with open(p, "w") as fh:
        yaml.safe_dump(raw, fh)
    return p


@pytest.fixture()
def safs_raw() -> dict:
    return load_raw("safs_alt")


# --------------------------------------------------------------------------- happy path
@pytest.mark.parametrize("stem", ["safs_alt", "safs_preferred", "demo_planar"])
def test_every_shipped_descriptor_loads(stem):
    """Schema validity of all shipped descriptors (files are not required to exist)."""
    cfg = Project.load(PROJECTS / f"{stem}.yaml", require_files=False)
    assert cfg.name == stem
    assert cfg.default_mesh in cfg.meshes


def test_safs_alt_values():
    cfg = Project.load(PROJECTS / "safs_alt.yaml", require_files=False)
    assert cfg.crs.epsg == "EPSG:32611"
    assert cfg.strike.azimuth_deg == 314.0
    assert cfg.strike.origin_xy == (606971.0, 3707270.0)
    assert cfg.strike.hint_deg == 314.0          # defaults to azimuth
    box = cfg.stress_box
    assert (box.xmin, box.xmax) == (350000.0, 630000.0)
    assert (box.ymin, box.ymax) == (3680000.0, 3850000.0)
    assert (box.zmin, box.zmax) == (-17000.0, 0.0)
    assert (box.dx, box.dz) == (1000.0, 250.0)
    assert cfg.meshes["alt"].tag_to_bc == {101: 3, 102: 1, 103: 5, 104: 5}
    assert cfg.meshes["alt"].daylights is False
    assert [b.name for b in cfg.gate_bands] == ["San Gorgonio"]
    assert (cfg.gate_bands[0].s_start_km, cfg.gate_bands[0].s_end_km) == (28.0, 72.0)
    ph = cfg.physics
    assert ph.mu_shear_pa == 23.5e9 and ph.w_energy_m == 10.0e3
    assert ph.dc_m == 1.2 and ph.kappa_c == 0.9
    assert ph.seis_band_km == (3.0, 12.0) and ph.active_band_km == (0.3, 15.0)
    assert ph.l_coast_km == 6.0 and ph.gap_max_km == 4.0 and ph.sn_floor_mpa == 30.0
    assert ph.rs_b == 0.019 and ph.rs_sl0 == 0.10


def test_safs_preferred_daylights_and_has_two_bands():
    cfg = Project.load(PROJECTS / "safs_preferred.yaml", require_files=False)
    assert cfg.meshes["preferred"].daylights is True
    assert cfg.meshes["preferred"].daylight_min_depth_m == 100.0
    assert cfg.stress_box.zmax == 3000.0          # the daylight z-extension
    assert [b.name for b in cfg.gate_bands] == ["San Gorgonio", "NW Big Bend"]


def test_demo_is_not_safs():
    """The demo must share no SAFS constant, or a leak would go unnoticed."""
    cfg = Project.load(PROJECTS / "demo_planar.yaml", require_files=False)
    assert cfg.crs.epsg != "EPSG:32611"
    assert cfg.strike.azimuth_deg != 314.0
    assert cfg.strike.origin_xy != (606971.0, 3707270.0)
    assert cfg.raw.orientation.kind == "constant"
    assert cfg.raw.velocity.kind == "layered_1d"
    assert cfg.raw.thermal.kind == "depth_profile"


def test_physics_defaults_apply_when_block_is_empty():
    cfg = Project.load(PROJECTS / "demo_planar.yaml", require_files=False)
    assert cfg.physics == PhysicsDefaults()


def test_accessors():
    cfg = Project.load(PROJECTS / "safs_alt.yaml", require_files=False)
    assert cfg.mesh() is cfg.meshes["alt"]
    assert cfg.mesh("alt") is cfg.meshes["alt"]
    assert cfg.hypocenter().label.startswith("v4_0_0")
    with pytest.raises(ConfigError, match="unknown mesh"):
        cfg.mesh("nope")
    with pytest.raises(ConfigError, match="no hypocenter"):
        cfg.hypocenter("nope")


def test_summary_and_as_dict_roundtrip():
    cfg = Project.load(PROJECTS / "safs_alt.yaml", require_files=False)
    s = cfg.summary()
    assert "safs_alt" in s and "EPSG:32611" in s
    d = as_dict(cfg)
    assert d["strike"]["azimuth_deg"] == 314.0
    assert "source_path" not in d and "data_dir" not in d


def test_sha256_is_stable_and_nonempty():
    cfg = Project.load(PROJECTS / "safs_alt.yaml", require_files=False)
    assert len(cfg.sha256()) == 64
    assert cfg.sha256() == cfg.sha256()


# ------------------------------------------------- acceptance: three rejection tests
def test_unknown_top_level_key_raises_with_the_key_named(tmp_path, safs_raw):
    """A YAML with an unknown top-level key must RAISE, not ignore.

    Silent key-drop is how a retargeted project keeps a SAFS default -- e.g. a misspelt
    'strike_frame:' would leave the real 'strike:' untouched.
    """
    safs_raw["strike_frame"] = {"azimuth_deg": 1.0}
    with pytest.raises(ConfigError) as exc:
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)
    assert "strike_frame" in str(exc.value)


def test_bad_reader_kind_raises_with_the_kind_named(tmp_path, safs_raw):
    safs_raw["raw"]["velocity"]["kind"] = "cvm_slice"   # typo: missing the s
    with pytest.raises(ConfigError) as exc:
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)
    msg = str(exc.value)
    assert "cvm_slice" in msg and "raw.velocity" in msg


def test_reversed_gate_band_raises_with_the_band_named(tmp_path, safs_raw):
    safs_raw["gate_bands"] = [{"s_start_km": 72.0, "s_end_km": 28.0, "name": "San Gorgonio"}]
    with pytest.raises(ConfigError) as exc:
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)
    msg = str(exc.value)
    assert "San Gorgonio" in msg and "gate_bands[0]" in msg


# ------------------------------------------------- other validation rules
def test_unknown_nested_key_raises(tmp_path, safs_raw):
    safs_raw["crs"]["epsq"] = "EPSG:1"
    with pytest.raises(ConfigError, match="epsq"):
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)


def test_default_mesh_not_in_meshes_raises(tmp_path, safs_raw):
    safs_raw["default_mesh"] = "preferred"
    with pytest.raises(ConfigError, match="default_mesh"):
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)


def test_hypocenter_without_matching_mesh_raises(tmp_path, safs_raw):
    safs_raw["hypocenters"]["ghost"] = {"x": 1.0, "y": 2.0, "z": -3.0}
    with pytest.raises(ConfigError, match="ghost"):
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)


def test_hypocenter_must_be_projected_xor_geographic(tmp_path, safs_raw):
    # BOTH complete -> caught by the mixing check (R-003), which fires first and names
    # exactly which keys clash.
    safs_raw["hypocenters"]["alt"] = {"x": 1.0, "y": 2.0, "z": -3.0,
                                      "lon": -116.0, "lat": 33.0, "depth_m": 10.0}
    with pytest.raises(ConfigError, match="mixes projected"):
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)

    safs_raw["hypocenters"]["alt"] = {"x": 1.0}          # incomplete, one system only
    with pytest.raises(ConfigError, match="EITHER"):
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)

    safs_raw["hypocenters"]["alt"] = {}                  # neither
    with pytest.raises(ConfigError, match="EITHER"):
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)


def test_R003_partial_coordinate_mixing_is_rejected(tmp_path, safs_raw):
    """A complete projected point PLUS a stray lon must not silently drop the lon.

    This is the shape of a half-finished retarget: the user fills in lon/lat and forgets to
    delete x/y/z.  The old complete-set XOR passed it, and the run would have used the
    ORIGINAL project's hypocentre with a ~0 m snap distance -- undetectable downstream.
    """
    safs_raw["hypocenters"]["alt"] = {
        "x": 604446.944, "y": 3704576.3853, "z": -10067.9819, "lon": -121.0}
    with pytest.raises(ConfigError) as exc:
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)
    msg = str(exc.value)
    assert "mixes projected" in msg and "'lon'" in msg


def test_geographic_hypocenter_is_accepted():
    cfg = Project.load(PROJECTS / "demo_planar.yaml", require_files=False)
    h = cfg.hypocenter()
    assert h.is_geographic and not h.is_projected


@pytest.mark.parametrize("mutate,needle", [
    ({"dx": 0.0}, "dx"),
    ({"dz": -1.0}, "dz"),
    ({"xmax": 1.0}, "xmax"),
    ({"ymax": 1.0}, "ymax"),
    ({"zmax": -1e9}, "zmax"),
])
def test_bad_grid_box_raises(tmp_path, safs_raw, mutate, needle):
    safs_raw["stress_box"].update(mutate)
    with pytest.raises(ConfigError, match=needle):
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)


def test_fault_bc_absent_from_tag_map_raises(tmp_path, safs_raw):
    safs_raw["meshes"]["alt"]["fault_bc"] = 7
    with pytest.raises(ConfigError, match="fault_bc"):
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)


def test_missing_required_key_raises(tmp_path, safs_raw):
    del safs_raw["strike"]
    with pytest.raises(ConfigError, match="strike"):
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)


def test_missing_descriptor_file_raises():
    with pytest.raises(ConfigError, match="not found"):
        Project.load(PROJECTS / "does_not_exist.yaml", require_files=False)


def test_snap_tol_must_be_positive(tmp_path, safs_raw):
    safs_raw["hypocenters"]["alt"]["snap_tol_m"] = 0.0
    with pytest.raises(ConfigError, match="snap_tol_m"):
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)


# ------------------------------------------------- numeric typing (the YAML 1.1 trap)
def test_unsigned_exponent_is_rejected_with_a_fix_hint(tmp_path, safs_raw):
    """PyYAML is YAML 1.1: `23.5e9` parses as a STRING, not a float.

    Caught for real during Phase 0: the descriptors were written `23.5e9` and the shear
    modulus silently arrived as the string '23.5e9'.  A string flowing into physics is a
    silent-wrong-value bug, so the loader must reject it and say how to fix it.
    """
    safs_raw["physics"]["mu_shear_pa"] = "23.5e9"       # exactly what YAML hands us
    with pytest.raises(ConfigError) as exc:
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)
    msg = str(exc.value)
    assert "mu_shear_pa" in msg
    assert "must be a number" in msg
    assert "23.5e+9" in msg          # the suggested fix, with the explicit sign


def test_shipped_descriptors_have_numeric_physics():
    """Guards the descriptors themselves against regressing to the unsigned form."""
    for stem in ("safs_alt", "safs_preferred", "demo_planar"):
        cfg = Project.load(PROJECTS / f"{stem}.yaml", require_files=False)
        for fname in ("mu_shear_pa", "w_energy_m", "dc_m", "kappa_c", "rs_b", "rs_sl0"):
            val = getattr(cfg.physics, fname)
            assert isinstance(val, float), f"{stem}.{fname} is {type(val).__name__}"
    alt = Project.load(PROJECTS / "safs_alt.yaml", require_files=False)
    assert alt.physics.mu_shear_pa == 23500000000.0
    assert alt.physics.w_energy_m == 10000.0


@pytest.mark.parametrize("block,key,bad", [
    ("stress_box", "dx", "1.0e3"),
    ("strike", "azimuth_deg", "3.14e2"),
])
def test_numeric_guard_covers_other_blocks(tmp_path, safs_raw, block, key, bad):
    safs_raw[block][key] = bad
    with pytest.raises(ConfigError, match="must be a number"):
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)


def test_bool_is_not_accepted_as_a_number(tmp_path, safs_raw):
    safs_raw["physics"]["dc_m"] = True
    with pytest.raises(ConfigError, match="must be a number"):
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)


def test_origin_xy_elements_must_be_numeric(tmp_path, safs_raw):
    safs_raw["strike"]["origin_xy"] = [606971.0, "3.70727e6"]
    with pytest.raises(ConfigError, match=r"origin_xy\[1\]"):
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)


# ------------------------------------------------- require_files semantics
def test_require_files_true_reports_the_absolute_path_and_a_hint(tmp_path, safs_raw):
    p = write_tmp(tmp_path, safs_raw)
    with pytest.raises(ConfigError) as exc:
        Project.load(p, require_files=True, data_dir=tmp_path / "data")
    msg = str(exc.value)
    assert "mesh_alt.puml.h5" in msg
    assert str(tmp_path) in msg              # the RESOLVED absolute path
    assert "download_data.py" in msg         # the hint


def test_require_files_false_skips_only_existence_not_schema(tmp_path, safs_raw):
    """require_files=False must not weaken schema validation."""
    Project.load(write_tmp(tmp_path, safs_raw), require_files=False)   # ok
    safs_raw["default_mesh"] = "ghost"
    with pytest.raises(ConfigError, match="default_mesh"):
        Project.load(write_tmp(tmp_path, safs_raw, "q.yaml"), require_files=False)


def test_require_files_true_passes_when_files_exist(tmp_path, safs_raw):
    data = tmp_path / "data"
    (data / "raw" / "cvm").mkdir(parents=True)
    (data / "raw" / "ctm").mkdir(parents=True)
    (data / "mesh_alt.puml.h5").write_bytes(b"stub")
    (data / "raw" / "CSM_orientation.csv").write_text("x")
    (data / "raw" / "cvm" / "CVM_1_h_data.csv").write_text("x")     # matches the glob
    (data / "raw" / "ctm" / "CTM_1_h_data_final.csv").write_text("x")
    cfg = Project.load(write_tmp(tmp_path, safs_raw), require_files=True, data_dir=data)
    assert cfg.resolve_path("mesh_alt.puml.h5") == data / "mesh_alt.puml.h5"


# ------------------------------------------------- acceptance: legacy field-by-field
LEGACY = os.environ.get(
    "DECKBUILD_LEGACY_COMMON",
    str(Path.home() / "projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs/"
        "seisol_quakeworx/v3_under_construction/toolbox/combined_workflow/lib/_common.py"),
)


@pytest.mark.skipif(not Path(LEGACY).is_file(),
                    reason=f"legacy _common.py not available at {LEGACY}")
def test_safs_alt_matches_legacy_common_field_by_field():
    """Acceptance criterion 2: every descriptor field equals the legacy constant.

    The legacy module lives in an untracked tree outside this repo, so this test skips
    when it is absent.  Point DECKBUILD_LEGACY_COMMON at it to run.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location("_legacy_common", LEGACY)
    legacy = importlib.util.module_from_spec(spec)
    sys.modules["_legacy_common"] = legacy
    spec.loader.exec_module(legacy)          # may need numpy; that is a real dependency

    cfg = Project.load(PROJECTS / "safs_alt.yaml", require_files=False)

    assert cfg.strike.azimuth_deg == legacy.STRIKE_HINT_AZ
    assert cfg.physics.mu_shear_pa == legacy.MU_SHEAR
    assert cfg.physics.w_energy_m == legacy.W_ENERGY
    assert cfg.physics.dc_m == legacy.DC
    assert cfg.physics.kappa_c == legacy.KAPPA_C
    assert cfg.physics.seis_band_km == tuple(legacy.SEIS_BAND_KM)
    assert cfg.physics.active_band_km == tuple(legacy.ACTIVE_BAND_KM)
    assert cfg.physics.l_coast_km == legacy.L_COAST_KM
    assert cfg.physics.gap_max_km == legacy.GAP_MAX_KM
    assert cfg.physics.sn_floor_mpa == legacy.SN_FLOOR_MPA

    for key in ("xmin", "xmax", "ymin", "ymax", "dx", "dz"):
        assert getattr(cfg.stress_box, key) == legacy.ALT_BOX[key], key

    legacy_bands = {(a, b, n) for a, b, n in legacy.GATE_BANDS}
    ours = {(b.s_start_km, b.s_end_km, b.name) for b in cfg.gate_bands}
    assert ours <= legacy_bands, f"gate bands not a subset of legacy: {ours - legacy_bands}"

    pref = Project.load(PROJECTS / "safs_preferred.yaml", require_files=False)
    for key in ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax", "dx", "dz"):
        assert getattr(pref.stress_box, key) == legacy.PREF_BOX[key], key
    pref_bands = {(b.s_start_km, b.s_end_km, b.name) for b in pref.gate_bands}
    assert pref_bands == legacy_bands


# ------------------------------------------------- review findings R-002 / R-004 / R-006
def test_R002_outputs_placeholder_is_committable():
    """The plan lists outputs/.gitkeep; `outputs/` in .gitignore would swallow it.

    git never descends into an excluded DIRECTORY, so `outputs/` + `!outputs/.gitkeep`
    does not work -- the contents form `outputs/*` is required.
    """
    import subprocess
    assert (ROOT / "outputs" / ".gitkeep").is_file()
    r = subprocess.run(["git", "check-ignore", "-v", "outputs/.gitkeep"],
                       cwd=ROOT, capture_output=True, text=True)
    # Either no pattern matches, or the matching pattern is the negation.
    assert r.returncode != 0 or r.stdout.lstrip().split(":")[-1].strip().startswith(
        "!outputs/.gitkeep") or "!outputs/.gitkeep" in r.stdout, \
        f"outputs/.gitkeep is ignored: {r.stdout!r}"
    # ...and build products must still be ignored.
    r2 = subprocess.run(["git", "check-ignore", "outputs/scratch.nc"],
                        cwd=ROOT, capture_output=True, text=True)
    assert r2.returncode == 0, "outputs/*.nc must stay ignored"


def test_R004_descriptor_mappings_are_read_only():
    """frozen=True blocks rebinding only; the dict fields must be proxies too."""
    cfg = Project.load(PROJECTS / "safs_alt.yaml", require_files=False)
    with pytest.raises(TypeError):
        cfg.meshes["alt"].tag_to_bc[999] = 1
    with pytest.raises(TypeError):
        cfg.raw.velocity.params["grid_dx"] = 1.0
    with pytest.raises(TypeError):
        cfg.meshes["nope"] = None
    with pytest.raises(TypeError):
        cfg.hypocenters["nope"] = None
    with pytest.raises(TypeError):
        cfg.provenance["x"] = "y"
    with pytest.raises((TypeError, AttributeError)):
        cfg.gate_bands.append(None)


def test_R004_frozen_mappings_still_compare_and_serialise():
    """MappingProxyType must not break equality checks or as_dict/json."""
    import json
    cfg = Project.load(PROJECTS / "safs_alt.yaml", require_files=False)
    assert cfg.meshes["alt"].tag_to_bc == {101: 3, 102: 1, 103: 5, 104: 5}
    d = as_dict(cfg)
    assert isinstance(d["meshes"]["alt"]["tag_to_bc"], dict)
    json.dumps(d)          # must not raise on a mappingproxy


def test_R006_hint_covers_a_leading_dot_exponent(tmp_path, safs_raw):
    safs_raw["physics"]["rs_b"] = ".19e2"      # YAML 1.1 reads this as a string
    with pytest.raises(ConfigError) as exc:
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)
    assert "e+" in str(exc.value)              # the fix hint fired
