"""Phase 6 acceptance tests: deck assembly and the P1-P8 pre-flight battery."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from deckbuild.config import Project  # noqa: E402
from deckbuild.contract import Manifest  # noqa: E402
from deckbuild.deck import DeckError, DeckSpec, DeckStage, resolve_paths  # noqa: E402
from deckbuild.friction import FrictionStage, FwDesign  # noqa: E402
from deckbuild.material import MaterialStage, PlasticitySpec  # noqa: E402
from deckbuild.mesh import MeshStage  # noqa: E402
from deckbuild.stress import KDesign, StressStage  # noqa: E402

PROJECTS = ROOT / "projects"


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """One full build, reused across the deck tests (it is the expensive part)."""
    out = tmp_path_factory.mktemp("build")
    cfg = Project.load(PROJECTS / "demo_planar.yaml")
    mat = MaterialStage().build(cfg, out, plasticity=PlasticitySpec())
    st = StressStage().build(cfg, out, sv_profile=mat.sv_profile,
                             design=KDesign(k_values=(1.7,)), freeze_above_depth_m=500.0)
    fr = FrictionStage().build(cfg, out, case=1)
    fw = FrictionStage().build_fw_map(
        cfg, out, FwDesign(fw_values=(0.0, 0.05), boundaries_s_km=((10.0, 20.0),)))
    me = MeshStage().build(cfg, out)
    arts = {"material": mat.material, "plasticity": mat.plasticity,
            "stress": st, "friction": fr, "mesh": me}
    return cfg, arts, fw, mat


def _spec(cfg, fw, notes=""):
    rx = np.array([[cfg.stress_box.xmin + 5000.0 * i,
                    cfg.stress_box.ymin + 5000.0 * i, -1.0] for i in range(1, 6)])
    return DeckSpec(prefix=f"{cfg.name}_", plasticity=True, receivers=rx,
                    rs_muw_lua=Path(fw.path).read_text(), mu_s=0.6, notes=notes)


# --------------------------------------------------------------------------- paths
def test_resolve_paths_writes_nothing_and_encodes_the_design(built, tmp_path):
    cfg, arts, fw, _ = built
    m = resolve_paths(tmp_path / "d", arts, _spec(cfg, fw))
    assert not (tmp_path / "d").exists()               # predicted, not created
    assert m["stress"]["filename"] == "demo_planar_stress_k1.7_freeze500m.nc"
    assert m["plasticity"]["filename"] == "demo_planar_plasticity_phi35_45.nc"
    assert m["friction"]["filename"] == "demo_planar_friction_case1.nc"
    assert m["stress"]["referenced_by"] == "demo_planar_initial_stress.yaml"


def test_prefix_comes_from_the_descriptor_not_safs(built, tmp_path):
    cfg, arts, fw, _ = built
    for e in resolve_paths(tmp_path / "d", arts, _spec(cfg, fw)).values():
        assert not Path(e["filename"]).name.startswith("safs_")
    assert all(e["filename"].startswith("demo_planar_") or e["filename"].endswith(".puml.h5")
               for e in resolve_paths(tmp_path / "d", arts, _spec(cfg, fw)).values())


def test_predicted_map_equals_the_assembled_map(built, tmp_path):
    cfg, arts, fw, _ = built
    d = tmp_path / "deck"
    predicted = {k: v["filename"] for k, v in resolve_paths(d, arts, _spec(cfg, fw)).items()}
    DeckStage().assemble(cfg, d, arts, _spec(cfg, fw))
    for name in predicted.values():
        assert (d / name).is_file(), name


# --------------------------------------------------------------------------- assemble
def test_assemble_copies_and_verifies(built, tmp_path):
    cfg, arts, fw, _ = built
    d = DeckStage().assemble(cfg, tmp_path / "deck", arts, _spec(cfg, fw))
    for f in d.glob("*.nc"):
        assert not f.is_symlink(), "a deck must COPY, never symlink -- it gets rsynced"
    assert (d / Manifest.FILENAME).is_file()
    assert (d / "parameters.par").is_file()
    assert len(list(d.glob("*.yaml"))) == 3


def test_assemble_refuses_to_overwrite_by_default(built, tmp_path):
    cfg, arts, fw, _ = built
    d = tmp_path / "deck"
    DeckStage().assemble(cfg, d, arts, _spec(cfg, fw))
    with pytest.raises(DeckError, match="overwrite=True"):
        DeckStage().assemble(cfg, d, arts, _spec(cfg, fw))
    DeckStage().assemble(cfg, d, arts, _spec(cfg, fw), overwrite=True)


def test_hand_written_notes_survive_regeneration(built, tmp_path):
    """The shipped decks' value is largely their provenance prose."""
    cfg, arts, fw, _ = built
    d = tmp_path / "deck"
    DeckStage().assemble(cfg, d, arts, _spec(cfg, fw, notes="# why k=1.7: see run 52538589"))
    assert "why k=1.7" in (d / f"{cfg.name}_fault.yaml").read_text()
    DeckStage().assemble(cfg, d, arts, _spec(cfg, fw), overwrite=True)   # no notes passed
    assert "why k=1.7" in (d / f"{cfg.name}_fault.yaml").read_text()


def test_assemble_refuses_mixed_descriptors(built, tmp_path):
    """Artifacts from two descriptors would produce plausible, meaningless numbers."""
    import copy
    cfg, arts, fw, _ = built
    bad = dict(arts)
    other = copy.deepcopy(arts["stress"])
    other.provenance = {"descriptor_sha256": "0" * 64}
    bad["stress"] = other
    with pytest.raises(DeckError, match="different descriptors|refusing to assemble"):
        DeckStage().assemble(cfg, tmp_path / "d2", bad, _spec(cfg, fw))


def test_yaml_file_field_matches_the_copied_filename(built, tmp_path):
    """The two are written from the same variable; P1 then re-checks on disk."""
    cfg, arts, fw, _ = built
    d = DeckStage().assemble(cfg, tmp_path / "deck", arts, _spec(cfg, fw))
    for y in d.glob("*.yaml"):
        for line in y.read_text().splitlines():
            s = line.strip()
            if s.startswith("file:"):
                assert (d / s.split(":", 1)[1].strip()).is_file(), s


# --------------------------------------------------------------------------- preflight
def test_preflight_all_pass(built, tmp_path):
    cfg, arts, fw, _ = built
    spec = _spec(cfg, fw)
    d = DeckStage().assemble(cfg, tmp_path / "deck", arts, spec)
    rep = DeckStage().preflight(cfg, d, spec=spec)
    assert rep.ok, "\n".join(rep.lines())
    assert {g.name for g in rep.gates} == {"P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"}


def test_P1_flags_an_unaccounted_extra_nc(built, tmp_path):
    """An unexplained nc in a deck is how a stale field gets shipped."""
    cfg, arts, fw, _ = built
    spec = _spec(cfg, fw)
    d = DeckStage().assemble(cfg, tmp_path / "deck", arts, spec)
    (d / "leftover_stress_k9.nc").write_bytes(b"stale")
    rep = DeckStage().preflight(cfg, d, spec=spec)
    p1 = [g for g in rep.gates if g.name == "P1"][0]
    assert not p1.passed and "leftover_stress_k9.nc" in p1.detail


def test_P1_flags_a_missing_target(built, tmp_path):
    cfg, arts, fw, _ = built
    spec = _spec(cfg, fw)
    d = DeckStage().assemble(cfg, tmp_path / "deck", arts, spec)
    next(d.glob("*stress*.nc")).unlink()
    rep = DeckStage().preflight(cfg, d, spec=spec)
    assert not [g for g in rep.gates if g.name == "P1"][0].passed


def test_P1_flags_a_tampered_file(built, tmp_path):
    """sha256 against the manifest, so an edited nc cannot pass."""
    cfg, arts, fw, _ = built
    spec = _spec(cfg, fw)
    d = DeckStage().assemble(cfg, tmp_path / "deck", arts, spec)
    f = next(d.glob("*friction*.nc"))
    with open(f, "r+b") as fh:
        fh.seek(0)
        fh.write(b"\x00" * 8)
    rep = DeckStage().preflight(cfg, d, spec=spec)
    p1 = [g for g in rep.gates if g.name == "P1"][0]
    assert not p1.passed and "sha256" in p1.detail


def test_P1_accounts_for_the_mesh_named_in_parameters_par(built, tmp_path):
    """The mesh is referenced by MeshFile, not a yaml `file:`."""
    cfg, arts, fw, _ = built
    spec = _spec(cfg, fw)
    d = DeckStage().assemble(cfg, tmp_path / "deck", arts, spec)
    rep = DeckStage().preflight(cfg, d, spec=spec)
    p1 = [g for g in rep.gates if g.name == "P1"][0]
    assert p1.passed and "extras none" in p1.detail


def test_P5_interpolates_across_independent_grids(built, tmp_path):
    """The four grids are independent by design; P5 must sample, not skip."""
    cfg, arts, fw, _ = built
    spec = _spec(cfg, fw)
    d = DeckStage().assemble(cfg, tmp_path / "deck", arts, spec)
    rep = DeckStage().preflight(cfg, d, spec=spec)
    p5 = [g for g in rep.gates if g.name == "P5"][0]
    assert p5.severity != "skip", p5.detail
    assert p5.passed and "yield" in p5.detail


def test_P6_is_a_warning_not_a_prediction(built, tmp_path):
    """kappa is necessary, NOT sufficient -- a design passed it and still arrested."""
    cfg, arts, fw, _ = built
    spec = _spec(cfg, fw)
    d = DeckStage().assemble(cfg, tmp_path / "deck", arts, spec)
    rep = DeckStage().preflight(cfg, d, spec=spec)
    p6 = [g for g in rep.gates if g.name == "P6"][0]
    assert p6.severity == "warn"
    assert "NECESSARY NOT SUFFICIENT" in p6.detail


def test_P8_catches_a_surface_receiver(built, tmp_path):
    """SeisSol v1.1.3 silently DROPS a receiver at z = 0."""
    cfg, arts, fw, _ = built
    spec = _spec(cfg, fw)
    spec.receivers = np.array([[cfg.stress_box.xmin, cfg.stress_box.ymin, 0.0]])
    d = DeckStage().assemble(cfg, tmp_path / "deck", arts, spec)
    rep = DeckStage().preflight(cfg, d, spec=spec)
    p8 = [g for g in rep.gates if g.name == "P8"][0]
    assert not p8.passed and "silently drops" in p8.detail


def test_P8_catches_a_literal_zero_in_output_bounds(built, tmp_path):
    """A literal 0.0 in OutputRegionBounds writes the WHOLE volume."""
    cfg, arts, fw, _ = built
    spec = _spec(cfg, fw)
    d = DeckStage().assemble(cfg, tmp_path / "deck", arts, spec)
    par = d / "parameters.par"
    par.write_text(par.read_text().replace(
        "wavefieldoutput = 1",
        "wavefieldoutput = 1\nOutputRegionBounds = 0.0 1.0 0.0 1.0 0.0 1.0"))
    rep = DeckStage().preflight(cfg, d, spec=spec)
    p8 = [g for g in rep.gates if g.name == "P8"][0]
    assert not p8.passed and "WHOLE volume" in p8.detail


def test_P8_requires_wavefieldoutput(built, tmp_path):
    """0 hangs in teardown -- 0 of 25 runs completed without it."""
    cfg, arts, fw, _ = built
    spec = _spec(cfg, fw)
    d = DeckStage().assemble(cfg, tmp_path / "deck", arts, spec)
    par = d / "parameters.par"
    par.write_text(par.read_text().replace("wavefieldoutput = 1", "wavefieldoutput = 0"))
    rep = DeckStage().preflight(cfg, d, spec=spec)
    assert not [g for g in rep.gates if g.name == "P8"][0].passed


def test_P7_uses_the_snapped_hypocentre(built, tmp_path):
    cfg, arts, fw, _ = built
    spec = _spec(cfg, fw)
    d = DeckStage().assemble(cfg, tmp_path / "deck", arts, spec)
    rep = DeckStage().preflight(cfg, d, spec=spec)
    p7 = [g for g in rep.gates if g.name == "P7"][0]
    assert p7.passed and "snap" in p7.detail


# --------------------------------------------------------------------------- diff
def test_diff_reports_exactly_one_changed_file(built, tmp_path):
    """An A/B pair must be shown to differ in exactly one input."""
    cfg, arts, fw, mat = built
    spec = _spec(cfg, fw)
    a = DeckStage().assemble(cfg, tmp_path / "a", arts, spec)
    b = DeckStage().assemble(cfg, tmp_path / "b", arts, spec)
    d0 = DeckStage().diff(a, b)
    # deck_manifest.json embeds nothing time-varying, so the decks are identical.
    assert d0["differ"] == [] and not d0["only_in_a"] and not d0["only_in_b"]

    f = next(b.glob("*friction*.nc"))
    with open(f, "r+b") as fh:
        fh.seek(0)
        fh.write(b"\x01" * 8)
    d1 = DeckStage().diff(a, b)
    assert d1["differ"] == [f.name]
