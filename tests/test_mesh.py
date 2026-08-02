"""Phase 5 acceptance tests: mesh ingest, gating, and Stage F."""
from __future__ import annotations

import difflib
import re
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from deckbuild.config import MeshSpec, Project  # noqa: E402
from deckbuild.geometry import FACE_VERTS, load_fault  # noqa: E402
from deckbuild.mesh import (  # noqa: E402
    MeshError, MeshStage, count_face_multiplicity, msh_to_puml, orient_tets_positive,
    pack_boundary, read_puml, signed_volumes, unpack_boundary, write_puml,
)
from deckbuild.stage_f import verify_mesh_against_deck  # noqa: E402

PROJECTS = ROOT / "projects"
DEMO_MESH = ROOT / "data" / "demo_planar" / "demo_planar.puml.h5"


@pytest.fixture()
def demo():
    return Project.load(PROJECTS / "demo_planar.yaml", require_files=False)


# --------------------------------------------------------------------------- packing
def test_boundary_packing_roundtrips():
    """byte i of `boundary` holds the BC of face i -- the SeisSol contract."""
    bc = np.array([[3, 1, 5, 0], [0, 0, 0, 0], [5, 5, 5, 5]], np.int64)
    assert np.array_equal(unpack_boundary(pack_boundary(bc)), bc)


def test_boundary_packing_is_little_endian_by_face():
    packed = pack_boundary(np.array([[3, 1, 5, 0]], np.int64))
    assert int(packed[0]) == 3 | (1 << 8) | (5 << 16) | (0 << 24)


def test_packing_rejects_an_oversized_code():
    with pytest.raises(MeshError, match="does not fit"):
        pack_boundary(np.array([[300, 0, 0, 0]], np.int64))


def test_packing_rejects_a_bad_shape():
    with pytest.raises(MeshError, match=r"\(Ntet, 4\)"):
        pack_boundary(np.zeros((3, 3), np.int64))


def test_unpack_rejects_an_unknown_format():
    with pytest.raises(MeshError, match="boundary-format"):
        unpack_boundary(np.zeros(3, np.int64), fmt="i8")


# --------------------------------------------------------------------------- orientation
def test_orient_flips_only_negative_tets():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    good = np.array([[0, 1, 2, 3]], np.int64)
    bad = np.array([[0, 2, 1, 3]], np.int64)
    _, nf = orient_tets_positive(pts, good)
    assert nf == 0
    fixed, nf = orient_tets_positive(pts, bad)
    assert nf == 1 and signed_volumes(pts, fixed)[0] > 0


def test_orientation_is_idempotent():
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    t1, _ = orient_tets_positive(pts, np.array([[0, 2, 1, 3]], np.int64))
    t2, nf = orient_tets_positive(pts, t1)
    assert nf == 0 and np.array_equal(t1, t2)


# --------------------------------------------------------------------------- face keys
def test_face_multiplicity_uses_searchsorted_not_a_packed_key():
    """The legacy packed int64 key overflows above ~2.1M nodes; this path must not.

    Node indices near 2^21 would overflow a 3x21-bit packed key.  A structured-view
    lexsort has no such ceiling.
    """
    big = 2_200_000
    pts_idx = np.array([[0, 1, 2, 3], [1, 2, 3, big]], np.int64)
    mult, n_uniq = count_face_multiplicity(pts_idx)
    assert mult.shape == (8,)
    # faces {1,2,3} is shared by both tets; everything else is unique.
    assert int((mult == 2).sum()) == 2
    assert n_uniq == 7


def test_face_multiplicity_on_the_demo_mesh():
    geom, conn, bc, _ = read_puml(DEMO_MESH)
    mult, _ = count_face_multiplicity(conn)
    mult = mult.reshape(4, len(conn)).T
    # A fault face is INTERIOR: shared by exactly two tets.
    assert np.all(mult[bc == 3] == 2)
    # A free-surface/absorbing face is on the hull: exactly one.
    assert np.all(mult[(bc > 0) & (bc != 3)] == 1)


# --------------------------------------------------------------------------- puml io
def test_puml_roundtrip(tmp_path):
    geom, conn, bc, grp = read_puml(DEMO_MESH)
    p = write_puml(tmp_path / "rt.puml.h5", geom, conn, bc, grp)
    g2, c2, b2, _ = read_puml(p)
    assert np.array_equal(geom, g2) and np.array_equal(conn, c2)
    assert np.array_equal(bc, b2)


def test_puml_attrs_are_the_seissol_contract(tmp_path):
    import h5py
    geom, conn, bc, _ = read_puml(DEMO_MESH)
    p = write_puml(tmp_path / "a.puml.h5", geom, conn, bc)
    with h5py.File(p) as f:
        assert f.attrs["boundary-format"] == "i32"
        assert f.attrs["topology-format"] == "geometric"
        assert f["geometry"].dtype == np.float64
        assert f["connect"].dtype == np.uint64


def test_read_puml_missing_dataset_raises(tmp_path):
    import h5py
    p = tmp_path / "bad.puml.h5"
    with h5py.File(p, "w") as f:
        f.create_dataset("geometry", data=np.zeros((3, 3)))
    with pytest.raises(MeshError, match="missing dataset"):
        read_puml(p)


def test_read_puml_missing_file_raises(tmp_path):
    with pytest.raises(MeshError, match="not found"):
        read_puml(tmp_path / "nope.puml.h5")


# --------------------------------------------------------------------------- the stage
def test_ingest_an_existing_puml(demo):
    art = MeshStage().build(demo, ROOT / "outputs")
    assert Path(art.path) == DEMO_MESH
    assert art.params["ingested"] is True


def test_build_without_a_mesh_says_where_to_get_one(tmp_path):
    cfg = Project.load(PROJECTS / "demo_planar.yaml", require_files=False,
                       data_dir=tmp_path)
    with pytest.raises(MeshError) as exc:
        MeshStage().build(cfg, tmp_path)
    msg = str(exc.value)
    assert "INGESTS a mesh" in msg and "skills/" in msg and "MESHING.md" in msg


def test_gates_pass_on_the_demo_mesh(demo):
    rep = MeshStage().verify(demo, MeshStage().build(demo, ROOT / "outputs"))
    assert rep.ok, "\n".join(rep.lines())
    assert {g.name for g in rep.gates} >= {"A", "B", "C", "D", "E", "G"}


def test_gate_thresholds_come_from_the_descriptor(demo):
    """500 m / 0.6667 Hz are SAFS values; hardcoding them would be a constant leak."""
    assert demo.mesh().gates.fault_edge_max_m == 4000.0        # the demo's own
    safs = Project.load(PROJECTS / "safs_alt.yaml", require_files=False)
    assert safs.mesh().gates.fault_edge_max_m == 500.0         # the SAFS default
    assert safs.mesh().gates.f_gate_hz == pytest.approx(0.6667)


def test_gate_B_fails_at_a_tighter_limit(demo):
    art = MeshStage().build(demo, ROOT / "outputs")
    rep = MeshStage().verify(demo, art, fault_edge_max_m=100.0)
    b = [g for g in rep.gates if g.name == "B"][0]
    assert not b.passed and "3535" in b.detail


def test_gate_E_skips_without_a_material_nc_and_a_skip_is_not_a_pass(demo):
    rep = MeshStage().verify(demo, MeshStage().build(demo, ROOT / "outputs"))
    e = [g for g in rep.gates if g.name == "E"][0]
    assert e.severity == "skip" and e.passed is False


def test_gate_C_catches_a_fault_face_on_the_hull(tmp_path, demo):
    """A fault face MUST be interior (x2).  Tagging a hull face as fault must fail C."""
    geom, conn, bc, _ = read_puml(DEMO_MESH)
    mult, _ = count_face_multiplicity(conn)
    mult = mult.reshape(4, len(conn)).T
    t, f = np.argwhere(mult == 1)[0]
    bc2 = bc.copy()
    bc2[t, f] = 3                       # mislabel a hull face as the fault
    p = write_puml(tmp_path / "bad.puml.h5", geom, conn, bc2)
    cfg = Project.load(PROJECTS / "demo_planar.yaml", require_files=False,
                       data_dir=tmp_path)
    import dataclasses
    spec = dataclasses.replace(demo.mesh(), path=str(p))
    object.__setattr__(cfg, "meshes", {"planar": spec})
    rep = MeshStage().verify(cfg, MeshStage().build(cfg, tmp_path))
    c = [g for g in rep.gates if g.name == "C"][0]
    assert not c.passed and "not interior" in c.detail


def test_gate_D_catches_an_inverted_tet(tmp_path, demo):
    geom, conn, bc, _ = read_puml(DEMO_MESH)
    conn2 = conn.copy()
    conn2[0] = conn2[0][[0, 2, 1, 3]]              # invert one tet
    p = write_puml(tmp_path / "inv.puml.h5", geom, conn2, bc)
    import dataclasses
    cfg = Project.load(PROJECTS / "demo_planar.yaml", require_files=False,
                       data_dir=tmp_path)
    object.__setattr__(cfg, "meshes",
                       {"planar": dataclasses.replace(demo.mesh(), path=str(p))})
    rep = MeshStage().verify(cfg, MeshStage().build(cfg, tmp_path))
    d = [g for g in rep.gates if g.name == "D"][0]
    assert not d.passed and "inverted tets: 1" in d.detail


# --------------------------------------------------------------------------- Stage F
def test_stage_f_refuses_a_descriptor_mismatch(demo):
    """Checking a mesh against fields from a DIFFERENT descriptor is worse than not."""
    rep = verify_mesh_against_deck(demo, descriptor_sha256="0" * 64)
    assert not rep.ok
    assert "refusing to cross-check" in rep.gates[0].detail


def test_stage_f_catches_an_off_fault_hypocentre(tmp_path):
    """The gate that caught a real bug in the demo descriptor itself."""
    import yaml
    raw = yaml.safe_load((PROJECTS / "demo_planar.yaml").read_text())
    raw["hypocenters"]["planar"] = {"lon": -121.0, "lat": 36.0, "depth_m": 8000.0,
                                    "snap_tol_m": 2000.0}
    p = tmp_path / "off.yaml"
    p.write_text(yaml.safe_dump(raw))
    cfg = Project.load(p, require_files=False,
                       data_dir=ROOT / "data" / "demo_planar")
    rep = verify_mesh_against_deck(cfg)
    f5 = [g for g in rep.gates if g.name.startswith("F5")]
    assert f5 and not f5[0].passed
    assert "does not land on this mesh" in f5[0].detail


def test_stage_f_passes_on_the_demo(demo):
    rep = verify_mesh_against_deck(demo)
    assert rep.ok, "\n".join(rep.lines())


def test_receiver_surface_check():
    from deckbuild.contract import GateReport
    from deckbuild.stage_f import check_receivers_under_the_surface
    rep = GateReport("t")
    check_receivers_under_the_surface(np.array([[0, 0, 0.0], [0, 0, -1.0]]), rep, "P8")
    assert not rep.ok and "silently drops" in rep.gates[0].detail
    rep2 = GateReport("t")
    check_receivers_under_the_surface(np.array([[0, 0, -1.0]]), rep2, "P8")
    assert rep2.ok


# --------------------------------------------------------------------------- the skill
def test_vendored_skill_differs_only_in_the_two_allowed_edits():
    """Only two edits are allowed when vendoring: the reference-path relabel and the
    tag-contract pointer.  Everything else ships verbatim."""
    src = Path.home() / ".claude/skills/code-mesh-build-improve/SKILL.md"
    dst = ROOT / "skills" / "code-mesh-build-improve" / "SKILL.md"
    assert dst.is_file()
    if not src.is_file():
        pytest.skip("source skill not available on this machine")
    diff = [l for l in difflib.unified_diff(src.read_text().split("\n"),
                                            dst.read_text().split("\n"),
                                            lineterm="", n=0)
            if l.startswith(("+", "-")) and not l.startswith(("+++", "---"))]
    added = "\n".join(l[1:] for l in diff if l[0] == "+")
    assert "WORKED EXAMPLES" in added, "the reference-path relabel is missing"
    assert "MeshSpec.tag_to_bc" in added, "the tag-contract pointer is missing"
    assert len(diff) <= 12, f"more than the two allowed hunks changed:\n" + "\n".join(diff)


def test_meshing_doc_names_the_stages_and_the_rules():
    t = (ROOT / "MESHING.md").read_text()
    for stage in "ABCDEF":
        assert f"\n{stage}  " in t or f"| {stage} |" in t, f"stage {stage} missing"
    assert "v2.2" in t and "Never move the fault surface" in t
    assert "meshing_<variant>" in t
    assert "does not build your mesh" in t.lower() or "not BUILD" in t


# ------------------------------------------------- review findings R-501..R-503
def _puml_to_msh22(src, dst):
    """Round-trip helper: write a .msh v2.2 carrying the same tagged surfaces."""
    import meshio
    geom, conn, bc, _ = read_puml(src)
    bc2tag = {3: 101, 1: 102, 5: 104}
    tris, tags = [], []
    for f, fv in enumerate(FACE_VERTS):
        m = bc[:, f] > 0
        if m.any():
            tris.append(conn[m][:, list(fv)])
            tags.append([bc2tag[int(v)] for v in bc[m, f]])
    tris = np.concatenate(tris)
    tags = np.concatenate(tags).astype(np.int32)
    meshio.write(str(dst), meshio.Mesh(
        points=geom, cells=[("tetra", conn), ("triangle", tris)],
        cell_data={"gmsh:physical": [np.ones(len(conn), np.int32), tags]}),
        file_format="gmsh22", binary=False)
    return geom, conn, bc


def test_R501_msh_roundtrip_tags_both_sides_of_the_fault(tmp_path):
    """A fault triangle is INTERIOR: both tet-face slots must carry the BC.

    Taking only the first searchsorted match left one side at 0 -- 192 of 9216 face slots
    on the demo mesh, exactly the fault-triangle count.
    """
    msh = tmp_path / "rt.msh"
    geom, conn, bc = _puml_to_msh22(DEMO_MESH, msh)
    out, info = msh_to_puml(msh, tmp_path / "rt.puml.h5", {101: 3, 102: 1, 103: 5, 104: 5})
    g2, c2, b2, _ = read_puml(out)
    assert np.array_equal(geom, g2)
    assert np.array_equal(bc, b2), f"{int((bc != b2).sum())} face slots differ"
    # and the fault is still interior on both sides
    mult, _ = count_face_multiplicity(c2)
    assert np.all(mult.reshape(4, len(c2)).T[b2 == 3] == 2)


def test_R501_converted_mesh_passes_gate_C(tmp_path, demo):
    """The end-to-end consequence: a converted mesh must pass the BC round-trip gate."""
    import dataclasses
    msh = tmp_path / "rt.msh"
    _puml_to_msh22(DEMO_MESH, msh)
    out, _ = msh_to_puml(msh, tmp_path / "rt.puml.h5", {101: 3, 102: 1, 103: 5, 104: 5})
    cfg = Project.load(PROJECTS / "demo_planar.yaml", require_files=False,
                       data_dir=tmp_path)
    object.__setattr__(cfg, "meshes",
                       {"planar": dataclasses.replace(demo.mesh(), path=str(out))})
    rep = MeshStage().verify(cfg, MeshStage().build(cfg, tmp_path))
    assert rep.ok, "\n".join(rep.lines())


def test_R502_no_removed_numpy_core_api():
    """np.core was removed in NumPy 2; np.rec is the replacement.

    Scans CODE only -- the module comment legitimately names the removed API to explain
    why it is not used.
    """
    code = [ln for ln in (ROOT / "deckbuild" / "mesh.py").read_text().splitlines()
            if not ln.lstrip().startswith("#")]
    code = re.sub(r'""".*?"""', "", "\n".join(code), flags=re.S)
    assert "np.core." not in code


def test_msh_unmapped_tag_raises(tmp_path):
    """Silently dropping a tag would lose a boundary condition."""
    msh = tmp_path / "rt.msh"
    _puml_to_msh22(DEMO_MESH, msh)
    with pytest.raises(MeshError, match="not in tag_to_bc"):
        msh_to_puml(msh, tmp_path / "x.puml.h5", {101: 3})       # 102/104 unmapped


def test_msh_missing_file_raises(tmp_path):
    with pytest.raises(MeshError, match="not found"):
        msh_to_puml(tmp_path / "nope.msh", tmp_path / "o.puml.h5", {101: 3})


def test_msh_v4_gives_an_actionable_message(tmp_path):
    """The v4 symptom is misleading; the message must be searchable."""
    bad = tmp_path / "b.msh"
    bad.write_text("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
    with pytest.raises(MeshError) as exc:
        msh_to_puml(bad, tmp_path / "o.puml.h5", {101: 3})
    m = str(exc.value)
    assert "msh22" in m and "vertices indices are not unique" in m
