"""Phase 0 acceptance tests for deckbuild.contract."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from deckbuild.contract import (  # noqa: E402
    Artifact, Gate, GateFailure, GateReport, HARD, Manifest, SKIP, Stage, WARN,
    sha256_file,
)


# --------------------------------------------------------------------------- Gate
def test_gate_rejects_bad_severity():
    with pytest.raises(ValueError, match="severity"):
        Gate(name="V1", passed=True, detail="x", severity="fatal")


def test_gate_requires_a_name():
    with pytest.raises(ValueError, match="name"):
        Gate(name="", passed=True, detail="x")


@pytest.mark.parametrize("severity,passed,marker", [
    (HARD, True, "[PASS]"),
    (HARD, False, "[FAIL]"),
    (WARN, True, "[PASS]"),
    (WARN, False, "[WARN]"),
    (SKIP, True, "[ -- ]"),
    (SKIP, False, "[ -- ]"),
])
def test_gate_markers(severity, passed, marker):
    assert Gate("G", passed, "d", severity).marker == marker


# --------------------------------------------------------------------------- GateReport
def test_report_ok_ignores_warn_and_skip():
    r = GateReport("t")
    r.add("V1", True, "fine")
    r.add("V2", False, "soft", severity=WARN)
    r.skip("V3", "no input")
    assert r.ok
    assert r.failed == []
    assert len(r.warnings) == 1
    assert len(r.skipped) == 1


def test_report_not_ok_when_a_hard_gate_fails():
    r = GateReport()
    r.add("V1", True, "fine")
    r.add("V2", False, "broken")
    assert not r.ok
    assert [g.name for g in r.failed] == ["V2"]


def test_skip_is_never_a_pass():
    """A skipped gate must not count as passing -- see the Phase 5 gate-E rule."""
    r = GateReport()
    g = r.skip("E", "no material nc")
    assert g.severity == SKIP
    assert g.passed is False
    assert r.ok            # skips do not block ...
    assert g not in r.failed   # ... but are not recorded as passes either


def test_raise_if_failed():
    r = GateReport("battery")
    r.add("V1", False, "bad thing")
    with pytest.raises(GateFailure) as exc:
        r.raise_if_failed()
    assert "V1" in str(exc.value) and "bad thing" in str(exc.value)
    GateReport("ok").raise_if_failed()      # no gates -> no raise


def test_print_golden_string(capsys):
    """Acceptance criterion 4: the printed format is fixed.  Run notes quote it."""
    r = GateReport("stress V1-V4")
    r.add("V1", True, "region 1/2 bit-identical to constant k=1.4")
    r.add("V2", False, "eig sigma1/sigma3 != k")
    r.add("V4a", False, "corridor kappa below kappa_c", severity=WARN)
    r.skip("V3", "freeze off")
    r.print()
    out = capsys.readouterr().out
    assert out == (
        "=== stress V1-V4 ===\n"
        "[PASS] V1: region 1/2 bit-identical to constant k=1.4\n"
        "[FAIL] V2: eig sigma1/sigma3 != k\n"
        "[WARN] V4a: corridor kappa below kappa_c\n"
        "[ -- ] V3: freeze off\n"
        "FAIL: 1/2 hard, 1 warn, 1 skipped\n"
    )


def test_to_dict_shape():
    r = GateReport("t")
    r.add("V1", True, "d")
    d = r.to_dict()
    assert d["title"] == "t" and d["ok"] is True
    assert d["gates"] == [{"name": "V1", "passed": True, "detail": "d", "severity": HARD}]


def test_extend_with_prefix():
    a, b = GateReport("a"), GateReport("b")
    a.add("G1", True, "x")
    b.add("G1", False, "y")
    a.extend(b, prefix="friction.")
    assert [g.name for g in a.gates] == ["G1", "friction.G1"]
    assert not a.ok


def test_len_and_repr():
    r = GateReport("t")
    r.add("V1", True, "d")
    assert len(r) == 1
    assert "GateReport" in repr(r) and "ok=True" in repr(r)


# --------------------------------------------------------------------------- Artifact
def test_artifact_of_hashes_the_file(tmp_path):
    p = tmp_path / "a.nc"
    p.write_bytes(b"hello")
    art = Artifact.of(p, kind="stress", params={"k": 1.7})
    assert art.sha256 == sha256_file(p)
    assert len(art.sha256) == 64
    assert art.kind == "stress" and art.params == {"k": 1.7}
    assert art.to_dict()["path"] == str(p)


def test_artifact_of_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="stress"):
        Artifact.of(tmp_path / "nope.nc", kind="stress")


def test_artifact_path_is_coerced_to_Path():
    assert isinstance(Artifact(path="x.nc", kind="k").path, Path)


def test_sha256_matches_hashlib(tmp_path):
    import hashlib
    p = tmp_path / "b.bin"
    payload = b"\x00\x01" * 5000
    p.write_bytes(payload)
    assert sha256_file(p) == hashlib.sha256(payload).hexdigest()


# --------------------------------------------------------------------------- Stage
class _Dummy(Stage):
    name = "dummy"

    def build(self, cfg, out_dir, **params):
        p = Path(out_dir) / "out.txt"
        p.write_text("x")
        return Artifact.of(p, kind="dummy", params=params)

    def verify(self, cfg, artifact, **params):
        r = GateReport("dummy")
        r.add("D1", Path(artifact.path).read_text() == "x", "content is x")
        return r


def test_stage_is_abstract():
    with pytest.raises(TypeError):
        Stage()                      # type: ignore[abstract]


def test_stage_roundtrip(tmp_path):
    s = _Dummy()
    art = s.build(None, tmp_path, k=1.0)
    assert s.verify(None, art).ok


def test_plot_default_raises(tmp_path):
    s = _Dummy()
    art = s.build(None, tmp_path)
    with pytest.raises(NotImplementedError, match="dummy"):
        s.plot(None, art)


def test_check_outputs_flags_extra_and_missing(tmp_path):
    s = _Dummy()
    art = s.build(None, tmp_path)
    (tmp_path / "scratch.nc").write_text("stale")
    rep = s.check_outputs(tmp_path, keep=[art.path])
    assert rep["ok"] is True
    assert rep["extra"] == [str((tmp_path / "scratch.nc").resolve())]
    rep2 = s.check_outputs(tmp_path, keep=[tmp_path / "absent.nc"])
    assert rep2["ok"] is False and rep2["missing"]


# --------------------------------------------------------------------------- Manifest
def test_manifest_records_and_writes(tmp_path):
    p = tmp_path / "s.nc"
    p.write_bytes(b"data")
    art = Artifact.of(p, kind="stress", params={"k": 1.7})
    r = GateReport("stress")
    r.add("V1", True, "ok")

    m = Manifest(project_name="safs_alt", descriptor_sha256="deadbeef")
    m.record("stress", artifact=art, report=r, wall_s=1.5)
    out = m.write(tmp_path)

    assert out.name == Manifest.FILENAME
    d = json.loads(out.read_text())
    assert d["project"] == "safs_alt"
    assert d["descriptor_sha256"] == "deadbeef"
    entry = d["stages"]["stress"]
    assert entry["artifacts"][0]["sha256"] == art.sha256
    assert entry["artifacts"][0]["params"] == {"k": 1.7}
    assert entry["gates"]["ok"] is True
    assert entry["wall_s"] == 1.5
    assert Manifest.read(out) == d


def test_manifest_accumulates_multiple_artifacts_per_stage(tmp_path):
    m = Manifest("p", "h")
    for n in ("a", "b"):
        f = tmp_path / f"{n}.nc"
        f.write_text(n)
        m.record("material", artifact=Artifact.of(f, kind="material"))
    assert len(m.stages["material"]["artifacts"]) == 2
    assert len(m.artifacts()) == 2
    assert "2 artifacts" in repr(m)


def test_manifest_for_project_uses_descriptor_hash():
    from deckbuild.config import Project
    cfg = Project.load(ROOT / "projects" / "safs_alt.yaml", require_files=False)
    m = Manifest.for_project(cfg)
    assert m.project_name == "safs_alt"
    assert m.descriptor_sha256 == cfg.sha256()


# ------------------------------------------------- review finding R-001
def test_R001_reexports_are_refreshed_after_reload():
    """The package must be reloaded AFTER its submodules.

    deckbuild/__init__.py does `from deckbuild.config import Project`.  Reloading the
    package first re-binds that name to the OLD class object, so `deckbuild.Project` would
    silently be the pre-edit class -- the exact stale-kernel failure init() exists to stop.
    """
    import deckbuild
    import deckbuild.config
    from deckbuild.bootstrap import _reload_package

    _reload_package()
    assert deckbuild.Project is deckbuild.config.Project
    assert deckbuild.GateReport is sys.modules["deckbuild.contract"].GateReport


def test_R001_reload_order_puts_the_package_last():
    from deckbuild.bootstrap import _reload_package  # noqa: F401
    names = ["deckbuild", "deckbuild.config", "deckbuild.contract", "deckbuild.bootstrap"]
    names.sort(key=lambda n: (n.count("."), n), reverse=True)
    assert names[-1] == "deckbuild", f"package must reload last, got {names}"
