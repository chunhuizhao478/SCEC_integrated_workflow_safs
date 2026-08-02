"""Phase 7 acceptance tests: the integrated notebook and the demo project."""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

NB = ROOT / "deck_workflow.ipynb"
PY_EXE = sys.executable


def _cells(nb=NB):
    return json.loads(nb.read_text())["cells"]


def _src(c):
    return "".join(c["source"])


# --------------------------------------------------------------------------- structure
def test_notebook_exists_and_is_valid_nbformat():
    import nbformat
    nbformat.validate(nbformat.read(str(NB), as_version=4))


def test_all_eight_sections_are_present():
    md = "\n".join(_src(c) for c in _cells() if c["cell_type"] == "markdown")
    for i, name in enumerate(["Setup and project", "What we are building", "Mesh",
                              "Material", "Stress", "Friction", "Assemble the deck",
                              "Design read-outs"]):
        assert f"## [{i}]" in md, f"section [{i}] missing"
        assert name in md, f"section [{i}] should mention {name!r}"


def test_every_parameters_cell_is_assignments_only():
    """A PARAMETERS cell is what the user edits; it must contain no logic."""
    import ast
    n = 0
    for c in _cells():
        s = _src(c)
        if "# ---- PARAMETERS" not in s:
            continue
        n += 1
        # everything between the PARAMETERS banner and the closing dashed rule
        after = s.split("# ---- PARAMETERS", 1)[1]
        lines = []
        for ln in after.splitlines()[1:]:
            if ln.startswith("# ---") and set(ln.strip("# -")) <= set(""):
                break
            if ln.lstrip().startswith("# ---"):
                break
            if ln.strip() and not ln.lstrip().startswith("#"):
                lines.append(ln)
        body = ast.parse("\n".join(lines))
        for node in body.body:
            assert isinstance(node, ast.Assign), \
                f"a PARAMETERS cell must be assignments only, found {type(node).__name__}"
    assert n >= 5, f"expected a PARAMETERS cell per editable section, found {n}"


def test_section_2_says_it_ingests_rather_than_builds():
    md = "\n".join(_src(c) for c in _cells() if c["cell_type"] == "markdown")
    seg = md.split("## [2]", 1)[1].split("## [3]", 1)[0]
    assert "does **not** build a mesh" in seg or "not** build a mesh" in seg
    assert "skills/" in seg and "MESHING.md" in seg


def test_material_runs_before_stress_and_friction():
    """The dependency order is load-bearing, not cosmetic."""
    md = "\n".join(_src(c) for c in _cells() if c["cell_type"] == "markdown")
    assert md.index("## [3]") < md.index("## [4]") < md.index("## [5]")
    assert "material runs first" in md.lower() or "Material** runs first" in md
    intro = md.split("## [0]", 1)[0]
    assert "Order is not cosmetic" in intro


def test_the_gotchas_a_user_would_otherwise_hit_are_stated():
    all_src = "\n".join(_src(c) for c in _cells())
    for needle in ("silently DROPS a receiver", "silently ignores it",
                   "COEFFICIENT", "necessary, not sufficient", "sea-level",
                   "strike-slip assumption"):
        assert needle.lower() in all_src.lower(), f"missing: {needle}"


def test_no_absolute_paths_in_the_notebook():
    all_src = "\n".join(_src(c) for c in _cells())
    assert "/Users/" not in all_src and "/home/" not in all_src


# --------------------------------------------------------------------------- execution
@pytest.mark.slow
def test_notebook_runs_end_to_end_under_five_minutes(tmp_path):
    """Acceptance: Run-All on the demo, zero downloads, under 5 minutes."""
    env = dict(os.environ, MPLBACKEND="Agg")
    t0 = time.time()
    r = subprocess.run(
        [PY_EXE, "-m", "nbconvert", "--to", "notebook", "--execute",
         "--ExecutePreprocessor.timeout=600",
         "--output", str(tmp_path / "out.ipynb"), str(NB)],
        cwd=ROOT, capture_output=True, text=True, env=env)
    dt = time.time() - t0
    assert r.returncode == 0, r.stderr[-4000:]
    assert dt < 300, f"took {dt:.0f}s, budget is 300s"

    out = json.loads((tmp_path / "out.ipynb").read_text())
    text = "\n".join(
        "".join(o.get("text", "")) for c in out["cells"] for o in c.get("outputs", [])
        if o.get("output_type") == "stream")
    # every gate battery printed, and none of them failed
    for battery in ("M1", "V2", "G2a", "F0", "P1", "P8"):
        assert battery in text, f"gate {battery} never printed"
    assert "[FAIL]" not in text, "a gate failed during the notebook run"
    assert text.count("PASS:") >= 5, "expected a verdict line per stage"

    # ...and it must actually have WRITTEN a deck, not just printed about one.
    decks = sorted((ROOT / "decks").glob("*_nb"))
    assert decks, "the notebook printed gates but produced no deck"
    d = decks[-1]
    assert (d / "deck_manifest.json").is_file()
    assert list(d.glob("*.puml.h5")) and list(d.glob("*stress*.nc"))
    assert (d / "parameters.par").is_file()


@pytest.mark.slow
def test_notebook_needs_no_downloads():
    """The demo must be fully self-contained: no http, no fetch, no download_data."""
    all_src = "\n".join(_src(c) for c in _cells())
    for bad in ("http://", "https://", "urlretrieve", "download_data", "wget", "curl "):
        assert bad not in all_src, f"the demo must not need {bad}"


def test_demo_inputs_are_committed_or_synthetic():
    """Zero downloads means the mesh is in the repo and the fields are synthetic."""
    from deckbuild.config import Project
    cfg = Project.load(ROOT / "projects" / "demo_planar.yaml", require_files=True)
    assert Path(cfg.resolve_path(cfg.mesh().path)).is_file()
    assert cfg.raw.velocity.kind == "layered_1d"      # no community model
    assert cfg.raw.orientation.kind == "constant"
    assert cfg.raw.thermal.kind == "depth_profile"


def test_readme_retarget_checklist_covers_every_required_field():
    from deckbuild.config import Project
    readme = (ROOT / "README.md").read_text().lower()
    for field in ("mesh", "velocity", "orientation", "friction", "hypocentre",
                  "strike", "crs"):
        assert field in readme, f"retarget checklist missing {field}"
    assert "tectonic regime" in readme
