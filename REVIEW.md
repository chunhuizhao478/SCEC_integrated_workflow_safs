# Code Review: Phase 0 — scaffold, descriptor, stage contract (2026-08-01)

## Review Scope
- Plan: `docs/PLAN_integrated_workflow_notebook_2026-08-01.md`, Phase 0
- Files reviewed: `deckbuild/{__init__,config,contract,bootstrap}.py`,
  `projects/{safs_alt,safs_preferred,demo_planar}.yaml`,
  `tests/{test_config,test_contract,test_no_hardcoded_constants}.py`,
  `.gitignore`, `environment.yml`, `README.md`
- Domain context: `docs/EXPLORE_deck_components_2026-08-01.md` section 10 (the constant table),
  and the Phase 0 acceptance criteria.
- Baseline: 71/71 tests pass. The bugs below are **not** caught by the current suite — that
  is itself a finding about test coverage.

## Findings

### [R-001] CRITICAL [bootstrap.py:_reload_package] — reload order is inverted, so the package re-exports stale objects

**Category:** BUG

**Description:**
The sort intends "submodules before the package" (its own comment says so), but
`reverse=True` combined with the key `(n == "deckbuild", n)` sorts `True` first, putting the
**package before its submodules**. Verified:

```
actual : ['deckbuild', 'deckbuild.contract', 'deckbuild.config', 'deckbuild.bootstrap']
wanted : submodules first, 'deckbuild' last
```

`deckbuild/__init__.py` does `from deckbuild.config import Project, ...`. Reloading the
package first re-binds those names to the **old** module objects; the submodules are then
reloaded, leaving `deckbuild.Project` pointing at a stale class while
`deckbuild.config.Project` is fresh. This defeats the single purpose of the function — the
plan's requirement 6 calls out "a stale Jupyter kernel silently running old code is the
single most common way to 'fix' something and see nothing change."

**Trigger:** Edit `deckbuild/config.py` in a live kernel, re-run `init()`, then use the
top-level re-export `deckbuild.Project`.

**Actual behavior:** `deckbuild.Project is not deckbuild.config.Project`; the top-level name
is the pre-edit class.

**Expected behavior:** After `init()`, every `deckbuild.*` name refers to the reloaded object.

**Suggested fix:**
```diff
-    # Reload submodules before the package itself so re-exports pick up new objects.
-    names.sort(key=lambda n: (n == "deckbuild", n), reverse=True)
+    # Reload submodules BEFORE the package itself, so the package's `from .config import ...`
+    # re-exports bind to the freshly reloaded objects.  Deepest module first.
+    names.sort(key=lambda n: (n.count("."), n), reverse=True)
```

**Test case:**
```python
def test_R001_reexports_are_refreshed_after_reload(tmp_path, monkeypatch):
    import deckbuild, deckbuild.config
    from deckbuild.bootstrap import _reload_package
    _reload_package()
    assert deckbuild.Project is deckbuild.config.Project, (
        "top-level re-export is stale: the package was reloaded before its submodules")
```

---

### [R-002] MODERATE [.gitignore] — `outputs/.gitkeep` is ignored, so the directory the plan requires is never committed

**Category:** DEVIATION

**Description:**
Phase 0's "Files to Create" list ends with `outputs/.gitkeep`. The file was created, but
`.gitignore` line 9 is `outputs/`, which ignores the directory **and everything in it**.
Verified: `git check-ignore -v outputs/.gitkeep` → `.gitignore:9:outputs/`, and
`git status --porcelain outputs/.gitkeep` returns nothing.

A fresh clone therefore has no `outputs/`. `bootstrap.init()` happens to `mkdir` it, but the
plan asked for the committed placeholder and anything that reads `outputs/` before `init()`
runs will fail.

**Trigger:** `git clone` the repo; `ls outputs`.

**Actual behavior:** No `outputs/` directory exists after clone.

**Expected behavior:** `outputs/.gitkeep` is tracked; the directory exists on clone; its
build products remain ignored.

**Suggested fix:**
```diff
 # workflow artifacts -- rebuilt, never committed
 outputs/
+!outputs/.gitkeep
 decks/
```
Then `git add -f outputs/.gitkeep`.

**Test case:**
```python
def test_R002_outputs_placeholder_is_tracked():
    import subprocess
    from pathlib import Path
    root = Path(__file__).resolve().parent.parent
    out = subprocess.run(["git", "check-ignore", "outputs/.gitkeep"],
                         cwd=root, capture_output=True)
    assert out.returncode != 0, "outputs/.gitkeep is gitignored; it must be committable"
    assert (root / "outputs" / ".gitkeep").is_file()
```

---

### [R-003] MODERATE [config.py:Project.validate] — hypocenter XOR check accepts a silently-ignored partial coordinate

**Category:** EDGE_CASE

**Description:**
The check is `if hypo.is_projected == hypo.is_geographic: raise`. `is_projected` requires all
of x/y/z; `is_geographic` requires all of lon/lat/depth_m. A descriptor giving a **complete
projected point plus a stray `lon:`** satisfies `is_projected=True, is_geographic=False`, so
validation passes and the `lon` is silently discarded.

This is precisely the class of error the strict loader exists to prevent: a user editing a
SAFS descriptor toward their own fault will typically fill in `lon`/`lat` and forget to
delete `x`/`y`/`z`. They then get a run centred on the **San Andreas** hypocentre with no
warning. Phase 1's `snap_hypocenter` would not catch it either — the projected point is on
the SAFS fault, so the snap distance is ~0.

**Trigger:**
```yaml
hypocenters:
  alt: {x: 604446.944, y: 3704576.3853, z: -10067.9819, lon: -121.0}
```

**Actual behavior:** Loads cleanly; `lon` ignored.

**Expected behavior:** `ConfigError` naming the mixed keys.

**Suggested fix:**
```diff
-            if hypo.is_projected == hypo.is_geographic:
+            proj_set = [k for k in ("x", "y", "z") if getattr(hypo, k) is not None]
+            geo_set = [k for k in ("lon", "lat", "depth_m") if getattr(hypo, k) is not None]
+            if proj_set and geo_set:
+                raise ConfigError(
+                    f"{where}: hypocenters.{mname} mixes projected {proj_set} with "
+                    f"geographic {geo_set}; give one coordinate system only")
+            if hypo.is_projected == hypo.is_geographic:
                 raise ConfigError(
                     f"{where}: hypocenters.{mname} must give EITHER a complete projected "
                     f"(x, y, z) OR a complete geographic (lon, lat, depth_m) point, "
                     f"not both and not neither")
```

**Test case:**
```python
def test_R003_partial_coordinate_mixing_is_rejected(tmp_path, safs_raw):
    safs_raw["hypocenters"]["alt"] = {
        "x": 604446.944, "y": 3704576.3853, "z": -10067.9819, "lon": -121.0}
    with pytest.raises(ConfigError, match="mixes projected"):
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)
```

---

### [R-004] MODERATE [config.py] — mutable dicts inside frozen dataclasses defeat immutability and break hashing

**Category:** ASSUMPTION

**Description:**
`MeshSpec.tag_to_bc`, `SourceSpec.params` and `Project.provenance` are plain `dict` fields on
`@dataclass(frozen=True)`. `frozen=True` blocks attribute *rebinding* only — the dict
contents are freely mutable, so `cfg.meshes["alt"].tag_to_bc[999] = 1` silently rewrites the
shared config. `Project.meshes` / `hypocenters` are likewise mutable dicts.

Separately, `frozen=True` generates `__hash__` from the fields, so `hash(mesh_spec)` raises
`TypeError: unhashable type: 'dict'`. Nothing hashes these today, but the manifest/provenance
work in Phase 6 plausibly will (e.g. deduplicating specs in a set).

The descriptor is meant to be the single immutable source of truth for a run; a stage that
mutates it in place would produce artifacts whose recorded provenance no longer matches what
was used.

**Trigger:** `cfg.meshes["alt"].tag_to_bc.clear()` then `cfg.validate()` — passes silently on
a now-empty map, because validation already ran at load.

**Actual behavior:** Mutation succeeds and is invisible.

**Expected behavior:** Mapping fields are read-only.

**Suggested fix:** Freeze the mappings at construction with `MappingProxyType`.
```diff
+from types import MappingProxyType
...
     if dc is MeshSpec and "tag_to_bc" in kwargs and kwargs["tag_to_bc"] is not None:
         try:
-            kwargs["tag_to_bc"] = {int(k): int(v) for k, v in kwargs["tag_to_bc"].items()}
+            kwargs["tag_to_bc"] = MappingProxyType(
+                {int(k): int(v) for k, v in kwargs["tag_to_bc"].items()})
```
and wrap `params`, `provenance`, `meshes`, `hypocenters` the same way in `_build_source` /
`from_dict`. Note `_default_tag_to_bc()` must return a proxy too, and `as_dict` already
handles `Mapping` via its `dict` branch — confirm it still does (`MappingProxyType` is not a
`dict` subclass, so `as_dict` needs `isinstance(obj, Mapping)`).

**Test case:**
```python
def test_R004_descriptor_mappings_are_read_only():
    cfg = Project.load(PROJECTS / "safs_alt.yaml", require_files=False)
    with pytest.raises(TypeError):
        cfg.meshes["alt"].tag_to_bc[999] = 1
    with pytest.raises(TypeError):
        cfg.raw.velocity.params["grid_dx"] = 1.0
```

---

### [R-005] LOW [config.py] — `__import__("re")` inline instead of a module-level import

**Category:** QUALITY

**Description:**
`_SCI_NO_SIGN = __import__("re").compile(...)` sits in the middle of the module. It works,
but it hides a dependency from the import block, defeats static analysis, and is the kind of
thing a reader stops on. Every other import in the package is conventional.

**Trigger:** n/a (readability).

**Suggested fix:**
```diff
 import hashlib
+import re
 from dataclasses import dataclass, field, fields, is_dataclass
...
-_SCI_NO_SIGN = __import__("re").compile(r"^[+-]?\d+(\.\d*)?[eE]\d+$")
+_SCI_NO_SIGN = re.compile(r"^[+-]?\d+(\.\d*)?[eE]\d+$")
```

---

### [R-006] LOW [config.py:_numeric] — the "did you mean" hint misses unsigned exponents with no integer part

**Category:** EDGE_CASE

**Description:**
`_SCI_NO_SIGN` is `^[+-]?\d+(\.\d*)?[eE]\d+$`, which requires at least one digit before the
decimal point. YAML `.5e9` also parses as a string, but takes the generic
"must be a number" path with no fix suggestion. The error still fires (so this is not a
correctness bug), but the most useful part of the message is missing for a plausible input.

**Suggested fix:**
```diff
-_SCI_NO_SIGN = re.compile(r"^[+-]?\d+(\.\d*)?[eE]\d+$")
+_SCI_NO_SIGN = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)[eE]\d+$")
```

**Test case:**
```python
def test_R006_hint_covers_leading_dot(tmp_path, safs_raw):
    safs_raw["physics"]["rs_b"] = ".19e-1".replace("-", "")   # '.19e1' -> a YAML string
    with pytest.raises(ConfigError, match=r"e\+"):
        Project.load(write_tmp(tmp_path, safs_raw), require_files=False)
```

---

### [R-007] LOW [environment.yml] — pins `pytest=7.*` but the suite was only ever run on pytest 9

**Category:** DEVIATION

**Description:**
`environment.yml` pins `pytest=7.*`. The 71 passing tests were run under the `pythonenv`
environment, which has pytest 9.0.2 and Python 3.13 — neither matches this file
(`python=3.11`). Nobody has executed the suite in the environment this file describes, so
the pin is an untested claim. Phase 8 depends on `scipy` being pinned for reproducibility,
so this file needs to be real, not aspirational.

**Suggested fix:** Either relax the pin and record what was tested, or create the env and
run the suite in it before claiming the pin.
```diff
-  - pytest=7.*
+  - pytest>=7,<10          # validated on 9.0.2
```
and add to README: "validated on Python 3.13 / pytest 9.0.2 via the `pythonenv` conda env;
`environment.yml` describes the intended pinned environment and has not yet been built."

---

## Summary
- Critical issues: 1  (R-001)
- Moderate issues: 3  (R-002, R-003, R-004)
- Low issues: 3       (R-005, R-006, R-007)
- Plan compliance: **PARTIAL** — every listed file exists and all five acceptance criteria
  have tests that pass, but `outputs/.gitkeep` is not actually committable (R-002) and the
  reload guarantee in requirement 6 does not hold (R-001).
- Verdict: **PASS WITH FIXES** — R-001 and R-002 must be fixed before Phase 1, since Phase 1
  onward is developed in live notebooks that depend on the reload behaviour.

## Unreviewed Areas
- `docs/` — copied verbatim from the planning session, not re-reviewed here.
- `skills/SKILL.md` — staged only; Phase 5 vendors it properly with the two required edits.
- The Colab branch of `bootstrap.init()` (`_colab_setup`, `_find_on_drive`) cannot be
  exercised on this machine and has no test. Flagged for Phase 7, which is when a Colab run
  first matters.
