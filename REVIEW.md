# Code Review: Phase 1 — geometry and ASAGI I/O (2026-08-01)

Phase 0's review is in the git history (commit `9fd3394`); all seven of its findings were
fixed before this phase started.

## Review Scope
- Plan: `docs/PLAN_integrated_workflow_notebook_2026-08-01.md`, Phase 1
- Files reviewed: `deckbuild/geometry.py`, `deckbuild/asagi.py`, `tests/nc_compare.py`,
  `tests/test_geometry.py`, `tests/test_asagi.py`, `deckbuild/__init__.py`
- Domain context: the legacy `generate_stress_nc_from_raw.py` (the source of every ported
  expression), and the exploration doc's convention and gotcha sections.
- Baseline: 139/139 tests pass, including three 0-ULP comparisons against the real legacy
  module. The findings below are not caught by that suite.

## Findings

### [R-101] MODERATE [geometry.py:Fault] — comparing two Fault objects raises ValueError

**Category:** BUG

**Description:**
`Fault` is `@dataclass(frozen=True)` holding numpy arrays. The generated `__eq__` compares
field tuples, which for arrays yields elementwise arrays and then an ambiguous truth value.
Verified:

```
a, b = mk(), mk()          # two distinct, element-wise identical Faults
bool(a == b)  ->  ValueError: The truth value of an array with more than one element is ambiguous
hash(a)       ->  TypeError: unhashable type: 'numpy.ndarray'
```

`a == a` happens to work (tuple comparison short-circuits on identity), which is why no
current test caught it. Any code that compares two independently loaded faults — an obvious
thing for Phase 5's "fault identity vs the deployed mesh" check, and for Phase 8's
reproduction diffing — hits this.

**Trigger:** `load_fault(m) == load_fault(m)`.

**Actual behavior:** `ValueError`.

**Expected behavior:** Either a correct elementwise-aware comparison, or no `__eq__` at all
(identity), which is the honest choice for a bulk array container.

**Suggested fix:** Disable the generated `__eq__`/`__hash__` and provide an explicit,
named comparison so the intent is unambiguous at the call site.
```diff
-@dataclass(frozen=True)
+@dataclass(frozen=True, eq=False)
 class Fault:
```
and add:
```python
    def same_geometry_as(self, other: "Fault", tol: float = 0.0) -> bool:
        """True when both faults have identical facet centroids and normals."""
        if len(self) != len(other):
            return False
        if tol == 0.0:
            return (np.array_equal(self.cent, other.cent)
                    and np.array_equal(self.normals, other.normals))
        return (np.allclose(self.cent, other.cent, rtol=0, atol=tol)
                and np.allclose(self.normals, other.normals, rtol=0, atol=tol))
```

**Test case:**
```python
def test_R101_fault_comparison_does_not_raise(planar_mesh):
    a = load_fault(planar_mesh, strike=SAFS_STRIKE)
    b = load_fault(planar_mesh, strike=SAFS_STRIKE)
    assert a is not b
    assert (a == b) is False          # identity semantics, no exception
    assert a.same_geometry_as(b)
```

---

### [R-102] MODERATE [asagi.py:roundtrip_selfcheck] — the max-error gate can never fail, deviating from the legacy G3 contract

**Category:** DEVIATION

**Description:**
The legacy G3 guard is: median within `1e-5` (hard); max is **WARN-only at `5e-4`
PROVIDED every excursion sits in a grid cell whose corners straddle a kink in the source
field**. A `>5e-4` excursion in a kink-free cell **is a hard fail** — the exploration doc
records this explicitly as "a real bug".

The implementation drops the proviso and marks the max gate `WARN` unconditionally, so the
hard case is unreachable. A genuine interpolation bug in a smooth region would report as a
warning and pass.

**Trigger:** Any grid whose stored values disagree with the evaluator by more than
`max_tol` in a smooth region.

**Actual behavior:** `[WARN]`, and `report.ok` is True.

**Expected behavior:** hard fail unless the caller can attest the excursion is kink-related.

**Suggested fix:** Take an optional predicate; when it is not supplied, the max gate is
HARD (the safe default), and the caller opts into leniency explicitly.
```diff
 def roundtrip_selfcheck(path, evaluator, field, n=1000, seed=12345,
                         median_tol=1e-5, max_tol=5e-4,
-                        gate="G3", report=None, bbox=None):
+                        gate="G3", report=None, bbox=None,
+                        kink_straddling=None):
```
```diff
-    rep.add(f"{gate}max", mx <= max_tol,
-            f"...",
-            severity=WARN)
+    if mx <= max_tol:
+        rep.add(f"{gate}max", True, f"{field}: max |nc - direct| = {mx:.3e} "
+                f"(tol {max_tol:.0e})")
+    else:
+        bad = np.flatnonzero(err > max_tol)
+        excused = (kink_straddling is not None
+                   and bool(np.all(kink_straddling(pts[bad]))))
+        rep.add(f"{gate}max", excused,
+                f"{field}: max |nc - direct| = {mx:.3e} > tol {max_tol:.0e} at "
+                f"{bad.size} point(s); "
+                + ("all sit in cells straddling a kink in the source field, which is "
+                   "expected for a piecewise-linear field"
+                   if excused else
+                   "no kink_straddling predicate was supplied, so these are treated as "
+                   "real interpolation errors"),
+                severity=WARN if excused else HARD)
```

**Test case:**
```python
def test_R102_large_smooth_excursion_is_a_hard_fail(nc):
    p = nc[0]
    rep = roundtrip_selfcheck(p, lambda pts: np.zeros(len(pts)), field="f",
                              n=50, median_tol=1e9, max_tol=1e-6)
    assert not rep.ok, "a large excursion with no kink excuse must fail HARD"

def test_R102_kink_excuse_downgrades_to_warn(nc):
    p = nc[0]
    rep = roundtrip_selfcheck(p, lambda pts: np.zeros(len(pts)), field="f",
                              n=50, median_tol=1e9, max_tol=1e-6,
                              kink_straddling=lambda pts: np.ones(len(pts), bool))
    assert rep.ok
```

---

### [R-103] MODERATE [tests/test_asagi.py] — acceptance criterion 5 is not tested against a shipped nc

**Category:** DEVIATION

**Description:**
Phase 1's fifth acceptance criterion is: *"`trilinear_sample` on the shipped
`safs_stress_andersonian_k1.7.nc` matches the legacy
`generate_stress_nc_from_raw.trilinear_sample` to 0 ULP at 1000 fixed-seed points."*

The suite tests trilinear behaviour only on synthetic grids. I confirmed 0-ULP agreement
against the legacy function manually on a synthetic compound nc, but the criterion names a
real shipped file, whose axes are non-trivial (69 z-levels, 250 m spacing, a daylight
extension) and whose values are float32 at ~1e7 Pa — a regime the synthetic test does not
cover.

**Trigger:** n/a — a missing test.

**Suggested fix:** Add a skipif-gated test that uses a shipped deck when present.
```python
DECK = Path(os.environ.get("DECKBUILD_DECKS", Path.home() / "Downloads/seisol_quakeworx"))
SHIPPED_STRESS = next(DECK.glob("*/safs_stress_andersonian_k1.7.nc"), None) if DECK.is_dir() else None

@pytest.mark.skipif(SHIPPED_STRESS is None, reason="shipped deck not available")
def test_R103_trilinear_matches_legacy_on_the_shipped_nc():
    legacy = _load_legacy()
    x, y, z = asagi_axes(SHIPPED_STRESS)
    rng = np.random.default_rng(12345)
    qx = rng.uniform(x[0], x[-1], 1000)
    qy = rng.uniform(y[0], y[-1], 1000)
    qz = rng.uniform(z[0], z[-1], 1000)
    mine = trilinear_sample(SHIPPED_STRESS, qx, qy, qz)
    theirs = legacy.trilinear_sample(str(SHIPPED_STRESS), qx, qy, qz)
    for f in theirs:
        assert np.array_equal(mine[f], theirs[f]), f
```

---

### [R-104] LOW [geometry.py:project_point] — documented deviation from the plan's topo-relative depth

**Category:** DEVIATION

**Description:**
The plan's `Hypocenter.depth_m` comment says "positive down; converted to z with the local
topo". `project_point` uses `z = -depth_m`, i.e. sea-level referencing, and says so in its
docstring.

This is the right call — the local topographic elevation lives on the mesh free surface,
which this function does not receive, and the SAFS descriptors are already sea-level
referenced. But the deviation is currently only visible to someone reading the function.
Anyone giving a catalogue depth for a fault under significant topography will be off by the
local elevation (up to ~3 km for the PREFERRED domain).

**Suggested fix:** Surface it where a user will see it — in the `Hypocenter` docstring in
`config.py`, and in the snap `GateReport` detail when the geographic path is used.
```diff
     depth_m: float | None = None
```
```diff
-    Give EITHER projected coords (x, y, z) OR geographic (lon, lat, depth_m); the loader
-    converts via `Project.crs`.
+    Give EITHER projected coords (x, y, z) OR geographic (lon, lat, depth_m); the loader
+    converts via `Project.crs`.  `depth_m` is referenced to SEA LEVEL (z = -depth_m), not
+    to the local ground surface -- under topography those differ by the local elevation.
```

---

### [R-105] LOW [asagi.py:write_asagi] — an undocumented extra check that could reject a legitimate field

**Category:** ASSUMPTION

**Description:**
`write_asagi` rejects any field containing a non-finite value. The plan does not ask for
this. It is defensive and almost certainly right (a NaN in an ASAGI grid propagates into
every element that samples it), but it is an added constraint: a stage that legitimately
wants a NaN sentinel — for instance a masked region outside a model's coverage — cannot
write one, and will discover this only at write time.

**Suggested fix:** Keep the check (it is the safer default) but make it opt-out, so the
constraint is a decision rather than an accident.
```diff
 def write_asagi(path, x, y, z, fields, attrs=None, dtype=np.float32):
+                # allow_nonfinite: ASAGI propagates NaN into every element that samples
+                # the grid, so this defaults to False.
```
Add `allow_nonfinite: bool = False` and gate the check on it.

---

## Summary
- Critical issues: 0
- Moderate issues: 3  (R-101, R-102, R-103)
- Low issues: 2       (R-104, R-105)
- Plan compliance: **PARTIAL** — every Phase 1 interface exists and four of five acceptance
  criteria are tested and pass; criterion 5 (0 ULP on the *shipped* nc) is untested
  (R-103), and the G3 max-error semantics deviate from the legacy contract (R-102).
- Verdict: **PASS WITH FIXES** — none of these block Phase 2 from starting, but R-102
  should be fixed before any stage relies on `roundtrip_selfcheck` as a real gate, which
  Phases 3 and 4 both will.

## Notes on things that were checked and found correct
- `strike_s_km` is 0 ULP against the legacy `strike_distance_km` over 10,000 random points,
  and reproduces the two Lua coefficients that appear in every shipped `rs_muw` map.
- `load_fault` reproduces the legacy centroids, normals, areas, strikes and dips exactly on
  a synthetic mesh.
- `trilinear_sample` is 0 ULP against the legacy function on a synthetic compound nc
  (see R-103 for the shipped-file gap).
- The `_reload_package` fix from Phase 0 holds, and the test-isolation problem it caused
  (reload replaces class objects, breaking `except` for other modules' bound references)
  is fixed by running that probe in a subprocess and is documented in `bootstrap.py`.

## Unreviewed Areas
- The Colab branch of `bootstrap.init()` — still unexercisable on this machine.
- `nc_compare._axes_and_fields` reads non-compound variables into the field set. Correct
  for comparison, but no shipped file exercises that path yet.
