# Code Review: Phase 2 — the stress stage (2026-08-01)

Phases 0 and 1 reviews are in the git history (`9fd3394`, `3866c65`); all twelve of their
findings were fixed.

## Review Scope
- Plan: `docs/PLAN_...md`, Phase 2
- Files: `deckbuild/stress.py`, `deckbuild/orientation.py`, `tests/test_stress.py`
- Baseline: 182/182 tests pass.

## Findings

### [R-201] MODERATE [stress.py:_apply_daylight_patch] — the daylight patch is a constant slab, not the per-column topographic overburden

**Category:** DEVIATION

**Description:**
The legacy `step6_write_nc._apply_daylight_patch` takes a `daylight_mesh` and reads the
BC-1 free surface to get `z_topo(x, y)`, then fills each column's `z >= 0` levels with that
column's own overburden. This implementation fills every `z >= 0` level with a **single
constant slab** evaluated at one depth.

For the SAFS PREFERRED fault — which daylights to +2209 m under +3 km of topography — those
are not the same field. The shipped `safs_stress_andersonian_*.nc` for PREFERRED was built
the legacy way, so Phase 8's E4 rung **cannot** reproduce it data-identically with this
implementation.

It is correct for a flat-top domain (ALT), which is why every current test passes.

**Trigger:** Build a PREFERRED-style project with `daylight_patch=True` and compare against
a shipped PREFERRED stress nc.

**Actual behavior:** One constant slab above `z = 0`.

**Expected behavior:** Per-column fill at the local topographic depth.

**Suggested fix:** Take the mesh, extract the free surface, and evaluate per column.
```diff
-def _apply_daylight_patch(comps, gz, depth_grid, sv_grid, az_col, R_col, k_col,
-                          nx, ny, min_depth_m):
+def _apply_daylight_patch(comps, gz, depth_grid, sv_grid, az_col, R_col, k_col,
+                          nx, ny, min_depth_m, colx=None, coly=None, topo=None):
```
where `topo(colx, coly) -> z_topo` comes from the mesh's BC-1 faces, and the per-column
depth is `max(min_depth_m, z_topo - z_level)`. Until that lands, the limitation must be
stated in the module docstring and the build must **refuse** `daylight_patch=True` when no
topography source is supplied, rather than silently writing the flat-top approximation.

**Test case:**
```python
def test_R201_daylight_patch_requires_a_topography_source(demo, sv, tmp_path):
    with pytest.raises(StressError, match="topograph"):
        StressStage().build(demo, tmp_path, sv_profile=sv, k=1.7,
                            daylight_patch=True)
```

---

### [R-202] MODERATE [stress.py:StressStage.verify] — a corrupt mesh is reported as a skipped gate

**Category:** BUG

**Description:**
```python
try:
    fault = load_fault(...)
except Exception as exc:
    rep.skip("V4", f"fault mesh unavailable ({type(exc).__name__})")
```
The bare `except Exception` cannot tell "no mesh file here" from "the mesh is corrupt", "the
BC map is wrong", or "there are no BC-3 faces". All four become a skip, and `report.ok`
stays True. The contract everywhere else in this package is that a skip is for a genuinely
absent optional input — a broken input must fail.

**Trigger:** `verify()` against a project whose mesh exists but has the wrong `fault_bc`.

**Actual behavior:** `[ -- ] V4: fault mesh unavailable (GeometryError)`, `ok == True`.

**Expected behavior:** absent mesh -> skip; present-but-unusable mesh -> HARD fail.

**Suggested fix:**
```diff
-        try:
-            fault = load_fault(cfg.mesh(), cfg.data_dir, strike=cfg.strike)
-        except Exception as exc:
-            rep.skip("V4", f"fault mesh unavailable ({type(exc).__name__}); "
-                           f"the on-fault screen did not run")
-            return rep
+        mesh_path = cfg.resolve_path(cfg.mesh(mesh_name).path)
+        if not Path(mesh_path).is_file():
+            rep.skip("V4", f"no mesh at {mesh_path}; the on-fault screen did not run")
+            return rep
+        try:
+            fault = load_fault(cfg.mesh(mesh_name), cfg.data_dir, strike=cfg.strike)
+        except Exception as exc:                     # present but unusable -> FAIL
+            rep.add("V4", False,
+                    f"mesh {mesh_path} exists but could not be read: "
+                    f"{type(exc).__name__}: {exc}")
+            return rep
```

**Test case:**
```python
def test_R202_corrupt_mesh_fails_v4_rather_than_skipping(demo, sv, tmp_path):
    bad = tmp_path / "data" / "demo_planar.puml.h5"
    bad.parent.mkdir(parents=True, exist_ok=True)
    bad.write_bytes(b"not an hdf5 file")
    art = StressStage().build(demo, tmp_path, sv_profile=sv, k=1.7,
                              freeze_above_depth_m=0.0)
    rep = StressStage().verify(demo, art)
    v4 = [g for g in rep.gates if g.name == "V4"][0]
    assert v4.severity == "hard" and not v4.passed
```

---

### [R-203] MODERATE [tests/test_stress.py] — the demo's parameters make `s_xx` blind to k, so a broken grading would pass

**Category:** EDGE_CASE

**Description:**
The demo project has `R = 0.5` and SHmax `az = 45`. For that combination the xx component
is **algebraically independent of k**:

```
sigma_xx = sig1 sin^2(45) + sig3 cos^2(45) = 0.5 (sig1 + sig3) = 0.5 sig3 (1 + k)
sig3     = Sv / ((1-R)k + R) = 2 Sv / (k + 1)      when R = 0.5
=> sigma_xx = Sv      for every k
```

Measured: a graded build with `k = 1.2 -> 2.4` gives `s_xx = -166.77 MPa` at every column.
So a test that probes `s_xx` for lateral variation proves nothing, and a regression that
silently dropped the grading entirely would still pass the current suite. The eigenvalue
test (`test_built_field_reproduces_the_closure_ratio`) does constrain k, but only for the
uniform case.

`s_xy` is not degenerate: `sigma_xy = Sv (k - 1) / (k + 1)`, which I verified matches the
build to 3 decimal places at both plateaus.

**Suggested fix:** Add an analytic grading test on `s_xy`.

**Test case:**
```python
def test_R203_graded_k_varies_the_field_analytically(demo, sv, tmp_path):
    """R=0.5, az=45 makes s_xx k-INVARIANT; s_xy is the component that moves."""
    d = KDesign(k_values=(1.2, 2.4), boundaries_s_km=((20.0, 30.0),))
    art = StressStage().build(demo, tmp_path, sv_profile=sv, design=d,
                              freeze_above_depth_m=0.0)
    _, _, z, f, _ = read_asagi(art.path)
    kz = int(np.argmin(np.abs(z - (-10000.0))))
    row = f["s_xy"][kz].mean(axis=1) / 1e6
    Sv = 1700.0 * 9.81 * 10000.0 / 1e6
    assert row[0] == pytest.approx(-Sv * (1.2 - 1) / (1.2 + 1), abs=0.05)
    assert row[-1] == pytest.approx(-Sv * (2.4 - 1) / (2.4 + 1), abs=0.05)
```

---

### [R-204] LOW [stress.py:_project_onto_fault] — an obfuscated dynamic import

**Category:** QUALITY

**Description:**
```python
_, _, zg = (np.asarray(a) for a in __import__(
    "deckbuild.asagi", fromlist=["asagi_axes"]).asagi_axes(nc_path))
```
`asagi_axes` is already imported at the top of the function's module neighbours; this
construct hides the dependency and is hard to read for no benefit.

**Suggested fix:**
```diff
-    from deckbuild.asagi import trilinear_sample
+    from deckbuild.asagi import asagi_axes, trilinear_sample
     from deckbuild.geometry import resolve_tractions
-    _, _, zg = (np.asarray(a) for a in __import__(
-        "deckbuild.asagi", fromlist=["asagi_axes"]).asagi_axes(nc_path))
+    _, _, zg = asagi_axes(nc_path)
```

---

### [R-205] LOW [stress.py:StressStage.verify] — V1 is always reported as skipped

**Category:** QUALITY

**Description:**
`verify()` unconditionally emits `rep.skip("V1", ...)`. The real V1 lives in
`verify_regions()`, which the tests exercise. That is a defensible split (V1 needs the Sv
profile and a scratch directory, which `verify` does not take), but a reader of a
`GateReport` sees `V1 skipped` with no indication that a real V1 ran elsewhere.

**Suggested fix:** Make the skip detail point at the method that does run it, and have
`verify_regions` return a report that can be merged, so a caller can present one battery.
The detail already names `verify_regions`; add a note that the merged form is the complete
V1–V4 and reference it from the stage docstring.

---

## Summary
- Critical: 0 | Moderate: 3 (R-201, R-202, R-203) | Low: 2 (R-204, R-205)
- Plan compliance: **PARTIAL** — the closure, grading, freeze and V1–V4 names are all in
  place and correct, but the daylight patch is a flat-top approximation of the legacy
  per-column treatment (R-201), which blocks PREFERRED reproduction in Phase 8.
- Verdict: **PASS WITH FIXES** — R-201 and R-202 must be fixed before Phase 8; R-203 before
  anyone trusts the grading tests.

## Checked and found correct
- `magnitudes_C1` and `build_tensor_andersonian` are verbatim ports; the built tensor's
  horizontal eigenvalue ratio equals k to 2e-3 (float32 storage), sigma2 is exactly
  vertical, and the eigenvalues are exactly (sig3, sig2, sig1).
- Array orientation: fields are written `(nz, ny, nx)` and a graded design varies along the
  correct axis, verified against the closed form `sigma_xy = Sv (k-1)/(k+1)`.
- The freeze touches only the shallow band: below the freeze depth a frozen and an unfrozen
  build are bit-identical, and `mu_app` is preserved inside the band.
- SHmax is interpolated as a doubled-angle unit vector, so the 179/1 wrap gives 0, not 90.
