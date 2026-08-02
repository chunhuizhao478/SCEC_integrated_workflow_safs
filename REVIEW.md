# Code Review: Phases 3 + 4 — friction and material (2026-08-01)

Reviews for Phases 0–2 are in the git history (`9fd3394`, `3866c65`, `7e09f26`); all
seventeen of their findings were fixed.

## Review Scope
- Plan: `docs/PLAN_...md`, Phases 3 and 4
- Files: `deckbuild/friction.py`, `deckbuild/material.py`, `run_workflow.py`,
  `tests/test_friction_material.py`
- Baseline: 219 tests pass and `run_workflow.py --project demo_planar` completes with
  13 hard gates passing. The findings below are not caught by either.

## Findings

### [R-301] CRITICAL [material.py:sv_profile_from_material] — the Sv profile's depth axis is grid-top referenced, but every consumer reads it as sea-level depth

**Category:** BUG

**Description:**
`sv_profile_from_material` builds `depth = ztop[0] - ztop`, i.e. depth **measured from the
top of the material grid**. `StressStage._assemble` then asks for `sv_eff_at(depth=-z)` —
depth **from sea level**. The two agree only when the material grid's `z_max` is exactly 0.

It is not. The whole point of `extend_z_top` is to push the material grid above sea level
so it covers topography; the shipped SAFS CVM tops out at **+3250 m**. For any such grid
the stress field is built from an Sv profile shifted by `z_max`.

Measured on a demo grid with `z_max = 500`:

```
stress asks Sv at depth=    0.0 m  ->    0.000 MPa   (should be ~ 0)
stress asks Sv at depth=10000.0 m  ->  151.656 MPa   (should be ~157.6)
```

That is ~4 % low at 10 km and proportionally much worse near the surface — where the
shallow-freeze band lives and where `sigma_n` decides whether the trace runs away. The
whole stress field is wrong, silently, with every gate passing: V2 only checks the
eigenvalue *ratio* `k`, which is invariant under an overall scale of `Sv`.

**Trigger:** Any project whose velocity grid has `z_max != 0` — which is every real one.

**Actual behavior:** `Sv(depth)` is offset by `z_max`.

**Expected behavior:** The profile is sea-level referenced (`z = 0` -> `depth = 0`), or it
declares its reference and the consumer honours it.

**Suggested fix:** Reference to sea level explicitly, and record it in the npz so a
consumer cannot guess wrong.
```diff
-    depth = ztop[0] - ztop                                      # 0 at the grid top
+    # SEA-LEVEL referenced: depth = -z, so depth 0 is z = 0, matching what the stress
+    # closure asks for.  Referencing to the grid top would offset the whole field by
+    # z_max, which is +3250 m for the shipped SAFS CVM.
+    depth = -ztop
```
and clamp the integration to start at the shallowest node while keeping the mapping, plus:
```diff
-    np.savez(spath, depth_m=depth, sv_eff_mpa=sv_eff)
+    np.savez(spath, depth_m=depth, sv_eff_mpa=sv_eff, reference="sea_level")
```
with `_sv_arrays` asserting `reference == "sea_level"` when present.

**Test case:**
```python
def test_R301_sv_profile_is_sea_level_referenced(tmp_path):
    """depth 0 must mean z = 0, whatever the material grid's z_max is."""
    raw = yaml.safe_load((PROJECTS / "demo_planar.yaml").read_text())
    raw["raw"]["velocity"]["params"]["z_max"] = 500.0
    p = tmp_path / "t.yaml"; p.write_text(yaml.safe_dump(raw))
    cfg = Project.load(p, require_files=False, data_dir=tmp_path)
    out = MaterialStage().build(cfg, tmp_path)
    d = np.load(out.sv_profile.path)
    i = int(np.argmin(np.abs(d["depth_m"] - 0.0)))
    assert d["sv_eff_mpa"][i] == pytest.approx(0.0, abs=1e-9)
    j = int(np.argmin(np.abs(d["depth_m"] - 10000.0)))
    assert d["sv_eff_mpa"][j] == pytest.approx(157.6, rel=0.02)
```

---

### [R-302] MODERATE [friction.py:build] — the depth-profile path silently borrows the STRESS grid

**Category:** ASSUMPTION

**Description:**
```python
gx, gy, gz = build_grid(cfg.stress_box)
```
The friction nc is written on `cfg.stress_box`. The four grids are deliberately
independent — the exploration doc measures the shipped ones at
`material 1500/250`, `stress 1000/250`, `friction 1500/200` — and the descriptor has no
friction grid, so this path picks the stress one by default without saying so.

It is not wrong today (any ASAGI grid that contains the fault works), but it means a user
who tunes the stress grid silently re-grids their friction field, and it makes the friction
nc's resolution un-settable. Phase 8 needs `dz = 200` for friction and `250` for stress,
which this cannot express.

**Suggested fix:** Read the grid from the thermal source's params, falling back to the
stress box with a printed note.
```diff
-            gx, gy, gz = build_grid(cfg.stress_box)
+            box = cfg.stress_box
+            gdx = float(p.get("grid_dx", box.dx)); gdz = float(p.get("dz", box.dz))
+            zlo = float(p.get("z_min", box.zmin)); zhi = float(p.get("z_max", box.zmax))
+            gx = np.arange(box.xmin, box.xmax + 0.5 * gdx, gdx)
+            gy = np.arange(box.ymin, box.ymax + 0.5 * gdx, gdx)
+            gz = np.arange(zlo, zhi + 0.5 * gdz, gdz)
```

**Test case:**
```python
def test_R302_friction_grid_is_settable(demo_with_friction_dz_200, tmp_path):
    art = FrictionStage().build(cfg, tmp_path)
    _, _, z, _, _ = read_asagi(art.path)
    assert np.diff(z)[0] == pytest.approx(200.0)
```

---

### [R-303] LOW [friction.py:_eval_lua_text] — the FL gate's regex can silently find no transitions

**Category:** EDGE_CASE

**Description:**
`_eval_lua_text` re-parses the emitted Lua with a multi-line regex. If the emitter's
formatting changes, the regex matches nothing and the function returns the base value.
For a graded design that makes FL **fail** (fail-safe, which is right), but for a UNIFORM
design it would pass while proving nothing.

**Suggested fix:** Assert the parse found exactly the expected number of transitions.
```diff
+    found = pat.findall(text)
+    n_want = text.count("-- transition ")
+    if len(found) != n_want:
+        raise FrictionError(
+            f"FL: parsed {len(found)} transition(s) from the emitted Lua but the text "
+            f"declares {n_want}; the emitter and this parser have drifted")
-    for lo, hi, _lo2, hv, av in pat.findall(text):
+    for lo, hi, _lo2, hv, av in found:
```

---

### [R-304] LOW [material.py] — `M4` is always skipped, so the round-trip guard never runs

**Category:** DEVIATION

**Description:**
The plan's `M4` is the fixed-seed round-trip self-check, the material's equivalent of `G3`.
`verify` unconditionally skips it. The machinery exists (`asagi.roundtrip_selfcheck`) and
the `layered_1d` reader has a trivially available evaluator (the 1-D interpolant), so this
one could actually run for the implemented path.

**Suggested fix:** Run it when the source is `layered_1d`, skip only for readers with no
cheap evaluator.

---

## Summary
- Critical: 1 (R-301) | Moderate: 1 (R-302) | Low: 2 (R-303, R-304)
- Plan compliance: **PARTIAL** — the friction profiles, graded f_w, Lua emission, plasticity
  and Sv derivation are all in place and gated; the `cvm_slices` and `ctm_slices` raw
  readers raise "not implemented" rather than being silently wrong, which is the right
  interim state but leaves Phase 8's E3 rung blocked.
- Verdict: **FAIL — must fix before proceeding.** R-301 makes every stress field built
  through the real pipeline quantitatively wrong, and no existing gate detects it.

## Checked and found correct
- The a-b and V_w profiles match the legacy hand values at all nine kink temperatures for
  both cases; the self-test runs before any baking.
- `bulkFriction` is `tan(phi)`, not degrees; cohesion is linear in mu (verified by a
  midpoint test), which is what keeps ASAGI-interpolated cohesion consistent with mu.
- The emitted Lua reproduces its design to 0 over s in [-50, 500] km, and carries the
  descriptor's strike frame — now in `%.17g`, matching the shipped maps' formatting.
- `check_seissol_supports_spatial_muw` correctly rejects 1.1.3 and 1.3.2, accepts 1.4.0.
