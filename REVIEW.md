# Code Review: Phase 8 — the SAFS reproduction exercise (2026-08-02)

Reviews for Phases 0–7 are in the git history; all 37 of their findings were fixed.

## Review Scope
- Plan: `docs/PLAN_...md`, Phase 8
- Files: `deckbuild/introspect.py`, `exercise_safs_reproduction.py`,
  `tests/test_introspect.py`, `docs/EXERCISE_safs_reproduction.md`
- Run against the real shipped deck
  `safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_intermediate_..._deep40km`
  (683 MB mesh, 362 MB CVM) with the real raw tree.
- Baseline: 303 tests pass.

## Result of the exercise itself

| Rung | Verdict | Evidence |
|:--|:--|:--|
| E0 mesh | **PASS** | 3,210,006 nodes / 15,979,903 tets / 320,560 BC-3 faces, all interior |
| E1 provenance | **PASS** | 37 CVM slices, 105 CTM slices, CSM present |
| E2 descriptor | **PASS** | 37 fields recovered, **0 unknown** |
| E2c freeze | **PASS** | measured on the FIELD: freeze is OFF, agreeing with the filename |
| E3 material | **BLOCKED** | the `cvm_slices` reader is not implemented |
| E4–E6 | **BLOCKED** | inherit E3 |

**The reproduction claim is therefore UNPROVEN.** E0–E2 establish that the deck is readable
and its design fully recoverable — the prerequisite for E3–E6, not a substitute.

## Findings

### [R-801] MODERATE [material.py] — `cvm_slices` / `ctm_slices` block the whole exercise

**Category:** DEVIATION

**Description:**
`VELOCITY_READERS` contains only `layered_1d`; `cvm_slices` raises "not implemented yet",
and the `ctm_slices` thermal path likewise. The plan's Phase 4 requirement 1 specifies the
two-stage raw pipeline (slice read → reproject → inscribed UTM grid → per-slice
`LinearNDInterpolator` → depth-to-elevation flip → moduli at source nodes → z-resample),
and Phase 8's E3 depends on it.

Raising is the right *interim* behaviour — it is strictly better than silently producing a
wrong field — but it means the plan's central acceptance test cannot run.

**Suggested fix:** Implement `_cvm_slices` and `_ctm_slices` in `material.py`, porting the
legacy `generate_velocity_nc_from_raw.py` / `generate_thermal_nc_from_raw.py` verbatim.
Then re-run the ladder; expect E3 to be the rung that reveals whether scipy's Delaunay
matches the version the shipped nc was built with (the plan's "might not match" category).

**Estimated scope:** ~200 lines and one expensive test (452×361×194 grid from 37 slices).

---

### [R-802] LOW [exercise_safs_reproduction.py] — E1 does not actually prove correspondence

**Category:** ASSUMPTION

**Description:**
E1 counts the raw slices and prints the shipped nc's shape, then reports PASS. It does not
compare the slices' lon/lat extent or depth levels against the nc's axes, so it establishes
*availability*, not *correspondence* — which is what the rung is named for.

The honest reading of the current PASS is "the raw inputs are present and plausibly the
right products", and the evidence string says the slice count and the z-level count need
not match. It should say what it did NOT check.

**Suggested fix:** Either compare the reprojected slice hull against the nc's x/y extent,
or rename the rung's verdict to `PARTIAL` with an explicit "correspondence not verified".

---

### [R-803] LOW [introspect.py:_parse_tnuc] — the z-sign parse is fragile

**Category:** EDGE_CASE

**Description:**
The hypocentre's z is recovered from `dz = x["z"] + 10067.9819` by string-inspecting the
sign. It happens to be right for the shipped deck (verified: `-10067.9819`), but it depends
on the emitter writing `+` for a negative z. A generated Lua that wrote
`x["z"] - -10067.98` would parse wrong.

**Suggested fix:** Evaluate the arithmetic rather than inspecting the sign character, or
require the emitter's exact form and assert it.

---

## Summary
- Critical: 0 | Moderate: 1 (R-801) | Low: 2 (R-802, R-803)
- Plan compliance: **PARTIAL** — the introspector, the ladder harness, the
  declared-before-comparison categories, and the results document all exist and work; the
  ladder itself stops at E3 for a named, recorded reason.
- Verdict: **PASS WITH A BLOCKING GAP** — nothing here is wrong, but Phase 8's purpose is
  not yet served. R-801 must be closed before the workflow can claim reproduction.

## Checked and found correct
- The introspector recovers **37 fields with zero unknowns** from a real deck, including
  three that exist only as literals inside emitted Lua: the strike azimuth (314.0) and
  origin (606971, 3707270) from the `rs_muw` coefficients, and the hypocentre
  (604446.944, 3704576.3853, −10067.9819) from `Tnuc_s`.
- The 7-region `f_w` design read out of the Lua matches the shipped yaml's own filename,
  value for value: `[0, 0, 0.045, 0.03, 0.06, 0.0175, 0.05]`.
- `phi` is read from `bulkFriction = tan(phi)` in the FIELD, not from the yaml prose —
  giving 30/40, the production values, not Roten's 35/45.
- The shallow freeze is verified **against the field**, not the filename, exactly as the
  plan requires. Both agree: it is off for this deck.
- The four grids are confirmed independent on the real deck: material 1500/250,
  stress 1000/250, friction 1500/200.
