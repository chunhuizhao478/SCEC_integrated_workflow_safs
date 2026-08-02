# Code Review: Phase 5 — mesh ingest, gating, Stage F (2026-08-02)

Reviews for Phases 0–4 are in the git history; all 25 of their findings were fixed.

## Review Scope
- Plan: `docs/PLAN_...md`, Phase 5
- Files: `deckbuild/mesh.py`, `deckbuild/stage_f.py`, `MESHING.md`,
  `skills/code-mesh-build-improve/SKILL.md`, `tools/make_demo_mesh.py`,
  `data/demo_planar/demo_planar.puml.h5`, `tests/test_mesh.py`
- Baseline: 253 tests pass and the full workflow is green on demo_planar.

## Findings

### [R-501] CRITICAL [mesh.py:msh_to_puml] — an interior fault face is tagged on only ONE side

**Category:** BUG

**Description:**
`_searchsorted_rows` returned a single index per tagged triangle. A fault triangle is
**interior**: it is a face of exactly two tets, and BOTH tet-face slots must carry the BC.
Taking the first match leaves the other side at 0.

Measured by round-tripping the demo mesh through `.msh` and back:

```
BC identical : False   mismatches: 192 of 9216
```

192 is exactly the fault-triangle count. Every fault face lost one of its two sides.

Gate C *does* catch the result ("fault face not interior"), which is the gate working —
but the converter should not produce it, and a user converting a real `.msh` would see a
confusing C failure on a mesh that is geometrically fine.

**Trigger:** `msh_to_puml` on any mesh with an embedded (interior) fault.

**Suggested fix:** Return the half-open `[lo, hi)` range and tag every match.
```diff
-    idx = _searchsorted_rows(sorted_faces, key)
+    lo, hi = _match_rows(sorted_faces, key)
...
-        flat = order[i]
-        bc[flat % len(tets), flat // len(tets)] = tag_to_bc[int(tri_tags[j])]
+        for i in range(lo[j], hi[j]):
+            flat = order[i]
+            bc[flat % len(tets), flat // len(tets)] = tag_to_bc[int(tri_tags[j])]
```

**Test case:** `test_R501_msh_roundtrip_tags_both_sides_of_the_fault` (added).

---

### [R-502] MODERATE [mesh.py] — `np.core.records` was removed in NumPy 2

**Category:** BUG

**Description:**
`_searchsorted_rows` used `np.core.records.fromarrays`. `np.core` is removed in NumPy 2.x
(this environment runs 2.3.3); it currently resolves only through a deprecation shim and
will stop. The replacement is `np.rec.fromarrays`. Fixed as part of R-501.

---

### [R-503] MODERATE [tests] — the `.msh` conversion path had no test at all

**Category:** DEVIATION

**Description:**
Phase 5's first acceptance criterion is about `msh_to_puml`, and 27 tests covered
everything *except* it. R-501 lived there undetected. A round-trip test is now added, and
it is the test that found the bug.

---

### [R-504] LOW [mesh.py:MeshStage.verify] — gate G is a permanent skip

**Category:** DEVIATION

**Description:**
Gate G (the pickpoint/partition-boundary screen) always reports skip, because it needs a
receiver list that only exists once `DeckStage` runs. That is honest, but it means the
gate is currently decorative. It should move to Phase 6 and run there, with `verify`
noting that it lives downstream rather than implying it could run here.

---

## Summary
- Critical: 1 (R-501) | Moderate: 2 (R-502, R-503) | Low: 1 (R-504)
- Plan compliance: **FULL** for ingest, gating, Stage F, the vendored skill and MESHING.md;
  gate G is deferred to Phase 6 by design (R-504).
- Verdict: **PASS WITH FIXES** — R-501 and R-502 fixed in this round.

## Checked and found correct
- The PUML contract round-trips exactly: dataset names, dtypes (`f8`/`u8`/`i4`), the
  `boundary = sum_i code_i << (8*i)` packing and both file attributes.
- Face multiplicity uses a structured lexsort, not a packed int64 key: verified against a
  node index of 2.2M, which would overflow a 3x21-bit key.
- Gate C correctly rejects a hull face mislabelled as fault; gate D correctly counts a
  planted inverted tet; gate E is a skip, never a pass, without a material nc.
- Stage F refuses to cross-check when the fields came from a different descriptor.
- The vendored skill differs from the source in exactly the two allowed hunks.
- Stage F earned its keep immediately: it caught an off-fault hypocentre in the **demo
  descriptor I had just written** — `lon = -121` projects 180 km east of the demo fault.

---

### [R-505] MODERATE [mesh.py:msh_to_puml] — a malformed `.msh` kills the interpreter

**Category:** BUG

**Description:**
`meshio.read` calls `sys.exit()` on some malformed headers. `SystemExit` derives from
`BaseException`, so `except Exception` does not catch it: a v4 or truncated `.msh` would
terminate the caller's process instead of raising the actionable "re-export as msh22"
message. Found by the v4 test, which failed with `SystemExit: 1`.

Same rule already applied to `geometry._load_puml`: a library must not exit.

**Suggested fix:**
```diff
-    except Exception as exc:                                  # noqa: BLE001
+    except (Exception, SystemExit) as exc:
```

**Test case:** `test_msh_v4_gives_an_actionable_message` (now passing).
