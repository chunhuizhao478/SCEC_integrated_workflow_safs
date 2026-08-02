# Code Review: Phase 6 — deck assembly and the P1–P8 pre-flight (2026-08-02)

Reviews for Phases 0–5 are in the git history; all 30 of their findings were fixed.

## Review Scope
- Plan: `docs/PLAN_...md`, Phase 6
- Files: `deckbuild/deck.py`, `run_workflow.py`, `tests/test_deck.py`
- Baseline: 282 tests pass; the full workflow is green end to end on `demo_planar`,
  producing a runnable deck with all eight pre-flight gates passing.

## Findings

### [R-601] CRITICAL [deck.py:preflight] — P1's sha256 check could never fire

**Category:** BUG

**Description:**
P1 looked the deck's files up in the manifest keyed on `Path(artifact["path"]).name` —
the **source** filename in `outputs/` (e.g. `friction_depth_profile.nc`). But the deck
copy is deliberately **renamed** so the filename encodes the design
(`demo_planar_friction_case1.nc`). The names never match, so `r in by_name` was always
False and the hash comparison was dead code.

The tamper test proved it: corrupting eight bytes of a shipped nc left P1 reporting
"bad none". The gate advertised integrity checking and delivered none.

**Trigger:** Edit any nc inside an assembled deck and run `preflight`.

**Suggested fix:** Record the DECK filenames at assemble time.
```diff
+        man.stages["deck_files"] = {
+            e["filename"]: sha256_file(Path(e["dst"])) for e in paths.values()}
```
```diff
-        by_name = {Path(a["path"]).name: a for st in man.get("stages", {}).values()
-                   for a in st.get("artifacts", [])}
+        by_name = man.get("stages", {}).get("deck_files", {})
```

**Test case:** `test_P1_flags_a_tampered_file` (now passing; it failed with a traceback
before, which is how R-602 surfaced).

---

### [R-602] MODERATE [deck.py:preflight] — a corrupt deck produced a traceback, not a report

**Category:** BUG

**Description:**
With a corrupted nc in the deck, `preflight` crashed with
`OSError: [Errno -101] NetCDF: HDF error` from a downstream gate before any report was
returned. The entire point of a pre-flight is to tell you what is wrong *before* you
queue a job; a traceback with no gate output is the opposite.

**Suggested fix:** Short-circuit after P1 when a referenced file is missing or unreadable.
```diff
+        if missing or bad_hash:
+            for g in ("P2", "P3", "P4", "P5", "P6", "P7", "P8"):
+                rep.skip(g, "not attempted: P1 found a missing or unreadable referenced "
+                            "file, so the downstream checks cannot be trusted")
+            return rep
```

---

### [R-603] MODERATE [deck.py:_collect_refs] — the mesh was reported as an unaccounted extra

**Category:** BUG

**Description:**
`_collect_refs` scanned only the YAMLs for `file:` targets. The mesh is named by
`MeshFile` in `parameters.par` (without its `.puml.h5` suffix), so P1 saw the mesh sitting
in the deck, could not match it to any reference, and reported it as an unexplained extra
— failing a deck that was in fact correct. A gate that cries wolf on every valid deck gets
ignored, which then hides the real case it exists for.

**Suggested fix:** Parse `MeshFile` from `parameters.par` and resolve the suffix.

---

### [R-604] MODERATE [deck.py:_check_yield] — P5 skipped whenever the grids differ, i.e. always

**Category:** DEVIATION

**Description:**
P5 compared the stress and plasticity fields elementwise and skipped if their shapes
differed. The four grids are **independent by design** — the shipped SAFS ones are
material 1500/250, stress 1000/250 — so they essentially always differ and P5 never ran.

**Suggested fix:** Sample the plasticity onto the stress grid with `trilinear_sample`,
which is what ASAGI does at run time anyway. P5 now runs and reports
`0 of 350,811 grid points yield`.

---

## Summary
- Critical: 1 (R-601) | Moderate: 3 (R-602, R-603, R-604) | Low: 0
- Plan compliance: **FULL** — assemble (copy-not-symlink, verified by hash), the filename
  contract, `resolve_paths` as a dry run, preserved hand-written notes, `diff`, and all
  eight pre-flight gates.
- Verdict: **PASS WITH FIXES** — all four fixed in this round.

## Checked and found correct
- Copy, never symlink; every copy is sha256-verified against its source.
- The prefix comes from the descriptor; no filename starts with a hardcoded `safs_`.
- The YAML `file:` field and the copied filename come from the same variable, and P1
  re-checks that on disk.
- Regenerating a deck preserves the hand-written notes block verbatim.
- Assembly refuses artifacts built from two different descriptors.
- P6 is a WARN and says so in its own detail: kappa is necessary, not sufficient.
- P8 catches a `z = 0` receiver, a literal `0.0` in `OutputRegionBounds`, and
  `wavefieldoutput = 0` — three defects that each cost real runs.
- `diff` reports exactly one changed file for a one-byte edit.
