# Code Review: Phase 7 — the integrated notebook and the demo project (2026-08-02)

Reviews for Phases 0–6 are in the git history; all 34 of their findings were fixed.

## Review Scope
- Plan: `docs/PLAN_...md`, Phase 7
- Files: `deck_workflow.ipynb`, `tests/test_notebook.py`, `pytest.ini`,
  `data/demo_planar/demo_planar.puml.h5`, `README.md`
- Baseline: 293 tests pass; the notebook executes end to end in ~2 s with every gate
  battery printing and no `[FAIL]`.

## Findings

### [R-701] MODERATE [deck_workflow.ipynb / run_workflow.py] — both write to the SAME deck directory

**Category:** BUG

**Description:**
The notebook derives `DECK = decks/{name}_{RUN_TAG}` and `run_workflow.py` derives the
identical path from the same tag. They do **not** produce the same deck — the notebook
adds a `Tnuc_s` nucleation LuaMap that the CLI runner does not — yet both call
`assemble(..., overwrite=True)`. Whichever ran last silently wins, and a user who runs the
notebook, then the CLI, then inspects `decks/…` is looking at a deck that is not the one
their notebook reported on.

This is the same class as the legacy "shared `outputs/graded_fw/`" defect the plan cites,
where a stale pre-edit file was copied into a deck.

**Trigger:** Run the notebook, then `python run_workflow.py`, then read the deck.

**Suggested fix:** Make the producer part of the path, and record it.
```diff
-DECK = B.root / "decks" / f"{cfg.name}_{RUN_TAG}"
+DECK = B.root / "decks" / f"{cfg.name}_{RUN_TAG}_nb"     # _nb: the notebook's own deck
```
and in `run_workflow.py`:
```diff
-    deck = ROOT / "decks" / f"{cfg.name}_{tag}"
+    deck = ROOT / "decks" / f"{cfg.name}_{tag}_cli"
```

---

### [R-702] LOW [deck_workflow.ipynb] — the notebook re-derives paths the runner already computes

**Category:** QUALITY

**Description:**
Section [1] hand-builds `OUT` and `DECK` and `mkdir`s the subdirectories, duplicating
`run_workflow.py`. The two can drift — and R-701 is exactly that drift. A shared helper
(`deckbuild.layout.run_paths(cfg, tag, producer)`) would make one the definition and the
other a caller.

Not fixed now: it is a refactor, and the plan's rule is to implement the phase, not
restructure. Recorded for Phase 9's cleanup.

---

### [R-703] LOW [tests/test_notebook.py] — the execution test does not assert the deck it produced

**Category:** DEVIATION

**Description:**
`test_notebook_runs_end_to_end_under_five_minutes` asserts the gate text and the wall
clock, but not that a deck folder now exists with the expected files. A notebook that
printed every gate and then failed to write anything would pass.

**Suggested fix:** Assert the deck directory and its manifest exist after the run.

---

## Summary
- Critical: 0 | Moderate: 1 (R-701) | Low: 2 (R-702, R-703)
- Plan compliance: **FULL** — eight sections, one assignments-only PARAMETERS cell each,
  physics stated before the code, the demo running with zero downloads and no gmsh under
  the 5-minute budget, section [2] explicit that it ingests rather than builds, and the
  retarget checklist in the README.
- Verdict: **PASS WITH FIXES** — R-701 and R-703 fixed in this round.

## Checked and found correct
- All 20 PARAMETERS values are actually consumed downstream (probed by counting
  references); none is decorative.
- Every PARAMETERS cell parses to assignments only — no logic hidden where a user edits.
- The notebook contains no absolute path and no network access of any kind.
- Material precedes stress and friction, and the intro says why the order is load-bearing.
- The six gotchas that cost real runs are stated in the prose the user reads: the dropped
  `z = 0` receiver, the silently-ignored spatial `rs_muw`, `bulkFriction` being a
  coefficient, kappa being necessary-not-sufficient, sea-level depth referencing, and the
  strike-slip assumption behind the Andersonian closure.
