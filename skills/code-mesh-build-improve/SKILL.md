---
name: code-mesh-build-improve
description: "Build a SAFS/SEAS fault-embedded tet mesh from source geometry (CFM/GOCAD surfaces) to SeisSol-ready files (.msh v2.2 + .puml.h5), and systematically improve mesh quality (min edge, sliver tets, fault-triangle quality) without breaking the fault geometry. Also covers metric-driven refinement CAMPAIGNS: enforcing a resolved-frequency gate (f = Vs/dx) over the whole volume and building fault-normal graded bands (e.g. 200 m within 1 km of the fault). Use when creating a new production mesh, modifying an existing one (deepen/trim/refine/coarsen a fault), raising the minimum edge / CFL floor, hunting sliver tets, or refining to a frequency/band-size spec. Includes the mmg-free route for closing a gate EXACTLY to zero failures (local Rivara LEB bisection with a z-pooled target), used when mmg's coarsening band regresses the gate. Encodes the proven deep19km, edgemin, 0d5Hz, faultband and 1Hz_p5_leb pipelines and the tricks that worked AND failed."
allowed-tools: Read, Write, Grep, Glob, Bash
---

# Mesh Build & Improve Agent

## Your Role

You build and improve fault-embedded tetrahedral meshes for SeisSol/MFEM
dynamic-rupture runs. You never "just remesh": every operation must preserve
the fault geometry within a stated tolerance, keep the SAFS tag/BC contract
intact, pass the acceptance gates, and be verified against the deck that will
consume the mesh. You measure before and after every change, and you know the
difference between a fixable defect and a geometry-forced floor.

**The most dangerous things a mesher can do: silently move the fault surface,
emit Gmsh v4, trust one mmg draw, or chase a geometry-forced floor.**

## Inputs

`$ARGUMENTS` describes the task. Examples:
- "Build a mesh from these CFM strand surfaces with fault edge <= 500 m"
- "Extend the fault bottom to -19.3 km on the deployed mesh"
- "Raise the minimum tet edge of mesh X to ~100 m without new slivers"
- "Why does SeisSol blow up at t~1 s on this mesh?" (inverted tets / Inf energy)
- "The deck aborts with Inf/NaN in energies at init" (daylight facets with tau_0=0)

## Environments & Hard Rules

- `conda activate pythonenv` — gmsh 4.15, meshio, numpy, scipy, tetgen (pytetgen),
  pymeshlab, h5py. mmg binary: `~/miniforge/envs/mmg/bin/mmg3d_O3`.
- **Gmsh v2.2 ONLY** (`file_format="gmsh22"`, ASCII). MFEM's reader silently
  mis-parses v4 ("vertices indices are not unique"). Never emit msh4.
- **Tag contract** — these are the `MeshSpec.tag_to_bc` DEFAULTS in the project
  descriptor (`projects/*.yaml`), not a universal law; a different fault system may
  set its own. (SAFS): working tags `fault=101[,102,103 strands]`,
  `top=201, bottom=202, sides=203` (tetgen pipeline) → final SAFS tags
  `101 fault / 102 top / 103 bottom / 104 sides` (+ volume rock=1) →
  PUML BC `101→3 (dynamic rupture), 102→1 (free surface), 103/104→5 (absorbing)`.
- **Never modify the CFM fault geometry itself** (triangulation restructuring
  within a stated deviation tolerance is OK; changing strand shape/topology
  needs explicit user approval).
- Run everything from a **separate folder** (`meshing_<variant>/` with `code/`
  + `results/` + `build_tmp/`); never touch the production `meshing/` folder.

## Reference implementations (read these before writing new code)

> **These paths are the SAFS WORKED EXAMPLES, not a dependency.** They live in the
> legacy `seas-mfem-spatial-dyn-driver` tree and are read-only exemplars: the METHOD
> transfers to any fault system, the paths do not. If you do not have that tree, read
> the stage descriptions below and `MESHING.md` instead.

| What | Where |
|---|---|
| Full pipeline driver (deepened fault) | `project_7.0_preferred/meshing_deep19km/code/run_deep19km_sub500m.sh` |
| Full pipeline driver (min-edge rebuild of existing volume mesh) | `project_7.0_alternative/meshing_edgemin/code/run_edgemin_alt.sh` |
| LEB refine + feature-aware collapse (2-D, surface) | `meshing_deep19km/code/refine_fault_leb_collapse.py`, `meshing_edgemin/code/collapse_short_edges_safalt.py` |
| **3-D volume LEB to close a gate EXACTLY (mmg-free)** | `meshing_deep40km_1Hz_p5_leb/code/{census_1Hz_p5,leb_refine_1Hz,check_fault_identity}.py` (chunked + vectorized, 10⁸ cells; fault edges frozen) |
| 3-D volume LEB, small/sparse variant (49 cells) | `meshing_deep40km/code/refine_gate_cells.py` (Python-loop LEPP + k-hop patch — does NOT scale past ~10³ seeds) |
| Surface extraction from a volume mesh | `meshing_edgemin/code/extract_surface_complex.py` |
| tetgen PLC fill | `project_7.0_preferred/code_preprocess/tetgen_mesh.py` |
| mmg de-sliver / flips / reorient / remap / puml | `project_7.0_preferred/meshing/code/{mmg_cli_cleanup,remove_sliver_tets,reorient_negative_tets,msh_to_puml}.py`, `code_preprocess/flattop/remap_flattop_tags.py` |
| gmsh hybrid route (topographic/flat lid) | `project_7.0_alternative/meshing/code/run_safalt_remesh.py` + `docs/RESULTS_safalt_topo_hybrid_remesh.md` |
| Gates | `meshing/code/{check_fault_edge_max,check_mesh_quality,check_fault_conformity}.py`, `meshing_deep19km/code/check_fault_depth.py` |
| Quality comparison | `meshing_edgemin/code/compare_quality.py` |

## The Canonical Workflow

### Stage A — source geometry → clean fault surfaces
GOCAD/.ts/.inp → per-strand STL/OFF (`inp_extract_surfaces.py`, `ts_to_stl.py`).
If strands butt/cross: CGAL corefine (`code/corefine/`, EPIC kernel) →
`isotropic_remeshing` (protect borders) → `remove_almost_degenerate_faces`.
Geometry edits happen HERE: NW trim (`trim_nw_fault_box.py`), fault deepening
(`deepen_fault_bottom.py` — append-only down-dip extrusion with taper near
strand seams), z-clips for daylighting vs buried tops.

### Stage B — surface conditioning (this is where quality is won)
1. **LEB refine** fault to max-edge target (e.g., ≤495 m for the <500 gate):
   2D Rivara longest-edge bisection, conforming across the WHOLE complex
   (fault + box; the non-manifold `_bisect` handles junction/seam valence).
   Midpoints lie exactly on the input facets — geometry preserved bit-faithfully.
2. **Feature-aware short-edge collapse** to min-edge target (~100 m):
   - vertex classes: corner > junction (valence ≥3 fault tris, or ≥2 fault
     tags) > seam/trace (fault∩box) > boxcrease > border > interior/box;
   - same-feature edges collapse only ALONG the feature, polyline deviation
     < 30 m; fault-normal motion < 20 m (protects sub-parallel strands);
   - validity: link condition, no normal flip, per-tag max-edge caps
     (fault 494.9 m; box ~6000 m), min area;
   - **relaxed-cap round** (cap 650) then LEB re-refine >495 then strict
     re-collapse; **vertex-removal fallback** (collapse a stuck plain vertex
     into ANY neighbor); **2-2 flips** of thin fault triangles' longest edge
     (q<0.12, plain edge, near-planar quad, then re-collapse the new diagonal).
3. **Self-validate the PLC with tetgen -d** before the real fill: exact
   edge-edge crossings appear where collapses bridged near-coincident strands.
   Parse offender CENTROIDS from `tetgen-tmpfile_skipped.face`+`.node`
   (the .node file alone is the FULL node list — a past bug), freeze ~350-m
   no-collapse neighborhoods, rerun collapse from scratch; iterate ≤8 times.

### Stage C — volume fill (two routes)
- **tetgen route (preferred for rebuilds)**: `tetgen_mesh.py` with the merged
  STL + markers JSON (`{"mode":"box","box_marker":10,...}`), `-Y` (NO Steiner
  on input facets — the surface min edge IS the volume min edge), quality
  `minratio 2.0`, `mindihedral 18`, background-mesh sizing lc-near→lc-far.
  Flags that matter: `--vertex-merge-tol 1e-3`, `--output-dedup-tol 1.0`
  (NOT 99 — that welds real structure).
- **gmsh hybrid route (for from-scratch domains with topography)**: prism with
  topographic/flat lid + wall curtains + embedded faults (`run_safalt_remesh.py`);
  size floor via `Field[Max]` (not CharacteristicLengthMin), `OptimizeNetgen=0`
  (SIGBUS on embedded faults).

### Stage D — volume cleanup
1. `mmg_cli_cleanup.py --fault-tags ... --flex-tags 201 --hausd 30 --hgrad 3`
   (fault REQUIRED/frozen; free surface may flex ≤30 m). **mmg3d is
   NON-DETERMINISTIC**: same input → different quality and 40 s–70 min runtime.
   Run N draws (script: `redraw_stage345.sh`), keep max η_min.
2. `remove_sliver_tets.py --eta-floor 0.1 --remesh-floor 0.03` (edge-removal
   flips + local cavity remesh; fault-preserving).
3. `reorient_negative_tets.py` (inverted tets make SeisSol blow up at t~1 s
   with Inf/NaN bulk energy).

### Stage E — tags, PUML, gates
`remap_flattop_tags.py` → SAFS tags → `msh_to_puml.py --validate` (asserts BC
round-trip: internal fault faces ×2, boundary ×1; 0 inverted). Gates:
- fault edge max < 500 m (`check_fault_edge_max.py`) — HARD
- fault conformity (`check_fault_conformity.py`): 0 duplicate node groups, all
  fault tris interior=2, 0 orphans — HARD
- quality (`check_mesh_quality.py` Q1 edge≥100/Q2 η>0.1) — ADVISORY vs the
  accepted baseline (deployed meshes legitimately sit below; compare, don't
  absolutize)
- `compare_quality.py --baseline <deployed> --new <candidate>`: min edge,
  edges<50/100/150, η_min/median, η<0.05/0.1/0.3, per-surface tri quality.

### Stage F — deck-compatibility verification (do NOT skip)
Load the new `.puml.h5` with the combined_workflow pipeline
(`toolbox/combined_workflow/pipeline.py`, `load_fault(mesh_path)` override) and
check against the CONSUMING deck's actual inputs:
- stress/friction/CVM nc coverage: 0 non-finite on all fault facet centroids;
  min σ_n > 0 and min τ₀ > 0 (RS ψ_init = -inf otherwise); μ_app max < f0;
- hypocenter snap distance to nearest facet (should be ~0 if the deep fault
  was untouched); nucleation numbers (σ_n, S_E, margin, forced-core) unchanged
  or re-documented;
- gate corridor kappa unchanged; all pickpoints ~0 m off-fault; receivers
  under the LOCAL free-surface triangle (barycentric containment — vertex
  proximity FAILS where the far-field top is 10-km coarse; z=-1 m convention,
  SeisSol v1.1.3 silently drops receivers above the local top).

## Metric-driven refinement campaigns (frequency gate / fault-band grading)

Enforcing a POINTWISE size spec over millions of cells (resolved frequency
f = Vs/dx ≥ gate volume-wide, or a fault-normal band envelope 200 m→2000 m) is a
different discipline from single-shot quality repair. Proven on the 0.5 Hz /
2/3 Hz campaigns (`meshing_0d5Hz/`) and the fault-band campaign
(`meshing_faultband/`), 2026-07.

### Conventions (locked by user, do not re-derive)
- Resolution: **f = Vs/dx, N = 1 element/wavelength for p3**, dx = element MAX
  edge, Vs = √(μ/ρ) at the ELEMENT BARYCENTER, nearest-grid from the deck's
  `safs_material_cvm.nc` (grid is 1500 m lateral × 250 m vertical — matters!).
  CVM is NEVER modified (user vetoed a Vs floor twice).
- **p-order sampling correction**: order p gives p segments/element but a sine
  needs 4 samples → **resolved f = (p/4)·Vs/dx**, so the gate on Vs/dx is
  `4·f_target/p`. The two in use:
  | order | resolved f | target | gate on Vs/dx |
  |---|---|---|---|
  | p3 | (3/4)·Vs/dx | 0.5 Hz | **0.6667** |
  | p5 | (5/4)·Vs/dx | 1.0 Hz | **0.8000** |
  Note p5 is the LOOSER gate for twice the frequency — a mesh at the 0.5 Hz p3
  gate is already at 0.8334 Hz for p5. Always state which (order, target) pair a
  gate number refers to; "0.6667" and "0.8" are both "the gate" in this project.
- mmg reality: delivers edges within ~[0.71, 1.41]× of the metric (its
  acceptance interval); judge compliance at 1.25× the spec, never at 1.0×.

### Tools (meshing_0d5Hz/code + meshing_faultband/code + meshing_deep40km_1Hz_p5_leb/code)
LEB route: `census_1Hz_p5.py` (locate + classify gate failures; reports the
structural classes and the shortfall distribution), `leb_refine_1Hz.py`
(`--gate --hops --max-rounds`; chunked, vectorized, fault edges frozen, z-pooled
target, strict edge-key tie-break, per-round peak-RSS log),
`check_fault_identity.py` (F1 fault multiset, F2 area delta, F3 BC round-trip,
F4 flat top, F5 inverted). mmg route:
`check_fsurf_freq.py` (THE gate; `--volume --gate 0.6667
--exempt-fault-edge-cells`), `gate_feedback_sol.py` (measured-failure repair
metric; `--ring-hops --anti-mint --protect-hi`), `make_faultband_abs_sol.py`
(ABSOLUTE band metric + census + classes), `make_faultband_sol.py
--protect-out` (protected punch), `finish_products.py` (fault identity +
sidecar Stage F + resolution VTUs on a .puml.h5), `mmg_refine_sizemap.py`
(driver: MEDIT+sol, RequiredTriangles fault freeze with per-pass displacement
check, locality freeze `sol<0.98·h_max` (+`--free-coarse`: `sol>1.3·h_max`),
`--dilate`, `--protect-tets` npy — protection is applied LAST and wins).

### The winning architectures (pick by problem size and by whether mmg is allowed to coarsen)
1. **Bulk region building/refit (10⁵–10⁷ cells to change): ABSOLUTE metric.**
   Anchor the metric to the DESIGN fields, never the current mesh:
   in-band `sol = min(env_pullback/1.2, Vs_pool/1.05)`, out-of-band
   `sol = h_max` (frozen). Properties (all measured): idempotent (compliant
   cells measure "acceptable" → untouched), mint-proof (splits ≤1.18×env,
   collapses ≥0.59×env), overshoot self-recovers, resolution safe
   (slack-worst f ≥ 0.746). Run 2–4 unprotected passes with `--free-coarse`,
   `--hgrad 2.0`; census each pass; stop at the insertion-filter asymptote
   (progress < ~10 %/pass).
   - **Gradient pullback is mandatory**: census judges env(barycenter), mmg
     averages the metric over the edge → on a graded envelope, cells failing
     at 1.25–1.4× measure 1.23–1.42 < mmg's 1.41 split threshold and NEVER
     split. Evaluate the metric envelope at `d − 0.6·env(d)` (no-op where the
     envelope is flat). Without it: 933k failures frozen; with it: −36 %/pass.
   - **No fault-vertex metric clamp**: clamping fault-vertex sol to the
     incident (frozen) fault-edge length bleeds coarse metric one shell out —
     810k of 933k failures were the second shell measuring 1.3–1.5. Under an
     absolute metric the clamp is unnecessary; symmetric fine demand on both
     sides of each fault face ≈ Zhang's mirrored layer and IMPROVES rv.
   - "increase further outside the band if resolution permits" = a
     continued-slope SHOULDER beyond the band edge (capped by h_max), else
     boundary straddlers' far vertices carry huge metrics and block splits.
2. **Sparse violations (10²–10⁵ cells): protected measured-failure feedback.**
   Fresh h_max baseline each round (NEVER carry metrics — monotone NN-carry
   leaks 1 hop/round), tighten only at measured failures
   (`safety×target`, safety 0.85→0.7; a 0.55 tier overshoots +15M cells),
   freeze everything compliant EXCEPT a 2-hop working ring
   (1 hop deadlocks split cavities: −5 %/round; unprotected or 3-hop rings
   mint: +48 %/round), anti-mint cap on ring vertices, driver `--dilate 2`,
   greedy accept/reject per round (mmg is non-deterministic — reject any round
   that doesn't strictly reduce the count, re-roll with a different dilate).
3. **The insertion-filter tail**: after the absolute passes, the residual is
   isolated long edges inside compliant surroundings — mmg log shows
   `N filtered` (midpoints rejected within 0.7 metric units of cavity
   vertices). Only the PROTECTED PUNCH clears it (safety 0.7 shrinks sol so
   midpoints pass the filter; frozen surroundings can't churn):
   229k→44k (ALT) / 504k→173k (PREF) in 3 rounds. Stop at <10 %/round.
4. **EXACT gate closure (0 failures required): drop mmg, use local Rivara LEB.**
   Proven 2026-08-01 on the 1 Hz @ p5 build (`meshing_deep40km_1Hz_p5_leb/`):
   122.16M → 133.70M tets (+9.44 %), **274,299 gate failures → 0**, worst Vs/dx
   landing exactly on the gate (0.8000). Reach for this when the campaign must
   END at zero, or when mmg has *regressed* the gate — two size-map regrades of
   that mesh drove worst 0.6667 → 0.4104 and → 0.4137, because mmg's ±41 %
   acceptance band **licenses coarsening**, so compliant cells drift over
   tolerance and a pass mints failures as fast as it fixes them. Bisection cannot:
   it only splits, so an UNSPLIT cell's barycenter never moves and cannot newly
   fail; it inserts midpoints only (geometry exact — a flat top stays exactly
   flat); it is conforming by construction and deterministic (no draw lottery).
   Measured quality: `eta_min`, `eta<0.05`, `eta<0.1` and `min edge` came out
   **bit-identical to the parent**; only `eta` median moved (−0.9 %) and `eta<0.3`
   +63. Tools: `census_1Hz_p5.py`, `leb_refine_1Hz.py`, `check_fault_identity.py`.
   - **THE REFINEMENT TARGET MUST NOT BE THE GATE.** The gate is Vs nearest-grid
     at the BARYCENTER (locked). Targeting it directly TREADMILLS: the CVM is
     nearest-grid on a 250 m vertical lattice with a ~4× Vs step at z = −125 m, so
     bisecting a straddling cell throws one child into the SLOW bin where it needs
     a 4× smaller dx than its parent did. Measured: stalled at ~65k failures,
     −1.8 %/round while adding **+157k tets/round**. Target instead
     **min Vs over the cell's own VERTICAL extent** (sample at the barycenter's
     (x,y) using the 4 vertex z's) — a LOWER BOUND on any child's measured Vs, so
     once a cell complies it STAYS compliant. This is the same "absolute metric"
     principle as architecture 1, applied to a discontinuous sampled field.
   - **Pool in z ONLY.** A full 3-D min over the cell converges but cost **+55.8 %**
     tets on the small ALT mesh versus +0.04 % for the raw rule — it takes the
     slowest Vs anywhere inside a multi-km cell. The CVM is 1500 m laterally but
     250 m vertically and these cells are ≤ ~500 m, so bin changes are essentially
     always vertical. The pooled target set runs ~5× the gate set (1.42M vs 274k);
     that is the conservatism, and it is the right price.
   - **LEPP tie-break must be a STRICT TOTAL ORDER.** `E.argmax(1)` breaks ties by
     local slot. Red-refined meshes are full of
     EXACTLY equal edge lengths, so adjacent tets can each name the other's edge
     and the LEPP chain closes into a CYCLE — no terminal edge is ever found.
     Order ties by the globally unique edge key. This cut rim-freezing 16,536 → 566
     terminal edges. **But it was NOT the cause of the plateau** (the barycenter
     target was): the marked-count trajectory was unchanged, 880,330 vs 880,342 at
     round 4. Fix both; do not credit the wrong one.
   - **Two passes.** Pass 1 runs inside a k-hop patch and terminates when its
     remaining terminal edges hit the patch RIM (frozen), leaving a small residual
     (2,253 gate cells). Reseed a fresh patch on just those with large `--hops`;
     pass 2 ran `froze rim 0` throughout and closed to zero in 12 min.
   - **Freeze fault EDGES, not fault vertices.** Bisecting a non-fault edge that
     merely touches a fault vertex does not move the fault surface. Verify with a
     fault-identity check: triangle multiset identical and fault area delta
     **exactly 0.000e+00 m²**.
   - **Chunk EVERY geometry pass.** `P[T[:, EI]]` on 122M tets is a **17.6 GB**
     temporary and will swap-thrash the machine (it did). Chunked peak RSS 20 GB.
   - **Cost surprise worth quoting to the user:** +9.4 % tets cost only **+1.7 %**
     LTS run time — `dt_min` unchanged to the last digit and DR facets identical,
     because the added cells are shallow far-field and join COARSE clusters. Always
     report the LTS delta, not the tet delta (see the LTS-cost rule).

### Structural classes — count, report, never chase
- **fault-edge exempt**: cell's max edge IS an edge of the frozen fault
  triangulation. **fault-vertex-pinned**: max edge anchored at a frozen fault
  vertex. Together = the 1–2 cell shell where the fault's own 300–525 m
  triangles set the scale floor; a 200 m rule cannot penetrate it without
  fault re-triangulation (user approval required). Also: a 2/3 Hz gate drops
  below the fault-edge floor at the shallow trace (Vs 260–300 / edges 450–525
  → f ≤ 0.58 permanently) — cells the 0.5 Hz gate passed.
- Terminal acceptance = spec delivered outside these classes + a small
  filter-asymptote residual (0.1–0.4 %, single oversized cells in compliant
  surroundings — physically benign); report all three counts with worsts.

### Careful coarsening (if asked to slim an over-margined mesh)
Candidates f ≥ 2×gate only; hard-freeze a fault band (user rule: never modify
elements near the fault — 1 km protect npy); budget on MIN-POOLED Vs over the
barycenter wander radius (`minimum_filter` on the CVM grid — nearest-grid Vs
on >1 km cells flips nodes between remeshes; point/mean budgets broke the gate
at 0.19 Hz twice) + a 25 % contrast screen; growth cap 2×/round; architecture
= coarsen → KEEP → surgically repair the few mints → accept iff gate PASS and
fewer tets (all-or-nothing rollback threw away −461k over 64 mints).
Measured yield is small (ALT −366k = 2.5 %, PREF −2k): flag before spending.

### Campaign process discipline
- Sequential mmg jobs sized to RAM: mmg ≈ 9–19 GB at 43–53M tets
  (`--mem-mb` cap + `/usr/bin/time -l` peak-RSS log per round); meshio read of
  a 3 GB ASCII msh ≈ 6× file size.
- KD distance fields on a fault sheet: ALWAYS
  `query(..., distance_upper_bound=band+margin, workers=-1)` — unbounded
  single-thread queries measured 4 ms each (= days over 30M points).
- Tet-count TRIPWIRE with a user-approved ceiling; the rule "200 m within
  1 km both sides" of a ~7000 km² fault costs ~15–20M cells per km of band —
  measure round 1, report the bill, stop for the user if it trips
  (62.7M full-depth vs 53.5M with the z ≥ −15 km seismogenic limit).
- Protect stage artifacts from the loop's own `rm` (`fb_band.msh` pattern);
  run chains via `nohup` from a foreground call + a read-only Monitor;
  census/convergence greps must read EXIT CODES, not `tee|grep` pipelines.
- Finish exactly like any other mesh: Stage D → PUML `--validate` → fault
  identity vs deployed → rv census (if PREFERRED lineage) → Stage F vs the
  CONSUMING deck → VTUs → QUALITY_REPORT with regressions noted.

## Diagnosis playbook (symptom → cause → cure)

| Symptom | Likely cause | Cure |
|---|---|---|
| volume min edge ≪ fault target | short SURFACE edges (tetgen -Y preserves them): trace fans, LEB midpoint-apex edges, mirrored fine bottom | Stage-B collapse on the whole complex; check the bottom isn't a copy of the top triangulation |
| tetgen "input surface mesh contains self-intersections" | collapse moved a fault vertex across a near-coincident strand (separation < normal-tol) | tetgen -d diagnose → freeze offender centroids (350 m) → re-collapse; do NOT raise dedup tol to weld it |
| η<0.1 slivers all touching fault, mmg/flips can't fix | faces pinned on required fault (splice-locked) | fix the SURFACE: thin fault tris (q<0.1) force sliver tets — collapse/flip them in Stage B |
| slivers at z≈−15 m along the daylight trace | flat tets between flat top and shallowest fault row — the trace-band family | usually a FLOOR (see below); do NOT refine the shallow band (fills the wedge with more, equally flat tets) |
| one stubborn ~15 m edge, both endpoints on fault | strand-gap chord (strands pass that close) | geometry-forced; document, don't chase |
| fault edges > 500 after collapse | growth cap missing/too loose | cap 494.9 during collapse; relaxed round must be followed by LEB re-refine |
| "Inf/NaN in energies" at SeisSol init | daylight facets with σ_n=0/τ₀=0 (stress nc capped at z=0) or inverted tets | z≥0 stress fill from the z=−250 slice; reorient tets; `--validate` |
| MFEM "vertices indices are not unique" | msh v4 emitted | re-export msh22 |
| refinement loop mints failures as fast as it fixes (treadmill/oscillation) | metric anchored to CURRENT sizes (`min(h_max, target)`) — mmg's ±41 % acceptance band lets compliant cells drift over tolerance | ABSOLUTE design-anchored metric (see campaigns section); protection/rings only for sparse tails |
| mmg log: `N filtered` huge, `0 splitted`, mesh barely changes | insertion filter: split midpoints land within 0.7 metric units of existing cavity vertices (isolated long edges in fine surroundings) | protected punch at safety 0.7 (shrinks sol → midpoints clear the filter; freeze surroundings) |
| graded-envelope cells stuck at 1.25–1.4× spec forever | gradient shadow: census uses env(bary), mmg averages metric over the edge (~17 % looser on the slope) | evaluate metric envelope at `d − 0.6·env(d)` (pullback) |
| band/gate failures hug the fault despite fault-edge exemption | fault-vertex metric clamp bleeding coarse metric one shell out, or the pinned shell itself | drop the clamp under an absolute metric; count `fault-vertex-pinned` as a structural class |
| coarsening pass breaks a Vs-based gate at basin edges (worst f collapses) | nearest-grid Vs flips nodes as >1 km cells' barycenters wander between remeshes | budget on `minimum_filter`-pooled Vs over the wander radius + 25 % contrast screen; coarsen→keep→repair, never all-or-nothing rollback |
| an mmg regrade makes the gate WORSE (worst f drops), not better | mmg's ±41 % band licenses COARSENING of compliant cells | stop using mmg for this; switch to local LEB bisection (architecture 4) — it can only split |
| REFINEMENT plateaus: count falls a few %/round while tets climb fast | target anchored to nearest-grid Vs at the BARYCENTER; bisected children flip across a CVM bin edge into slow Vs | target min-Vs over the cell's own VERTICAL extent (lower bound on any child) — architecture 4 |
| LEB/LEPP finds few or no terminal edges; refinement barely moves | equal-length edge ties broken by LOCAL SLOT → chains cycle. Endemic on RED-refined meshes (red: TRIANGLE→4 similar children, TET→8 with only the 4 corners similar) | break ties by the globally unique edge key (strict total order) |
| LEB stops with "all terminal edges frozen (rim …)" and a small residual | the k-hop working patch is too tight for the last chains | reseed a fresh patch on just the residual with large `--hops` and run a second pass |
| a geometry pass swap-thrashes / OOMs the whole machine on a 10⁸-cell mesh | `P[T[:, EI]]` materialises an (n,6,3) float64 temporary — 17.6 GB at 122M tets | chunk every geometry pass (~2–4M tets/chunk); log peak RSS per round |

## Tricks that WORK (proven)

1. **Fix the surface, not the volume** — tetgen/mmg/flips never create edges
   shorter than the surface minimum; every volume-quality battle is won or
   lost in Stage B.
2. **Freeze-and-retry with tetgen -d as the oracle** (pymeshlab's
   self-intersection filter false-positives on benign non-manifold welds).
3. **Relaxed-cap → re-refine → strict re-collapse** (deep19km v2) unlocks
   collapses blocked by the max-edge cap.
4. **Vertex-removal fallback**: when the edge collapse flips a fan neighbor,
   remove the plain vertex into ANY neighbor instead.
5. **Thin-tri 2-2 flips before giving up on q<0.1 fault triangles**; accept a
   short new diagonal, then re-collapse it.
6. **mmg draw lottery**: 3–4 draws of stages mmg→flips→reorient, keep max
   η_min (score OUTSIDE bash — a broken $() parse once picked the wrong draw).
7. **Locate before fixing**: cluster bad tets/edges by position and by
   fault-vertex count; 2-fault-vertex flat tets ≠ 4-fault-vertex gap tets ≠
   trace fans — different cures.
8. **Byte-compare / structurally compare regenerated PUML** against the
   shipped one; conversions must be reproducible.
9. **When a spec must be met EXACTLY, use an operation that cannot undo itself.**
   LEB bisection only splits, so compliance is monotone and the loop provably
   converges; mmg's acceptance band can always re-mint. Reach for bisection the
   moment "0 failures" is the acceptance criterion.
10. **Separate the ACCEPTANCE metric from the DRIVING metric.** Measure the gate
    exactly as the locked convention says, but drive refinement with a
    conservative lower bound that is stable under the operation you are applying.
    Report both (the pooled target set ran ~5× the gate set — that is expected,
    not a defect).
11. **Smoke-test the refiner on a small mesh of the same family first** (60
    failing cells, 2 s). It caught a BC-propagation and orientation bug before
    burning 75 min on 122M cells, and after every rewrite it re-verified
    bit-identical output (+560 tets) in seconds.
12. **Keep the logs of the runs that FAILED** (`*_STALLED_*.log`). The plateau
    trajectories are what let you tell "wrong metric" from "wrong tie-break"
    apart — and in that case they proved the tie-break was NOT the cause.

## Tricks that DO NOT work (verified failures — don't repeat)

- **Shallow-band LEB refine to fix trace flat tets**: η_min got WORSE
  (0.029→0.019, η<0.05 38→154). The wedge is a thin REGION; subdividing fills
  it with more, equally flat tets.
  *Scope note (2026-08-01):* this is about using LEB to fix SHAPE in a thin
  wedge — it does not generalise. Using LEB to reduce EDGE LENGTH for a
  frequency gate left `η_min`, `η<0.05`, `η<0.1` and `min edge` **bit-identical**
  to the parent over +11.5M new tets. Refine to fix size, not to fix slivers.
- **Targeting a discontinuous sampled field directly** (nearest-grid Vs at the
  barycenter) in a refinement loop: bisected children flip across the bin edge
  into slow Vs and the loop treadmills (−1.8 %/round, +157k tets/round). Target
  a min-pooled LOWER BOUND instead.
- **Full 3-D min-pooling of Vs over a cell** as that lower bound: correct but
  wildly over-conservative on coarse cells (+55.8 % tets on the small ALT mesh
  vs +0.04 % for the raw rule). Pool along the SHORT grid axis only (z).
- **Judging a refinement campaign by the tet delta**: +9.4 % tets was +1.7 % run
  cost, because the new cells joined COARSE LTS clusters and `dt_min` did not
  move. Always quote the LTS-cost delta.
- **mmg size-map without -nosurf to fix fault-adjacent slivers**: overshoots
  the on-fault max edge (±40 % metric, 600 m > gate) and bloats 7×.
- **CGAL isotropic_remeshing at strand junctions**: moves vertices → PLC
  errors / thousands of sub-100 m edges at the triple junction. Use LEB
  bisection (never moves vertices) instead.
- **Snap-welding sub-parallel strands** (knot-tol): surfaces are sub-tet
  distance apart over km² — welding creates overlapping facets.
- **Raising vertex dedup tol to ~99 m** to "fix" crossings: welds real fault
  structure (input floor is ~100 m).
- **Chasing η_min below the structural floor**: if N independent mmg draws,
  a stricter tetgen (mindihedral 25 / minratio 1.7), and cavity remesh all
  land within ~5 %, it is geometry-forced. Report it with locations and stop.

## Acceptance & reporting

Always produce a `QUALITY_REPORT_<variant>.md` with: baseline-vs-new table
(min edge, edges<50/100/150, fault edge min/med/max, fault tri q, η_min/median,
η<0.05/0.1/0.3, tets, inverted), gate results, the deck-compat verification,
explicit *notes for every metric that got worse*, and the identified
geometry-forced floors with coordinates. Improvement claims without the
regression notes are not acceptable. Update the folder README with reproduce
commands and knobs, and leave `build_tmp/` logs + the final PLC for provenance.
