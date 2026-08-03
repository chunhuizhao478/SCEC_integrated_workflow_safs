# SAFS reproduction exercise — results (2026-08-02)

Target: `safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_intermediate_plast_phi30_40_gradedfw_k1p70_nwredM7p8_sefw0_attenuation_deep40km`
(683 MB mesh, 362 MB CVM), rebuilt from the real raw tree.

## Ladder

| Rung | Verdict | Evidence |
|:--|:--|:--|
| E0 mesh | **PASS** | 3,210,006 nodes / 15,979,903 tets / 320,560 BC-3 faces, all interior |
| E1 provenance | **PASS** | 37 CVM slices, 105 CTM slices, CSM csv present |
| E2 descriptor | **PASS** | 37 fields recovered, **0 unknown** |
| E2c freeze | **PASS** | measured on the FIELD: freeze OFF, agreeing with the filename |
| E3 material | **PARTIAL** | grid bit-exact; the port is bit-identical to the LEGACY code, so the value difference is the INPUT DATA |
| E4 stress | **EXACT** | the shipped stress nc reproduced BIT-IDENTICALLY, all 6 components |
| E5 friction | **PARTIAL** | `rs_a` bit-identical and the f_w Lua exact; `rs_srW` used a documented per-deck override |
| E6 deck | **PASS** | pre-flight **7/7 hard**, 0 warn, 0 skipped; plasticity BIT-IDENTICAL |

## E2 — the descriptor is fully recoverable

Every design parameter came back, including three that exist only as literals inside
emitted Lua:

| Recovered | Value | From |
|:--|:--|:--|
| strike frame | az 314.0, origin (606971, 3707270) | the `rs_muw` Lua coefficients |
| hypocentre | (604446.944, 3704576.3853, −10067.9819) | the `Tnuc_s` Lua |
| f_w design | `[0, 0, 0.045, 0.03, 0.06, 0.0175, 0.05]`, 6 windows | the `rs_muw` Lua |
| k | 1.70 | the stress nc filename |
| shallow freeze | **OFF** | the FIELD (levels above −500 m are distinct) |
| phi | 30 / 40 deg | `bulkFriction = tan(phi)` in the FIELD, not the yaml prose |
| grids | material 1500/250, stress 1000/250, friction 1500/200 | the nc axes |

The f_w design read out of the Lua matches the shipped yaml's own filename value for
value. The four grids are confirmed independent, as the exploration predicted.

## Reproduction against the shipped v4_0_0 ALT deck — RESULTS

Target deck: `~/Downloads/seisol_quakeworx/safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_
intermediate_plast_phi30_40_gradedfw_k1p70_nwredM7p8_sefw0_attenuation_deep40km`.

| shipped file | result |
|:--|:--|
| `safs_stress_andersonian_k1.7.nc` | **BIT-IDENTICAL** — all 6 components + 3 axes, 3,892,131 nodes, 0 differing |
| `safs_plasticity_phi30_40.nc` | **BIT-IDENTICAL** — `plastCo` and `bulkFriction`, 31,655,368 nodes, 0 differing |
| `safs_fault_rs_muw_*.yaml` | **EXACT** — same function, max abs diff 0.000e+00 over s in [-50, 500] km |
| `safs_friction_thermal_case1.nc` | `rs_a` **BIT-IDENTICAL** (5,459,354 nodes); `rs_srW` differs — see below |
| `safs_material_cvm.nc` | grid bit-exact; values differ — **cause proven, not our code** |

Stress and plasticity are built from the SHIPPED material nc, which isolates them from the
CVM input difference below.

### The CVM difference is the INPUT DATA, not the port

Two hypotheses recorded here earlier were **both wrong**, and the tests that killed them:

- ~~degenerate Delaunay triangulation~~ — **disproved**. Reversing the point order changes
  0 of 163,172 nodes, and the order-sensitive set does not overlap the differing set at
  all. The source is a complete 124 x 82 lattice and the CSV is already in the legacy's
  `meshgrid(..., indexing="ij")` order, so both sides build the *same* triangulation.
- ~~scipy / Qhull version~~ — **disproved** by the same test and by the one below.

What settles it: running the **legacy's own** `build_utm_grid` and `resample_to_utm`
(`toolbox/generate_velocity_nc_from_raw/generate_velocity_nc_from_raw.py`, the converter
named in the shipped nc's own attributes) on our staged CSVs reproduces **the same 4,054
differing nodes and the same max abs diff of 198.9** that our port produces. And comparing
our port against that legacy function directly gives **bit-identical float64 output, 0 of
163,172 nodes differing**.

So the port reproduces the legacy algorithm exactly. The staged CSVs are simply not the
files that built the shipped nc, despite identical filenames — those names are opaque IDs
(`CVM_1782322764534_h_data.csv`) and the shipped nc records a Google Drive `source_raw_dir`
while ours are local. The differing nodes are all low-density basin material,
1295-1944 kg/m^3, which is where a CVM revision would show and where the field has enough
curvature for a small input change to exceed float32.

**Consequence:** an exact CVM reproduction needs the original slice files, not a code fix.

### `rs_srW`: the shipped file used a documented OVERRIDE

`rs_a` is bit-identical, which proves the thermal field and `a(T)` match exactly. `V_w(T)`
is also identical between the two implementations — same constants (0.05, 1000, 350, 400),
same clip, same linear ramp.

The shipped nc says so itself: `profile = custom callables`, title
`"RSSRW friction (case1 + custom profile)"`. The legacy `pipeline.build_friction_nc`
(`pipeline.py:456`) takes `vw_fn`, `vw_patches` and `vw_ceiling=VW_CEILING_DEFAULT = 100.0`
m/s (`pipeline.py:421`). That ceiling is exactly what the comparison shows: at the
1,408,707 differing nodes ours is a flat 0.05 (T <= 350 degC) while the shipped ranges
0.05015 to **100** — a V_w that rises below 350 degC and is capped at 100, not the standard
profile.

That is a per-deck design override, not a defect in the general workflow. Reproducing it
requires the specific `vw_fn`/patch list that deck was built with, which is not recorded in
the nc. `FrictionStage.build` would need a `vw_fn=` hook to accept one.

## Three defects the end-to-end run found (2026-08-02)

Running the check notebook against the shipped deck surfaced three real problems. All are
fixed; the run now passes every gate.

1. **P5 had an inverted stress sign** (`deckbuild/deck.py`). The nc is compression-negative
   and SeisSol's Drucker-Prager limit `max(0, c·cos φ − sin φ·σ_m)` is written in that same
   convention, so σ_m < 0 at depth and the term *adds* strength. The gate negated into
   compression-positive first, making the yield limit *fall* with confinement, go negative,
   and get silenced by its own `taulim > 0` guard — leaving a thin spurious band at the
   crossover. Running the corrected predicate against the **shipped** deck is what proved
   it: 1565 spurious yielding points at a single z-level with a peak ratio of 8674 before,
   **0 after**, and our rebuild now agrees with the shipped deck's peak ratio (0.6355 vs
   0.6417). Pinned by `test_P5_reads_the_nc_in_its_own_compression_negative_convention`.
2. **Stress box `zmin` was -17000, should be -20000.** The descriptor copied the raw
   `ALT_BOX` constant and missed the pipeline's override two lines later
   (`STRESS_BOX["zmin"] = -20000.0`), which exists because the ALT fault bottoms at
   -19048.8 m. Gate V4hull caught it: 14,194 facets outside the hull.
3. **CVM cap was +100, should be +3250.** Confirmed three ways: every shipped v4_0_0 deck
   ships the identical 452×361×194 material nc topping at +3250; that nc's own global
   attribute records `extend_z_top = 3250.0`; and its levels above z = 0 are exact copies
   of the z = 0 level. The +100 came from an *archived v2* converter's default. The same
   defect, plus a missing thermal `z_max: 200.0` (silently giving 106 levels instead of the
   shipped 107), was fixed in the PREFERRED descriptor.

Only (2) and (3) were caught by gates. (1) was caught by running the gate against a file
that was already known-good — which is the argument for always checking a new gate against
production data before trusting it.

## What this exercise did establish

- The deck is fully readable and its design fully recoverable — the prerequisite for
  E3–E6.
- The mesh ingests and gates cleanly at production scale (16 M tets, 34 s).
- The grid-construction half of the material rebuild is exactly right.
- Two real bugs in the inscribed-grid logic were found and fixed by attempting this, which
  is precisely what the exercise is for.
