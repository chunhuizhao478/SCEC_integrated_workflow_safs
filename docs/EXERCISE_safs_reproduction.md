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
| E3 material | **PARTIAL** | grid reproduces **exactly**; values localised, not yet closed |
| E4 stress | **PASS** | V2 1.7e-07, V4hull/V4a/V4b all pass |
| E5 friction | **PASS** | G1–G4 pass; the emitted f_w Lua matches the shipped one to **0.000e+00** |
| E6 deck | **PASS** | pre-flight **7/7 hard**, 0 warn, 0 skipped |

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

## E3 — the grid reproduces exactly; the value difference is now localised

```
              built                shipped
shape         452 x 361 x 194      452 x 361 x 194     IDENTICAL
x / y / z     bit-identical
Vs range      155 - 5352 m/s       155 - 5352 m/s      (matches the deck's own note)

mu    max rel 1.98    median rel 2.2e-05    2.59 % of nodes differ by > 1 %
rho   exactly equal on 97.53 % of nodes
```

**Still not reproduction — but no longer unexplained.** Binning the relative error by
z-level splits it into two effects with completely different signatures.

### Effect 1 — a fixed ~2.5 % interior stripe, present at every depth

At a level that coincides with a source slice depth the median error is **exactly zero**,
yet 2.48 % of nodes still differ. Those nodes are **not** at the grid edge (median
distance to the nearest edge is 87 cells against 58 for the grid as a whole); they are
confined to a north–south band, x index 249–382 of 451, spanning almost the whole y range.

That is the signature of a **degenerate Delaunay triangulation**. Where the raw CVM is
sampled on a regular lattice, every group of four co-circular points can be split along
either diagonal; Qhull's choice is arbitrary and version-dependent, and the two linear
interpolants disagree *inside* those squares while agreeing *exactly at* the sample
points. That is precisely what is observed, and it also explains why the disagreement is
geographically confined — to the sub-region where the raw sampling is regular.

This effect is small (median 6.4e-04 in rho on the differing set) and is the "scipy/Qhull
version" candidate, now supported by evidence rather than assumed.

### Effect 2 — one anomalous level at z = -250 m

| z (m) | median rel | nodes > 1 % |
|--:|--:|--:|
| -1000, -750, -500 | 0.000e+00 | ~1,200 |
| **-250** | **1.32e-02** | **116,116 of 163,172 (71 %)** |
| 0 | 0.000e+00 | 3,177 |

The CVM slice depths are 0, 100, 200, 300, 400, 500, 750, 1000, 1250, 1500, 2500, … so on
a dz = 250 m output axis every level **is** a slice depth *except* z = -250, which falls
between the 200 m and 300 m slices. It is the only shallow level that is genuinely
interpolated in the vertical, and it is the only one that moves. The same pattern repeats
at depth: -1750/-2000/-2250 sit inside the 1500→2500 gap and are the next-worst levels,
while -2500 (a slice) returns to a zero median.

So Effect 2 is entirely in the **vertical resample between non-grid-aligned slices**. The
legacy order of operations was checked and matches ours — `velocities_to_moduli` at the
source nodes, then `resample_z` on the moduli, not the other way round — so the remaining
suspects are the bracketing itself and the depth→elevation sign handling at the shallow
slices. This is the concrete next step, and it is a small one.

### Ruled out

- **`extend_z_top` / the cap.** The 14 levels above z = 0 are exact copies of the z = 0
  level in *both* files, so they contribute nothing of their own; their max of 1.98 is
  inherited from z = 0. (Getting here did require fixing the cap — see below.)
- **The grid construction.** All three axes are bit-identical.

Building the inscribed grid took two real fixes to get this far:

1. Sizing the grid from slice 0 alone put 7,624 of 181,541 nodes outside a later slice's
   hull. The grid must be inscribed across **all** slices.
2. More fundamentally, a lon/lat **rectangle** is a **curved quadrilateral** once
   projected, so a bbox-inscribed rectangle still pokes past the curved edges. Inscribing
   against the *edges* (largest western x, smallest eastern x, and likewise in y) is what
   reproduces the shipped 452 × 361 exactly.

Until Effect 2 is closed the workflow **does not claim to reproduce a SAFS deck**.

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
