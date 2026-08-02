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
| E3 material | **PARTIAL** | grid reproduces **exactly**; values do **not** — see below |
| E4–E6 | **NOT RUN** | gated on closing E3 |

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

## E3 — the grid reproduces exactly; the values do not

```
              built                shipped
shape         452 x 361 x 194      452 x 361 x 194     IDENTICAL
x axis        bit-identical
y axis        bit-identical
z axis        bit-identical
Vs range      155 - 5352 m/s       155 - 5352 m/s      (matches the deck's own note)

mu    max |d| 1.31e+10 Pa   max rel 1.98    median rel 2.2e-05
rho   max |d| 199 kg/m^3    exactly equal on 0.53 % of nodes
```

**This is not yet reproduction, and the difference is not yet explained.**

Building the inscribed grid took two real fixes to get this far:

1. Sizing the grid from slice 0 alone put 7,624 of 181,541 nodes outside a later slice's
   hull. The grid must be inscribed across **all** slices.
2. More fundamentally, a lon/lat **rectangle** is a **curved quadrilateral** once
   projected, so a bbox-inscribed rectangle still pokes past the curved edges. Inscribing
   against the *edges* (largest western x, smallest eastern x, and likewise in y) is what
   reproduces the shipped 452 × 361 exactly.

That the axes now match bit-for-bit says the grid construction is right. The value
difference is a separate question with at least three candidates, none yet eliminated:

- **scipy/Qhull version.** The plan anticipates this ("might not match"). A median
  relative difference of 2.2e-5 is plausible for a different Delaunay triangulation of the
  same scattered points — but a **max of 1.98 (198 %) is not**, so this cannot be the
  whole story.
- **A different vertical resample.** The legacy converts to moduli at the source nodes and
  then z-resamples; the large excursions may be concentrated in the `extend_z_top` cap or
  at the deepest levels, which would point at the cap handling.
- **A different source revision.** The raw slices in the tree may post-date the shipped nc.
  E1 established the slices are *present and plausible*; it did not prove correspondence.

**Next step:** localise the excursions — are they at the cap, the edges, or scattered? A
histogram of relative error by z-level would separate the three candidates in one pass.
Until that is done, the workflow **cannot claim to reproduce a SAFS deck**.

## What this exercise did establish

- The deck is fully readable and its design fully recoverable — the prerequisite for
  E3–E6.
- The mesh ingests and gates cleanly at production scale (16 M tets, 34 s).
- The grid-construction half of the material rebuild is exactly right.
- Two real bugs in the inscribed-grid logic were found and fixed by attempting this, which
  is precisely what the exercise is for.
