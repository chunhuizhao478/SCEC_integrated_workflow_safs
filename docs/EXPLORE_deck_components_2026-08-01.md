# SeisSol dynamic-rupture deck components — mesh / material / stress / friction

**Explored 2026-08-01.** Sources: `Downloads/seisol_quakeworx/` (30 shipped v4_0_0 decks),
`miniapps/seas/safs/seisol_quakeworx/` (v3 tree, `toolbox/`, `thermal/`, `document/`,
`v3_under_construction/toolbox/{stress_build_workflow,combined_workflow}`), and the 18
`project_7.0_{alternative,preferred}/meshing*` mesh pipelines.

Purpose: establish exactly what a deck **is**, which code produces each piece, and what is
generic physics vs. SAFS-specific hardcoding — so a standalone, reusable workflow can be
designed on top (companion: `PLAN_integrated_workflow_notebook_2026-08-01.md`).

---

## 1. Overview — what a deck is

A shipped, runnable SeisSol dynamic-rupture deck is **one folder** with five kinds of thing.
Measured on
`safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_intermediate_plast_phi30_40_gradedfw_k1p70_nwredM7p8_sefw0_attenuation_deep40km`:

| Kind | File | Size | Producer |
|---|---|---|---|
| **mesh** | `safalt_0d5Hz_p3_deep40km.puml.h5` | 683 MB | `meshing*/code/` chain → `msh_to_puml.py` |
| **material** | `safs_material_cvm.nc` | 362 MB | `toolbox/generate_velocity_nc_from_raw/` |
| ↳ derived | `safs_plasticity_phi30_40.nc` | 242 MB | `lib/build_plasticity_roten2014.py` |
| **stress** | `safs_stress_andersonian_k1.7.nc` | 89 MB | `pipeline.build_stress_nc` / `graded_k.build` |
| **friction** | `safs_friction_thermal_case1.nc` | 42 MB | `lib/build_friction_nc_thermal.py` + `friction_build.py` |
| **wiring** | `safs_material_cvm.yaml`, `safs_initial_stress.yaml`, `safs_fault.yaml`, `parameters.par` | KB | hand-written + generated `!LuaMap` blocks |
| **sampling** | `safs_pgv_receivers_100k.dat`, `safs_pickpoints.dat` | 3 MB | `lib/pgv_receiver_map.py`, `lib/onfault_receivers.py` |
| **launch** | `run_*.sbatch` | KB | hand-written per machine |

Every gridded field is an **ASAGI NetCDF**: COARDS NETCDF4, `x`/`y`/`z` float64 axes
(strictly increasing, **equidistant** — ASAGI requires this), one compound variable
`data(z,y,x)` whose members are the field names easi asks for. `z` is **elevation**
(positive up), CRS UTM 11N (EPSG:32611), metres.

The four grids are **independent** — they do not share axes:

```
material CVM   452 x 361 x 194   dx=dy=1500 m, dz=250 m   x[28.5,705] km  y[3543,4083] km  z[-45000,+3250] m
plasticity     same grid as CVM (derived pointwise)
stress         dx=dy=1000 m, dz=250 m   ALT box or PREF box (see below)
friction/therm dx=dy=1500 m, dz=200 m   z[-21000,+200] m
mesh           unstructured tets, its own bbox
```

Consistency between them is enforced only by **hull-containment guards**, not by construction.

---

## 2. Architecture — the two toolbox generations

Two generations of builder code; the newer one **vendors** the older:

```
toolbox/                                  GEN-1: per-task folders, CLI scripts
  generate_velocity_nc_from_raw/          CVM raw -> material nc
  generate_stress_nc_from_raw/            CSM raw -> stress nc  (the verified physics)
  on_fault_stress_projection[_csm]/       projection + gate analysis
  h5_to_vtu/, nc_to_vtu/, ground_motion/, hypocenter_facts/, generate_period/

v3_under_construction/toolbox/
  stress_build_workflow/                  GEN-2a: 8 numbered steps, self-contained, Colab-ready
    functions/{_common,step1..step8}.py   + vendored GEN-1 physics
    step*.ipynb + workflow.ipynb
    data/ (109 MB), artifacts/
  combined_workflow/                      GEN-2b: THE current production toolbox (3.5 GB, UNTRACKED)
    pipeline.py                           driver: build_stress_nc / build_friction_nc /
                                          check_gate / check_nucleation / run_pipeline
    lib/  (38 modules)                    every reused function, self-contained
    data/ (mesh_alt, mesh_preferred, safs_material_cvm.nc, safs_thermal_T.nc,
           CSM_orientation.csv, density_*.npz, pgv/)
    *.ipynb (14 notebooks)                one per task, each Run-All
```

`combined_workflow` is the right foundation. Its design contract is already what a reusable
workflow needs:

- **self-contained** — no absolute paths; `lib/nb_bootstrap.py:init()` walks up from cwd to
  find the folder containing `pipeline.py`, and on Colab mounts Drive and searches it;
- **fail-loud** — `pipeline.require_inputs()` hard-errors on any missing `data/` input; there
  are no silent fallbacks;
- **reload-on-run** — `init()` re-imports every already-imported `lib` module, so a stale
  Jupyter kernel cannot silently run old code;
- **build / verify / plot triple** — every design module exposes `build(...)`, `verify(...)`
  (a *named* gate battery), `plot_design(...)`, `plot_stress(...)`, `check_outputs(...)`.

That contract is the single most reusable thing in the tree; the plan preserves it verbatim.

---

## 3. Component 1 — MESH

### Data flow

```
CFM/GOCAD source           CFM_data_step/*.step  (19 fault strands, 500/1000/2000 m LOD)
  -> inp_extract_surfaces.py / ts_to_stl.py      -> per-strand STL
  -> corefine/ (CGAL)                            -> mutually conforming fault sheets
  -> clean / clip / nw-cut                       -> results/stl_*/
  -> run_z0cut_meshing.py (gmsh Python API)      -> tet .msh v2.2, tags 101/102/103/104
       * OCC box + embedded fault surfaces
       * size field (near-fault lc, far lc)
       * fault z=0 boundary edges embedded into the box top face
  -> mmg_optimize.py / refine_fault_edges.py     -> quality + metric refinement
  -> msh_to_puml.py                              -> *.puml.h5   (SeisSol PUML/HDF5)
  -> verify_*.py                                 -> acceptance gates A..E
```

`run_z0cut_meshing.py` exists because the naive path PLC-errors: the OCC box top face is
triangulated **without knowledge of the fault trace**, so fault edge segments slice through
top-face triangles → invalid PLC → empty volume. The fix is to build a discrete curve from
the fault's `z = 0` boundary edges (reusing the fault's own merged node tags) and
`gmsh.model.mesh.embed` it into the box top face.

### `msh_to_puml.py` — the format contract

Reverse-engineered from a PUMGen reference; no PUMGen needed on macOS.

- datasets: `geometry (Nnode,3) f8`, `connect (Ntet,4) u8` 0-based, `boundary (Ntet,) i4`,
  `group (Ntet,) i4 = 1`
- file attrs: `boundary-format='i32'`, `topology-format='geometric'`
- **BC packing**: byte *i* of `boundary` holds the BC code of tet face *i*:
  `boundary = sum_i code_i << (8*i)`
- **face order**: `f0={0,2,1} f1={0,1,3} f2={1,2,3} f3={0,3,2}`
- **tag → BC**: `101 fault→3` (dynamic rupture), `102 top→1` (free surface),
  `103 bottom→5`, `104 sides→5` (absorbing)
- each tet's local vertex order is normalised to **positive signed volume**
  (`orient_tets_positive`); gmsh does not guarantee orientation and inverted tets make
  SeisSol diverge (bulk energy → Inf/NaN at t ≈ 1 s)

### Acceptance gates (`verify_deep40km.py` and siblings)

```
A  bbox / node & tet counts
B  fault edge max < 500 m                                       HARD
C  BC round-trip: tagged faces == geometric hull + fault x2     HARD
D  inverted tets == 0
E  volume resolution gate f = Vs/dx >= --gate (0.6667 Hz)       HARD when --check-gate
   for a derived mesh the meaningful form is: introduce NO NEW failures vs the frozen parent
```

### State of the art

There is **no single mesh driver**. There are **18** `meshing*` folders, each a fork of the
previous for one campaign (`deep40km`, `deep40km_refine2/3`, `faultband`,
`faultband_cut30_refine1/2`, `faultband_cut_s30`, `0d5Hz`, `edgemin`, `deep19km`,
`faultband_refine1/2/3`, ...). They share `puml_io.py` / `msh_to_puml.py` / `puml_to_msh.py`
by **copy, not import**. This is the least reusable component by a wide margin.

---

## 4. Component 2 — MATERIAL (CVM) and its derivatives

### `generate_velocity_nc_from_raw.py` — raw ASCII slices → ASAGI nc

Two stages in one file.

**Stage 1**
- read every `CVM_*_h_data.csv` horizontal slice (one per depth; lon/lat regular grid;
  columns `lon,lat,vp,vs,density`)
- reproject EPSG:4326 → EPSG:32611 with pyproj
- build the **inscribed** axis-aligned UTM grid at `--grid-dx`, rounded *inward* so every
  node is inside the source hull
- per slice, per field: `scipy.interpolate.LinearNDInterpolator` onto the UTM grid;
  **NaN cells fail loud**
- flip depth → elevation and clone the surface slice `+--extend-z-top` metres up (covers
  topography / mesh pad)

**Stage 2**
- moduli **at the source nodes**: `mu = rho*Vs^2`, `lambda = rho*(Vp^2 - 2*Vs^2)`
  (fails if `lambda <= 0` anywhere)
- linear **z**-resample onto the uniform output axis — converting *before* resampling keeps
  the stored grid the exact piecewise-linear interpolant of node moduli
- write COARDS NETCDF4 `data(z,y,x){rho, mu, lambda}` float32
- fixed-seed round-trip self-check (trilinear NetCDF samples vs an independent evaluation of
  the source stack); non-zero exit on FAIL

Production parameters: `--grid-dx 1500 --extend-z-top 100 --z-min -45000 --z-max 100
--dz 250 --dtype float32`.

`generate_thermal_nc_from_raw.py` is the same recipe applied to the SCEC Community Thermal
Model (`CTM_*_h_data_final.csv`, columns `Lon,Lat,Temperature`), producing `data(z,y,x){T}`
at `dx=1500, dz=200, z[-21000,+200]`.

### `build_plasticity_roten2014.py` — derived pointwise on the CVM grid

Roten, Olsen, Day, Cui & Faeh (2014, GRL, doi:10.1002/2014GL059411) Drucker-Prager:

```
phi = 35 deg  where Vs <  2500 m/s        (production decks lower this to 30)
phi = 45 deg  where Vs >= 2500 m/s        (production decks lower this to 40)
c   = 1e-4 * mu                            (their cohesion model 3, eq. 5)
```

converted to what SeisSol actually queries (`src/Model/Plasticity.h` computes
`angularFriction = atan(bulkFriction)`):

```
bulkFriction = tan(phi)      -- the COEFFICIENT, not degrees
plastCo      = c  [Pa]
```

Because `c = 1e-4*mu` is **linear in mu**, ASAGI's linear interpolation of `plastCo` equals
`1e-4 x` the linearly interpolated `mu` SeisSol uses for the elastic material — the cohesion
an element sees is exactly consistent with its shear modulus. `bulkFriction` is the only
thresholded field.

Fluid pressure: the paper's yield includes `P_f`; SeisSol's has no explicit `P_f` term — but
the SAFS bulk loading stresses are **effective** (`Sv_eff = lithostat - hydrostatic Pp`), so
`sigma_m(effective) = tau_m + P_f` and SeisSol reproduces the paper's yield exactly.

### Q (attenuation)

Not a separate file — `Qs = 0.05*Vs`, `Qp = 2*Qs`, derived inside `safs_material_cvm.yaml`
from the CVM's own `mu`/`rho`. Requires a **viscoelastic build** of SeisSol
(`EQUATIONS=viscoelastic2, NUMBER_OF_MECHANISMS=3` — build-time, not a runtime flag).
`FreqCentral`/`FreqRatio` set the band over which Q is held ~constant; SeisSol fits the 3
mechanisms by an Emmerich & Korn (1987) least-squares solve at 5 log-spaced collocation
frequencies. **Outside the band Q → infinity, i.e. the medium goes elastic.**
`fitAttenuation` OVERWRITES `mu`/`lambda` with the unrelaxed moduli, so the CVM values are
interpreted as the moduli *at* `FreqCentral`.

---

## 5. Component 3 — STRESS

### Pipeline (`stress_build_workflow` steps 1–7, wired by `pipeline.build_stress_nc`)

```
step1_orientation   read CSM_orientation.csv (Yang & Hauksson YHSM-2013):
                    per-cell stress tensor -> SHmax azimuth (author col 10) + shape ratio R
                    (col 13).  NOT eigendecomposition -- the author azimuth is used directly.
step2_grid          build the (gx, gy, gz) box; interpolate SHmax(x,y) and R(x,y) onto it
step3_vertical_st.  Sv_eff(z) = lithostat - HYDROSTATIC Pp, from the lateral-mean density
                    profile (precomputed into a 6 KB density_profile.npz so the 360 MB CVM
                    need not ship)
step4_closure       magnitudes_C1(Sv_eff, R, k):
                      sig2 = Sv_eff                      (Anderson strike-slip: sigma2 vertical)
                      sig3 = Sv_eff / ((1-R)*k + R)
                      sig1 = k * sig3
step5_tensor        build_tensor_andersonian(az, s1, s2, s3): sigma1 -> e_H(SHmax),
                    sigma3 -> e_h, sigma2 -> vertical
step6_write_nc      assemble all z levels -> compound data{s_xx..s_xz}(z,y,x) float32,
                    compression-NEGATIVE Pa.  Options: daylight patch, shallow freeze, k_map
step7_project       trilinear-sample the nc at fault-facet centroids -> (sigma_n, tau, mu_app)
                    in the Tandem (strike, dip) basis, compression-POSITIVE MPa
```

`mu_app = tau/sigma_n` is **Sv-invariant** (a function of orientation, `R` and `k` only), so
`_common.onfault_mu_app` builds the Andersonian tensor with `Sv_eff = 1` for design scans.

### The two shallow treatments in `step6` (both are bug fixes, both matter)

1. **daylight patch** (`daylight_patch`, PREFERRED only) — the fault daylights to +2209 m
   under +3 km topography. The grid is extended to `z=+3000` and every `z>=0` level filled
   with the shallowest sub-surface slice (`z=-250 m`). Without it, 7146 daylighting facets
   had `sigma_n=0, tau_0=0` → `psi_ini = a*ln((2*sr0/V_ini)*sinh(tau_0/(a*sigma_n))) = -inf`
   → "Inf/NaN in energies" abort. The seismogenic band is byte-identical (`max|d| = 0 Pa`).
2. **shallow freeze** (`freeze_above_depth_m`, default 500 m) — every node shallower than
   this depth carries **its own column's** stress evaluated **at** that depth. `sigma_n` and
   `tau` freeze together (the C1 tensor scales with `Sv_eff` at fixed `k`), so `mu_app` and
   the lateral orientation/`k` structure are unchanged and **below the depth the field is
   bit-identical**. Columnwise analogue of the depth-independent `sigma_n` of TPV101–104 /
   BP5. Motivation: unfrozen, the trace band `z ∈ [-120, 0] m` sits at `sigma_n < 1 MPa`; in
   `PRE_k_2_5_CASE2` those cells slid at ~4 m/s for the whole run (446 m final slip — mesh
   quality was exonerated). At 500 m the trace carries `sigma_n ≈ 5.4–7.6 MPa`, touching only
   ~2.9 % of fault area.

### Graded `k` (`lib/graded_k.py`)

`k` is a **column property** like `R` and SHmax, so a per-column `k_map` broadcasts exactly
like `R` through `magnitudes_C1`. `graded_k.build(design=...)` supports any number of
along-strike regions joined by smoothstep windows in the along-strike coordinate `s`.

Acceptance battery `graded_k.verify` → **V1–V4**:

```
V1  regional bit-identity  graded nc == constant-k_i nc at every grid column of region i
V2  closure eigen-check    sigma1/sigma3 of the WRITTEN tensor == k_i
V3  freeze                 every level above the freeze depth == the freeze-depth slice
V4a gate crossing          windowed corridor kappa_bar over the gate band clears CROSSED
V4b NW starvation          >= min_starve_km below the BLOCKED level
V4c nucleation             hypocenter overstress passes
```

V1 needs one constant-`k` nc per region value as scratch; `check_outputs` quarantines those
so exactly ONE stress nc + ONE friction nc are left as deliverables.

---

## 6. Component 4 — FRICTION

Split across four layers.

### (a) The thermal nc — `build_friction_nc_thermal.py`

Reads the CTM grid `safs_thermal_T.nc` and applies piecewise-linear-in-T profiles nodewise
(`b = 0.019` constant, slope `m = 8.0e-5 /degC`):

```
CASE1 (3 zones)  a-b(T) = +0.004                for T <= 50
                        = +0.004 - m*(T-50)     for 50..150     (0 @ 100)
                        = -0.004                for 150..300
                        = -0.004 + m*(T-300)    for T >= 300    (0 @ 350)
CASE2 (2 zones)  a-b(T) = -0.004                for T <= 300
                        = -0.004 + m*(T-300)    for T >= 300
V_w(T) both      = 0.05 m/s for T <= 350; linear to 1000 m/s over 350..400; 1000 above
```

Output `data{rs_a, rs_srW}`. CASE2's trace is therefore **fully velocity-weakening** with
`V_w = 0.05` — which is why CASE1 stress/friction calibration does not transfer to CASE2.

Guards **G1–G4** (hard-fail unless noted):
```
G1  every fault vertex AND centroid inside the grid hull; per-axis margins reported
G2  rs_a > 0; rs_srW in [0.05, 1000]; re-checked after the float32 round-trip
G3  fixed-seed self-check at 1000 fault-bbox points, median tol 1e-5; max is WARN-only at
    5e-4 PROVIDED the excursion sits in a cell whose corner temperatures straddle an a(T)
    kink (a is piecewise linear in T, so trilinear-of-a == a-of-trilinear exactly on any
    single branch).  A >5e-4 excursion in a kink-free cell IS a hard fail.
G4  a at the reference hypocenter == 0.015 and V_w == 0.05, both cases
```
A unit self-test of the closed-form profiles at hand-computed T values runs **before** any
baking; assertion failure aborts.

### (b) Zone surgery — `lib/friction_build.py`

Post-processes the baked fields without touching the underlying physics:
- **`CASE2_SPECIAL`** — raises `rs_srW` to a ceiling (default 100 m/s) inside cos-tapered
  cylindrical patches near the free surface on the branch, so those cells stay on the `f_LV`
  branch (`f_ss(3.7 m/s) ≈ 0.54` → ~4 MPa strength) and **arrest**, instead of
  `f_ss ≈ 0.007` → ~0.06 MPa and 404 m of runaway slip. `a-b` is untouched (the trace stays
  velocity-weakening) and `V_w` is only ever *raised*, so the deep VS zone is unaffected.
- **`apply_vs_barrier`** — from a given `s_km` to the NW end, over full depth, make the fault
  velocity-**strengthening** with a cosine lead-in on the SE side, to arrest the rupture
  before it reaches the fault edge.

### (c) Graded `f_w` — `lib/graded_fw.py` (163 KB, the largest design module)

`f_w` (SeisSol `rs_muw`) is the strong-rate-weakening floor. Grading it limits the **dynamic
stress drop** past the gate and starves the front, while the stress nc stays at ONE constant
`k`. Ships as an `!LuaMap` in `safs_fault.yaml`: an along-strike `s` computed inline in Lua
from a fixed azimuth and origin, then N−1 smoothstep increments.

Acceptance battery `graded_fw.verify` → **F0–F6 + FL**:

```
F0  admissibility     0 <= f_w < f0 in every region; windows ordered, non-overlapping
F1  regional identity f_w == f_w_1 EXACTLY over s <= a_1 (hypocenter + nucleation patch)
F2  gate-side kappa   corridor kappa SE of the first transition bit-identical to the
                      constant-f_w_gate reference
F3  gate crossing     gate-band kappa_bar_min clears the CROSSED level of THIS mesh
F4  NW starvation     >= min_starve_km of windowed kappa below the BLOCKED level
F5  stress drop       dtau_dyn REDUCED over the raised ground (s >= b_1) vs the constant
                      baseline f_w = PRODUCTION_FW = 0.0
F6  nucleation        unchanged: f_w does not enter S_E = a*sigma_n*ln(V_dyn/V_init)
FL  the emitted Lua evaluates to the design values
```

**Version gate**: spatial `rs_muw` needs SeisSol **> v1.3.2** (post-master). v1.1.3 accepts
the yaml and **silently ignores** it. Section [0] of the notebook checks the source tree
before anything is built.

### (d) Scalars

`rs_b = 0.019`, `rs_sl0 = 0.10`, `f0`, `V_init` are `!ConstantMap` / `parameters.par`
scalars. SeisSol v1.1.3 FL=103 **never reads `rs_b` spatially** — that is why `b` is constant
in the thermal build. Nucleation ships as a second `!LuaMap` (`Tnuc_s`) using the SCEC
compact bell `F(r) = exp(r^2/(r^2 - R^2))` for `r < R`, 0 otherwise.

---

## 7. Cross-component checks (what makes a deck shippable)

These live in `pipeline.py` and read the **baked files**, not in-memory arrays — the right
design, because it tests what SeisSol will actually see.

```
sample_stress_on_fault(stress_nc, fault)     trilinear -> (sigma_n, tau, mu_app) MPa, comp-POSITIVE
sample_friction_on_fault(friction_nc, fault) trilinear -> (a, V_w)
check_gate(...)     kappa = G/Gc over the VW seismogenic band; corridor kappa_bar(s) with a
                    L_COAST_KM = 6 km coasting moving average vs kappa_c; longest contiguous
                    sub-kappa_c run > GAP_MAX_KM = 4 km  =>  BLOCKING gate.
                    Named restraining bends: GATE_BANDS = [(28,72,'San Gorgonio'),
                    (334,380,'NW Big Bend')]
check_nucleation(...) barrier S_E = a*sigma_n*ln(V_dyn/V_init) at the hypocenter vs the chosen
                    overstress amplitude; L_b, L_nuc, and whether the mesh RESOLVES L_b
                    (>= 4-5 elements).
                      LSW:   L_b = mu*Dc/((mu_s-mu_d)*sigma_n)
                             L_nuc = 1.158*mu/(1-nu) * Dc/((mu_s-mu_d)*sigma_n)
                      RSSRW: L_b = mu*Dc/(b*sigma_n),  L_nuc = mu*Dc/((b-a)*sigma_n)
check_yield.py      (deck-local) 0 of N stress-grid points may yield at t=0 under the
                    plasticity fields; reports grid-wide max tau/taulim and the binding point
initial_state.py    psi_ini = a*ln[(2*RS_sr0/V_ini)*sinh(tau_0/(a*sigma_n))], reading every
                    input from the deck's live `file:` targets; checks the identity
                    mu_ini == tau_0/sigma_n, negative/-inf state, and exp overflow
mw_estimate.py      design-time Mw by rescaling a completed run's measured slip field with
                    this design's stress-drop ratio (calibrate="auto")
```

---

## 8. Conventions

| Thing | Convention |
|---|---|
| CRS | UTM 11N EPSG:32611, metres, `x`=Easting `y`=Northing `z`=Up |
| Elevation | `z` is elevation, positive up; depth = `-z` |
| Stress in the nc | compression-**NEGATIVE** Pa (SeisSol/easi) |
| Stress on the fault | compression-**POSITIVE** MPa (projection modules) |
| Fault basis | Tandem `(strike, dip)`; normals harmonised into the SW half-space via a strike hint |
| Along-strike `s` | `s_km` from a fixed origin along azimuth `STRIKE_AZ` |
| ASAGI axes | float64, strictly increasing, **equidistant** (hard requirement) |
| ASAGI data | one compound variable `data(z,y,x)`, float32 members |
| Gmsh | **v2.2 ASCII only** (`-format msh22`); MFEM's reader has no v4 branch |
| Failure mode | fail loud; no silent fallbacks; `!ConstantMap` only as far-field graceful degradation |
| Artifacts | `build()` writes to `outputs/`, never into a shipped deck folder, unless told to |

---

## 9. Gotchas

These are the non-obvious things that cost real runs. Each is a design constraint on the new
workflow.

1. **`z = 0` receivers are dropped** by SeisSol v1.1.3 — place them at `z = -1 m`. Dev
   checkouts falsely pass.
2. **A literal `0.0` in `OutputRegionBounds`** writes the whole volume. Bounds are
   mesh-specific (ALT and PREF differ).
3. **`wavefieldoutput=1` is required** to avoid a teardown deadlock (0 hangs in 25 runs
   without it); a healthy teardown is 0 s.
4. **Pickpoint rename race**: a point exactly on a partition boundary is claimed by two ranks
   and aborts the run. Delete that point.
5. **Fault output cell order is partition-dependent** — lexsort-align before ANY run-vs-run
   differencing. Frames compare along-strike, not by XY.
6. **`printtimeinterval` counts DR slots, not seconds** under LTS.
7. **Spatial `rs_muw` needs SeisSol > v1.3.2.** v1.1.3 ignores it silently — the run
   completes and is wrong.
8. **v1.1.3 FL=103 spatial fields are `a` / `srW` / `sl0` only.** `rs_b` is never spatial.
9. **`bulkFriction` is `tan(phi)`, not degrees.** SeisSol atans it internally.
10. **The unfrozen shallow stress field is a runaway-slip generator** (446 m of slip in
    `PRE_k_2_5_CASE2`). Use `freeze_above_depth_m >= 500`.
11. **`psi_ini` is `-inf` when `tau_0 = 0`** — hence the daylight patch and the RSSRW
    `min_depth = 100 m`.
12. **`mu_ini` is ALWAYS `tau/sigma_n`**, whatever the state variable says.
13. **`GRAPH_PARTITIONING_LIBS=ParMETIS`** (mixed case) silently degrades to *none* — an 11x
    slowdown. Use a `-pm` build.
14. **`srun --mpi=pmi2` pins all ranks to one core set** on Frontera (9.7x slower). Use
    `ibrun` with `module load impi`.
15. **The `Vs = 2500 m/s` threshold field** is the one discontinuous input; ASAGI smooths it
    over a single cell (1500 m horizontal / 250 m vertical). Benign — but it means the
    plasticity nc is *not* reproducible by thresholding an interpolated `Vs`.
16. **Elements are barycenter-sampled** for plasticity — v1.1.3 logs "Material Averaging is
    not implemented for plastic materials".
17. **The resolution gate is non-monotone under refinement** — a 4x `Vs` step at CVM
    `z = -250 m` plus nearest-grid barycenters means refining can *create* gate failures.
18. **`msh_to_puml`'s packed int64 triangle key overflows above ~2.1M nodes.** Use a
    structured-view `searchsorted`, not a `Counter` (which needs >12 GB).
19. **Inverted tets** must be normalised to positive signed volume or SeisSol diverges at
    t ≈ 1 s.
20. **kappa is necessary, not sufficient.** An `f_w = 0.10` RSSRW design passed the kappa
    screen and still arrested at `s ≈ 60 km` on Frontera. Only the run decides crossing.
21. **`build()` and `verify()` must read the same file.** Several past defects were
    verify-on-memory / ship-a-different-array. Every gate above re-opens the nc.

---

## 10. What is SAFS-specific (the generalisation surface)

Everything below is hardcoded today and must become configuration for a second fault system.

| Constant | Value | Where |
|---|---|---|
| CRS | `EPSG:32611` | `bbp_export.py:82`, `build_friction_nc_thermal.py:480`, `build_plasticity_roten2014.py:138`, both raw generators |
| Strike azimuth | `314.0 deg` (N46W) | `_common.py:126 STRIKE_HINT_AZ`, `pipeline.STRIKE_AZ`, `bbp_export.py:88` |
| `s` origin | `(606971.0, 3707270.0)` | `pipeline.GATE_ORIGIN_XY`, `bbp_export.py:89`, inline in emitted Lua |
| Grid boxes | `ALT_BOX`, `PREF_BOX` literals | `_common.py:59-62` |
| Named gates | `[(28,72,'San Gorgonio'), (334,380,'NW Big Bend')]` | `_common.py:140 GATE_BANDS` |
| Mesh registry | `mesh_alt` / `mesh_preferred_rvfix4` | `_common.py:53-54`, `pipeline.MESHES` |
| Hypocenter | per-mesh literals + Lua-embedded coordinates | `pipeline.mesh_config`, `safs_fault.yaml` |
| Raw formats | CSM `orientation.csv` cols 10/13; CVM `lon,lat,vp,vs,density`; CTM `Lon,Lat,Temperature` | the three raw readers |
| BC tags | `101/102/103/104 -> 3/1/5/5` | `msh_to_puml.py` |
| Friction constants | `b=0.019`, `sl0=0.10`, `m=8.0e-5`, T kinks 50/100/150/300/350/400 | `build_friction_nc_thermal.py` |
| Physics defaults | `MU_SHEAR=23.5e9`, `W_ENERGY=10e3`, `DC=1.2`, `KAPPA_C=0.9`, `SEIS_BAND_KM=(3,12)`, `ACTIVE_BAND_KM=(0.3,15)`, `L_COAST_KM=6`, `GAP_MAX_KM=4`, `SN_FLOOR_MPA=30` | `_common.py:124-140` |

By contrast, these are **already generic** and need no change: the ASAGI writer, the
trilinear sampler, `magnitudes_C1`, `build_tensor_andersonian`, `resolve_tractions`,
`tandem_basis`, `harmonise_normals`, the `L_b` / `L_nuc` / `kappa` / `S_E` formulas, the
Roten-2014 derivation, `msh_to_puml`'s format layer, and the whole
build/verify/plot/check_outputs contract.

---

## 11. Open questions

1. **Mesh generalisation depth.** Should the new workflow own mesh *generation* (gmsh box +
   embedded fault from an arbitrary surface), or only mesh *ingestion + gating* (accept a
   `.msh`/`.puml.h5` built elsewhere, verify it, report `f = Vs/dx`)? The 18 forked
   `meshing*` folders suggest generation is campaign-specific and resists a single API.
   **Recommendation in the plan: ingest + gate + one reference generator, not a universal one.**
2. **Which stress recipe is "the" recipe** for a new fault system? The Andersonian C1 closure
   assumes Anderson strike-slip (`sigma2` vertical). A normal- or thrust-regime fault needs a
   different axis assignment, and the CSM orientation csv is a California product.
   **Partially answered by the tree itself:** `seisol_quakeworx/tpv13/` is a complete SCEC
   TPV12/13 deck — a 60° dipping *normal* fault — and it does not use the Andersonian recipe
   at all. Its stress is an analytic depth-dependent `!LuaMap`
   (`s_zz = -9.8*(2700-1000)*depth`, with `s_xx`/`s_yy` scaled by `s3_to_s1 = 0.3496` above
   11951.15 m and isotropic below), and its material and plasticity are `!ConstantMap`. So a
   second, non-strike-slip stress recipe already ships here, with its own mesh
   (`tpv13_training.puml.h5`), and it is the natural generality test for any reusable
   workflow.
3. **Demo dataset.** `combined_workflow/data/` is 562 MB and untracked (the whole folder is 3.5 GB). What ships with a
   shareable notebook — a downsampled SAFS subset, a synthetic planar-fault example, or a
   download script?
4. **Does `stress_build_workflow` stay?** Its 8 step-notebooks are pedagogically the clearest
   artifact in the tree, but its `functions/` are an older copy of what `combined_workflow/lib/`
   now has (`step6_write_nc.py` 6.5 KB vs 9.3 KB — the newer one has `k_map` and the freeze).
   Merging them is a prerequisite for a single integrated notebook.
5. **Thermal dependency.** The friction build is CTM-driven. For a fault system with no
   community thermal model, is there a depth-parameterised fallback for `a(z)` / `V_w(z)`?
   The spatial MFEM driver already has a `depth_profile` mode; the SeisSol side does not.
