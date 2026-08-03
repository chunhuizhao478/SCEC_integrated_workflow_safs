# deckbuild — codebase guide

A SeisSol dynamic-rupture deck needs four ingredients. Three are a well-defined
transformation from raw data to a NetCDF file and live in this package; the fourth — the
mesh — is a judgement-heavy craft and is handed to a skill. See `MESHING.md`.

## Module map

```
config.py       the project descriptor: every fault-specific number, one YAML per system.
                Strict loading -- an unknown key, an unknown reader kind, a reversed band
                or a missing file RAISES.  Mappings are read-only proxies.
contract.py     Stage / Gate / GateReport / Artifact / Manifest.  The three-verb contract:
                build writes, verify RE-OPENS THAT FILE and gates it, plot draws it.
bootstrap.py    init(): locate the folder, path, RELOAD every deckbuild.* module, load the
                descriptor.  Reloading replaces class objects -- see its docstring.
geometry.py     the fault frame: strike_s_km, load_fault, resolve_tractions,
                snap_hypocenter.  Ported verbatim; only the signatures changed.
asagi.py        the ASAGI NetCDF layer: write/read/trilinear_sample + the G1/G3 gates.
                Rejects a non-equidistant axis, which ASAGI does not check.
rawslices.py    horizontal-slice ASCII -> a uniform grid (the CVM/CTM two-stage recipe).
orientation.py  SHmax + R from a community model, a constant, or a callable.
material.py     velocity -> {rho, mu, lambda}; plasticity; the Sv(z) profile.  M1-M7.
stress.py       orientation x Sv x closure k -> {s_xx..s_xz}.  V1-V4.
friction.py     a(T)/V_w(T) or a depth profile -> {rs_a, rs_srW}; graded f_w -> LuaMap.
                F0/F1/FL, G1-G4.
mesh.py         .msh <-> .puml.h5 and gates A-G.  Does NOT build a mesh.
stage_f.py      the mesh<->deck seam.  Its predicates are SHARED with deck.py's P1-P8 so
                the two batteries cannot drift.
deck.py         assemble a runnable deck + the P1-P8 pre-flight + diff.
introspect.py   recover a descriptor from a SHIPPED deck.  Marks what it cannot infer.
plots.py        the two check figures: material slices + a fault-following section, and
                the on-fault (s, depth) maps.  Both read FILES already written, so what
                you look at is what SeisSol reads.  Three things it must do and any
                replacement must too: take the fault trace AT each slice depth (a 3-D
                surface moves with depth), scale the depth window and the section binning
                to the MESH (constants tuned on a 160k-facet fault silently blank out a
                192-facet one), and use PERCENTILE colour limits.
```

## The dependency chain (order is load-bearing)

```
material -> Sv(z) -----> stress
         -> thermal ---> friction
mesh ----------------->  every projection
                          \
                           -> deck (assemble + P1-P8)
```

Material runs **first**. The stress closure needs its `Sv(z)` and the friction zoning needs
its thermal field. Running out of order is the likeliest way to pair a new velocity model
with a stale stress field.

## Conventions

| | |
|:--|:--|
| `z` | elevation, positive up; depth is `-z`, **sea-level** referenced |
| stress in a nc | compression-**negative** Pa |
| stress on the fault | compression-**positive** MPa |
| ASAGI axes | float64, strictly increasing, **equidistant** |
| ASAGI data | one compound variable `data(z,y,x)`, float32 members |
| Gmsh | **v2.2 ASCII** only |
| failure | fail loud; no silent fallbacks |
| `verify()` | re-opens the written file; never checks an in-memory array |
| a skipped gate | is **never** a pass |

## Gate batteries, by stage

`M1–M7` material · `V1–V4` stress · `G1–G4` + `F0/F1/FL` friction · `A–G` mesh ·
`Stage F` mesh-vs-deck · `P1–P8` deck pre-flight. Names are preserved from the legacy
verifiers so existing run notes stay readable.

## Gotchas the code encodes

Each of these cost a real run, and each is asserted somewhere:

- a receiver at `z = 0` is **silently dropped** by SeisSol v1.1.3 — use `z = -1 m`
- a spatial `rs_muw` is **silently ignored** below v1.3.2 — the run completes and is wrong
- `bulkFriction` is `tan(phi)`, the **coefficient**; SeisSol `atan`s it
- a literal `0.0` in `OutputRegionBounds` writes the **whole volume**
- `wavefieldoutput = 0` **hangs** in teardown
- an **inverted tet** diverges at t ≈ 1 s with Inf/NaN bulk energy
- `sigma_n` or `tau_0` = 0 gives `psi_ini = -inf` and aborts at init
- kappa is **necessary, not sufficient** — a design passed the screen and still arrested
- the resolution gate is **non-monotone** under refinement
- `23.5e9` in YAML is a **string** (YAML 1.1 needs `23.5e+9`)

## Legacy tree decision (Phase 9)

- **`combined_workflow/` — KEEP, frozen.** Its post-processing notebooks (PGV maps, BBP
  export, on-fault comparison, Mw estimation) are out of scope here and still in use.
  `integrated_workflow/` is the maintained path for *building* a deck.
- **`stress_build_workflow/` — retire** once its eight step-notebooks are represented as
  the stress section's deep dive. Its `functions/` are already a stale copy of
  `combined_workflow/lib/`.

No file in this package imports from either tree. The byte-exact tests import the legacy
modules only as a *reference*, behind `skipif`.
