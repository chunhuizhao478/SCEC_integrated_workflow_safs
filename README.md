# SCEC integrated workflow — SeisSol dynamic-rupture deck builder

Build every input a SeisSol dynamic-rupture simulation needs, from **raw data**, for **any
fault system** — then use those freshly built files to determine the physical parameters,
and write a folder you can run.

> **Status: under construction.** Phases 0–6 of 9 are implemented, and the chain runs end to
> end: `python run_workflow.py --project demo_planar` builds material, plasticity, Sv,
> stress, friction and the `rs_muw` LuaMap from raw inputs, gating each stage. `run_workflow.py` now assembles a runnable deck and
> runs the P1-P8 pre-flight.  The notebook is Phase 7 and the SAFS reproduction
> exercise Phase 8. Full plan: `docs/PLAN_integrated_workflow_notebook_2026-08-01.md`.

## The idea

A deck needs four ingredients. Three are a well-defined transformation from raw data to a
NetCDF file, so they live in code here. The fourth — the mesh — is a judgement-heavy craft,
so it is handed to a **skill document** you give to Claude.

| Ingredient | Built from | Owner |
|:--|:--|:--|
| material (+ plasticity, + Q) | raw CVM velocity slices | `deckbuild/material.py` **(done)** |
| stress | raw stress-orientation data + Sv from the material | `deckbuild/stress.py` **(done)** |
| friction | raw temperature slices, or a depth profile | `deckbuild/friction.py` **(done)** |
| **mesh** | your fault geometry | **`skills/code-mesh-build-improve`** — ingest/gating in `deckbuild/mesh.py` **(done)** |

Everything fault-specific — map projection, strike direction, grid boxes, named restraining
bends, hypocentre, which raw readers to use — lives in **one YAML per fault system** under
`projects/`. Retargeting the workflow means editing that file, not the code. A test
(`tests/test_no_hardcoded_constants.py`) enforces it.

## Layout

```
deckbuild/          the package
  config.py           the project descriptor: schema, loader, strict validation
  contract.py         Stage / Gate / GateReport / Artifact / Manifest
  bootstrap.py        one call that makes a notebook runnable anywhere
projects/           one descriptor per fault system
  safs_alt.yaml       San Andreas, ALT geometry
  safs_preferred.yaml San Andreas, PREFERRED geometry
  demo_planar.yaml    a synthetic example that runs with zero downloads
docs/               the plan and the component exploration
skills/             the vendored mesh build-and-improve skill (see MESHING.md)
tools/              make_demo_mesh.py
tests/              pytest suite
outputs/            build artifacts (gitignored)
data/<project>/     raw inputs and caches (gitignored, except the demo mesh)
```

## Quick start

```bash
conda env create -f environment.yml
conda activate deckbuild
pytest -q
python run_workflow.py --project demo_planar    # the chain, end to end
```

> **Validated on:** Python 3.13 / pytest 9.0.2 / PyYAML 6.0.3, via an existing conda env.
> `environment.yml` describes the *intended* pinned environment (Python 3.11) and has not
> itself been built yet — build it before relying on the pins, especially `scipy`, which
> Phase 8 depends on for reproducibility.

In a notebook or REPL:

```python
from deckbuild.bootstrap import init
B = init(project="demo_planar", require_files=False)
print(B.cfg.summary())
```

`init()` locates this folder by walking up for `deckbuild/__init__.py`, puts it on the path,
**reloads every already-imported `deckbuild.*` module** (so editing a module and re-running
the cell always takes effect), ensures `outputs/`, and loads the descriptor.

## Retargeting to your fault system

Copy `projects/demo_planar.yaml` and edit it. You must supply:

- **fault surface geometry** → a mesh, built with the skill in `skills/`
- **a velocity model** — raw slices (`kind: cvm_slices`) or a 1-D profile (`layered_1d`)
- **a stress orientation** — raw data (`csm_csv`) or a constant (`constant`)
- **a friction parameterisation** — a thermal model (`ctm_slices`) or a depth profile
- **a hypocentre** — in lon/lat/depth or projected coordinates; it is snapped onto the real
  fault, and the snap distance is reported
- **the strike frame and CRS**

You may keep every value under `physics:` — those are generic. You must decide one thing
yourself: **the tectonic regime.** The Andersonian closure assumes `sigma2` is vertical,
which is a strike-slip assumption. A normal or thrust fault needs a different stress recipe
(Phase 9 adds an `analytic` kind for exactly this).

## Conventions

- `z` is **elevation**, positive up; depth is `-z`.
- Stress inside a NetCDF is compression-**negative** Pa; on the fault it is
  compression-**positive** MPa.
- ASAGI axes must be float64, strictly increasing and **equidistant**.
- Gmsh output must be **v2.2 ASCII**.
- **Fail loud.** A missing input raises. There are no silent fallbacks.
- `verify()` re-opens the file that was written; it never checks an in-memory array.
