# SCEC integrated workflow — SeisSol dynamic-rupture deck builder

Build every input a SeisSol dynamic-rupture simulation needs, from **raw data**, for **any
fault system** — then use those freshly built files to determine the physical parameters,
and write a folder you can run.

> **Status: under construction.** Phases 0–7 complete; 8 partial, 9 partial, and the chain runs end to
> end: `python run_workflow.py --project demo_planar` builds material, plasticity, Sv,
> stress, friction and the `rs_muw` LuaMap from raw inputs, gating each stage. Open **`deck_workflow.ipynb`** and Run All, or use
> `run_workflow.py` headlessly.  Phase 8 (SAFS reproduction) is PARTIAL -- see
> `docs/EXERCISE_safs_reproduction.md`: the grid reproduces exactly, the values do not
> yet.  `CODEBASE_GUIDE.md` is the module map. Full plan: `docs/PLAN_integrated_workflow_notebook_2026-08-01.md`.

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
  tpv13.yaml          a 60-deg dipping NORMAL fault (the generality test)
docs/               the plan and the component exploration
skills/             the vendored mesh build-and-improve skill (see MESHING.md)
tools/              make_demo_mesh.py
tests/              pytest suite
outputs/            build artifacts (gitignored)
data/<project>/     raw inputs and caches (gitignored, except the demo mesh)
```

## Setup

You need `conda` or `mamba` ([miniforge](https://github.com/conda-forge/miniforge) is the
easiest way to get either). `mamba` resolves this environment in seconds where `conda` can
take minutes — use it if you have it.

```bash
git clone https://github.com/chunhuizhao478/SCEC_integrated_workflow_safs.git
cd SCEC_integrated_workflow_safs

mamba env create -f environment.yml     # or: conda env create -f environment.yml
conda activate deckbuild

# register this env as a Jupyter kernel -- REQUIRED, see below
python -m ipykernel install --user --name deckbuild --display-name "deckbuild"

pytest -q                               # expect: 312 passed
```

Then open the notebook and **select the kernel**:

```bash
jupyter lab deck_workflow.ipynb
```

> **Kernel → Change kernel → deckbuild.**
>
> Do not skip this. Creating the environment is not enough: JupyterLab only offers kernels
> that have been *registered*, and if you launch it from a different environment the
> notebook silently runs against that environment's packages instead — which looks exactly
> like "missing packages" no matter how many times you reinstall.

Prefer to stay in the terminal? The same chain runs headless, no kernel needed:

```bash
python run_workflow.py --project demo_planar
```

`demo_planar` is a synthetic example that needs **zero downloads** — it builds its own mesh
and runs the full seven-stage chain, so it is the fastest way to confirm a working install.

### Validated environment

Built from scratch and tested on macOS (darwin, arm64) on 2026-08-02:

| | |
|:--|:--|
| python | 3.11 |
| numpy / scipy | 1.26 / 1.11 |
| netCDF4 / h5py | 1.6 / 3.10 |
| meshio / pyproj | 5.3 / 3.6 |
| matplotlib | 3.10 |
| pytest | 9.x — 312 passed |

Two pins are deliberate and should not be relaxed casually:

- **`scipy=1.11.*`** — `scipy.spatial.Delaunay` decides the CVM interpolation, and a
  different triangulation of the same scattered points changes the rebuilt material field.
  See `docs/EXERCISE_safs_reproduction.md`.
- **`matplotlib=3.10.*`** — matplotlib 3.8 calls `IPython.core.pylabtools.backend2gui`,
  which IPython **removed in 9.16.0**. That pair raises `ImportError` the moment `pyplot`
  selects a backend, killing every cell that plots. 3.10 never calls it, so the IPython
  version stops mattering.

### Troubleshooting

| symptom | cause | fix |
|:--|:--|:--|
| `ModuleNotFoundError` for `netCDF4`, `pyproj`, `yaml`, … in the notebook | the notebook is on a different kernel | Kernel → Change kernel → **deckbuild** |
| `deckbuild` is not offered as a kernel | `ipykernel install` was never run | re-run the register step above |
| `ImportError: cannot import name 'backend2gui'` | `matplotlib` 3.8 with IPython ≥ 9.16 | you are on an old `environment.yml`; `git pull` and rebuild |
| a stage cannot find its raw inputs | descriptor paths are relative to `data/<project>/` | check `Project.load(..., require_files=True)` output — it names the missing file |

To rebuild the environment from scratch:

```bash
conda deactivate && mamba env remove -n deckbuild -y && mamba env create -f environment.yml
```

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
