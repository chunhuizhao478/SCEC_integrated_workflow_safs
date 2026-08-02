"""bootstrap.py -- one call that makes a notebook runnable anywhere.

    from deckbuild.bootstrap import init
    B = init(project="safs_alt")
    B.cfg, B.root, B.out, B.data

Generalises the legacy `lib/nb_bootstrap.py:init()`.  Same behaviour, one marker change
(`deckbuild/__init__.py` instead of `pipeline.py`) and one addition (it loads the project
descriptor).  In order it:

  1. LOCATES the workflow folder by walking up from cwd for the marker; on Colab it also
     mounts Drive and searches it.  No absolute path is ever hardcoded.
  2. pip-installs the light dependencies on Colab.
  3. chdirs into the folder and puts it on sys.path.
  4. RELOADS every already-imported `deckbuild.*` module, so editing a module and re-running
     the cell always takes effect.  A stale Jupyter kernel silently running old code is the
     single most common way to "fix" something and see nothing change.
  5. ensures outputs/ exists.
  6. loads and validates the descriptor, failing loudly on a missing input.
  7. prints a one-screen summary.

Calling init() twice in one kernel is idempotent and still reloads.

CONSEQUENCE OF RELOADING, worth knowing: `importlib.reload` creates NEW class objects.
Any name bound before the reload -- including exception classes captured in an earlier
notebook cell -- no longer satisfies `is` or `except` against the new ones.  In a notebook
that is the accepted price of picking up edits (re-run the cell that imported them).  In a
test suite it causes cross-test contamination, which is why the reload test runs in a
subprocess.
"""
from __future__ import annotations

import importlib
import os
import sys
from dataclasses import dataclass
from pathlib import Path

__all__ = ["init", "Bootstrap", "in_colab", "find_workflow_dir"]

MARKER = Path("deckbuild") / "__init__.py"
_COLAB_DEPS = ("pyyaml", "netCDF4", "h5py", "pyproj", "meshio")


def in_colab() -> bool:
    try:
        import google.colab  # noqa: F401
        return True
    except ImportError:
        return False


def find_workflow_dir(start: str | Path | None = None) -> Path:
    """Walk up from `start` looking for the folder that contains deckbuild/__init__.py."""
    here = Path(start or os.getcwd()).resolve()
    for cand in (here, *here.parents):
        if (cand / MARKER).is_file():
            return cand
    # Also try the directory this file lives in -- covers `python -m` from elsewhere.
    pkg_parent = Path(__file__).resolve().parent.parent
    if (pkg_parent / MARKER).is_file():
        return pkg_parent
    raise FileNotFoundError(
        f"could not locate the workflow folder (no {MARKER} found walking up from {here}). "
        f"Run from inside the integrated_workflow folder, or pass root= explicitly.")


def _find_on_drive() -> Path | None:
    """On Colab, search the mounted Drive for the workflow folder."""
    import glob
    hits = sorted(glob.glob(f"/content/drive/**/{MARKER}", recursive=True))
    if not hits:
        return None
    if len(hits) > 1:
        print(f"  ! {len(hits)} candidate folders on Drive; using the first:")
        for h in hits:
            print(f"      {Path(h).parent.parent}")
    return Path(hits[0]).parent.parent


def _colab_setup() -> Path | None:
    import subprocess
    print("  - Colab detected: installing dependencies")
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", *_COLAB_DEPS], check=False)
    try:
        from google.colab import drive  # type: ignore
        if not Path("/content/drive").exists():
            drive.mount("/content/drive")
    except Exception as exc:  # pragma: no cover - Colab only
        print(f"  ! could not mount Drive: {exc}")
    return _find_on_drive()


def _reload_package() -> int:
    """Re-import every already-loaded deckbuild.* module so edits take effect."""
    names = [n for n in list(sys.modules) if n == "deckbuild" or n.startswith("deckbuild.")]
    # Reload the DEEPEST submodules first and the package itself LAST, so that when
    # `deckbuild/__init__.py` re-runs its `from deckbuild.config import ...` it binds the
    # freshly reloaded objects.  Reloading the package first would leave the top-level
    # re-exports pointing at the pre-edit classes -- exactly the stale-kernel failure this
    # function exists to prevent.
    names.sort(key=lambda n: (n.count("."), n), reverse=True)
    reloaded = 0
    for name in names:
        mod = sys.modules.get(name)
        if mod is None or getattr(mod, "__file__", None) is None:
            continue
        try:
            importlib.reload(mod)
            reloaded += 1
        except Exception as exc:
            print(f"  ! could not reload {name}: {exc}")
    return reloaded


@dataclass
class Bootstrap:
    """What a notebook needs: the config and the three directories."""

    cfg: object | None
    root: Path
    out: Path
    data: Path
    project: str = ""

    def __repr__(self) -> str:
        return f"Bootstrap(project={self.project!r}, root={self.root})"


def init(project: str = "safs_alt", require_files: bool = True,
         root: str | Path | None = None, chdir: bool = True,
         verbose: bool = True) -> Bootstrap:
    """Locate the workflow, put it on the path, reload it, and load `project`.

    project       -- descriptor stem under projects/ (or a path to a YAML).
    require_files -- forwarded to Project.load; False skips only existence checks.
    root          -- override the auto-located workflow folder.
    chdir         -- cd into the workflow folder (notebooks expect relative paths).
    """
    if root is not None:
        wf = Path(root).resolve()
        if not (wf / MARKER).is_file():
            raise FileNotFoundError(f"{wf} does not contain {MARKER}")
    else:
        wf = None
        if in_colab():
            wf = _colab_setup()
        if wf is None:
            wf = find_workflow_dir()

    if chdir:
        os.chdir(wf)
    if str(wf) not in sys.path:
        sys.path.insert(0, str(wf))

    n_reloaded = _reload_package()

    out = wf / "outputs"
    out.mkdir(parents=True, exist_ok=True)
    data = wf / "data"

    # Import AFTER the path is set up and the package reloaded.
    from deckbuild.config import Project

    cand = Path(project)
    desc = cand if cand.suffix in (".yaml", ".yml") else wf / "projects" / f"{project}.yaml"
    if not desc.is_file():
        available = sorted(p.stem for p in (wf / "projects").glob("*.yaml"))
        raise FileNotFoundError(
            f"no descriptor for project {project!r} at {desc}; available: {available}")

    cfg = Project.load(desc, require_files=require_files)

    b = Bootstrap(cfg=cfg, root=wf, out=out, data=data, project=cfg.name)
    if verbose:
        print(f"deckbuild bootstrap  (reloaded {n_reloaded} module(s))")
        print(cfg.summary())
        print(f"outputs      : {out}")
        if not require_files:
            print("  ! require_files=False -- input existence was NOT checked")
    return b
