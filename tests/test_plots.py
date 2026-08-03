"""test_plots.py -- the figure helpers, and the two mesh-density bugs they shipped with.

Both regressions here were invisible on the 160,280-facet production fault and obvious on
the 192-facet demo.  A figure module tested against one mesh density is not tested.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from deckbuild.config import GridBox, StrikeFrame            # noqa: E402
from deckbuild.geometry import Fault, fault_trace            # noqa: E402
from deckbuild.plots import (                                # noqa: E402
    _depth_tol_km, _is_effectively_constant, _section_along_fault, plot_onfault,
)

STRIKE = StrikeFrame(azimuth_deg=0.0, origin_xy=(0.0, 0.0))


def _synthetic_fault(n_along: int, n_down: int, depth_max_m: float = 14000.0) -> Fault:
    """A vertical N-S planar fault as a regular (n_along x n_down) grid of facets."""
    sy = np.linspace(-20_000.0, 20_000.0, n_along)
    sz = np.linspace(-depth_max_m, -500.0, n_down)
    yy, zz = np.meshgrid(sy, sz, indexing="ij")
    cent = np.column_stack([np.zeros(yy.size), yy.ravel(), zz.ravel()])
    n = len(cent)
    return Fault(cent=cent,
                 normals=np.tile([1.0, 0.0, 0.0], (n, 1)),
                 strikes=np.tile([0.0, 1.0, 0.0], (n, 1)),
                 dips=np.tile([0.0, 0.0, -1.0], (n, 1)),
                 areas=np.full(n, 1.0e6), mesh_name="synthetic")


class _Cfg:
    """The two attributes the plot helpers actually read."""
    strike = STRIKE
    gate_bands = ()


# --------------------------------------------------------------- the depth tolerance
def test_depth_tolerance_adapts_to_mesh_density():
    """A fixed 0.5 km window is right for a fine mesh and wrong for a coarse one.

    On the 192-facet demo it caught no facet at all, so every slice reported "fault does
    not reach 5 km" for a fault spanning 0.83-14.17 km.
    """
    coarse = _synthetic_fault(16, 12)      # ~1.2 km between facet depths
    fine = _synthetic_fault(400, 400)      # ~34 m
    assert _depth_tol_km(coarse) > 0.5, "a coarse mesh must widen the window"
    assert _depth_tol_km(fine) == pytest.approx(0.5), "a fine mesh collapses to the floor"


@pytest.mark.parametrize("n_along,n_down", [(16, 12), (60, 40), (300, 120)])
def test_every_requested_depth_finds_facets_at_any_mesh_density(n_along, n_down):
    """The regression itself: a slice inside the fault's depth range must find facets."""
    f = _synthetic_fault(n_along, n_down)
    tol = _depth_tol_km(f)
    for dkm in (2.0, 5.0, 10.0):
        fx, _ = fault_trace(f, depth_km=dkm, tol_km=tol)
        assert len(fx) > 0, (
            f"no facets within {tol:.2f} km of {dkm} km on a {len(f)}-facet fault, whose "
            f"depths span {f.depth_m.min()/1e3:.2f}-{f.depth_m.max()/1e3:.2f} km")


# ------------------------------------------------------------------- the section bins
def test_section_bins_are_not_mostly_empty_on_a_coarse_fault():
    """240 fixed bins over 192 facets left 239 of them NaN -- a near-white section."""
    f = _synthetic_fault(16, 12)
    box = GridBox(xmin=-30_000., xmax=30_000., ymin=-30_000., ymax=30_000.,
                  zmin=-20_000., zmax=0., dx=2000., dz=1000.)
    gx = np.arange(box.xmin, box.xmax + box.dx, box.dx)
    gy = np.arange(box.ymin, box.ymax + box.dx, box.dx)
    gz = np.arange(box.zmin, box.zmax + box.dz, box.dz)
    vol = np.ones((len(gz), len(gy), len(gx)))

    from deckbuild.plots import _section_nbin
    nbin = _section_nbin(f, _Cfg())                    # the adaptive rule in plots.py
    assert nbin <= 16, ("the bin count must not exceed the fault's 16 distinct "
                        f"along-strike columns; got {nbin}")
    ctr, prof, ztop, zbot = _section_along_fault(f, _Cfg(), gx, gy, gz, vol, nbin)
    filled = np.mean(np.isfinite(ztop))
    assert filled > 0.8, f"only {filled:.0%} of section bins carry a facet"
    assert prof.shape == (len(gz), nbin)
    assert np.all(ztop[np.isfinite(ztop)] >= zbot[np.isfinite(zbot)])


# ----------------------------------------------------------------- constant detection
def test_constant_field_is_recognised_and_a_varying_one_is_not():
    """mu_app on a uniform planar fault really is constant; matplotlib renders that as
    an offset like `1e-8 + 2.592592e-1`, which reads as structure and is not."""
    assert _is_effectively_constant(np.full(500, 0.2592592) + 1e-12 * np.arange(500))
    assert not _is_effectively_constant(np.linspace(1.0, 2.0, 500))
    assert not _is_effectively_constant(np.full(10, np.nan))


# ------------------------------------------------------------------------- the figure
def test_plot_onfault_adds_the_derived_panels_only_when_it_can():
    """dtau_dyn needs f_w; S needs f0.  Without them the panels must be ABSENT, not
    silently computed from a guessed default."""
    f = _synthetic_fault(40, 20)
    sn = 10.0 + f.depth_m / 1e3
    fields = {"sigma_n (MPa)": sn, "tau_0 (MPa)": 0.3 * sn, "mu_app": np.full(len(f), .3)}

    bare = plot_onfault(f, _Cfg(), fields, map_view=False)
    with_f0 = plot_onfault(f, _Cfg(), fields, f0=0.6, map_view=False)
    with_fw = plot_onfault(f, _Cfg(), fields, f0=0.6, fw=np.full(len(f), .1),
                           map_view=False)
    n_bare = len([a for a in bare.axes if a.get_ylabel() == "depth (km)"])
    n_f0 = len([a for a in with_f0.axes if a.get_ylabel() == "depth (km)"])
    n_fw = len([a for a in with_fw.axes if a.get_ylabel() == "depth (km)"])
    assert n_bare == 3, "no derived panel is possible without f0 or f_w"
    assert n_f0 == 4, "f0 alone buys the S panel"
    assert n_fw == 5, "f0 + f_w buys both S and dtau_dyn"


def test_plot_onfault_map_view_is_optional_and_uses_real_coordinates():
    f = _synthetic_fault(40, 20)
    fields = {"sigma_n (MPa)": 10.0 + f.depth_m / 1e3}
    fig = plot_onfault(f, _Cfg(), fields, map_view=True)
    maps = [a for a in fig.axes if a.get_xlabel() == "Easting (km)"]
    assert len(maps) == 1
    xs, _ = maps[0].collections[0].get_offsets().T
    assert np.allclose(xs, f.cent[:, 0] / 1e3), "the map view must plot real coordinates"
