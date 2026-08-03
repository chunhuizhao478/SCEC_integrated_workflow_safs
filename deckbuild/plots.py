"""plots.py -- the figures that let you SEE whether a deck is right.

Two families, both driven by files that were already written to disk, never by in-memory
state.  What you look at is therefore exactly what SeisSol will read.

    plot_material_slices()  horizontal Vs maps at chosen depths + a vertical section,
                            with the REAL fault outline (from the mesh) drawn on every one
    plot_onfault()          the on-fault maps in the fault's own (s, depth) frame

Design notes, learned from the legacy `graded_fw.plot_stress`:

* colour limits come from PERCENTILES, not min/max.  One bad facet at a domain corner
  otherwise flattens the whole field to a single colour -- the failure visible in the
  before-figures this module replaces.
* the fault is drawn from the MESH at the depth of the slice.  A fault is a 3-D surface,
  so its map-view trace moves with depth; drawing the z=0 trace on a 10 km slice quietly
  misplaces it by kilometres on a dipping or curved fault.
* named gate bands and the hypocentre are annotated on every panel, because "is the
  nucleation inside the band I think it is" is the question these plots exist to answer.
"""
from __future__ import annotations

import numpy as np

__all__ = ["plot_material_slices", "plot_onfault", "vs_from_material"]

# The depth-aware map trace already exists -- do not reimplement it here.
from deckbuild.geometry import fault_trace           # noqa: E402


# --------------------------------------------------------------------------- helpers


def vs_from_material(material_nc, fields=("mu", "rho")):
    """(x, y, z, Vs) from a material nc.  Vs = sqrt(mu/rho), the quantity that sets the
    resolvable frequency and the one worth looking at."""
    from deckbuild.asagi import read_asagi
    x, y, z, f, _ = read_asagi(material_nc, fields=list(fields))
    with np.errstate(divide="ignore", invalid="ignore"):
        vs = np.sqrt(f["mu"] / f["rho"])
    return np.asarray(x), np.asarray(y), np.asarray(z), vs      # vs is (nz, ny, nx)


def _clim(v, lo=2.0, hi=98.0):
    """Percentile colour limits.  Guards the degenerate all-equal case."""
    v = np.asarray(v)
    good = np.isfinite(v)
    if not good.any():
        return 0.0, 1.0
    a, b = np.nanpercentile(v[good], [lo, hi])
    if not np.isfinite(a) or not np.isfinite(b) or b <= a:
        a, b = float(np.nanmin(v[good])), float(np.nanmax(v[good]))
    return (a, b) if b > a else (a, a + 1.0)


def _is_effectively_constant(v, rtol=1e-6):
    """True when a field's spread is negligible against its own magnitude."""
    v = np.asarray(v, float)
    good = np.isfinite(v)
    if not good.any():
        return False
    m = np.abs(np.nanmean(v[good]))
    return bool(m > 0 and (np.nanmax(v[good]) - np.nanmin(v[good])) / m < rtol)


def _nearest_k(z, target_z):
    return int(np.argmin(np.abs(np.asarray(z) - target_z)))


def _section_nbin(fault, cfg, cap=240, floor=8):
    """How many along-strike bins the section can actually fill.

    Must follow the number of DISTINCT along-strike positions, not the facet count: the
    demo fault has 192 facets but only 16 columns, so a facet-count rule (192//8 = 24)
    asks for more bins than there are columns and leaves a third of them empty -- visible
    as white stripes through the section.
    """
    if cfg is None or getattr(cfg, "strike", None) is None:
        return cap
    s = fault.s_km(cfg.strike)
    n_col = int(np.unique(np.round(s, 1)).size)      # distinct s to the nearest 100 m
    return int(np.clip(n_col, floor, cap))


def _depth_tol_km(fault, floor_km=0.5):
    """A depth window wide enough to catch facets on THIS mesh.

    Half the median spacing of the facet depths, floored.  On a fine mesh this collapses
    to the floor; on a coarse one it opens up enough that a slice still finds the fault.
    """
    d = np.sort(np.unique(np.round(fault.depth_m / 1000.0, 3)))
    if len(d) < 2:
        return floor_km
    return float(max(floor_km, 0.5 * np.median(np.diff(d)) * 2.0))


def _annotate_bands(ax, cfg, ymax=None, label=True):
    """Shade the named gate bands along s and name them.

    The label is placed in AXES coordinates near the top, not at a data y -- a data-space
    y has to be recomputed per panel and lands outside the view (invisible) the moment a
    panel's y-limits differ.
    """
    if cfg is None or not getattr(cfg, "gate_bands", None):
        return
    for b in cfg.gate_bands:
        ax.axvspan(b.s_start_km, b.s_end_km, alpha=0.16, color="tab:orange", zorder=0)
        if label:
            ax.text(0.5 * (b.s_start_km + b.s_end_km), 0.94, b.name,
                    transform=ax.get_xaxis_transform(), ha="center", va="top",
                    fontsize=7.5, color="saddlebrown", zorder=7,
                    bbox=dict(fc="white", ec="none", alpha=.6, pad=1.2))


# ------------------------------------------------------------------ material figure
def plot_material_slices(material_nc, fault=None, *, depths_km=(0.0, 5.0, 10.0),
                         cfg=None, sv_profile=None, tol_km=None, cmap="viridis",
                         nbin=None, zlim_km=None, figsize=None):
    """Horizontal Vs slices at `depths_km` + a FAULT-FOLLOWING vertical section.

    The section is cut along the fault's own curved trace, not a straight chord, and its
    x-axis is the same along-strike `s` the on-fault figures use -- so the two read
    together: a low-Vs column here lines up with the facets it weakens there.

    A straight chord would be wrong on a curved fault.  It wanders off the trace, so the
    facets it captures form disconnected wedges rather than the continuous fault surface,
    which reads as a hole in the fault that is not there.

    Parameters
    ----------
    material_nc : path to the material nc (the file, not an in-memory field)
    fault       : a `geometry.Fault`, or None to omit the fault entirely
    depths_km   : the horizontal slices to cut, positive-down
    cfg         : the project descriptor -- supplies the strike frame and band names
    sv_profile  : optional .npz from MaterialStage; adds the Sv_eff(z) panel

    Every panel is a check you can fail by eye: a fault curve running off the coloured
    region means the grid does not cover the fault there -- what M3/V4hull test
    numerically.
    """
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    x, y, z, vs = vs_from_material(material_nc)
    depths_km = list(depths_km)
    # Both of these MUST scale with the mesh.  A fixed 0.5 km window and 240 bins are
    # right for a 160k-facet production fault and catastrophically wrong for a 192-facet
    # demo: the window captures no facet at all ("fault does not reach 5 km" on a fault
    # that plainly does), and 239 of 240 bins come back empty, leaving a section that is
    # mostly white.  Derive both from the fault itself.
    if fault is not None:
        if tol_km is None:
            tol_km = _depth_tol_km(fault)
        if nbin is None:
            nbin = _section_nbin(fault, cfg)
    else:
        tol_km, nbin = (tol_km or 0.5), (nbin or 240)
    ncol = len(depths_km)
    figsize = figsize or (4.9 * ncol, 9.0)
    # constrained_layout, NOT tight_layout: tight_layout cannot place colorbars added to
    # GridSpec axes and warns that results "might be incorrect" -- and they were.
    fig = plt.figure(figsize=figsize, layout="constrained")
    gs = GridSpec(2, ncol + 1, figure=fig, height_ratios=[1.0, 0.92],
                  width_ratios=[1.0] * ncol + [0.55])

    vmin, vmax = _clim(vs)          # ONE scale, so the slices are comparable

    # ---- row 1: horizontal slices, sharing a single colorbar ---------------------
    ims, row1 = None, []
    for j, dkm in enumerate(depths_km):
        ax = fig.add_subplot(gs[0, j])
        k = _nearest_k(z, -dkm * 1000.0)
        ims = ax.pcolormesh(x / 1e3, y / 1e3, vs[k], cmap=cmap, vmin=vmin, vmax=vmax,
                            shading="auto", rasterized=True)
        if fault is not None:
            fx, fy = fault_trace(fault, depth_km=dkm, tol_km=tol_km)
            if len(fx):
                ax.plot(fx / 1e3, fy / 1e3, ".", ms=1.0, color="k", alpha=.9)
                ax.plot([], [], "-", color="k", lw=2, label=f"fault @ {dkm:g} km")
                ax.legend(loc="upper right", fontsize=7, framealpha=.85)
            else:
                ax.text(.5, .03, f"fault does not reach {dkm:g} km",
                        transform=ax.transAxes, ha="center", fontsize=8, color="crimson")
        ax.set_title(f"Vs at {dkm:g} km depth   (z = {z[k]:+.0f} m)", fontsize=10)
        ax.set_xlabel("Easting (km)")
        ax.set_ylabel("Northing (km)" if j == 0 else "")
        if j:
            ax.set_yticklabels([])
        ax.set_aspect("equal")
        row1.append(ax)
    if ims is not None:
        # attach to the ROW, not to a spare grid cell: a colorbar given its own tall cell
        # floats detached at whatever height the cell happens to be.
        fig.colorbar(ims, ax=row1, label="Vs (m/s)", fraction=.040, pad=.012,
                     shrink=.92)

    # ---- row 2: the fault-following section --------------------------------------
    ax = fig.add_subplot(gs[1, :ncol])
    if fault is not None and cfg is not None:
        sbin, prof, ztop, zbot = _section_along_fault(fault, cfg, x, y, z, vs, nbin)
        im = ax.pcolormesh(sbin, np.asarray(z) / 1e3, prof, cmap=cmap, vmin=vmin,
                           vmax=vmax, shading="auto", rasterized=True)
        ax.plot(sbin, ztop, "-", color="k", lw=1.4, label="fault top / bottom")
        ax.plot(sbin, zbot, "-", color="k", lw=1.4)
        _annotate_bands(ax, cfg, ymax=float(np.nanmax(ztop)))
        ax.legend(loc="lower left", fontsize=8, framealpha=.85)
        ax.set_xlabel("along-strike s (km)  -- the same axis as the on-fault figures")
        ax.set_title("vertical section of Vs ALONG THE FAULT TRACE "
                     "(black = the fault's own top and bottom edges)", fontsize=10)
    else:
        j0 = len(y) // 2
        im = ax.pcolormesh(x / 1e3, np.asarray(z) / 1e3, vs[:, j0, :], cmap=cmap,
                           vmin=vmin, vmax=vmax, shading="auto", rasterized=True)
        ax.set_xlabel("Easting (km)")
        ax.set_title(f"vertical section of Vs at y = {y[j0] / 1e3:.0f} km "
                     f"(no fault supplied)", fontsize=10)
    for dkm in depths_km:
        ax.axhline(-dkm, color="w", lw=.9, ls="--", alpha=.75)
    # Crop to the fault plus a margin.  A 45 km-deep grid squeezes the 0-19 km the fault
    # actually occupies into the top third, which is the part you came here to read.
    if zlim_km is None and fault is not None:
        zbot_km = -float(fault.depth_m.max()) / 1e3
        zlim_km = (max(float(np.min(z)) / 1e3, 1.6 * zbot_km), float(np.max(z)) / 1e3)
    if zlim_km is not None:
        ax.set_ylim(*zlim_km)
    ax.set_ylabel("z (km, elevation)")
    fig.colorbar(im, ax=ax, label="Vs (m/s)", fraction=.030, pad=.012)

    if sv_profile is not None:
        axs = fig.add_subplot(gs[1, ncol])
        svp = np.load(sv_profile)
        axs.plot(svp["sv_eff_mpa"], svp["depth_m"] / 1e3, lw=1.8)
        axs.invert_yaxis()
        axs.set_xlabel("Sv_eff (MPa)")
        axs.set_ylabel("depth (km)")
        axs.set_title("vertical effective stress\n(SEA-LEVEL referenced)", fontsize=9)
        axs.grid(alpha=.3)

    fig.suptitle("MATERIAL CHECK -- the velocity model the deck will ship, "
                 "with the fault taken from the mesh", fontsize=12, fontweight="bold")
    return fig


def _section_along_fault(fault, cfg, x, y, z, vol, nbin=240):
    """Sample `vol` down the fault's own curved trace.

    Returns (s_km_bin_centres, (nz, nbin) values, z_top_km, z_bot_km).  Bins with no
    facets are NaN in every output, so a genuine gap in the fault stays a visible gap
    instead of being interpolated over.
    """
    s = fault.s_km(cfg.strike)
    edges = np.linspace(float(s.min()), float(s.max()), nbin + 1)
    ctr = 0.5 * (edges[:-1] + edges[1:])
    idx = np.clip(np.digitize(s, edges) - 1, 0, nbin - 1)

    prof = np.full((len(z), nbin), np.nan)
    ztop = np.full(nbin, np.nan)
    zbot = np.full(nbin, np.nan)
    for b in range(nbin):
        m = idx == b
        if not m.any():
            continue
        xc, yc = np.median(fault.cent[m, 0]), np.median(fault.cent[m, 1])
        ix = int(np.clip(np.searchsorted(x, xc), 0, len(x) - 1))
        iy = int(np.clip(np.searchsorted(y, yc), 0, len(y) - 1))
        prof[:, b] = vol[:, iy, ix]
        ztop[b] = fault.cent[m, 2].max() / 1e3
        zbot[b] = fault.cent[m, 2].min() / 1e3
    return ctr, prof, ztop, zbot


# ------------------------------------------------------------------ on-fault figure
def plot_onfault(fault, cfg, fields: dict, *, snap=None, f0=None, fw=None,
                 cmap="viridis", figsize=None, ncol=2, map_view=True):
    """On-fault maps in the fault's own (s, depth) frame, plus an optional map view.

    `fields` is {label: per-facet array}.  Beyond whatever you pass, two derived panels
    are added when the inputs allow it, because they are what actually decides a design:

        dtau_dyn = (mu_app - f_w) * sigma_n     the drop each patch can pay out
        S        = (f0*sigma_n - tau)/(tau - f_w*sigma_n)
                                                the strength parameter; higher = harder
                                                to break.  Capped at 10 for readability.

    Every panel gets percentile colour limits, the named gate bands, and the hypocentre.
    """
    import matplotlib.pyplot as plt

    s = fault.s_km(cfg.strike)
    d = fault.depth_m / 1000.0
    panels = [(k, np.asarray(v, float)) for k, v in fields.items()]

    sn = fields.get("sigma_n (MPa)")
    tau = fields.get("tau_0 (MPa)")
    mu = fields.get("mu_app")
    if sn is not None and mu is not None and fw is not None:
        panels.append((r"$\Delta\tau_{dyn}$ (MPa)",
                       (np.asarray(mu, float) - np.asarray(fw, float)) * np.asarray(sn, float)))
    if sn is not None and tau is not None and f0 is not None:
        fwv = np.zeros_like(np.asarray(sn, float)) if fw is None else np.asarray(fw, float)
        with np.errstate(divide="ignore", invalid="ignore"):
            S = (f0 * np.asarray(sn, float) - np.asarray(tau, float)) / \
                (np.asarray(tau, float) - fwv * np.asarray(sn, float))
        panels.append(("S (capped at 10)", np.clip(S, 0, 10)))

    # one (s, depth) panel per field, stacked full width -- a 270 km x 20 km fault is
    # extremely wide, so side-by-side columns squash it into an unreadable strip.
    nrow = len(panels) + (1 if map_view else 0)
    figsize = figsize or (13.0, 2.65 * nrow + 0.8)
    fig, ax = plt.subplots(nrow, 1, figsize=figsize, squeeze=False, layout="constrained")
    ax = ax.ravel()

    for a_, (label, v) in zip(ax, panels):
        vmin, vmax = _clim(v)
        sc = a_.scatter(s, d, c=v, s=3, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
        # label the colorbar, not the title: repeating the same string twice per panel
        # spends the reader's attention on nothing.
        cb = fig.colorbar(sc, ax=a_, label=label, fraction=.030, pad=.010)
        # A field that is constant to within roundoff gets rendered by matplotlib as an
        # offset like "1e-8 + 2.592592e-1", which reads as structure and is not.  Say
        # what it actually is instead -- on a uniform planar fault, mu_app really is
        # constant, and that is the useful statement.
        if _is_effectively_constant(v):
            cb.formatter.set_useOffset(False)
            cb.update_ticks()
            a_.text(.01, .06, f"{label} is CONSTANT to ~1 part in 1e6 "
                              f"({np.nanmean(v):.4g})", transform=a_.transAxes,
                    fontsize=7.5, color="0.25",
                    bbox=dict(fc="white", ec="0.7", alpha=.85, pad=1.5))
        _annotate_bands(a_, cfg)
        if snap is not None:
            a_.plot([snap.s_km], [snap.depth_m / 1000.0], "*", ms=14, mfc="crimson",
                    mec="w", mew=.9, zorder=6, label="hypocentre")
        a_.invert_yaxis()
        a_.set_ylabel("depth (km)")
        a_.margins(x=.01)
    ax[len(panels) - 1 if not map_view else len(panels) - 1].set_xlabel(
        "along-strike s (km)")
    if snap is not None:
        ax[0].legend(loc="lower right", fontsize=8, framealpha=.85)

    if map_view:
        a_ = ax[len(panels)]
        key, v = panels[0]
        vmin, vmax = _clim(v)
        sc = a_.scatter(fault.cent[:, 0] / 1e3, fault.cent[:, 1] / 1e3, c=v, s=3,
                        cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
        fig.colorbar(sc, ax=a_, label=key, fraction=.030, pad=.010)
        if snap is not None:
            a_.plot([snap.xyz[0] / 1e3], [snap.xyz[1] / 1e3], "*", ms=14, mfc="crimson",
                    mec="w", mew=.9, zorder=6)
        a_.set_aspect("equal")
        a_.set_xlabel("Easting (km)")
        a_.set_ylabel("Northing (km)")
        a_.set_title(f"map view of {key} -- the fault's real trace, "
                     f"so a kink here is a kink in the mesh", fontsize=9)

    fig.suptitle("ON-FAULT CHECK -- sampled from the WRITTEN nc at the real facets",
                 fontsize=12, fontweight="bold")
    return fig
