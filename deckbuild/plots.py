"""plots.py -- the check figures, in the SAFS workflow house style.

Every figure here is driven by a FILE already written to disk, never by in-memory state,
so what you look at is exactly what SeisSol will read.

    plot_material_slices()    2 x 3 horizontal depth slices of the velocity model
    plot_plasticity_slices()  the same grid for phi and cohesion (opt-in physics)
    plot_material_profiles()  slice statistics vs depth + Sv_eff(z) + the fault section
    plot_onfault()            the (s, z) fault maps, lettered, with derived panels
    onfault_report()          the per-region text diagnostics that go under it

HOUSE STYLE -- taken from the legacy toolbox (`lib/graded_fw.py::plot_design` and
`plot_stress`, and `plasticity_roten2014.ipynb`) so these figures sit alongside those
without looking foreign:

  * on-fault panels plot **z (km), negative down** -- NOT an inverted positive-depth axis
  * the hypocentre is a **lime star with a black edge**: `"*", ms=14, mfc="lime", mec="k"`
  * every along-strike region boundary is `axvline(color="tab:red", ls="--", lw=1.0)`;
    named gate bands are `axvspan(alpha=0.15, color="tab:orange")`
  * depth slices are `imshow(origin="lower", extent=[km], aspect="equal",
    interpolation="nearest")` with the **fault footprint** over them as `'.', ms=0.4,
    color='0.35', alpha=0.5`
  * each slice title carries a **statistic**, not just its depth
  * one **shared** colorbar per slice grid: `fig.colorbar(im, ax=axes, shrink=0.7)`
  * a **bold, wrapped** suptitle naming the design
  * panels are lettered (a), (b), (c)... and titled with a sentence saying what to look for
  * colour limits are PERCENTILES -- one bad facet at a domain corner otherwise flattens
    the whole field to a single colour
"""
from __future__ import annotations

import numpy as np

from deckbuild.geometry import fault_trace           # depth-aware; do not reimplement

__all__ = ["plot_fw_and_kappa", "plot_material_slices", "plot_plasticity_slices", "plot_material_profiles",
           "plot_onfault", "onfault_report", "vs_from_material"]

HYPO_KW = dict(marker="*", ms=14, mfc="lime", mec="k", mew=.8, ls="none", zorder=6)
BOUND_KW = dict(color="tab:red", ls="--", lw=1.0)
BAND_KW = dict(alpha=0.15, color="tab:orange")
FOOTPRINT_KW = dict(marker=".", ms=0.4, color="0.35", alpha=0.5, ls="none", zorder=3)
DEFAULT_SLICES_M = (0.0, -1000.0, -3000.0, -5000.0, -10000.0, -15000.0)


# --------------------------------------------------------------------------- helpers
def vs_from_material(material_nc, fields=("mu", "rho")):
    """(x, y, z, Vs) from a material nc.  Vs = sqrt(mu/rho) -- the quantity that sets the
    resolvable frequency, and the one worth looking at."""
    from deckbuild.asagi import read_asagi
    x, y, z, f, _ = read_asagi(material_nc, fields=list(fields))
    with np.errstate(divide="ignore", invalid="ignore"):
        vs = np.sqrt(f["mu"] / f["rho"])
    return np.asarray(x), np.asarray(y), np.asarray(z), vs      # vs is (nz, ny, nx)


def _clim(v, lo=2.0, hi=98.0):
    """Percentile colour limits, guarding the degenerate all-equal case."""
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


def _wrap_title(txt, fig_w_in, fontsize, min_chars=60):
    """Hard-wrap a long suptitle to the FIGURE width.  Returns (text, n_lines).

    Ported from the legacy toolbox because the failure mode is non-obvious: Jupyter's
    inline backend saves with `bbox_inches="tight"`, so a suptitle that overhangs the
    canvas EXPANDS the exported bitmap sideways and the panels end up as a thin strip down
    the middle of a very wide PNG -- the "why is this plot so small" report.
    """
    per_char = fontsize * 0.60 / 72.0
    budget = max(int(min_chars), int(float(fig_w_in) / per_char))
    lines, cur = [], ""
    for w in str(txt).split(" "):
        trial = f"{cur} {w}".strip()
        if len(trial) > budget and cur:
            lines.append(cur)
            cur = w
        else:
            cur = trial
    if cur:
        lines.append(cur)
    return "\n".join(lines), len(lines)


def _suptitle(fig, txt, figsize, fontsize=12):
    wrapped, n = _wrap_title(txt, figsize[0], fontsize)
    st = fig.suptitle(wrapped, fontsize=fontsize, fontweight="bold")
    st.set_wrap(True)
    return 1.0 - (n * fontsize * 1.45 / (figsize[1] * 72.0)) - 0.004


def _marker_area(ax, n, fill=0.55, lo=3.0, hi=600.0):
    """Scatter marker area (points^2) that TILES the axes for this many points.

    The legacy `s=3` is right for a 160,280-facet production fault and turns a 192-facet
    demo into a field of specks with white between them -- which reads as a broken plot,
    not as a coarse mesh.  Area per point is (axes area / n), scaled by `fill`.
    """
    fig = ax.get_figure()
    bb = ax.get_window_extent(renderer=fig.canvas.get_renderer()) \
        if fig.canvas.get_renderer is not None else None
    try:
        w_pt = bb.width * 72.0 / fig.dpi
        h_pt = bb.height * 72.0 / fig.dpi
    except Exception:                                   # no renderer yet -- estimate
        w_pt, h_pt = fig.get_size_inches()[0] * 72.0, fig.get_size_inches()[1] * 72.0 / 4
    return float(np.clip(fill * w_pt * h_pt / max(int(n), 1), lo, hi))


def _nearest_k(z, target_z):
    return int(np.argmin(np.abs(np.asarray(z) - target_z)))


def _depth_tol_km(fault, floor_km=0.5):
    """A depth window wide enough to catch facets on THIS mesh.

    A fixed 0.5 km is right for a 160k-facet production fault and catches nothing at all
    on a 192-facet demo -- which then reports "fault does not reach 5 km" for a fault that
    plainly does.
    """
    d = np.sort(np.unique(np.round(fault.depth_m / 1000.0, 3)))
    if len(d) < 2:
        return floor_km
    return float(max(floor_km, np.median(np.diff(d))))


def _section_nbin(fault, cfg, cap=240, floor=8):
    """How many along-strike bins the section can actually fill.

    Must follow the number of DISTINCT along-strike positions, not the facet count: the
    demo fault has 192 facets in only 16 columns, so a facet-count rule asks for more bins
    than there are columns and leaves a third of them empty -- white stripes.
    """
    if cfg is None or getattr(cfg, "strike", None) is None:
        return cap
    s = fault.s_km(cfg.strike)
    n_col = int(np.unique(np.round(s, 1)).size)
    # 3x the distinct columns: nearest-column sampling fills between them, so oversampling
    # only smooths, while undersampling throws real along-strike structure away.
    return int(np.clip(3 * n_col, floor, cap))


def _mark_bands(ax, cfg, label_first=True):
    """Shade named gate bands, one legend entry."""
    if cfg is None or not getattr(cfg, "gate_bands", None):
        return
    for i, b in enumerate(cfg.gate_bands):
        ax.axvspan(b.s_start_km, b.s_end_km,
                   label=f"{b.name} gate" if (i == 0 and label_first) else None, **BAND_KW)


def _design_bounds(fw_design):
    if fw_design is None:
        return []
    b = getattr(fw_design, "boundaries_s_km", None)
    return [tuple(w) for w in b] if b else []


def _mark_boundaries(ax, fw_design):
    """A tab:red dashed line at EVERY f_w transition edge."""
    for w0, w1 in _design_bounds(fw_design):
        for b in (w0, w1):
            ax.axvline(b, **BOUND_KW)


# ------------------------------------------------------ horizontal slices (2 x 3 grid)
def _slice_grid(x, y, z, vol, *, depths_m, fault, title_stat, cb_label, suptitle,
                cmap="viridis", norm=None, vmin=None, vmax=None, contour=None,
                contour_levels=(), figsize=(16.5, 9), tol_km=None):
    """The legacy 2 x 3 imshow slice grid: one shared colorbar, a statistic per title, the
    fault footprint over every panel."""
    import matplotlib.pyplot as plt

    depths_m = list(depths_m)
    nrow = 2 if len(depths_m) > 3 else 1
    ncol = int(np.ceil(len(depths_m) / nrow))
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize, constrained_layout=True,
                             squeeze=False)
    ext = [x[0] / 1e3, x[-1] / 1e3, y[0] / 1e3, y[-1] / 1e3]
    if fault is not None and tol_km is None:
        tol_km = _depth_tol_km(fault)

    im = None
    for ax, zm in zip(axes.ravel(), depths_m):
        k = _nearest_k(z, zm)
        kw = dict(cmap=cmap)
        if norm is None:
            kw.update(vmin=vmin, vmax=vmax)
        else:
            kw.update(norm=norm)
        im = ax.imshow(vol[k], origin="lower", extent=ext, interpolation="nearest",
                       aspect="equal", **kw)
        if (contour is not None and len(contour_levels)
                and not _is_effectively_constant(contour[k])):
            ax.contour(x / 1e3, y / 1e3, contour[k], levels=list(contour_levels),
                       colors="w", linewidths=0.7)
        if fault is not None:
            fx, fy = fault_trace(fault, depth_km=-z[k] / 1e3, tol_km=tol_km)
            if len(fx):
                # ms 0.4 is right for thousands of facets and invisible for tens
                kw = dict(FOOTPRINT_KW)
                kw["ms"] = float(np.clip(60.0 / max(np.sqrt(len(fx)), 1.0), 0.4, 4.0))
                ax.plot(fx / 1e3, fy / 1e3, **kw)
        if _is_effectively_constant(vol[k]):
            # A one-colour panel is indistinguishable from a broken one unless it says so.
            ax.text(.5, .5, "laterally uniform\n(1-D layered model)", ha="center",
                    va="center", transform=ax.transAxes, fontsize=11, color="w",
                    bbox=dict(fc="0.25", ec="none", alpha=.55, pad=4))
        ax.set_title(f"z = {z[k]:.0f} m   ({title_stat(vol[k])})", fontsize=10)
        ax.set_xlabel("x UTM [km]")
        ax.set_ylabel("y [km]")
    for ax in axes.ravel()[len(depths_m):]:
        ax.axis("off")
    cb = None
    if im is not None:
        cb = fig.colorbar(im, ax=axes, shrink=0.7)
        cb.set_label(cb_label)
    # wrap: an overhanging suptitle makes bbox_inches="tight" widen the export and squeeze
    # the panels into a strip down the middle (see _wrap_title)
    fig.suptitle(_wrap_title(suptitle, figsize[0], 12)[0], fontsize=12)
    fig._deckbuild_cb = cb          # so callers can retick a categorical colorbar
    return fig


def plot_material_slices(material_nc, fault=None, *, depths_m=DEFAULT_SLICES_M,
                         vs_threshold=2500.0, figsize=(16.5, 9), cmap="viridis"):
    """2 x 3 horizontal slices of Vs through the material nc, fault footprint on each.

    The white contour is `vs_threshold`, the Roten et al. (2014) soft/hard rock split --
    the same line that decides the plasticity friction angle, so this figure and the
    plasticity one overlay mentally.  Each title carries the slice median and the fraction
    below the threshold, because "how much soft rock is at this depth" is the question.
    """
    x, y, z, vs = vs_from_material(material_nc)
    vmin, vmax = _clim(vs)

    def stat(sl):
        return (f"median {np.nanmedian(sl):.0f} m/s, "
                f"{100.0 * np.nanmean(sl < vs_threshold):.1f}% < {vs_threshold:g}")

    lat_uniform = bool(np.nanmax(np.nanmax(vs, axis=(1, 2)) - np.nanmin(vs, axis=(1, 2)))
                       < 1e-9)
    sup = (f"MATERIAL nc -- $V_s$ on horizontal slices; white = {vs_threshold:g} m/s "
           f"contour (the Roten soft/hard split); grey = fault footprint at that depth")
    if lat_uniform:
        sup = ("MATERIAL nc -- this model is 1-D LAYERED, so every horizontal slice is one "
               "colour BY CONSTRUCTION.  Nothing is wrong with the figure; there is simply "
               "no lateral structure to show.  Use plot_material_profiles() instead, and a "
               "CVM reader (velocity.kind = cvm_slices) if you want 3-D structure.")
    return _slice_grid(
        x, y, z, vs, depths_m=depths_m, fault=fault, title_stat=stat,
        cb_label="$V_s$ [m/s]", cmap=cmap, vmin=vmin, vmax=vmax,
        contour=vs, contour_levels=(vs_threshold,), figsize=figsize, suptitle=sup)


def plot_plasticity_slices(plasticity_nc, fault=None, *, depths_m=DEFAULT_SLICES_M,
                           figsize=(16.5, 9)):
    """The same grid for the opt-in Drucker-Prager fields.  Returns [phi_fig, cohesion_fig].

    phi is drawn with a two-colour ListedColormap because it takes exactly two values --
    a continuous ramp would imply a gradient that is not there.  Cohesion is log-scaled
    because c = 1e-4*mu spans orders of magnitude between basin and basement.
    """
    from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm

    from deckbuild.asagi import read_asagi
    x, y, z, f, _ = read_asagi(plasticity_nc, fields=["bulkFriction", "plastCo"])
    phi = np.degrees(np.arctan(f["bulkFriction"]))
    coh = f["plastCo"]

    ang = np.unique(np.round(phi[np.isfinite(phi)], 4))
    lo, hi = float(ang.min()), float(ang.max())
    mid = 0.5 * (lo + hi)
    cmap = ListedColormap(["#e08214", "#542788"])
    norm = BoundaryNorm([lo - 1e-6, mid, hi + 1e-6], cmap.N)

    f_phi = _slice_grid(
        x, y, z, phi, depths_m=depths_m, fault=fault,
        title_stat=lambda sl: f"soft {100.0 * np.mean(sl < mid):.1f}%",
        cb_label=r"friction angle $\phi$ [deg]   (bulkFriction = $\tan\phi$)",
        cmap=cmap, norm=norm, figsize=figsize,
        suptitle=f"PLASTICITY nc -- Roten et al. (2014): $\\phi$ = {lo:g} deg (soft rock) "
                 f"/ {hi:g} deg (hard rock); grey = fault footprint at that depth")
    if getattr(f_phi, "_deckbuild_cb", None) is not None:
        f_phi._deckbuild_cb.set_ticks([0.5 * (lo + mid), 0.5 * (hi + mid)])
        f_phi._deckbuild_cb.set_ticklabels([f"{lo:g}", f"{hi:g}"])

    cmin = max(float(np.nanmin(coh)) / 1e6, 1e-3)
    f_coh = _slice_grid(
        x, y, z, coh / 1e6, depths_m=depths_m, fault=fault,
        title_stat=lambda sl: f"median {np.nanmedian(sl):.3f} MPa",
        cb_label="cohesion $c = 10^{-4}\\mu$   [MPa, log scale]",
        norm=LogNorm(vmin=cmin, vmax=float(np.nanmax(coh)) / 1e6), figsize=figsize,
        suptitle="PLASTICITY nc -- Roten et al. (2014) cohesion model 3 (their eq. 5): "
                 "$c = 10^{-4}\\mu$, so cohesion is LINEAR in the shear modulus")
    return [f_phi, f_coh]


# ------------------------------------------------------------------ profiles + section
def plot_material_profiles(material_nc, fault=None, *, cfg=None, sv_profile=None,
                           vs_threshold=2500.0, depths_m=DEFAULT_SLICES_M,
                           figsize=(17, 10)):
    """Slice statistics vs depth, Sv_eff(z), and the fault-following vertical section.

    The section is cut along the fault's OWN CURVED trace with the same along-strike `s`
    axis the on-fault figures use, so the two read together.  A straight chord would
    wander off a curved fault and render it as disconnected wedges -- a hole that is not
    there.
    """
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    x, y, z, vs = vs_from_material(material_nc)
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1.0, 1.15])

    a = fig.add_subplot(gs[0, 0])
    p10, p50, p90 = (np.nanpercentile(vs, q, axis=(1, 2)) for q in (10, 50, 90))
    a.plot(p50, z / 1e3, color="#2166ac", label="median")
    a.fill_betweenx(z / 1e3, p10, p90, alpha=.25, color="#2166ac", label="10-90%")
    a.axvline(vs_threshold, color="0.4", ls=":", lw=1.0, label=f"{vs_threshold:g} m/s")
    for zm in depths_m:
        a.axhline(z[_nearest_k(z, zm)] / 1e3, color="0.85", lw=.7, zorder=0)
    a.set_xlabel("$V_s$ [m/s]")
    a.set_ylabel("z [km]")
    a.grid(alpha=.3)
    a.legend(fontsize=8)
    a.set_title("(a) $V_s$ vs depth (grey = the slice depths)", fontsize=10)

    b = fig.add_subplot(gs[0, 1])
    b.plot(100.0 * np.nanmean(vs < vs_threshold, axis=(1, 2)), z / 1e3, color="#e08214")
    b.set_xlabel(f"area with $V_s$ < {vs_threshold:g} m/s  [%]")
    b.set_ylabel("z [km]")
    b.grid(alpha=.3)
    b.set_title("(b) soft-rock fraction vs depth", fontsize=10)

    c = fig.add_subplot(gs[0, 2])
    if sv_profile is not None:
        svp = np.load(sv_profile)
        c.plot(svp["sv_eff_mpa"], -svp["depth_m"] / 1e3, lw=1.8, color="#762a83")
        c.set_xlabel("$S_v^{eff}$ [MPa]")
        c.set_ylabel("z [km]")
        c.grid(alpha=.3)
        c.set_title("(c) vertical effective stress, SEA-LEVEL referenced", fontsize=10)
    else:
        c.axis("off")

    d = fig.add_subplot(gs[1, :])
    if fault is not None and cfg is not None:
        sbin, prof, ztop, zbot = _section_along_fault(
            fault, cfg, x, y, z, vs, _section_nbin(fault, cfg))
        vmin, vmax = _clim(vs)
        im = d.pcolormesh(sbin, z / 1e3, prof, cmap="viridis", vmin=vmin, vmax=vmax,
                          shading="auto", rasterized=True)
        d.plot(sbin, ztop, "-", color="k", lw=1.4, label="fault top / bottom")
        d.plot(sbin, zbot, "-", color="k", lw=1.4)
        _mark_bands(d, cfg)
        d.set_ylim(max(float(z.min()) / 1e3, 1.6 * float(np.nanmin(zbot))),
                   float(z.max()) / 1e3)
        d.set_xlabel("along-strike s (km)")
        d.set_ylabel("z (km)")
        d.legend(fontsize=8, loc="lower left")
        fig.colorbar(im, ax=d, label="$V_s$ [m/s]")
        d.set_title("(d) $V_s$ on a vertical section ALONG THE FAULT TRACE -- same s axis "
                    "as the on-fault figures, so a slow column lines up with the facets "
                    "it weakens", fontsize=10)
    else:
        d.axis("off")

    fig.suptitle("MATERIAL nc -- depth structure, and the section the fault actually sees",
                 fontsize=12, fontweight="bold")
    return fig


def _section_along_fault(fault, cfg, x, y, z, vol, nbin=240):
    """Sample `vol` down the fault's own curved trace at `nbin` evenly spaced s positions.

    Returns (s_km, (nz, nbin) values, z_top_km, z_bot_km).

    Each output column takes the fault's real (x, y) from the facets NEAREST it in s, not
    from "the facets that happen to fall inside this bin".  Binning leaves a column empty
    whenever the facet spacing is uneven -- which it always is on an unstructured mesh --
    and an empty column renders as a white stripe through the section.  Nearest-column
    assignment cannot produce one.

    A genuine gap in the fault is still reported: `gap_km` marks columns whose nearest
    facet is further away than the local column spacing, and the caller masks those.
    """
    s = fault.s_km(cfg.strike)
    s0, s1 = float(s.min()), float(s.max())
    sq = np.linspace(s0, s1, nbin)
    ds = (s1 - s0) / max(nbin - 1, 1)

    j = np.clip(np.rint((s - s0) / ds).astype(int), 0, nbin - 1)   # facet -> nearest column
    prof = np.full((len(z), nbin), np.nan)
    ztop = np.full(nbin, np.nan)
    zbot = np.full(nbin, np.nan)
    have = np.zeros(nbin, bool)

    order = np.argsort(j, kind="stable")
    js, starts = np.unique(j[order], return_index=True)
    groups = np.split(order, starts[1:])
    for jj, g in zip(js, groups):
        ix = int(np.clip(np.searchsorted(x, np.median(fault.cent[g, 0])), 0, len(x) - 1))
        iy = int(np.clip(np.searchsorted(y, np.median(fault.cent[g, 1])), 0, len(y) - 1))
        prof[:, jj] = vol[:, iy, ix]
        ztop[jj] = fault.cent[g, 2].max() / 1e3
        zbot[jj] = fault.cent[g, 2].min() / 1e3
        have[jj] = True

    if not have.all():
        # fill from the nearest populated column -- the material varies smoothly laterally,
        # so the nearest real fault column is the honest value, and it beats a white stripe.
        src = np.where(have)[0]
        if len(src):
            near = src[np.abs(np.arange(nbin)[:, None] - src[None, :]).argmin(1)]
            prof[:, ~have] = prof[:, near[~have]]
            ztop[~have] = ztop[near[~have]]
            zbot[~have] = zbot[near[~have]]
    return sq, prof, ztop, zbot


# ----------------------------------------------------------------- on-fault (s, z) maps
_WHAT = {
    "sigma_n (MPa)": "effective normal stress -- the confinement everything else scales "
                     "with",
    "tau_0 (MPa)": "initial shear traction resolved on each facet",
    "mu_app": "apparent friction $\\tau/\\sigma_n$ -- Sv-INVARIANT, so it reads the "
              "orientation and $k$ alone",
    r"$\Delta\tau_{dyn}$ (MPa)": "dynamic stress drop the run will actually pay out",
    "S (capped at 10)": "strength parameter: higher = harder to break",
    "$f_w$": "the dynamic friction floor the rs_muw LuaMap ships",
}


def plot_onfault(fault, cfg, fields: dict, *, snap=None, f0=None, fw=None,
                 fw_design=None, title=None, cmap="viridis", figsize=None,
                 map_view=True):
    """Lettered (s, z) panels on the real fault, in the legacy house style.

    `fields` is {label: per-facet array}.  Two DERIVED panels are appended when their
    inputs are supplied -- never from a guessed default, because a guessed f_w would make
    a wrong design look right:

        dtau_dyn = (mu_app - f_w) * sigma_n           what each patch can pay out
        S = (f0*sigma_n - tau)/(tau - f_w*sigma_n)    higher = harder to break (capped 10)

    Axis convention is the legacy one: **z (km), negative down**, not inverted depth.
    """
    import matplotlib.pyplot as plt

    s = fault.s_km(cfg.strike)
    zk = -fault.depth_m / 1000.0
    panels = [(k, np.asarray(v, float)) for k, v in fields.items()]

    sn = fields.get("sigma_n (MPa)")
    tau = fields.get("tau_0 (MPa)")
    mu = fields.get("mu_app")
    if sn is not None and mu is not None and fw is not None:
        panels.append((r"$\Delta\tau_{dyn}$ (MPa)",
                       (np.asarray(mu, float) - np.asarray(fw, float))
                       * np.asarray(sn, float)))
    if sn is not None and tau is not None and f0 is not None:
        fwv = np.zeros_like(np.asarray(sn, float)) if fw is None else np.asarray(fw, float)
        with np.errstate(divide="ignore", invalid="ignore"):
            S = ((f0 * np.asarray(sn, float) - np.asarray(tau, float))
                 / (np.asarray(tau, float) - fwv * np.asarray(sn, float)))
        panels.append(("S (capped at 10)", np.clip(S, 0, 10)))

    n = len(panels) + (1 if map_view else 0)
    figsize = tuple(figsize) if figsize else (15.0, 3.5 * n + 1.2)
    fig, ax = plt.subplots(n, 1, figsize=figsize, squeeze=False)
    ax = ax.ravel()
    s_lim = (float(s.min()) - 8.0, float(s.max()) + 8.0)
    letters = "abcdefghij"

    fig.canvas.draw()                    # realise the axes so marker sizing is honest
    ms_area = _marker_area(ax[0], len(fault))
    for i, (a_, (label, v)) in enumerate(zip(ax, panels)):
        vmin, vmax = _clim(v)
        sc = a_.scatter(s, zk, c=v, s=ms_area, cmap=cmap, vmin=vmin, vmax=vmax,
                        rasterized=True)
        fig.colorbar(sc, ax=a_, label=label)
        _mark_bands(a_, cfg, label_first=(i == 0))
        _mark_boundaries(a_, fw_design)
        if snap is not None:
            a_.plot([snap.s_km], [-snap.depth_m / 1000.0],
                    label="hypocenter" if i == 0 else None, **HYPO_KW)
        if _is_effectively_constant(v):
            a_.text(.01, .07, f"{label} is CONSTANT to ~1 part in $10^6$ "
                              f"({np.nanmean(v):.4g})", transform=a_.transAxes,
                    fontsize=8, color="0.25",
                    bbox=dict(fc="white", ec="0.7", alpha=.85, pad=1.5))
        a_.set_xlim(s_lim)
        a_.set_xlabel("along-strike s (km)")
        a_.set_ylabel("z (km)")
        a_.set_title(f"({letters[i]}) {_WHAT.get(label, label)}", fontsize=10)
        # legend() with nothing labelled warns and draws an empty box; only the first
        # panel carries labels, and only when there is a band or a hypocentre to name.
        if i == 0 and a_.get_legend_handles_labels()[0]:
            a_.legend(fontsize=8, loc="lower right")

    if map_view:
        a_ = ax[len(panels)]
        key, v = panels[0]
        vmin, vmax = _clim(v)
        sc = a_.scatter(fault.cent[:, 0] / 1e3, fault.cent[:, 1] / 1e3, c=v,
                        s=min(ms_area, 40.0), cmap=cmap, vmin=vmin, vmax=vmax,
                        rasterized=True)
        fig.colorbar(sc, ax=a_, label=key)
        if snap is not None:
            a_.plot([snap.xyz[0] / 1e3], [snap.xyz[1] / 1e3], label="hypocenter",
                    **HYPO_KW)
            a_.legend(fontsize=8, loc="lower left")

        a_.set_aspect("equal")
        a_.set_anchor("C")
        a_.set_xlabel("UTM easting (km)")
        a_.set_ylabel("UTM northing (km)")
        a_.grid(alpha=.3)
        a_.set_title(f"({letters[len(panels)]}) the same field on the REAL fault (map "
                     f"view) -- a kink here is a kink in the MESH", fontsize=10)

    top = _suptitle(fig, title or
                    f"ON-FAULT CHECK -- {len(fault):,} facets, sampled from the WRITTEN "
                    f"nc at the real facet centroids", figsize)
    fig.tight_layout(rect=(0, 0, 1, top))
    return fig


def onfault_report(fault, cfg, sigma_n, tau, mu_app, *, f0=None, fw=None,
                   fw_design=None, seis_band_km=None) -> str:
    """The text diagnostics that belong under `plot_onfault`.

    A figure shows a pattern; these are the numbers you quote and compare between designs.
    Mirrors the legacy per-region printout.
    """
    s = fault.s_km(cfg.strike)
    dkm = fault.depth_m / 1000.0
    sn, ta, mu = (np.asarray(a, float) for a in (sigma_n, tau, mu_app))
    band = seis_band_km or getattr(getattr(cfg, "physics", None), "seis_band_km", None)

    L = [f"on-fault summary over {len(fault):,} facets, "
         f"s {s.min():.1f} to {s.max():.1f} km, depth 0 to {dkm.max():.1f} km",
         f"  sigma_n  min {np.nanmin(sn):8.3f}  median {np.nanmedian(sn):8.3f}  "
         f"max {np.nanmax(sn):8.3f}  MPa",
         f"  tau_0    min {np.nanmin(ta):8.3f}  median {np.nanmedian(ta):8.3f}  "
         f"max {np.nanmax(ta):8.3f}  MPa",
         f"  mu_app   min {np.nanmin(mu):8.4f}  median {np.nanmedian(mu):8.4f}  "
         f"max {np.nanmax(mu):8.4f}"]

    if band is not None:
        m = (dkm >= band[0]) & (dkm <= band[1])
        if m.any():
            L.append(f"  seismogenic corridor {band[0]:g}-{band[1]:g} km: "
                     f"{int(m.sum()):,} facets, median mu_app {np.nanmedian(mu[m]):.4f}")

    fwv = None if fw is None else np.asarray(fw, float)
    if fwv is not None:
        dt = (mu - fwv) * sn
        L.append(f"  dtau_dyn median {np.nanmedian(dt):.2f} MPa "
                 f"({100.0 * np.nanmean(dt <= 0):.2f}% locked, i.e. f_w >= mu_app)")
        bounds = _design_bounds(fw_design)
        if bounds:
            edges = [-np.inf] + [w[0] for w in bounds] + [np.inf]
            for i in range(len(edges) - 1):
                m = (s >= edges[i]) & (s < edges[i + 1])
                if not m.any():
                    continue
                d_new = np.nanmedian((mu[m] - fwv[m]) * sn[m])
                d_ref = np.nanmedian(mu[m] * sn[m])
                rem = 100.0 * (1.0 - d_new / d_ref) if d_ref else 0.0
                L.append(f"    region {i + 1}: f_w={np.nanmedian(fwv[m]):<8.4g} "
                         f"dtau_dyn {d_ref:6.2f} -> {d_new:6.2f} MPa "
                         f"({rem:5.1f}% removed), "
                         f"{100.0 * np.nanmean((mu[m] - fwv[m]) * sn[m] <= 0):.2f}% locked")

    if f0 is not None:
        ok = np.nanmax(mu) < f0
        L.append(f"  no pre-slip requires max mu_app < f0: {np.nanmax(mu):.4f} < {f0:g} "
                 f" ->  {'OK' if ok else 'VIOLATED'}")
    return "\n".join(L)


# --------------------------------------------------------------- the design read-out
def _band_profile(s, v, mask, sbin_km, min_per_bin=5):
    """(centres, median, p10, p90) of `v` binned along strike over `mask`."""
    edges = np.arange(np.floor(np.nanmin(s)), np.ceil(np.nanmax(s)) + sbin_km, sbin_km)
    med = np.full(len(edges) - 1, np.nan)
    lo = np.full_like(med, np.nan)
    hi = np.full_like(med, np.nan)
    for i in range(len(edges) - 1):
        sel = mask & (s >= edges[i]) & (s < edges[i + 1])
        if int(sel.sum()) >= min_per_bin:
            med[i], lo[i], hi[i] = np.nanpercentile(v[sel], [50, 10, 90])
    return 0.5 * (edges[:-1] + edges[1:]), med, lo, hi


def plot_fw_and_kappa(fault, cfg, fw_design, *, mu_app, sigma_n_mpa, fw_on_fault,
                      tau_mpa=None, snap=None, f0=None, a_minus_b=None, figsize=(14, 13)):
    """The whole design on ONE along-strike axis, in causal order.

      (a) f_w(s)                      -- THE LEVER you set
      (b) sigma_n and tau_0           -- the stress it acts on (unchanged by f_w)
      (c) mu_app vs f_w               -- the MARGIN between them; the shaded gap IS the
                                         dynamic stress drop, in friction units
      (d) corridor kappa              -- what the lever costs

    They belong together because kappa ~ (mu_app - f_w)^2 / (f0 - f_w): a step in (a) that
    looks modest closes the gap in (c) and cuts (d) as roughly the SQUARE of the drop.
    Reading f_w alone hides that, which is why the legacy never plots one without the rest.

    Every panel is a binned median over the seismogenic band with a 10-90 percentile band,
    so a single ugly facet cannot move a curve.

    kappa_c is NECESSARY, NOT SUFFICIENT -- a design that cleared it still arrested at
    s ~ 60 km on the cluster.  Only the run decides.
    """
    import matplotlib.pyplot as plt

    from deckbuild.friction import fw_profile_1d, kappa_profile

    ph = cfg.physics
    sb = float(ph.kappa_sbin_km)
    s = fault.s_km(cfg.strike)
    d = fault.depth_m / 1000.0
    lo_km, hi_km = ph.seis_band_km
    band = (d >= lo_km) & (d <= hi_km) & np.isfinite(mu_app) & (np.asarray(sigma_n_mpa) > 0)
    if a_minus_b is not None:
        band &= np.asarray(a_minus_b) < 0.0
    s_lim = (float(s.min()) - 8.0, float(s.max()) + 8.0)

    fig, ax = plt.subplots(4, 1, figsize=figsize, sharex=True,
                           gridspec_kw=dict(height_ratios=[1.0, 1.15, 1.15, 1.3]))

    # (a) the lever
    ss = np.linspace(s_lim[0], s_lim[1], 1600)
    ax[0].plot(ss, fw_profile_1d(ss, fw_design), lw=2, color="tab:blue", label="$f_w(s)$")
    ax[0].set_ylabel("$f_w$")
    ax[0].set_title("(a) the graded strong-rate-weakening floor $f_w(s)$ -- THE LEVER",
                    fontsize=10)

    # (b) the stress it acts on
    c, m_sn, l_sn, h_sn = _band_profile(s, np.asarray(sigma_n_mpa, float), band, sb)
    ax[1].plot(c, m_sn, lw=1.8, color="tab:purple", label=r"$\sigma_n$ median")
    ax[1].fill_between(c, l_sn, h_sn, alpha=.22, color="tab:purple", label="10-90%")
    if tau_mpa is not None:
        _, m_t, l_t, h_t = _band_profile(s, np.asarray(tau_mpa, float), band, sb)
        ax[1].plot(c, m_t, lw=1.8, color="tab:green", label=r"$\tau_0$ median")
        ax[1].fill_between(c, l_t, h_t, alpha=.22, color="tab:green")
    ax[1].set_ylabel("stress (MPa)")
    ax[1].set_title(rf"(b) the on-fault stress over {lo_km:g}-{hi_km:g} km -- what the "
                    rf"design acts on.  $f_w$ does NOT change this", fontsize=10)

    # (c) the margin
    _, m_mu, l_mu, h_mu = _band_profile(s, np.asarray(mu_app, float), band, sb)
    _, m_fw, _, _ = _band_profile(s, np.asarray(fw_on_fault, float), band, sb)
    ax[2].plot(c, m_mu, lw=1.8, color="tab:orange", label=r"$\mu_{app}$ median")
    ax[2].fill_between(c, l_mu, h_mu, alpha=.22, color="tab:orange", label="10-90%")
    ax[2].plot(c, m_fw, lw=1.8, color="tab:blue", label="$f_w$")
    ax[2].fill_between(c, m_fw, m_mu, where=np.isfinite(m_mu) & np.isfinite(m_fw),
                       alpha=.28, color="tab:red",
                       label=r"the gap = $\Delta\tau_{dyn}/\sigma_n$")
    ax[2].set_ylabel("friction")
    ax[2].set_title(r"(c) apparent friction against the floor -- the SHADED GAP is what "
                    r"the rupture pays out", fontsize=10)

    # (d) the cost
    sc, med, coast = kappa_profile(fault, cfg, mu_app, sigma_n_mpa, fw_on_fault,
                                   f0=f0, a_minus_b=a_minus_b)
    ax[3].plot(sc, med, lw=1.0, color="0.6", label=rf"$\kappa$, {sb:g} km bins")
    ax[3].plot(sc, coast, lw=2.2, color="tab:red",
               label=rf"coasting $\kappa$ ({ph.l_coast_km:g} km mean) -- what the gate reads")
    kc = float(ph.kappa_c)
    ax[3].axhline(kc, color="k", ls="--", lw=1.2, label=rf"$\kappa_c$ = {kc:g}")
    ax[3].axhline(1.0, color="0.3", ls=":", lw=1.0, label=r"$\kappa$ = 1 break-even")
    if np.isfinite(coast).any():
        i = int(np.nanargmin(coast))
        ax[3].plot([sc[i]], [coast[i]], "v", ms=9, mfc="crimson", mec="k", zorder=6)
        ax[3].annotate(rf"min coasting $\kappa$ = {coast[i]:.3f} at s = {sc[i]:.0f} km",
                       xy=(sc[i], coast[i]), xytext=(6, 14), textcoords="offset points",
                       fontsize=8, color="crimson")
    ax[3].set_yscale("log")
    ax[3].set_ylabel(r"corridor $\kappa$")
    ax[3].set_xlabel("along-strike s (km)")
    ax[3].set_title(rf"(d) corridor $\kappa$ -- WHAT THE LEVER COSTS.  "
                    rf"$\kappa \sim (\mu_{{app}} - f_w)^2/(f_0 - f_w)$, so it falls as "
                    rf"the SQUARE of the drop", fontsize=10)

    # Per-panel legend placement: a legend that covers the curve it describes is worse
    # than no legend.  (d) rises to the top-left, so its box goes bottom-right.
    locs = ["upper left", "upper left", "upper right", "lower right"]
    for i, a_ in enumerate(ax):
        _mark_bands(a_, cfg, label_first=(i == 0))
        _mark_boundaries(a_, fw_design)
        if snap is not None:
            a_.axvline(snap.s_km, color="crimson", ls="--", lw=1.2,
                       label="hypocentre" if i == 0 else None)
        a_.set_xlim(s_lim)
        a_.grid(alpha=.3)
        a_.legend(fontsize=8, loc=locs[i], ncol=2, framealpha=.9)
    # headroom so the kappa curve never runs under its own title or legend
    fin = coast[np.isfinite(coast)]
    if fin.size:
        ax[3].set_ylim(min(0.5 * float(np.nanmin(med[np.isfinite(med)])), 0.6 * kc),
                       4.0 * float(fin.max()))

    top = _suptitle(fig, "GRADED $f_w$ DESIGN READ-OUT -- lever, stress, margin, cost.  "
                         "kappa is NECESSARY, NOT SUFFICIENT: a design that cleared this "
                         "screen still arrested on the cluster", figsize)
    fig.tight_layout(rect=(0, 0, 1, top))
    return fig


