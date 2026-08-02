"""rawslices.py -- horizontal-slice ASCII products -> a uniform ASAGI grid.

Ported from the legacy `generate_velocity_nc_from_raw.py` and
`generate_thermal_nc_from_raw.py`, which are the same two-stage recipe applied to
different scalars.  The ORDER of the stages is load-bearing:

  Stage 1  read each slice (one per depth; lon/lat grid) -> reproject to the project CRS
           -> build the INSCRIBED axis-aligned grid, rounded INWARD so every node is
           inside the source hull -> LinearNDInterpolator per slice per field; a NaN
           cell fails loud -> flip depth to elevation and clone the surface slice
           `extend_z_top` metres up (to cover topography).

  Stage 2  form the moduli AT THE SOURCE NODES (mu = rho Vs^2,
           lambda = rho(Vp^2 - 2Vs^2)) and only THEN resample onto the uniform z axis.
           Converting after resampling would store something that is not the
           piecewise-linear interpolant of the node moduli.

A missing depth level is a WARNING, not a failure: the output node is then the exact
linear interpolant of its neighbours (the CTM legitimately lacks one level).
"""
from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np

__all__ = ["read_slice", "slice_depth_m", "build_uniform_grid", "interp_slices",
           "SliceError"]


class SliceError(ValueError):
    """Raised when a slice stack is unusable."""


_DEPTH_RE = re.compile(r"Depth\(m\)\s*:\s*([-\d.eE+]+)", re.I)
_TITLE_RE = re.compile(r"Slice at\s+([-\d.eE+]+)\s*(m|km)", re.I)


def slice_depth_m(path) -> float:
    """The slice's depth, from its header.  Raises if it cannot be found."""
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            m = _DEPTH_RE.search(line)
            if m:
                return float(m.group(1))
            m = _TITLE_RE.search(line)
            if m:
                v = float(m.group(1))
                return v * (1000.0 if m.group(2).lower() == "km" else 1.0)
    raise SliceError(f"{path}: no depth in the header (want 'Depth(m):' or 'Slice at N m')")


def read_slice(path, n_fields: int):
    """(lon, lat, values[n_pts, n_fields]) from one slice.  '#' lines are the header."""
    lon, lat, vals = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split(",")
            if len(parts) < 2 + n_fields:
                raise SliceError(
                    f"{path}: expected >= {2 + n_fields} columns, got {len(parts)}")
            lon.append(float(parts[0]))
            lat.append(float(parts[1]))
            vals.append([float(v) for v in parts[2:2 + n_fields]])
    if not lon:
        raise SliceError(f"{path}: no data rows")
    return np.asarray(lon), np.asarray(lat), np.asarray(vals)


def build_uniform_grid(cx, cy, dx: float):
    """The INSCRIBED axis-aligned grid: rounded inward so every node is inside the hull.

    Rounding outward would put nodes beyond the source data, where the interpolator
    returns NaN -- which then fails loud, but only after the expensive pass.
    """
    x0 = np.ceil(cx.min() / dx) * dx
    x1 = np.floor(cx.max() / dx) * dx
    y0 = np.ceil(cy.min() / dx) * dx
    y1 = np.floor(cy.max() / dx) * dx
    if x1 <= x0 or y1 <= y0:
        raise SliceError(
            f"the inscribed grid at dx={dx} is empty; the source spans "
            f"x[{cx.min():.0f},{cx.max():.0f}] y[{cy.min():.0f},{cy.max():.0f}]")
    return np.arange(x0, x1 + 0.5 * dx, dx), np.arange(y0, y1 + 0.5 * dx, dx)


def interp_slices(paths, crs, field_names, grid_dx, z_min, z_max, dz,
                  extend_z_top=0.0, expect_spacing_m=None):
    """The full two-stage pipeline.  Returns (gx, gy, gz, {field: (nz,ny,nx)}, info)."""
    from pyproj import Transformer
    from scipy.interpolate import LinearNDInterpolator

    paths = sorted(paths, key=slice_depth_m)
    if len(paths) < 2:
        raise SliceError(f"need at least 2 slices, got {len(paths)}")
    depths = np.array([slice_depth_m(p) for p in paths])

    # A missing level is a WARNING: the output node becomes the exact linear interpolant
    # of its neighbours.  The CTM legitimately lacks one.
    gaps = []
    if expect_spacing_m:
        d = np.diff(depths)
        gaps = [(float(depths[i]), float(depths[i + 1]))
                for i in np.flatnonzero(d > 1.5 * expect_spacing_m)]
        if gaps:
            warnings.warn(f"{len(gaps)} missing depth level(s): {gaps}; those output "
                          f"nodes are the linear interpolant of their neighbours",
                          stacklevel=2)

    tf = Transformer.from_crs(crs.geographic_epsg, crs.epsg, always_xy=True)
    # A lon/lat RECTANGLE is a CURVED quadrilateral once projected, so a bbox-inscribed
    # rectangle still pokes out past the curved edges -- on the real CVM stack that put
    # 7,624 of 181,541 nodes outside the hull.  Inscribe against the EDGES instead:
    # take the largest western x, the smallest eastern x, and likewise in y.  Do it over
    # every slice, since they need not share a footprint.
    lo = np.array([-np.inf, -np.inf])
    hi = np.array([np.inf, np.inf])
    for pth in paths:
        lon, lat, _ = read_slice(pth, len(field_names))
        lon, lat = np.asarray(lon), np.asarray(lat)
        cx, cy = (np.asarray(v) for v in tf.transform(lon, lat))
        west, east = lon <= lon.min() + 1e-9, lon >= lon.max() - 1e-9
        south, north = lat <= lat.min() + 1e-9, lat >= lat.max() - 1e-9
        lo = np.maximum(lo, [cx[west].max(), cy[south].max()])
        hi = np.minimum(hi, [cx[east].min(), cy[north].min()])
    if not np.all(hi > lo):
        raise SliceError(
            f"the slices have no common inscribed footprint: x[{lo[0]:.0f},{hi[0]:.0f}] "
            f"y[{lo[1]:.0f},{hi[1]:.0f}]")
    gx, gy = build_uniform_grid(np.array([lo[0], hi[0]]), np.array([lo[1], hi[1]]),
                                grid_dx)
    XX, YY = np.meshgrid(gx, gy, indexing="ij")
    q = np.column_stack([XX.ravel(), YY.ravel()])

    # Stage 1: per slice, per field, onto the horizontal grid.
    stack = np.empty((len(paths), len(field_names), len(gx), len(gy)))
    for i, p in enumerate(paths):
        lon, lat, v = read_slice(p, len(field_names))
        cx, cy = tf.transform(lon, lat)
        interp = LinearNDInterpolator(np.column_stack([cx, cy]), v)
        got = interp(q)
        if not np.isfinite(got).all():
            n = int((~np.isfinite(got)).any(axis=1).sum())
            raise SliceError(
                f"{Path(p).name}: {n} of {len(q)} grid nodes fell outside the source "
                f"hull and interpolated to NaN.  The inscribed grid should prevent this; "
                f"a ragged slice is the usual cause.")
        for k in range(len(field_names)):
            stack[i, k] = got[:, k].reshape(len(gx), len(gy))

    # depth -> elevation, ascending, plus the cloned cap above the surface.
    z_src = -depths                                   # descending
    order = np.argsort(z_src)
    z_src, stack = z_src[order], stack[order]
    if extend_z_top > 0:
        z_src = np.concatenate([z_src, [z_src[-1] + extend_z_top]])
        stack = np.concatenate([stack, stack[-1:]], axis=0)

    # Stage 2: resample onto the uniform axis.  The CALLER converts to moduli BEFORE
    # calling this when that matters (see material._cvm_slices).
    gz = np.arange(z_min, z_max + 0.5 * dz, dz)
    out = {}
    for k, name in enumerate(field_names):
        src = stack[:, k]                              # (nz_src, nx, ny)
        flat = src.reshape(len(z_src), -1)
        res = np.empty((len(gz), flat.shape[1]))
        for j in range(flat.shape[1]):
            res[:, j] = np.interp(gz, z_src, flat[:, j])
        out[name] = np.transpose(res.reshape(len(gz), len(gx), len(gy)), (0, 2, 1))
    info = {"n_slices": len(paths), "depth_min": float(depths.min()),
            "depth_max": float(depths.max()), "gaps": gaps,
            "nx": len(gx), "ny": len(gy), "nz": len(gz)}
    return gx, gy, gz, out, info
