"""orientation.py -- where the horizontal stress orientation comes from.

The stress closure needs two per-column scalars:

    az_deg   the SHmax azimuth, degrees east of north
    R        the shape ratio (w1 - w0) / (w2 - w0), in [0, 1]

Three sources, dispatched on `SourceSpec.kind`.  This is the extension point for a fault
system with a different stress product -- or none at all.

  "csm_csv"   a community stress model csv (the SAFS path).  Uses the AUTHOR's SHmax
              azimuth column, not an eigendecomposition: the legacy recipe is explicit
              that the author azimuth is the intended input.
  "constant"  one azimuth and one R everywhere.  What a fault system with no community
              stress model uses, and what the demo does.
  "callable"  a dotted path to f(x, y) -> (az_deg, R), for anything else.
"""
from __future__ import annotations

import csv as csv_mod
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from deckbuild.config import CRS, ConfigError, SourceSpec

__all__ = ["OrientationField", "read_orientation", "OrientationError",
           "read_csm_csv", "csm_tensors_tension", "csm_axes_and_shape",
           "interpolate_csm_field"]

EPS = 1.0e-9


class OrientationError(ValueError):
    """Raised when an orientation source is unusable."""


@dataclass(frozen=True, eq=False)
class OrientationField:
    """Scattered orientation samples in PROJECTED coordinates.

    `cx`, `cy` are empty for a uniform source; `az_deg`/`R` are then scalars broadcast at
    interpolation time.  `is_uniform` tells the two apart.
    """

    cx: np.ndarray
    cy: np.ndarray
    az_deg: np.ndarray
    R: np.ndarray
    kind: str = ""
    n_dropped: int = 0

    @property
    def is_uniform(self) -> bool:
        return self.cx.size == 0

    def at(self, qx, qy) -> tuple[np.ndarray, np.ndarray, int]:
        """(az_deg, R, n_nearest_fallback) at query points."""
        qx = np.asarray(qx, float).ravel()
        qy = np.asarray(qy, float).ravel()
        if self.is_uniform:
            return (np.full(qx.shape, float(self.az_deg)),
                    np.full(qx.shape, float(self.R)), 0)
        # Interpolate the azimuth as a UNIT VECTOR doubled in angle, never as a raw
        # number: SHmax is an axis with period 180 deg, so averaging 179 and 1 must give
        # 0, not 90.  The legacy pipeline interpolated the tensor and re-derived the
        # azimuth, which has the same effect.
        a2 = np.radians(2.0 * self.az_deg)
        stack = np.stack([np.cos(a2), np.sin(a2), self.R], axis=1)
        vals, n_fb = interpolate_csm_field(self.cx, self.cy, stack, qx, qy)
        az = 0.5 * np.degrees(np.arctan2(vals[:, 1], vals[:, 0])) % 180.0
        return az, np.clip(vals[:, 2], 0.0, 1.0), n_fb


# --------------------------------------------------------------------------- csm csv
def read_csm_csv(path, col_lon: int, col_lat: int, col_s: slice | list,
                 col_shmax: int | None = None, col_R: int | None = None):
    """Parse a community-stress-model csv.

    VERBATIM parsing behaviour from the legacy `read_csm_csv`: '#' lines skipped, rows
    with a NaN tensor component dropped (orientation is undefined there).  The column
    indices are parameters here rather than module constants.
    """
    path = Path(path)
    if not path.is_file():
        raise OrientationError(f"orientation csv not found: {path.resolve()}")
    lon, lat, S, shmax, Rc = [], [], [], [], []
    idx_s = list(range(col_s.start, col_s.stop)) if isinstance(col_s, slice) else list(col_s)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            row = next(csv_mod.reader([line]))
            vals = [float(v) if v.strip() not in ("", "NaN", "nan") else np.nan
                    for v in row]
            lon.append(vals[col_lon])
            lat.append(vals[col_lat])
            S.append([vals[i] for i in idx_s])
            if col_shmax is not None:
                shmax.append(vals[col_shmax])
            if col_R is not None:
                Rc.append(vals[col_R])
    if not S:
        raise OrientationError(f"no data rows in {path}")
    lon, lat, S = np.asarray(lon), np.asarray(lat), np.asarray(S)
    bad = np.isnan(S).any(axis=1)
    keep = ~bad
    out = {"lon": lon[keep], "lat": lat[keep], "S": S[keep],
           "n_dropped": int(bad.sum())}
    if col_shmax is not None:
        out["shmax"] = np.asarray(shmax)[keep]
    if col_R is not None:
        out["R"] = np.asarray(Rc)[keep]
    return out


def csm_tensors_tension(S):
    """(N,6) See,Sen,Seu,Snn,Snu,Suu -> (N,3,3) TENSION-POSITIVE (e,n,u) tensors."""
    T = np.empty((len(S), 3, 3))
    T[:, 0, 0] = S[:, 0]
    T[:, 0, 1] = T[:, 1, 0] = S[:, 1]
    T[:, 0, 2] = T[:, 2, 0] = S[:, 2]
    T[:, 1, 1] = S[:, 3]
    T[:, 1, 2] = T[:, 2, 1] = S[:, 4]
    T[:, 2, 2] = S[:, 5]
    return T


def csm_axes_and_shape(T_tension):
    """Eigendecomposition -> (U, R).  legacy: csm_axes_and_shape."""
    w, U = np.linalg.eigh(T_tension)
    denom = w[:, 2] - w[:, 0]
    R_eig = np.where(denom > EPS, (w[:, 1] - w[:, 0]) / np.where(denom > EPS, denom, 1.0),
                     np.nan)
    return U, R_eig


def interpolate_csm_field(cx, cy, values, qx, qy):
    """Delaunay-linear interpolation with nearest-neighbour fallback outside the hull.

    VERBATIM from the legacy `interpolate_csm_field`.  Returns (interp, n_fallback); the
    fallback count matters, because a large one means the grid box reaches well outside
    the stress model's coverage.
    """
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
    pts = np.stack([cx, cy], axis=1)
    q = np.stack([qx, qy], axis=1)
    vals = np.asarray(values, dtype=float)
    flat = vals.reshape(len(vals), -1)
    out = LinearNDInterpolator(pts, flat)(q)
    bad = ~np.isfinite(out).all(axis=1)
    n_fallback = int(bad.sum())
    if n_fallback:
        out[bad] = NearestNDInterpolator(pts, flat)(q[bad])
    return out.reshape((len(q),) + vals.shape[1:]), n_fallback


# --------------------------------------------------------------------------- dispatch
def _read_csm(spec: SourceSpec, crs: CRS, data_dir) -> OrientationField:
    p = spec.params
    path = Path(spec.path)
    if not path.is_absolute() and data_dir is not None:
        path = Path(data_dir) / path
    col_shmax = int(p.get("shmax_col", 10))
    col_R = int(p.get("shape_ratio_col", 13))
    raw = read_csm_csv(path, int(p.get("lon_col", 0)), int(p.get("lat_col", 1)),
                       slice(int(p.get("tensor_col0", 4)), int(p.get("tensor_col0", 4)) + 6),
                       col_shmax=col_shmax, col_R=col_R)
    from pyproj import Transformer
    tf = Transformer.from_crs(crs.geographic_epsg, crs.epsg, always_xy=True)
    cx, cy = tf.transform(raw["lon"], raw["lat"])

    az = np.asarray(raw["shmax"], float)
    R = np.asarray(raw["R"], float)
    finite = np.isfinite(az) & np.isfinite(R)
    if not finite.any():
        raise OrientationError(
            f"{path}: no row has a finite SHmax (col {col_shmax}) and R (col {col_R})")
    n_dropped = raw["n_dropped"] + int((~finite).sum())
    return OrientationField(cx=np.asarray(cx)[finite], cy=np.asarray(cy)[finite],
                            az_deg=az[finite] % 180.0, R=np.clip(R[finite], 0.0, 1.0),
                            kind="csm_csv", n_dropped=n_dropped)


def _read_constant(spec: SourceSpec, crs: CRS, data_dir) -> OrientationField:
    p = spec.params
    if "azimuth_deg" not in p or "shape_ratio" not in p:
        raise OrientationError(
            "orientation kind 'constant' needs params.azimuth_deg and params.shape_ratio")
    R = float(p["shape_ratio"])
    if not 0.0 <= R <= 1.0:
        raise OrientationError(f"shape_ratio must be in [0, 1], got {R}")
    return OrientationField(cx=np.empty(0), cy=np.empty(0),
                            az_deg=np.float64(float(p["azimuth_deg"]) % 180.0),
                            R=np.float64(R), kind="constant")


def _read_callable(spec: SourceSpec, crs: CRS, data_dir) -> OrientationField:
    target = spec.params.get("target") or spec.path
    if not target or ":" not in str(target):
        raise OrientationError(
            "orientation kind 'callable' needs params.target as 'module.path:function'")
    mod_name, fn_name = str(target).rsplit(":", 1)
    import importlib
    fn = getattr(importlib.import_module(mod_name), fn_name)
    p = spec.params
    if "sample_xy" not in p:
        raise OrientationError(
            "orientation kind 'callable' needs params.sample_xy = [[x, y], ...] so the "
            "callable can be evaluated on a known point set")
    pts = np.asarray(p["sample_xy"], float)
    az, R = fn(pts[:, 0], pts[:, 1])
    return OrientationField(cx=pts[:, 0], cy=pts[:, 1],
                            az_deg=np.asarray(az, float) % 180.0,
                            R=np.clip(np.asarray(R, float), 0.0, 1.0), kind="callable")


READERS = {"csm_csv": _read_csm, "constant": _read_constant, "callable": _read_callable}


def read_orientation(spec: SourceSpec | None, crs: CRS,
                     data_dir=None) -> OrientationField:
    if spec is None:
        raise ConfigError("this project has no raw.orientation source")
    if spec.kind not in READERS:
        raise OrientationError(
            f"unknown orientation kind {spec.kind!r}; have {sorted(READERS)}")
    return READERS[spec.kind](spec, crs, data_dir)
