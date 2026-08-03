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
           "read_csm_csv", "csm_tensors_tension", "csm_axes_and_shape", "shmax_from_sigma",
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
    S: np.ndarray | None = None      # (N, 6) tension-positive tensor samples, csm_csv only

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
        if self.S is not None:
            # INTERPOLATE THE TENSOR, then derive az and R at the query points.  This is
            # what the legacy does (step2_grid.py:42-45) and the order is NOT
            # interchangeable with deriving first and interpolating the result: on the
            # shipped ALT grid the two differ by up to 1.0e8 Pa in s_xx/s_yy/s_xy.
            # Interpolating a doubled-angle unit vector is the right way to average
            # AZIMUTHS, but it is not what produced the shipped file.
            vals, n_fb = interpolate_csm_field(self.cx, self.cy, self.S, qx, qy)
            T = csm_tensors_tension(vals)
            return (shmax_from_sigma(-T), np.clip(csm_axes_and_shape(T)[1], 0.0, 1.0),
                    n_fb)
        # No tensor available (constant / callable readers): average the azimuth as a
        # UNIT VECTOR DOUBLED IN ANGLE, never as a raw number -- SHmax has period 180 deg,
        # so averaging 179 and 1 must give 0, not 90.
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


def shmax_from_sigma(sigma):
    """Azimuth (deg E of N, in [0,180)) of the most-compressive eigenvector of the
    horizontal 2x2 block of the COMPRESSION-POSITIVE tensor.

    VERBATIM from the legacy `project_csm_stress_to_vtu.shmax_from_sigma`.  Callers pass
    -T, where T is the tension-positive CSM tensor.
    """
    h = sigma[:, :2, :2]
    w, U = np.linalg.eigh(h)                  # ascending
    v = U[:, :, 1]                            # most compressive (e, n)
    az = np.degrees(np.arctan2(v[:, 0], v[:, 1]))   # E of N
    return np.mod(az, 180.0)


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
    # COLUMN INDICES ARE 0-BASED.  The CSM csv's own header numbers its columns 1-BASED
    # ("10) SHmax angle...", "13) R stress ratio...", "4) See..."), and using those
    # numbers as indices lands one column late on every one:
    #     10 -> SHmax_unc (an UNCERTAINTY in degrees, 0-21.75)
    #     13 -> Aphi      (Anderson's shape parameter, 0.055-2.732, clipped into [0,1])
    #      4 -> Sen       (the tensor read shifted by one component)
    # The defaults below are the 0-based indices of SHmax, R and See.  Gate O1 in
    # read_orientation() re-checks them against the data so this cannot recur silently.
    # DERIVE az and R FROM THE TENSOR -- exactly what the legacy does
    # (lib/step1_orientation.py:43-44, all three copies):
    #     az_pts = shmax_from_sigma(-T)
    #     R_pts  = clip(csm_axes_and_shape(T)[1], 0, 1)
    # and the legacy reads ONLY lon, lat and the six tensor components
    # (generate_stress_nc_from_raw.py:111-112, COL_S = slice(3, 9)).  It never reads an
    # author SHmax or R column at all.
    #
    # The author columns are still read, but ONLY to cross-check the derivation (gate O1).
    # They are not what the stress is built from.  Their indices are 0-BASED while the csv
    # header numbers its columns 1-BASED, which is how an earlier version of this reader
    # ended up on SHmax_unc and Aphi.
    col_shmax = int(p.get("shmax_col", 9))
    col_R = int(p.get("shape_ratio_col", 12))
    tcol0 = int(p.get("tensor_col0", 3))
    raw = read_csm_csv(path, int(p.get("lon_col", 0)), int(p.get("lat_col", 1)),
                       slice(tcol0, tcol0 + 6), col_shmax=col_shmax, col_R=col_R)
    from pyproj import Transformer
    tf = Transformer.from_crs(crs.geographic_epsg, crs.epsg, always_xy=True)
    cx, cy = tf.transform(raw["lon"], raw["lat"])

    S = np.asarray(raw["S"], float)
    T = csm_tensors_tension(S)                       # tension-positive, as stored
    az = shmax_from_sigma(-T)                        # -T -> compression-positive
    R = csm_axes_and_shape(T)[1]
    finite = np.isfinite(az) & np.isfinite(R)
    if not finite.any():
        raise OrientationError(
            f"{path}: no row yields a finite SHmax and R from the tensor columns "
            f"{tcol0}-{tcol0 + 5}")
    n_dropped = raw["n_dropped"] + int((~finite).sum())

    # GATE O1 -- are shmax_col and shape_ratio_col pointing at the right columns?
    # They sit next to look-alikes (SHmax_unc, phi, Aphi), and the csv header numbers its
    # columns 1-BASED, so an off-by-one is the expected mistake rather than an exotic one.
    # The tensor components are RIGHT THERE in the same rows, so the columns can be
    # checked against the physics instead of trusted:
    #   * R must be a ratio in [0, 1] before any clipping -- Aphi runs to 2.7 and
    #     SHmax_unc is in degrees, so either one fails this outright;
    #   * R must agree with the eigendecomposition of the tensor, which is what makes
    #     this a check on the COLUMN CHOICE and not merely on the range.
    _audit_csm_columns(raw, az, R, col_shmax, col_R, path)
    return OrientationField(cx=np.asarray(cx)[finite], cy=np.asarray(cy)[finite],
                            az_deg=az[finite] % 180.0, R=np.clip(R[finite], 0.0, 1.0),
                            kind="csm_csv", n_dropped=n_dropped, S=S[finite])


def _audit_csm_columns(raw, az_derived, R_derived, col_shmax, col_R, path,
                       az_tol_deg=2.0, r_tol=0.02):
    """GATE O1 -- cross-check the tensor derivation against the author columns.

    az and R are DERIVED from the tensor (the legacy path), so this is not what the
    stress is built from.  It exists because the author columns sit next to look-alikes
    -- SHmax_unc, phi, Aphi -- and the csv header numbers its columns 1-BASED while these
    indices are 0-BASED.  An earlier version of this reader took the header's numbers
    literally and read SHmax_unc as the azimuth and Aphi as R, which moved on-fault
    mu_app by 21% and passed every stress gate: V2 checks the eigenvalue RATIO k, which
    is blind to both the azimuth and R.

    Disagreement means the COLUMNS are mis-identified (or the tensor slice is), not that
    the physics is wrong -- so the message says which index to look at.
    """
    az_c = np.asarray(raw.get("shmax"), float)
    R_c = np.asarray(raw.get("R"), float)

    good = np.isfinite(R_c)
    if good.any():
        lo, hi = float(np.min(R_c[good])), float(np.max(R_c[good]))
        if lo < -r_tol or hi > 1.0 + r_tol:
            raise OrientationError(
                f"{path}: shape_ratio_col={col_R} has range [{lo:.3f}, {hi:.3f}], which "
                f"is not a ratio in [0, 1].  The csv header numbers its columns 1-BASED; "
                f"if you took the index from it, subtract 1 (R is 0-based 12; 13 is Aphi).")
        m = good & np.isfinite(R_derived)
        if m.sum() > 100:
            d = float(np.median(np.abs(R_derived[m] - R_c[m])))
            if d > r_tol:
                raise OrientationError(
                    f"{path}: R derived from the tensor columns disagrees with "
                    f"shape_ratio_col={col_R} by a median of {d:.4f}.  Either that index "
                    f"or tensor_col0 is wrong -- the header's numbers are 1-BASED, these "
                    f"are 0-BASED.")

    m = np.isfinite(az_c) & np.isfinite(az_derived)
    if m.any():
        span = float(np.max(az_c[m]) - np.min(az_c[m]))
        if span <= 25.0:
            raise OrientationError(
                f"{path}: shmax_col={col_shmax} spans only {span:.2f} deg, which looks "
                f"like an UNCERTAINTY column rather than an azimuth.  SHmax is 0-based "
                f"column 9; 10 is SHmax_unc.")
        if m.sum() > 100:
            dd = np.abs(np.mod(az_derived[m] - (az_c[m] % 180.0), 180.0))
            dd = np.minimum(dd, 180.0 - dd)
            med = float(np.median(dd))
            if med > az_tol_deg:
                raise OrientationError(
                    f"{path}: the azimuth derived from the tensor disagrees with "
                    f"shmax_col={col_shmax} by a median of {med:.2f} deg.  Either that "
                    f"index or tensor_col0 is wrong (0-BASED: SHmax 9, tensor 3-8).")


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
