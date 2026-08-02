"""asagi.py -- the ASAGI NetCDF format layer.

SeisSol's easi `!ASAGI` reader expects one exact layout, and it does not check most of it:

    COARDS NETCDF4
    x / y / z   float64, strictly increasing, EQUIDISTANT
    data        ONE compound variable, dims (z, y, x), float32 members

ASAGI silently produces wrong interpolation on a non-equidistant axis, so `write_asagi`
checks that itself.  `z` is elevation, positive up.

Ported from the legacy `generate_stress_nc_from_raw.py` and `build_friction_nc_thermal.py`,
which each owned a copy of this.  The expressions are unchanged; only the field set is
generalised (the legacy writer hardcoded the six stress components).
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Mapping, Sequence

import numpy as np

from deckbuild.contract import GateReport, HARD, WARN

__all__ = [
    "write_asagi", "read_asagi", "asagi_axes", "asagi_fields",
    "trilinear_sample", "hull_containment", "roundtrip_selfcheck",
    "AsagiError", "axis_spacing",
]

# Relative tolerance for "equidistant".  float64 axes built by np.arange drift by ~1e-16
# relative; 1e-9 is loose enough for that and tight enough to catch a real non-uniform axis.
EQUIDISTANT_RTOL = 1.0e-9


class AsagiError(ValueError):
    """Raised when a grid or file violates the ASAGI layout contract."""


# --------------------------------------------------------------------------- writing
def axis_spacing(a: np.ndarray, name: str) -> float:
    """Validate one axis and return its spacing.

    Raises AsagiError naming the axis and the worst deviation.  ASAGI does not check any
    of this: a non-equidistant axis is read as if it were uniform, which silently shifts
    every sample.
    """
    a = np.asarray(a)
    if a.dtype != np.float64:
        raise AsagiError(f"axis {name!r} must be float64, got {a.dtype}")
    if a.ndim != 1:
        raise AsagiError(f"axis {name!r} must be 1-D, got shape {a.shape}")
    if a.size < 2:
        raise AsagiError(f"axis {name!r} must have >= 2 nodes, got {a.size}")
    d = np.diff(a)
    if not np.all(d > 0):
        bad = int(np.argmin(d))
        raise AsagiError(
            f"axis {name!r} must be strictly increasing; "
            f"a[{bad}]={a[bad]!r} >= a[{bad + 1}]={a[bad + 1]!r}")
    dmin, dmax, dmean = float(d.min()), float(d.max()), float(d.mean())
    dev = (dmax - dmin) / abs(dmean)
    if dev > EQUIDISTANT_RTOL:
        worst = int(np.argmax(np.abs(d - dmean)))
        raise AsagiError(
            f"axis {name!r} is not equidistant (relative spread {dev:.3e} > "
            f"{EQUIDISTANT_RTOL:.0e}); spacing ranges [{dmin!r}, {dmax!r}], worst at "
            f"index {worst}.  ASAGI assumes a uniform axis and would mis-locate every "
            f"sample.")
    return dmean


def write_asagi(path: str | Path, x, y, z, fields: Mapping[str, np.ndarray],
                attrs: Mapping[str, object] | None = None,
                dtype=np.float32, allow_nonfinite: bool = False) -> Path:
    """Write a compound-variable ASAGI NetCDF.

    fields -- {name: array of shape (nz, ny, nx)}.  NOTE the ordering: the legacy
              `write_stress_asagi` took (nx, ny, nz) and transposed internally; this
              function takes the STORED order directly, so a caller porting legacy code
              must transpose (2, 1, 0) itself.  A wrong order is caught below, not
              silently written.
    allow_nonfinite -- default False: ASAGI propagates a NaN into every element that
              samples the grid, so a non-finite value is rejected at write time rather
              than discovered in a diverging run.  Set True only for a field where a
              sentinel is deliberate.
    """
    path = Path(path)
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    for a, n in ((x, "x"), (y, "y"), (z, "z")):
        axis_spacing(a, n)

    if not fields:
        raise AsagiError("write_asagi needs at least one field")
    nx, ny, nz = len(x), len(y), len(z)
    want = (nz, ny, nx)
    for name, arr in fields.items():
        arr = np.asarray(arr)
        if arr.shape != want:
            hint = ""
            if arr.shape == (nx, ny, nz):
                hint = ("  -- this is (nx, ny, nz); ASAGI stores (nz, ny, nx). "
                        "Transpose with np.transpose(arr, (2, 1, 0)).")
            raise AsagiError(
                f"field {name!r} has shape {arr.shape}, expected {want}{hint}")
        if not allow_nonfinite and not np.all(np.isfinite(arr)):
            n_bad = int((~np.isfinite(arr)).sum())
            raise AsagiError(
                f"field {name!r} has {n_bad} non-finite value(s); ASAGI would propagate "
                f"them into every element that samples this grid.  Pass "
                f"allow_nonfinite=True if the sentinel is deliberate.")

    from netCDF4 import Dataset

    names = list(fields)
    st = np.dtype([(f, dtype) for f in names])
    path.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(str(path), "w", format="NETCDF4") as ds:
        ds.createDimension("x", nx)
        ds.createDimension("y", ny)
        ds.createDimension("z", nz)
        ds.createVariable("x", "f8", ("x",))[:] = x
        ds.createVariable("y", "f8", ("y",))[:] = y
        ds.createVariable("z", "f8", ("z",))[:] = z
        mtype = ds.createCompoundType(st, "data_t")
        var = ds.createVariable("data", mtype, ("z", "y", "x"))
        buf = np.empty((nz, ny, nx), dtype=st)
        for f in names:
            buf[f] = np.asarray(fields[f]).astype(dtype)
        var[:] = buf
        for k, v in (attrs or {}).items():
            ds.setncattr(k, v)
    return path


# --------------------------------------------------------------------------- reading
def asagi_axes(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    from netCDF4 import Dataset
    with Dataset(str(path), "r") as ds:
        return (np.asarray(ds.variables["x"][:], float),
                np.asarray(ds.variables["y"][:], float),
                np.asarray(ds.variables["z"][:], float))


def asagi_fields(path: str | Path) -> list[str]:
    """The compound variable's member names, in file order."""
    from netCDF4 import Dataset
    with Dataset(str(path), "r") as ds:
        if "data" not in ds.variables:
            raise AsagiError(f"{path}: no variable named 'data'")
        dt = ds.variables["data"].datatype
        names = getattr(dt, "dtype", None)
        if names is None or names.names is None:
            raise AsagiError(f"{path}: 'data' is not a compound type")
        return list(names.names)


def read_asagi(path: str | Path, fields: Sequence[str] | None = None):
    """Return (x, y, z, {name: (nz,ny,nx) float64}, attrs)."""
    from netCDF4 import Dataset
    with Dataset(str(path), "r") as ds:
        x = np.asarray(ds.variables["x"][:], float)
        y = np.asarray(ds.variables["y"][:], float)
        z = np.asarray(ds.variables["z"][:], float)
        raw = ds.variables["data"][:]
        names = list(raw.dtype.names)
        if fields is not None:
            missing = [f for f in fields if f not in names]
            if missing:
                raise AsagiError(f"{path}: no such field(s) {missing}; have {names}")
            names = list(fields)
        out = {f: np.asarray(raw[f], float) for f in names}
        attrs = {k: ds.getncattr(k) for k in ds.ncattrs()}
    return x, y, z, out, attrs


def trilinear_sample(path: str | Path, qx, qy, qz,
                     fields: Sequence[str] | None = None,
                     return_clamped: bool = False):
    """Trilinear sample of an ASAGI nc at (qx, qy, qz), edge-clamped.

    The clamping is the legacy behaviour and is relied on (fault facets can sit a little
    outside the grid in z).  With return_clamped=True the boolean mask of points that were
    clamped on ANY axis is returned alongside, so callers can gate on it instead of
    silently accepting an extrapolated value.
    """
    x, y, z, flds, _ = read_asagi(path, fields)
    qx = np.atleast_1d(np.asarray(qx, float))
    qy = np.atleast_1d(np.asarray(qy, float))
    qz = np.atleast_1d(np.asarray(qz, float))
    if not (len(qx) == len(qy) == len(qz)):
        raise AsagiError(
            f"query arrays must be the same length, got {len(qx)}, {len(qy)}, {len(qz)}")

    def locate(g, q):
        i = np.clip(np.searchsorted(g, q) - 1, 0, len(g) - 2)
        w = (q - g[i]) / (g[i + 1] - g[i])
        return i, np.clip(w, 0.0, 1.0)

    ix, wx = locate(x, qx)
    iy, wy = locate(y, qy)
    iz, wz = locate(z, qz)

    out = {}
    for f, arr in flds.items():
        acc = np.zeros(len(qx))
        for dz_ in (0, 1):
            for dy_ in (0, 1):
                for dx_ in (0, 1):
                    w = ((wx if dx_ else 1 - wx) * (wy if dy_ else 1 - wy)
                         * (wz if dz_ else 1 - wz))
                    acc += w * arr[iz + dz_, iy + dy_, ix + dx_]
        out[f] = acc

    if not return_clamped:
        return out
    clamped = ((qx < x[0]) | (qx > x[-1]) | (qy < y[0]) | (qy > y[-1])
               | (qz < z[0]) | (qz > z[-1]))
    return out, clamped


# --------------------------------------------------------------------------- gates
def hull_containment(path: str | Path, pts: np.ndarray, label: str = "",
                     gate: str = "G1", severity: str = HARD,
                     report: GateReport | None = None) -> GateReport:
    """Gate G1, generalised: every point must lie inside the grid hull.

    Reports the per-axis margin, so a near-miss is visible before it becomes an
    edge-clamped (i.e. silently extrapolated) sample.
    """
    pts = np.atleast_2d(np.asarray(pts, float))
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise AsagiError(f"pts must be (N, 3), got {pts.shape}")
    x, y, z = asagi_axes(path)
    rep = report if report is not None else GateReport(f"hull containment {label}".strip())

    outside_total = 0
    margins = []
    for a, col, nm in ((x, 0, "x"), (y, 1, "y"), (z, 2, "z")):
        lo_gap = float(pts[:, col].min() - a[0])
        hi_gap = float(a[-1] - pts[:, col].max())
        n_out = int(((pts[:, col] < a[0]) | (pts[:, col] > a[-1])).sum())
        outside_total += n_out
        margins.append(f"{nm}: [{lo_gap:+.1f}, {hi_gap:+.1f}] m")
    detail = (f"{len(pts)} point(s) vs grid hull; per-axis margin (low, high) "
              f"{'; '.join(margins)}")
    if outside_total:
        detail = f"{outside_total} point(s) OUTSIDE the hull; {detail}"
    rep.add(gate, outside_total == 0, f"{label + ': ' if label else ''}{detail}",
            severity=severity)
    return rep


def roundtrip_selfcheck(path: str | Path, evaluator: Callable[[np.ndarray], np.ndarray],
                        field: str, n: int = 1000, seed: int = 12345,
                        median_tol: float = 1e-5, max_tol: float = 5e-4,
                        gate: str = "G3", report: GateReport | None = None,
                        bbox: tuple | None = None,
                        kink_straddling: Callable[[np.ndarray], np.ndarray] | None = None,
                        ) -> GateReport:
    """Gate G3, generalised: the written grid must reproduce an independent evaluation.

    Samples `n` fixed-seed points and compares `trilinear_sample(path)[field]` against
    `evaluator(pts)`.

    The legacy G3 contract, preserved here exactly:
      * median within `median_tol`  -- HARD;
      * max within `max_tol`        -- HARD **unless** every excursion sits in a grid cell
        whose corners straddle a kink in the source field.  A piecewise-linear field's
        trilinear interpolant is exact on any single branch, so only kink-straddling cells
        can legitimately disagree; an excursion in a smooth cell IS a real bug.

    `kink_straddling(pts) -> bool array` is how a caller attests to that.  Supplying
    nothing means "this field is smooth here", so a large excursion fails HARD -- the safe
    default.  Without this the max gate could never fail and a genuine interpolation error
    in a smooth region would be reported as a warning.
    """
    x, y, z = asagi_axes(path)
    lo = np.array([x[0], y[0], z[0]]) if bbox is None else np.asarray(bbox[0], float)
    hi = np.array([x[-1], y[-1], z[-1]]) if bbox is None else np.asarray(bbox[1], float)
    rng = np.random.default_rng(seed)
    pts = lo + (hi - lo) * rng.random((n, 3))

    got = trilinear_sample(path, pts[:, 0], pts[:, 1], pts[:, 2], fields=[field])[field]
    want = np.asarray(evaluator(pts), float)
    if want.shape != got.shape:
        raise AsagiError(
            f"evaluator returned shape {want.shape}, expected {got.shape}")
    err = np.abs(got - want)
    med, mx = float(np.median(err)), float(err.max())

    rep = report if report is not None else GateReport(f"round-trip {field}")
    rep.add(gate, med <= median_tol,
            f"{field}: median |nc - direct| = {med:.3e} (tol {median_tol:.0e}) "
            f"over {n} seed-{seed} points")

    if mx <= max_tol:
        rep.add(f"{gate}max", True,
                f"{field}: max |nc - direct| = {mx:.3e} (tol {max_tol:.0e})")
    else:
        bad = np.flatnonzero(err > max_tol)
        excused = False
        if kink_straddling is not None:
            excused = bool(np.all(np.asarray(kink_straddling(pts[bad]), bool)))
        rep.add(f"{gate}max", excused,
                f"{field}: max |nc - direct| = {mx:.3e} > tol {max_tol:.0e} at "
                f"{bad.size} of {n} point(s); "
                + ("every excursion sits in a cell straddling a kink in the source "
                   "field, which is expected for a piecewise-linear field"
                   if excused else
                   "no kink_straddling attestation was supplied, so these are treated "
                   "as real interpolation errors"),
                severity=WARN if excused else HARD)
    return rep
