"""nc_compare.py -- the two-tier NetCDF comparison harness (Phase 1 requirement 4).

`data_identical` is the acceptance criterion; `file_identical` is a bonus.

A NetCDF file's BYTES depend on the netCDF4/HDF5 library version -- chunking, compression
and header layout can all change without a single stored value changing.  Requiring byte
equality would make the suite fail on a library upgrade for a reason that affects no
simulation.  `data_identical` catches every real defect: a changed value, a changed axis,
a transposed field, a shifted grid.  No tolerance is involved -- every value must match
exactly.
"""
from __future__ import annotations

import filecmp
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

__all__ = ["NcDiff", "compare_nc", "assert_data_identical"]


@dataclass
class NcDiff:
    a: Path
    b: Path
    data_identical: bool = False
    file_identical: bool = False
    first_diff: tuple | None = None          # (variable, flat index, value_a, value_b)
    attr_diffs: dict = field(default_factory=dict)
    field_diffs: dict = field(default_factory=dict)   # name -> (max_abs, max_rel, n_diff)
    missing_in_a: list = field(default_factory=list)
    missing_in_b: list = field(default_factory=list)
    axis_diffs: list = field(default_factory=list)
    error: str = ""

    def summary(self) -> str:
        if self.error:
            return f"NcDiff ERROR: {self.error}"
        head = (f"data_identical={self.data_identical} "
                f"file_identical={self.file_identical}")
        if self.data_identical:
            return head
        parts = [head]
        if self.axis_diffs:
            parts.append(f"  axes differ: {self.axis_diffs}")
        if self.missing_in_a or self.missing_in_b:
            parts.append(f"  fields only in b: {self.missing_in_a}; "
                         f"only in a: {self.missing_in_b}")
        for name, (mabs, mrel, n) in sorted(self.field_diffs.items()):
            parts.append(f"  {name}: {n} differing value(s), "
                         f"max|d|={mabs:.6g}, max rel={mrel:.6g}")
        if self.first_diff:
            v, i, va, vb = self.first_diff
            parts.append(f"  first difference: {v}[flat {i}] {va!r} != {vb!r}")
        return "\n".join(parts)


def _axes_and_fields(path: Path):
    from netCDF4 import Dataset
    with Dataset(str(path), "r") as ds:
        axes = {n: np.asarray(ds.variables[n][:], float)
                for n in ("x", "y", "z") if n in ds.variables}
        flds = {}
        if "data" in ds.variables:
            raw = ds.variables["data"][:]
            names = raw.dtype.names or ()
            for n in names:
                flds[n] = np.asarray(raw[n])
        for n, v in ds.variables.items():
            if n in ("x", "y", "z", "data"):
                continue
            flds[n] = np.asarray(v[:])
        attrs = {k: ds.getncattr(k) for k in ds.ncattrs()}
    return axes, flds, attrs


def compare_nc(a: str | Path, b: str | Path,
               ignore_attrs: tuple[str, ...] = ()) -> NcDiff:
    """Compare two NetCDF files at both tiers.

    ignore_attrs -- attribute names expected to differ (timestamps, provenance strings).
                    Differences are still recorded in `attr_diffs`, just not fatal to
                    `data_identical`.
    """
    a, b = Path(a), Path(b)
    d = NcDiff(a=a, b=b)
    if not a.is_file() or not b.is_file():
        d.error = f"missing file: {'a=' + str(a) if not a.is_file() else ''} " \
                  f"{'b=' + str(b) if not b.is_file() else ''}".strip()
        return d

    d.file_identical = filecmp.cmp(str(a), str(b), shallow=False)

    try:
        ax_a, fl_a, at_a = _axes_and_fields(a)
        ax_b, fl_b, at_b = _axes_and_fields(b)
    except Exception as exc:                      # noqa: BLE001 - report, do not mask
        d.error = f"could not read: {type(exc).__name__}: {exc}"
        return d

    for name in sorted(set(ax_a) | set(ax_b)):
        va, vb = ax_a.get(name), ax_b.get(name)
        if va is None or vb is None:
            d.axis_diffs.append(f"{name}: present in only one file")
        elif va.shape != vb.shape:
            d.axis_diffs.append(f"{name}: shape {va.shape} vs {vb.shape}")
        elif not np.array_equal(va, vb):
            k = int(np.argmax(va != vb))
            d.axis_diffs.append(f"{name}: first differs at {k} ({va[k]!r} vs {vb[k]!r})")

    d.missing_in_a = sorted(set(fl_b) - set(fl_a))
    d.missing_in_b = sorted(set(fl_a) - set(fl_b))

    for name in sorted(set(fl_a) & set(fl_b)):
        va, vb = fl_a[name], fl_b[name]
        if va.shape != vb.shape:
            d.field_diffs[name] = (float("inf"), float("inf"), -1)
            if d.first_diff is None:
                d.first_diff = (name, -1, f"shape {va.shape}", f"shape {vb.shape}")
            continue
        if np.array_equal(va, vb):
            continue
        neq = va != vb
        n = int(neq.sum())
        fa, fb = va.astype(float).ravel(), vb.astype(float).ravel()
        diff = np.abs(fa - fb)
        denom = np.maximum(np.abs(fa), np.abs(fb))
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(denom > 0, diff / denom, 0.0)
        d.field_diffs[name] = (float(np.nanmax(diff)), float(np.nanmax(rel)), n)
        if d.first_diff is None:
            k = int(np.argmax(neq.ravel()))
            d.first_diff = (name, k, fa[k], fb[k])

    for k in sorted(set(at_a) | set(at_b)):
        va, vb = at_a.get(k, "<absent>"), at_b.get(k, "<absent>")
        same = bool(np.all(va == vb)) if not isinstance(va, str) else va == vb
        if not same:
            d.attr_diffs[k] = (va, vb)

    fatal_attrs = {k for k in d.attr_diffs if k not in ignore_attrs}
    d.data_identical = (not d.axis_diffs and not d.missing_in_a and not d.missing_in_b
                        and not d.field_diffs and not fatal_attrs)
    return d


def assert_data_identical(a, b, ignore_attrs: tuple[str, ...] = ()) -> NcDiff:
    d = compare_nc(a, b, ignore_attrs=ignore_attrs)
    assert d.data_identical, f"\n{d.summary()}"
    return d
