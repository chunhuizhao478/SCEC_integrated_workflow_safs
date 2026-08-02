#!/usr/bin/env python3
"""Generate data/demo/demo_planar.puml.h5 -- the committed demo mesh.

A vertical strike-slip fault embedded in a tet slab, sized to the demo_planar descriptor.
Deliberately coarse: it demonstrates the pipeline, not a production resolution.  Shipping
it (rather than generating at run time) keeps gmsh off the critical path and makes the
demo deterministic.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from deckbuild.mesh import orient_tets_positive, write_puml          # noqa: E402
from deckbuild.geometry import FACE_VERTS                            # noqa: E402

# The demo fault: x in [485, 515] km (30 km along strike is plenty), y = 500000 (the
# strike origin), z from 0 down to -15 km.  Strike azimuth is 0, so s runs with y --
# the fault plane is x = const?  No: a vertical plane containing the strike direction.
# Strike az = 0 means s increases northward, so the fault strikes N-S: the plane is
# x = x0, and the slab is a few cells thick in x.
X0, DY, DZ = 500000.0, 2500.0, 2500.0
Y0, Y1 = 3985000.0, 4025000.0          # 40 km along strike, inside the stress box
ZTOP, ZBOT = 0.0, -15000.0
DX = 2500.0


def build():
    ys = np.arange(Y0, Y1 + 0.5 * DY, DY)
    zs = np.arange(ZTOP, ZBOT - 0.5 * DZ, -DZ)
    xs = np.array([X0 - 2 * DX, X0 - DX, X0, X0 + DX, X0 + 2 * DX])

    pts, idx = [], {}
    for k, zz in enumerate(zs):
        for j, yy in enumerate(ys):
            for i, xx in enumerate(xs):
                idx[(i, j, k)] = len(pts)
                pts.append((xx, yy, zz))
    pts = np.asarray(pts, float)

    KUHN = [(0, 1, 3, 7), (0, 1, 7, 5), (0, 5, 7, 4),
            (0, 3, 2, 7), (0, 6, 4, 7), (0, 2, 6, 7)]
    tets = []
    for k in range(len(zs) - 1):
        for j in range(len(ys) - 1):
            for i in range(len(xs) - 1):
                c = [idx[(i + (n & 1), j + ((n >> 1) & 1), k + ((n >> 2) & 1))]
                     for n in range(8)]
                for t in KUHN:
                    tets.append([c[t[0]], c[t[1]], c[t[2]], c[t[3]]])
    tets = np.asarray(tets, np.int64)
    tets, n_flip = orient_tets_positive(pts, tets)

    on_fault = np.isclose(pts[:, 0], X0)
    xlo, xhi = np.isclose(pts[:, 0], xs[0]), np.isclose(pts[:, 0], xs[-1])
    ylo, yhi = np.isclose(pts[:, 1], ys[0]), np.isclose(pts[:, 1], ys[-1])
    top, bot = np.isclose(pts[:, 2], ZTOP), np.isclose(pts[:, 2], zs[-1])

    bc = np.zeros((len(tets), 4), np.int64)
    for f, fv in enumerate(FACE_VERTS):
        tri = tets[:, list(fv)]
        bc[tri_all(on_fault, tri), f] = 3                        # 101 fault -> 3
        for m in (xlo, xhi, ylo, yhi):
            bc[tri_all(m, tri), f] = 5                           # 104 sides -> 5
        bc[tri_all(bot, tri), f] = 5                             # 103 bottom -> 5
        bc[tri_all(top, tri), f] = 1                             # 102 top -> 1
    # The fault is INTERIOR; a free-surface/side face can never also be the fault.
    for f, fv in enumerate(FACE_VERTS):
        bc[tri_all(on_fault, tets[:, list(fv)]), f] = 3
    return pts, tets, bc, n_flip


def tri_all(mask, tri):
    return mask[tri].all(axis=1)


def main():
    pts, tets, bc, n_flip = build()
    out = ROOT / "data" / "demo_planar" / "demo_planar.puml.h5"
    write_puml(out, pts, tets, bc)
    n_fault = int((bc == 3).sum())
    print(f"wrote {out}")
    print(f"  {len(pts):,} nodes, {len(tets):,} tets, {n_fault:,} fault faces, "
          f"{n_flip} tets reoriented, {out.stat().st_size / 1024:.0f} KiB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
