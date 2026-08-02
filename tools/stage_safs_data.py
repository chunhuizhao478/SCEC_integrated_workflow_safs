#!/usr/bin/env python3
"""stage_safs_data.py -- collect the SAFS raw inputs into data/safs_alt/.

The raw products live scattered across the legacy tree.  This puts them where the
`safs_alt` descriptor expects them, so `safs_check.ipynb` can run end to end.

    python tools/stage_safs_data.py [--src <seisol_quakeworx>] [--link]

Nothing here is committed: data/*/raw/ is gitignored and the whole staged set is ~124 MB.
Use --link to symlink instead of copy (faster, but the deck you build then references a
link -- fine for a check run, never for a deck you ship).

WHAT GETS STAGED, and why each one:

  raw/CSM_orientation.csv          the Yang & Hauksson (YHSM-2013) community stress model.
                                   Byte-identical to the copy the legacy toolbox shipped
                                   (md5-verified), so the orientation path has no hidden
                                   preprocessing step.
  raw/cvm/CVM_*_h_data.csv         37 horizontal slices of the statewide "muscal" CVM
                                   (lon, lat, Vp, Vs, density), 0 - 70 km depth.
  raw/ctm/CTM_*_h_data_final.csv   105 slices of the SCEC Community Thermal Model
                                   (shinevar2024), 0 - 21 km, 200 m spacing.
  mesh_alt.puml.h5                 the ALT fault mesh -- the same file the legacy
                                   _common.ALT_MESH pointed at.
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_SRC = (Path.home() / "projects/seas-mfem-spatial-dyn-driver/miniapps/seas/safs"
               / "seisol_quakeworx")

# (destination relative to data/safs_alt/, source glob relative to --src, what it is)
PLAN = [
    ("raw/CSM_orientation.csv", "raw_data/yang_and_hauksson_orientation/CSM_data_*.csv",
     "CSM YHSM-2013 stress orientation"),
    ("raw/cvm/", "raw_data/multiscale_statewise_cvm/CVM_*_h_data.csv",
     "statewide CVM velocity slices"),
    ("raw/ctm/", "thermal/raw/CTM_*_h_data_final.csv",
     "SCEC Community Thermal Model slices"),
    ("mesh_alt.puml.h5", "safs_seisol_v3_2_0_RSSRW_ALTv2_freesurf/mesh_alt.puml.h5",
     "the ALT fault mesh"),
]


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default=str(DEFAULT_SRC))
    ap.add_argument("--project", default="safs_alt")
    ap.add_argument("--link", action="store_true", help="symlink instead of copy")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args(argv)

    src = Path(a.src).expanduser()
    if not src.is_dir():
        print(f"ERROR: source tree not found: {src}", file=sys.stderr)
        return 2
    dst_root = ROOT / "data" / a.project
    total = 0
    print(f"source : {src}\ntarget : {dst_root}\nmode   : "
          f"{'symlink' if a.link else 'copy'}\n")

    for dest, pattern, what in PLAN:
        hits = sorted(src.glob(pattern))
        if not hits:
            print(f"  MISSING  {what}\n           no match for {pattern}")
            return 1
        out = dst_root / dest
        if dest.endswith("/"):
            out.mkdir(parents=True, exist_ok=True)
            n = 0
            for h in hits:
                t = out / h.name
                if t.exists() and not a.force:
                    continue
                _place(h, t, a.link)
                n += 1
            size = sum(h.stat().st_size for h in hits)
            total += size
            print(f"  OK       {what}\n           {len(hits)} file(s), "
                  f"{size/1e6:.1f} MB -> {out.relative_to(ROOT)} ({n} new)")
        else:
            out.parent.mkdir(parents=True, exist_ok=True)
            h = hits[0]
            if len(hits) > 1:
                print(f"           ! {len(hits)} matches for {pattern}; using {h.name}")
            if not out.exists() or a.force:
                _place(h, out, a.link)
            total += h.stat().st_size
            print(f"  OK       {what}\n           {h.name} "
                  f"({h.stat().st_size/1e6:.1f} MB) -> {out.relative_to(ROOT)}")
            if dest.endswith(".csv"):
                print(f"           md5 {md5(out)}")

    print(f"\nstaged {total/1e6:.0f} MB into {dst_root.relative_to(ROOT)}")
    print("gitignored -- data/*/raw/ and data/**/*.puml.h5 are excluded.\n")
    print("Now run:  jupyter lab safs_check.ipynb")
    return 0


def _place(srcf: Path, dstf: Path, link: bool) -> None:
    if dstf.exists() or dstf.is_symlink():
        dstf.unlink()
    if link:
        dstf.symlink_to(srcf.resolve())
    else:
        shutil.copy2(srcf, dstf)


if __name__ == "__main__":
    raise SystemExit(main())
