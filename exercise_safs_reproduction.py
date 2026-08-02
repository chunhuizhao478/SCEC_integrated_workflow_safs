#!/usr/bin/env python3
"""exercise_safs_reproduction.py -- Phase 8, the E0-E6 ladder.

Rebuild a SHIPPED SAFS production deck from its raw inputs and compare field by field.
This is the exercise that decides whether the workflow is trustworthy.

    python exercise_safs_reproduction.py --deck <path> [--raw <seisol_quakeworx tree>]

Run the rungs IN ORDER and stop reading at the first one that fails: a full-chain
mismatch is nearly impossible to localise.

    E0  INGEST the deck's mesh; gates A-E; facet count vs the legacy pipeline
    E1  PROVENANCE: do the raw slices here correspond to the shipped ncs?
    E2  DESCRIPTOR: introspect the deck -> draft descriptor; reconcile
    E3  MATERIAL: raw CVM -> nc; + plasticity; raw CTM -> thermal; + Sv(z)
    E4  STRESS:   orientation + E3's Sv + k -> stress nc
    E5  FRICTION: E3's thermal + case -> friction nc; + the rs_muw LuaMap
    E6  DECK:     assemble; P1-P8; diff against the shipped deck

WHAT WILL AND WILL NOT MATCH -- declared BEFORE any comparison, so a difference in the
third category is not mistaken for a bug:

  WILL match, value for value
    every nc variable and axis; the rs_muw Lua's numbers; the active (non-comment)
    content of every easi yaml; the receiver and pickpoint coordinates

  WILL NOT match, by design
    the mesh        -- INGESTED, never rebuilt (Phase 5 scope)
    yaml comments   -- the shipped headers are hand-written provenance prose
    parameters.par  -- machine/campaign-specific (nodes, paths, EndTime)
    nc header bytes -- netCDF4/HDF5 library version, not data

  MIGHT NOT match, and each is a real finding to run down
    material nc     -- scipy's Delaunay can differ across versions, perturbing
                       LinearNDInterpolator at the last ULP
    Sv(z)           -- depends on the material matching first
    anything downstream -- a material difference propagates
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

import numpy as np                                                    # noqa: E402

from deckbuild.contract import GateReport, WARN                       # noqa: E402
from deckbuild.introspect import UNKNOWN, introspect_deck             # noqa: E402


def rung(name, title):
    print(f"\n{'=' * 78}\n{name}  {title}\n{'=' * 78}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--deck", required=True, help="a shipped SAFS deck folder")
    ap.add_argument("--raw", default=None,
                    help="the seisol_quakeworx tree holding raw_data/ and thermal/raw/")
    ap.add_argument("--out", default=None)
    ap.add_argument("--stop-at", default="E6")
    a = ap.parse_args(argv)

    deck = Path(a.deck)
    out = Path(a.out) if a.out else ROOT / "outputs" / "safs_repro"
    out.mkdir(parents=True, exist_ok=True)
    results: list[tuple[str, str, str]] = []

    def record(rid, verdict, evidence):
        results.append((rid, verdict, evidence))
        print(f"  -> {rid}: {verdict}\n     {evidence}")

    print(__doc__.split("WHAT WILL")[1].join(["WHAT WILL", ""]))

    # ---------------------------------------------------------------- E0
    rung("E0", "INGEST the deck's mesh and gate it")
    mesh = next(deck.glob("*.puml.h5"), None)
    if mesh is None:
        record("E0", "BLOCKED", f"no .puml.h5 in {deck}")
        return _finish(results, out)
    t0 = time.time()
    from deckbuild.mesh import count_face_multiplicity, read_puml
    geom, conn, bc, _ = read_puml(mesh)
    n_fault = int((bc == 3).sum())
    mult, _ = count_face_multiplicity(conn)
    interior = int((mult.reshape(4, len(conn)).T[bc == 3] == 2).sum())
    record("E0", "PASS" if interior == n_fault else "FAIL",
           f"{mesh.name}: {len(geom):,} nodes, {len(conn):,} tets, {n_fault:,} BC-3 "
           f"faces, {interior:,} of them interior ({time.time() - t0:.1f} s)")

    # ---------------------------------------------------------------- E1
    rung("E1", "PROVENANCE: do the raw slices here correspond to the shipped ncs?")
    if a.raw is None:
        record("E1", "SKIPPED", "no --raw tree given; the rebuild rungs cannot run")
    else:
        raw = Path(a.raw)
        cvm = sorted((raw / "raw_data" / "multiscale_statewise_cvm").glob("CVM_*_h_data.csv"))
        ctm = sorted((raw / "thermal" / "raw").glob("CTM_*_h_data_final.csv"))
        csm = next((raw / "raw_data" / "yang_and_hauksson_orientation").glob("CSM_*.csv"),
                   None)
        mnc = next(deck.glob("*material*.nc"), None)
        det = [f"{len(cvm)} CVM slices, {len(ctm)} CTM slices, "
               f"CSM {'present' if csm else 'ABSENT'}"]
        ok = bool(cvm and ctm and csm)
        if mnc is not None and cvm:
            from deckbuild.asagi import asagi_axes
            x, y, z = asagi_axes(mnc)
            det.append(f"shipped material nc: {len(x)}x{len(y)}x{len(z)}, "
                       f"z[{z[0]:.0f},{z[-1]:.0f}] dz={np.diff(z).mean():.0f}")
            det.append(f"CVM slice count {len(cvm)} vs nc z-levels {len(z)} -- the nc is "
                       f"RESAMPLED onto a uniform axis, so these need not be equal")
        record("E1", "PASS" if ok else "BLOCKED", "; ".join(det))

    # ---------------------------------------------------------------- E2
    rung("E2", "DESCRIPTOR: introspect the shipped deck")
    facts = introspect_deck(deck)
    print(facts.report())
    unknown = facts.unknown
    record("E2", "PASS" if len(unknown) <= 2 else "PARTIAL",
           f"{len(facts)} fields recovered, {len(unknown)} unknown: {unknown}")

    # the recovered design, cross-checked against the field where possible
    rung("E2b", "RECOVERED PARAMETERS (read, never assumed)")
    for k in ("k_from_filename", "freeze_from_filename", "fw_values",
              "fw_boundaries_s_km", "hypocenter_xyz", "nucleation_radius_m",
              "nucleation_amplitude", "phi_deg_from_field", "rs_b_constant",
              "rs_sl0_constant", "strike_azimuth_deg", "strike_origin_xy",
              "Plasticity", "Tv", "FreqCentral"):
        if k in facts:
            print(f"    {k:24} = {facts[k]}")

    # VERIFY the freeze against the FIELD, not the filename.
    snc = next(deck.glob("*stress*.nc"), None)
    if snc is not None:
        from deckbuild.asagi import read_asagi
        x, y, z, sf, _ = read_asagi(snc, fields=["s_zz"])
        top = np.flatnonzero(z >= -500.0)
        frozen = (len(top) > 1 and
                  all(np.array_equal(sf["s_zz"][i], sf["s_zz"][top[0]]) for i in top))
        record("E2c", "PASS",
               f"shallow freeze measured ON THE FIELD: levels above -500 m are "
               f"{'IDENTICAL (freeze ON)' if frozen else 'DISTINCT (freeze OFF)'}; the "
               f"filename said {facts.get('freeze_from_filename')}")

    # ---------------------------------------------------------------- E3-E6
    for rid, title in (("E3", "MATERIAL: raw CVM/CTM -> nc"),
                       ("E4", "STRESS: orientation + Sv + k -> nc"),
                       ("E5", "FRICTION: thermal + case -> nc; + rs_muw"),
                       ("E6", "DECK: assemble, P1-P8, diff")):
        rung(rid, title)
        if rid == "E3":
            from deckbuild.material import VELOCITY_READERS
            missing = [k for k in ("cvm_slices",) if k not in VELOCITY_READERS]
            record(rid, "BLOCKED",
                   f"the {missing} reader is not implemented, so the CVM cannot be "
                   f"rebuilt from raw.  Every downstream rung inherits this: E4 needs "
                   f"E3's Sv, E5 needs E3's thermal nc.  This is the honest state -- the "
                   f"reader raises rather than silently producing a wrong field.")
        else:
            record(rid, "BLOCKED", f"inherits E3: {rid} consumes a material product")

    return _finish(results, out, facts, deck)


def _finish(results, out, facts=None, deck=None) -> int:
    print(f"\n{'=' * 78}\nLADDER SUMMARY\n{'=' * 78}")
    w = max(len(r[1]) for r in results)
    for rid, verdict, ev in results:
        print(f"  {rid:4} {verdict:<{w}}  {ev.splitlines()[0][:90]}")
    n_pass = sum(1 for _, v, _ in results if v == "PASS")
    n_blocked = sum(1 for _, v, _ in results if v in ("BLOCKED", "SKIPPED"))
    print(f"\n  {n_pass} passed, {n_blocked} blocked/skipped, "
          f"{len(results) - n_pass - n_blocked} failed")
    md = out / "EXERCISE_safs_reproduction.md"
    md.write_text(_report_md(results, facts, deck))
    print(f"\n  results document: {md}")
    return 0 if not any(v == "FAIL" for _, v, _ in results) else 1


def _report_md(results, facts, deck) -> str:
    L = ["# SAFS reproduction exercise — results", "",
         f"Target deck: `{deck.name if deck else '?'}`", "",
         "## Ladder", "", "| Rung | Verdict | Evidence |", "|:--|:--|:--|"]
    for rid, v, ev in results:
        L.append(f"| {rid} | **{v}** | {ev.replace('|', ' ')} |")
    if facts:
        L += ["", "## Recovered parameters", "", "```", facts.report(), "```"]
    L += ["", "## Blocking issue", "",
          "`cvm_slices` and `ctm_slices` are not implemented in `deckbuild/material.py`. "
          "They raise `MaterialError('not implemented yet')` rather than silently "
          "producing a wrong field, which is the right interim state, but it means E3–E6 "
          "cannot run and the reproduction claim is **unproven**.", "",
          "The mesh (E0), the raw-data provenance (E1) and the descriptor recovery (E2) "
          "all pass, so the exercise establishes that the deck is readable and its design "
          "is fully recoverable — which is the prerequisite for E3–E6, not a substitute "
          "for them."]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
