#!/usr/bin/env python3
"""make_safs_check.py -- GENERATE safs_check.ipynb from deck_workflow.ipynb.

The SAFS check is not a second notebook that happens to resemble the first.  It IS the
general workflow, run on SAFS.  If the two are maintained by hand they drift, and then a
green check stops meaning "the workflow reproduces a production deck" and starts meaning
"a notebook that looks like the workflow reproduces a production deck" -- which is worth
nothing.

So: every structural cell is COPIED verbatim from deck_workflow.ipynb.  The only edits are

  1. the `---- PARAMETERS ----` blocks, swapped for the shipped SAFS deck's values, and
  2. a `[C]` section appended, which introspects the shipped deck, ASSERTS the parameters
     declared above match it, and diffs the rebuilt fields against the shipped ones.

Ordering matters and is deliberate: the parameters are DECLARED as literals and only then
checked against the shipped deck.  Reading them out of the shipped deck first and then
"reproducing" them would be circular -- it would reproduce perfectly and prove nothing.

    python tools/make_safs_check.py [--check]

`--check` regenerates in memory and exits non-zero if the file on disk differs, which is
what tests/test_notebook.py uses to keep the two in sync.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "deck_workflow.ipynb"
DST = ROOT / "safs_check.ipynb"

SHIPPED = ("~/Downloads/seisol_quakeworx/"
           "safs_seisol_v4_0_0_RSSRW_ALT_THERMAL_CASE1_intermediate"
           "_plast_phi30_40_gradedfw_k1p70_nwredM7p8_sefw0_attenuation_deep40km")

# ---------------------------------------------------------------- parameter overrides
# Each entry: (unique substring identifying the cell, old PARAMETERS text, new text).
# Matching on the STEP header makes a silent mis-hit impossible -- a header that stops
# existing raises rather than quietly leaving the demo values in a SAFS notebook.
SUBS = [
    ("STEP 0 · locate the workflow", '''PROJECT       = "demo_planar"    # a descriptor stem under projects/
RAW_DATA_DIR  = None             # where YOUR raw data lives.  None -> data/<PROJECT>/
                                 # e.g. "/data/my_fault/raw_exports"
REQUIRE_FILES = True             # False skips existence checks, never schema checks''',
     '''PROJECT       = "safs_alt"       # the San Andreas ALT descriptor
RAW_DATA_DIR  = None             # None -> data/safs_alt/, what tools/stage_safs_data.py
                                 # fills with the 37 CVM + 105 CTM slices and the CSM csv
REQUIRE_FILES = True             # every raw input must exist -- this is a CHECK'''),

    ("STEP 1 · declare the output layout",
     'RUN_TAG = "k1.7_case1_fz500"',
     'RUN_TAG = "k1.70_case1_gradedfw_phi30_40_Q"    # the shipped deck\'s design'),

    ("STEP 3 · MATERIAL", '''WITH_PLASTICITY  = False   # off-fault Drucker-Prager yielding (Roten et al. 2014).''',
     '''WITH_PLASTICITY  = True    # ON: the shipped deck has Plasticity = 1'''),
    ("STEP 3 · MATERIAL", '''WITH_ATTENUATION = False   # intrinsic (anelastic) Q, from the CVM's own Vs.''',
     '''WITH_ATTENUATION = True    # ON: the shipped deck is a VISCOELASTIC run'''),
    ("STEP 3 · MATERIAL", '''PHI_SOFT, PHI_HARD = 35.0, 45.0   # friction angle, deg: soft rock / hard rock.
                                  # 35/45 is Roten et al. (2014).  The SAFS production
                                  # decks use 30/40 -- neither is a silent default.''',
     '''PHI_SOFT, PHI_HARD = 30.0, 40.0   # the SAFS PRODUCTION values, not Roten's 35/45.
                                  # Read off the shipped plasticity nc's bulkFriction
                                  # field in [C], not off its yaml prose, which has
                                  # been wrong before.'''),

    ("STEP 4 · STRESS", '''K_VALUES        = [1.7]     # one k per along-strike region
K_BOUNDARIES    = []        # exactly len(K_VALUES)-1 smoothstep windows, e.g. [(150,160)]
FREEZE_ABOVE_M  = 500.0     # 0 disables the shallow freeze''',
     '''K_VALUES        = [1.70]    # the shipped deck's closure ratio (from its FILENAME,
                            # verified against the FIELD in [C])
K_BOUNDARIES    = []        # a single k region -- no graded-k on this deck
FREEZE_ABOVE_M  = 0.0       # OFF.  The shipped nc carries no _freeze tag, and [C]
                            # measures the freeze ON THE FIELD rather than trusting
                            # the filename.'''),

    ("STEP 5 · FRICTION", '''CASE          = 1     # which thermal zoning to use -- see the table above.''',
     '''CASE          = 1     # the shipped deck is THERMAL_CASE1'''),
    ("STEP 5 · FRICTION", '''FW_VALUES     = [0.0, 0.05]     # the strong-rate-weakening FLOOR f_w, one per region,
                                # ordered along increasing s.  Lower f_w = weaker at high
                                # slip rate = LARGER dynamic stress drop.
FW_BOUNDARIES = [(10.0, 20.0)]  # exactly len(FW_VALUES)-1 smoothstep windows, in s (km)''',
     '''FW_VALUES     = [0.0, 0.0, 0.045, 0.03, 0.06, 0.0175, 0.05]   # the shipped 7-region
                                # graded floor.  These literals are DECLARED here and
                                # checked against the shipped rs_muw LuaMap in [C];
                                # reading them out of the deck first would make the
                                # reproduction circular.
FW_BOUNDARIES = [(20.0, 28.0), (150.0, 160.0), (176.0, 182.0),
                 (188.0, 194.0), (214.0, 226.0), (244.0, 252.0)]'''),

    ("STEP 6 · ASSEMBLE", '''MU_S          = 0.6        # static strength, for the t=0 pre-slip check
END_TIME_S    = 20.0
NUC_RADIUS_M  = 2000.0
NUC_AMP_MPA   = 75.0
DECK_NOTES    = ""         # preserved verbatim across regeneration''',
     '''MU_S          = 0.6        # static strength, for the t=0 pre-slip check
END_TIME_S    = 20.0       # campaign-specific; not part of the reproduction claim
NUC_RADIUS_M  = 2000.0     # the shipped Tnuc_s bell, recovered in [C]
NUC_AMP_MPA   = 75.0
DECK_NOTES    = "SAFS reproduction check -- see docs/EXERCISE_safs_reproduction.md"'''),
]

TITLE_MD = """# SAFS check — the same workflow, run on the San Andreas

This notebook is **generated** from `deck_workflow.ipynb` by `tools/make_safs_check.py`.
Every structural cell is copied verbatim; only the `PARAMETERS` blocks differ, and a
`[C]` section is appended that compares the result against a **shipped production deck**.

Do not edit this file by hand — edit `deck_workflow.ipynb` (or the generator) and re-run

```bash
python tools/make_safs_check.py
```

`tests/test_notebook.py` fails if the two drift, because a check that has drifted from the
workflow it is supposed to check is worse than no check at all.

**Read the order carefully.** The design parameters are *declared as literals* in the
PARAMETERS blocks and only *then* compared with the shipped deck. Reading them out of the
shipped deck first and then "reproducing" them would be circular.
"""

COMPARE_MD = """## [C] SAFS only — compare against a shipped production deck

Everything above is the general workflow. What follows is specific to this check: read
the shipped deck's own design out of its files, assert the parameters declared above match
it, and diff the rebuilt fields field by field.
"""


def _cells(nb):
    return nb["cells"]


def _src(c):
    return "".join(c["source"])


def _lines(s):
    """nbformat source list: one entry per line, newline-terminated except the last.

    Strip a trailing newline FIRST -- otherwise the split leaves an empty final element,
    which is legal JSON but differs from what every other cell looks like, and made the
    round-trip check report permanent drift.
    """
    s = s.rstrip("\n")
    lines = s.split("\n")
    return [ln + "\n" for ln in lines[:-1]] + [lines[-1]]


def _set(c, s):
    c["source"] = _lines(s)


def _code(s, ident):
    return {"cell_type": "code", "execution_count": None, "id": ident, "metadata": {},
            "outputs": [], "source": _lines(s)}


def _md(s, ident):
    return {"cell_type": "markdown", "id": ident, "metadata": {}, "source": _lines(s)}


COMPARE_CELLS = [
    ("md", "c-head", COMPARE_MD),
    ("code", "c-intro", '''# ══ STEP C1 · read the SHIPPED deck's design out of its own files ══════════════
#   introspect_deck()  deckbuild/introspect.py
#     nc axes            -> the four grids
#     easi !ConstantMap  -> rs_b, rs_sl0
#     the rs_muw Lua     -> the strike frame AND the f_w design (they exist ONLY there)
#     the Tnuc_s Lua     -> the hypocentre, radius, amplitude
#     bulkFriction field -> phi (NOT the yaml prose, which has been wrong)
#     parameters.par     -> Plasticity, Tv, the attenuation band
#   Anything it cannot infer is marked UNKNOWN, never guessed.
# ---- PARAMETERS ----------------------------------------------------------------
SHIPPED_DECK = ("''' + SHIPPED + '''")
# ----------------------------------------------------------------------------------
from deckbuild.introspect import introspect_deck
shipped = Path(SHIPPED_DECK).expanduser()
if not shipped.is_dir():
    print(f"  ! shipped deck not found at {shipped}; [C] will be skipped")
    shipped = None
facts = introspect_deck(shipped) if shipped else {}
for k in ("k_from_filename", "freeze_from_filename", "fw_values", "fw_boundaries_s_km",
          "hypocenter_xyz", "nucleation_radius_m", "nucleation_amplitude",
          "phi_deg_from_field", "rs_b_constant", "strike_azimuth_deg", "Plasticity",
          "FreqCentral"):
    if k in facts:
        print(f"  {k:22} = {facts[k]}")
if facts:
    print(f"\\n  {len(facts)} fields recovered, {len(facts.unknown)} unknown")'''),

    ("code", "c-assert", '''# ══ STEP C2 · do the DECLARED parameters match the shipped deck? ═══════════════
# The literals in the PARAMETERS blocks above were written by hand.  This cell is what
# makes that honest: every one is checked against the shipped deck's own files, and a
# mismatch is reported as a FAILURE, not smoothed over.
from deckbuild.contract import GateReport
rep = GateReport("declared parameters vs the shipped deck")
if facts:
    def _close(a, b, tol=1e-9):
        import numpy as _np
        return _np.allclose(_np.asarray(a, float), _np.asarray(b, float), atol=tol)

    rep.add("C2k", facts.get("k_from_filename") == K_VALUES[0],
            f"k: declared {K_VALUES[0]} vs shipped {facts.get('k_from_filename')}")
    rep.add("C2freeze", float(facts.get("freeze_from_filename", 0.0)) == FREEZE_ABOVE_M,
            f"shallow freeze: declared {FREEZE_ABOVE_M:g} m vs shipped "
            f"{facts.get('freeze_from_filename')} m")
    rep.add("C2fw", _close(facts.get("fw_values", []), FW_VALUES),
            f"f_w values: declared {FW_VALUES} vs shipped {facts.get('fw_values')}")
    rep.add("C2fwb", _close([list(w) for w in facts.get("fw_boundaries_s_km", [])],
                            [list(w) for w in FW_BOUNDARIES]),
            f"f_w windows: declared {FW_BOUNDARIES} vs shipped "
            f"{facts.get('fw_boundaries_s_km')}")
    rep.add("C2phi", _close(sorted(facts.get("phi_deg_from_field", [])),
                            sorted([PHI_SOFT, PHI_HARD]), tol=1e-3),
            f"phi (read off the FIELD): declared {[PHI_SOFT, PHI_HARD]} vs shipped "
            f"{facts.get('phi_deg_from_field')}")
    rep.add("C2plast", bool(facts.get("Plasticity", 0)) == WITH_PLASTICITY,
            f"plasticity: declared {WITH_PLASTICITY} vs shipped "
            f"Plasticity = {facts.get('Plasticity')}")
    rep.add("C2nuc", (facts.get("nucleation_radius_m") == NUC_RADIUS_M
                      and facts.get("nucleation_amplitude") == NUC_AMP_MPA * 1e6),
            f"nucleation: declared R = {NUC_RADIUS_M:g} m, {NUC_AMP_MPA:g} MPa vs shipped "
            f"R = {facts.get('nucleation_radius_m')} m, "
            f"{facts.get('nucleation_amplitude', 0) / 1e6:g} MPa")
    # The freeze is measured ON THE FIELD, because the filename is only a naming rule.
    from deckbuild.asagi import read_asagi
    _snc = next(shipped.glob("*stress*.nc"), None)
    if _snc is not None:
        _, _, _z, _sf, _ = read_asagi(_snc, fields=["s_zz"])
        _top = np.flatnonzero(_z >= -500.0)
        _frozen = (len(_top) > 1 and
                   all(np.array_equal(_sf["s_zz"][i], _sf["s_zz"][_top[0]]) for i in _top))
        rep.add("C2field", _frozen == (FREEZE_ABOVE_M > 0),
                f"freeze measured ON THE FIELD: levels above -500 m are "
                f"{'IDENTICAL (freeze ON)' if _frozen else 'DISTINCT (freeze OFF)'}; "
                f"declared {FREEZE_ABOVE_M:g} m")
    rep.print()
else:
    print("  skipped: no shipped deck")'''),

    ("code", "c-fields", '''# ══ STEP C3 · diff the rebuilt fields against the shipped ones ═════════════════
#   asagi_axes(), read_asagi()  deckbuild/asagi.py
#   The two-tier harness is tests/nc_compare.py; this cell inlines it so the numbers show.
from deckbuild.asagi import asagi_axes, read_asagi
if facts:
    s_nc = next(shipped.glob("*material_cvm.nc"), None)
    if s_nc:
        gx, gy, gz = asagi_axes(mat.material.path)
        sx, sy, sz = asagi_axes(s_nc)
        print(f"  material grid  built {len(gx)}x{len(gy)}x{len(gz)}   "
              f"shipped {len(sx)}x{len(sy)}x{len(sz)}")
        for n_, a_, b_ in (("x", gx, sx), ("y", gy, sy), ("z", gz, sz)):
            print(f"    {n_} axis identical: {np.array_equal(a_, b_)}")
        if (len(gx), len(gy), len(gz)) == (len(sx), len(sy), len(sz)):
            _, _, _, bf, _ = read_asagi(mat.material.path, fields=["mu", "rho"])
            _, _, _, sf2, _ = read_asagi(s_nc, fields=["mu", "rho"])
            rel = np.abs(bf["mu"] - sf2["mu"]) / np.maximum(np.abs(sf2["mu"]), 1e-30)
            print(f"    mu   median rel {np.median(rel):.3e}   max rel {rel.max():.3e}   "
                  f"{100 * np.mean(rel > 0.01):.2f}% of nodes differ by > 1%")
            print(f"    rho  exactly equal on {100 * np.mean(bf['rho'] == sf2['rho']):.2f}%"
                  f" of nodes")
            print("\\n    OPEN: the grid reproduces bit-for-bit, the VALUES do not.")
            print("    Localised to two effects in docs/EXERCISE_safs_reproduction.md:")
            print("      - a fixed ~2.5% interior stripe at EVERY depth -- the signature")
            print("        of a degenerate Delaunay triangulation (agrees exactly AT the")
            print("        sample points, differs between them);")
            print("      - one anomalous level, z = -250 m, the only shallow output level")
            print("        that is not itself a CVM slice depth -> the vertical resample.")'''),

    ("code", "c-lua", '''# ══ STEP C4 · is our LuaMap the SAME FUNCTION as the shipped one? ══════════════
#   _eval_lua_text()  deckbuild/friction.py  evaluates the SHIPPED Lua's own arithmetic,
#                                            including the legacy `X + inc * g` form
#   fw_profile_1d()   deckbuild/friction.py  what we emit
# A text diff would fail on whitespace and prove nothing.  Comparing the two as FUNCTIONS
# over the fault's s-range is the question that matters.
from deckbuild.friction import _eval_lua_text, fw_profile_1d
if facts:
    ship_yaml = next((p for p in shipped.glob("*rs_muw*.yaml")), None)
    if ship_yaml is not None:
        ss_ = np.linspace(-50.0, 500.0, 4000)
        ours = fw_profile_1d(ss_, fwdes)
        theirs = _eval_lua_text(ship_yaml.read_text(), ss_)
        print(f"  shipped: {ship_yaml.name}")
        print(f"  f_w(s) max |ours - shipped| over s in [-50, 500] km: "
              f"{np.max(np.abs(ours - theirs)):.3e}")
        print("  (0 means our emitted LuaMap is the SAME FUNCTION as the shipped one)")'''),

    ("code", "c-verdict", '''# ══ STEP C5 · verdict ══════════════════════════════════════════════════════════
print(f"deck : {deck}")
for f_ in sorted(deck.iterdir()):
    print(f"  {f_.name:52} {f_.stat().st_size / 1024 / 1024:7.1f} MB")
print()
print("A file-by-file diff against the shipped deck is NOT run: the two use different")
print("meshes (ALT vs the deep40km variant) and different EndTime/output settings.")
print("The meaningful comparisons are the FIELD ones in [C2]-[C4].")
print()
print("See docs/EXERCISE_safs_reproduction.md for the E0-E6 ladder and what is still")
print("open.  The workflow does NOT claim to reproduce a SAFS deck until the CVM value")
print("difference closes.")'''),
]


def build() -> dict:
    nb = json.loads(SRC.read_text())
    cells = _cells(nb)

    # 1. the title cell
    _set(cells[0], TITLE_MD)

    # 2. the PARAMETERS substitutions -- every one must hit
    for header, old, new in SUBS:
        hits = [c for c in cells if c["cell_type"] == "code" and header in _src(c)]
        if not hits:
            raise SystemExit(f"generator is stale: no cell contains {header!r}")
        cell = hits[0]
        s = _src(cell)
        if old not in s:
            raise SystemExit(
                f"generator is stale: in the cell for {header!r}, this text is gone:\n"
                f"---\n{old}\n---\n"
                f"deck_workflow.ipynb changed; update SUBS in tools/make_safs_check.py.")
        _set(cell, s.replace(old, new, 1))

    # 3. drop the trailing "---" cell, append the [C] section, restore the rule
    tail = cells.pop() if _src(cells[-1]).strip().startswith("---") else None
    for kind, ident, text in COMPARE_CELLS:
        cells.append(_md(text, ident) if kind == "md" else _code(text, ident))
    if tail is not None:
        cells.append(tail)
    return nb


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if the file on disk differs from the regeneration")
    a = ap.parse_args(argv)

    nb = build()
    text = json.dumps(nb, indent=1, ensure_ascii=False)
    if a.check:
        # Compare the PARSED notebooks, not the serialized text.  Escaping (\u2014 vs a
        # literal em dash) and indentation are not drift, and treating them as drift makes
        # the guard cry wolf until someone disables it -- which is how a real drift then
        # gets through.  What matters is the cell sources and their order.
        def shape(o):
            return [(c["cell_type"], "".join(c["source"])) for c in o["cells"]]

        try:
            cur = json.loads(DST.read_text()) if DST.is_file() else {"cells": []}
        except json.JSONDecodeError as exc:
            print(f"safs_check.ipynb is not valid JSON: {exc}", file=sys.stderr)
            return 1
        a_, b_ = shape(cur), shape(nb)
        if a_ != b_:
            for i, (x, y) in enumerate(zip(a_, b_)):
                if x != y:
                    print(f"first divergence at cell {i}:\n"
                          f"--- on disk ---\n{x[1][:300]}\n"
                          f"--- regenerated ---\n{y[1][:300]}", file=sys.stderr)
                    break
            else:
                print(f"cell COUNT differs: {len(a_)} on disk vs {len(b_)} regenerated",
                      file=sys.stderr)
            print("\nsafs_check.ipynb is OUT OF SYNC with deck_workflow.ipynb.\n"
                  "Run:  python tools/make_safs_check.py", file=sys.stderr)
            return 1
        print("safs_check.ipynb is in sync with deck_workflow.ipynb")
        return 0
    DST.write_text(text)
    print(f"wrote {DST}  ({len(nb['cells'])} cells, generated from {SRC.name})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
