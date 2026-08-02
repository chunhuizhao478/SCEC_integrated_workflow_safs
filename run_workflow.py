#!/usr/bin/env python3
"""run_workflow.py -- the workflow end to end, for the stages that exist.

    python run_workflow.py [--project demo_planar] [--out outputs/<project>/<tag>]

Walks the real dependency chain:

    raw velocity  -> material nc -> plasticity nc
                                 -> Sv(z) profile ---.
    raw orientation --------------------------------+-> stress nc
    thermal / depth profile ----------------------------> friction nc
    along-strike design --------------------------------> rs_muw LuaMap
                                                          (+ Tnuc_s, given a mesh)

and gates every stage as it goes.  Order is not cosmetic: the stress closure needs Sv
from the material, and the friction zoning needs the thermal field, so material MUST run
first.  Running them out of order is the likeliest way to pair a new velocity model with
a stale stress field.

Not yet wired (later phases): mesh ingest (5), deck assembly (6), the notebook (7).
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from deckbuild.config import Project                                    # noqa: E402
from deckbuild.contract import GateReport, Manifest                     # noqa: E402
from deckbuild.friction import FrictionStage, FwDesign                  # noqa: E402
from deckbuild.material import (                                        # noqa: E402
    AttenuationSpec, MaterialStage, PlasticitySpec,
)
from deckbuild.deck import DeckSpec, DeckStage                          # noqa: E402
from deckbuild.mesh import MeshStage                                    # noqa: E402
from deckbuild.stage_f import verify_mesh_against_deck                  # noqa: E402
from deckbuild.stress import KDesign, StressStage                       # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", default="demo_planar")
    ap.add_argument("--out", default=None)
    ap.add_argument("--k", type=float, default=1.7)
    ap.add_argument("--case", type=int, default=1)
    ap.add_argument("--freeze-above-depth-m", type=float, default=500.0)
    ap.add_argument("--phi-soft", type=float, default=35.0)
    ap.add_argument("--phi-hard", type=float, default=45.0)
    ap.add_argument("--require-files", action="store_true",
                    help="fail if any descriptor-referenced input is missing")
    a = ap.parse_args(argv)

    desc = ROOT / "projects" / f"{a.project}.yaml"
    cfg = Project.load(desc, require_files=a.require_files)
    tag = f"k{a.k:g}_case{a.case}_fz{a.freeze_above_depth_m:g}"
    out = Path(a.out) if a.out else ROOT / "outputs" / cfg.name / tag
    for sub in ("material", "stress", "friction"):
        (out / sub).mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print(cfg.summary())
    print(f"run tag      : {tag}")
    print(f"outputs      : {out}")
    print("=" * 78)

    man = Manifest.for_project(cfg)
    all_ok = True

    def stage(label, fn):
        nonlocal all_ok
        t0 = time.time()
        print(f"\n--- {label} " + "-" * (72 - len(label)))
        res = fn()
        dt = time.time() - t0
        print(f"    ({dt:.2f} s)")
        return res, dt

    # -- 1. MATERIAL (first: stress needs its Sv, friction may need its thermal) -----
    def _material():
        st = MaterialStage()
        arts = st.build(cfg, out / "material",
                        plasticity=PlasticitySpec(phi_soft_deg=a.phi_soft,
                                                  phi_hard_deg=a.phi_hard),
                        attenuation=AttenuationSpec())
        rep = st.verify(cfg, arts, plasticity_artifact=arts.plasticity)
        rep.print()
        return arts, rep
    (mat, mrep), mdt = stage("[1] MATERIAL  raw velocity -> nc, plasticity, Sv(z)",
                             _material)
    all_ok &= mrep.ok
    for art in (mat.material, mat.plasticity, mat.sv_profile):
        man.record("material", artifact=art)
    man.record("material", report=mrep, wall_s=round(mdt, 3))

    # -- 2. STRESS ------------------------------------------------------------------
    def _stress():
        st = StressStage()
        art = st.build(cfg, out / "stress", sv_profile=mat.sv_profile,
                       design=KDesign(k_values=(a.k,)),
                       freeze_above_depth_m=a.freeze_above_depth_m)
        rep = st.verify(cfg, art)
        rep.print()
        return art, rep
    (sart, srep), sdt = stage("[2] STRESS    orientation x Sv x k -> nc", _stress)
    all_ok &= srep.ok
    man.record("stress", artifact=sart, report=srep, wall_s=round(sdt, 3))

    # -- 3. FRICTION ----------------------------------------------------------------
    def _friction():
        st = FrictionStage()
        art = st.build(cfg, out / "friction", case=a.case)
        rep = st.verify(cfg, art)
        rep.print()
        return art, rep
    (fart, frep), fdt = stage("[3] FRICTION  thermal/depth profile -> rs_a, rs_srW",
                              _friction)
    all_ok &= frep.ok
    man.record("friction", artifact=fart, report=frep, wall_s=round(fdt, 3))

    # -- 4. GRADED f_w --------------------------------------------------------------
    def _fw():
        st = FrictionStage()
        d = FwDesign(fw_values=(0.0, 0.05),
                     boundaries_s_km=((10.0, 20.0),), f0=0.6)
        art = st.build_fw_map(cfg, out / "friction", d)
        rep = st.verify_fw_map(cfg, art, d)
        rep.print()
        return art, rep
    (wart, wrep), wdt = stage("[4] rs_muw    along-strike design -> !LuaMap", _fw)
    all_ok &= wrep.ok
    man.record("friction", artifact=wart, report=wrep, wall_s=round(wdt, 3))

    # -- 5. MESH: ingest and gate ----------------------------------------------------
    mesh_path = cfg.resolve_path(cfg.mesh().path)
    have_mesh = Path(mesh_path).is_file()

    def _mesh():
        st = MeshStage()
        if not have_mesh:
            rep = GateReport("mesh")
            rep.skip("A-G", f"no mesh at {mesh_path}.  This workflow INGESTS a mesh; "
                            f"build yours with the skill in skills/ (MESHING.md).")
            rep.print()
            return None, rep
        art = st.build(cfg, out)
        rep = st.verify(cfg, art, material_nc=mat.material.path,
                        stress_nc=sart.path, friction_nc=fart.path)
        rep.print()
        return art, rep
    (meshart, mrep2), mdt2 = stage("[5] MESH      ingest -> convert / gate (A-G)", _mesh)
    all_ok &= mrep2.ok
    man.record("mesh", artifact=meshart, report=mrep2, wall_s=round(mdt2, 3))

    # -- 6. STAGE F: is the mesh compatible with these fields? -----------------------
    def _stage_f():
        if not have_mesh:
            rep = GateReport("Stage F")
            rep.skip("F", "no mesh; the mesh-vs-deck check did not run")
            rep.print()
            return rep
        rep = verify_mesh_against_deck(
            cfg, stress_nc=sart.path, friction_nc=fart.path,
            material_nc=mat.material.path, descriptor_sha256=cfg.sha256())
        rep.print()
        return rep
    frep2, fdt2 = stage("[6] STAGE F   mesh vs the deck's actual fields", _stage_f)
    all_ok &= frep2.ok
    man.record("stage_f", report=frep2, wall_s=round(fdt2, 3))

    # -- 7. ASSEMBLE the deck + pre-flight P1-P8 -------------------------------------
    def _deck():
        arts = {"material": mat.material, "plasticity": mat.plasticity,
                "stress": sart, "friction": fart, "mesh": meshart}
        arts = {k: v for k, v in arts.items() if v is not None}
        # z = -1 m, never 0: SeisSol v1.1.3 silently DROPS a receiver at the surface.
        rx = np.array([[cfg.stress_box.xmin + 5000.0 * i,
                        cfg.stress_box.ymin + 5000.0 * i, -1.0] for i in range(1, 6)])
        spec = DeckSpec(prefix=f"{cfg.name}_", plasticity=True,
                        attenuation=mat.attenuation, receivers=rx,
                        rs_muw_lua=Path(wart.path).read_text(), mu_s=0.6,
                        notes="")
        st = DeckStage()
        dpath = st.assemble(cfg, deck, arts, spec, overwrite=True)
        rep = st.preflight(cfg, dpath, spec=spec)
        rep.print()
        return dpath, rep
    deck = ROOT / "decks" / f"{cfg.name}_{tag}"
    (dpath, drep), ddt = stage("[7] DECK      assemble + pre-flight (P1-P8)", _deck)
    all_ok &= drep.ok
    man.record("deck", report=drep, wall_s=round(ddt, 3), path=str(dpath))

    mpath = man.write(out)
    print("\n" + "=" * 78)
    print(f"manifest     : {mpath}")
    print("artifacts:")
    for e in man.artifacts():
        print(f"  {e['kind']:<12} {Path(e['path']).name:<44} {e['sha256'][:12]}")
    print(f"deck         : {dpath}")
    print("=" * 78)
    print("WORKFLOW OK" if all_ok else "WORKFLOW: SOME HARD GATES FAILED")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
