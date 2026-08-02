"""stage_f.py -- the seam between meshing and this workflow.

The mesh skill's Stage F is "check a finished mesh against the CONSUMING deck's actual
inputs".  It is the one part of the mesh workflow that genuinely belongs in code, because
it is the only part that needs the other three ingredients.

These predicates are written ONCE here and composed by two callers:

  * `verify_mesh_against_deck` -- a CANDIDATE mesh with a deck's ncs, before a deck exists
    (what the skill calls at the end of a meshing campaign);
  * the Phase 6 pre-flight battery P1-P8 -- an ASSEMBLED deck.

Two batteries, one set of predicates.  If they were written twice they would drift and
disagree, which is worse than having only one.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from deckbuild.asagi import asagi_axes, hull_containment, trilinear_sample
from deckbuild.contract import GateReport, HARD, WARN
from deckbuild.geometry import load_fault, resolve_tractions, snap_hypocenter
from deckbuild.stress import COMP_IJ

__all__ = ["verify_mesh_against_deck", "project_stress_onto_fault",
           "check_finite_fields", "check_positive_normal_stress", "check_no_preslip",
           "check_receivers_under_the_surface"]


# --------------------------------------------------------------------------- predicates
def project_stress_onto_fault(stress_nc, fault):
    """(sigma_n, tau, mu_app) at facet centroids, compression-POSITIVE MPa."""
    _, _, zg = asagi_axes(stress_nc)
    smp = trilinear_sample(stress_nc, fault.cent[:, 0], fault.cent[:, 1],
                           np.clip(fault.cent[:, 2], zg.min(), zg.max()))
    sig = np.zeros((len(fault), 3, 3))
    for f, (a, b) in COMP_IJ.items():
        sig[:, a, b] = -smp[f] / 1.0e6
        sig[:, b, a] = sig[:, a, b]
    return resolve_tractions(sig, fault.strikes, fault.dips, fault.normals)


def check_finite_fields(nc, fault, label, rep: GateReport, gate: str) -> GateReport:
    """Every field must be finite at every facet centroid."""
    _, _, zg = asagi_axes(nc)
    smp = trilinear_sample(nc, fault.cent[:, 0], fault.cent[:, 1],
                           np.clip(fault.cent[:, 2], zg.min(), zg.max()))
    bad = {k: int((~np.isfinite(v)).sum()) for k, v in smp.items()}
    n = sum(bad.values())
    rep.add(gate, n == 0,
            f"{label}: {n} non-finite sample(s) at {len(fault):,} facet centroids"
            + (f" ({ {k: v for k, v in bad.items() if v} })" if n else ""))
    return rep


def check_positive_normal_stress(sn, tau, rep: GateReport, gate: str) -> GateReport:
    """sigma_n > 0 AND tau_0 > 0 on every facet.

    This is the psi_ini = -inf abort: the rate-and-state state init is
    a*ln[(2 sr0/V_ini) sinh(tau_0/(a sigma_n))], which is -inf at tau_0 = 0, and SeisSol
    dies with "Inf/NaN in energies" before the first step.
    """
    n_sn = int((sn <= 0).sum())
    n_tau = int((tau <= 0).sum())
    rep.add(gate, n_sn == 0 and n_tau == 0,
            f"sigma_n > 0 and tau_0 > 0 on every facet: min sigma_n = {np.nanmin(sn):.4f} "
            f"MPa ({n_sn} non-positive), min tau_0 = {np.nanmin(tau):.5f} MPa "
            f"({n_tau} non-positive).  A zero here is the psi_ini = -inf abort.")
    return rep


def check_no_preslip(mu_app, mu_s, rep: GateReport, gate: str) -> GateReport:
    """max(mu_app) < mu_s, or the fault slips at t = 0."""
    m = float(np.nanmax(mu_app))
    rep.add(gate, m < mu_s,
            f"no pre-slip at t=0: max mu_app = {m:.4f} vs strength {mu_s:.4f}")
    return rep


def check_receivers_under_the_surface(receivers, rep: GateReport, gate: str,
                                      z_surface: float = 0.0) -> GateReport:
    """No receiver exactly at z = 0.

    SeisSol v1.1.3 silently DROPS a receiver at the free surface; the convention is
    z = -1 m.  A dev checkout does not drop them, so this passes locally and loses
    stations on the cluster.
    """
    z = np.asarray(receivers, float)[:, 2]
    n = int(np.isclose(z, z_surface).sum())
    rep.add(gate, n == 0,
            f"{n} of {len(z)} receiver(s) sit exactly at z = {z_surface:g}; SeisSol "
            f"v1.1.3 silently drops those -- use z = -1 m")
    return rep


# --------------------------------------------------------------------------- Stage F
def verify_mesh_against_deck(cfg, mesh_path=None, *, stress_nc=None, friction_nc=None,
                             material_nc=None, receivers=None, mesh: str | None = None,
                             mu_s: float | None = None,
                             descriptor_sha256: str | None = None) -> GateReport:
    """Stage F: is this mesh compatible with the deck that will consume it?

    Every field is sampled at the real fault facets, so what is checked is exactly what
    SeisSol will read.
    """
    mspec = cfg.mesh(mesh or cfg.default_mesh)
    if mesh_path is not None:
        from dataclasses import replace
        mspec = replace(mspec, path=str(mesh_path))
    rep = GateReport(f"Stage F: mesh vs deck ({Path(mspec.path).name})")

    if descriptor_sha256 is not None and descriptor_sha256 != cfg.sha256():
        # Checking a mesh against fields built from a DIFFERENT descriptor is worse than
        # not checking it: every number would look plausible and mean nothing.
        rep.add("F0", False,
                f"the supplied fields were built from descriptor {descriptor_sha256[:12]} "
                f"but this config is {cfg.sha256()[:12]}; refusing to cross-check")
        return rep

    fault = load_fault(mspec, cfg.data_dir, strike=cfg.strike)
    rep.add("F1", len(fault) > 0,
            f"{len(fault):,} fault facets, {fault.n_degenerate} degenerate", severity=WARN)

    for label, nc, gate in (("stress", stress_nc, "F2s"), ("friction", friction_nc, "F2f"),
                            ("material", material_nc, "F2m")):
        if nc is None:
            rep.skip(gate, f"no {label} nc supplied")
            continue
        hull_containment(nc, fault.cent, label=f"fault vs {label}", gate=f"F-hull-{label}",
                         severity=HARD, report=rep)
        check_finite_fields(nc, fault, label, rep, gate)

    if stress_nc is not None:
        sn, tau, mu = project_stress_onto_fault(stress_nc, fault)
        check_positive_normal_stress(sn, tau, rep, "F3")
        if mu_s is not None:
            check_no_preslip(mu, mu_s, rep, "F4")
        else:
            rep.skip("F4", "no mu_s given; the t=0 pre-slip check did not run")
    else:
        rep.skip("F3", "no stress nc supplied")
        rep.skip("F4", "no stress nc supplied")

    # F5 -- the hypocentre must still land on this mesh.
    try:
        hypo = cfg.hypocenter(mesh or cfg.default_mesh)
    except Exception:                                           # noqa: BLE001
        rep.skip("F5", "no hypocentre defined for this mesh")
    else:
        try:
            snap = snap_hypocenter(hypo, fault, cfg.strike, cfg.crs,
                                   gate_bands=cfg.gate_bands)
            rep.extend(snap.report, prefix="F5-")
        except Exception as exc:                                # noqa: BLE001
            rep.add("F5", False, f"the hypocentre does not land on this mesh: {exc}")

    if receivers is not None:
        check_receivers_under_the_surface(receivers, rep, "F6")
    else:
        rep.skip("F6", "no receiver list supplied")
    return rep
