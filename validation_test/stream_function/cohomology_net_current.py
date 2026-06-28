"""Cohomology net-current for the stream-function coil designer (torus, b1=2).

A single-valued current potential psi gives K = n x grad_s psi, whose net
current around EVERY cycle of the winding surface is ZERO -- it can only shape
fields, never carry net current.  A real coil set DOES carry net current.  The
current potential gains the SECULAR (multivalued) term

    Psi = psi + sum_k c_k * (theta_k / 2pi)        (theta_k = the k-th angle)
    K   = n x grad_s psi + sum_k c_k h_k ,  h_k = n x grad(theta_k)

where the h_k are the harmonic 1-forms = first cohomology generators of the
winding surface (dim = first Betti number b1).  For the genus-1 torus b1 = 2:
the net-POLOIDAL current (TF solenoid -> toroidal B ~ 1/R) and the net-TOROIDAL
current (a plasma-current-like ring -> vertical B).

This 'finishes' the stream-function coil design by INTEGRATING the secular term
into the design solve (REGCOIL-style: the net current is a PRESCRIBED
engineering input -- its field is ~tangent to the target boundary, so it is not
fitted -- and psi is fit to the RESIDUAL target after the secular footprint is
removed).  It also VALIDATES the generators numerically:

  * div_s h_k ~ 0                         (div-free: a valid surface current)
  * h_k is NOT a single-valued gradient   (non-exact => genuine cohomology, not
                                           reproducible by any psi)
  * the TF generator's field obeys Ampere  B_tor * R = const inside, ~0 outside
  * the surface Euler characteristic (gmsh-free) confirms b1 = 2 (the COUNT)

Builds on demo_regcoil_fusion.py (analytic torus generators) + the surface-FE
stream function in calc_streamfunction.py.  Reuses radia.stream_function only
through the design solve.  Caller wraps NGSolve work in TaskManager.

Run:  python validation_test/stream_function/cohomology_net_current.py
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))
_EXAMPLES = os.path.join(_REPO, "examples", "stream_function")
sys.path.insert(0, _HERE)
sys.path.insert(0, _EXAMPLES)
sys.path.insert(0, os.path.join(_REPO, "src", "radia"))
sys.path.insert(0, os.path.join(_REPO, "src", "radia", "panels"))

MU0 = 4.0e-7 * math.pi


def _div_and_nonexactness(coil, n_cf, K):
    """Return (div_s K / |K|, non-exactness of n x K).

    div_s K via the weak surface divergence (H1 dual norm).  Non-exactness =
    || n x K - grad_s w* || / || n x K || where w* is the best single-valued H1
    potential: ~0 for an exact (psi-type) current, ~1 for a genuine cohomology
    generator (n x K = -grad(multivalued angle) is NOT a single-valued grad)."""
    from ngsolve import (H1, BilinearForm, LinearForm, GridFunction, grad, ds,
                         InnerProduct, Cross, Integrate)
    Kn = math.sqrt(Integrate(InnerProduct(K, K) * ds, coil))
    fes = H1(coil, order=2, definedon=coil.Boundaries(".*"))
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += grad(u).Trace() * grad(v).Trace() * ds + 1e-10 * u * v * ds
    a.Assemble()
    inv = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")

    lf = LinearForm(fes)
    lf += InnerProduct(K, grad(v).Trace()) * ds
    lf.Assemble()
    w = GridFunction(fes)
    w.vec.data = inv * lf.vec
    div_rel = math.sqrt(abs(InnerProduct(lf.vec, w.vec))) / (Kn + 1e-30)

    nxK = Cross(n_cf, K)
    lf2 = LinearForm(fes)
    lf2 += InnerProduct(nxK, grad(v).Trace()) * ds
    lf2.Assemble()
    w2 = GridFunction(fes)
    w2.vec.data = inv * lf2.vec
    resid = math.sqrt(abs(Integrate(
        InnerProduct(nxK - grad(w2).Trace(), nxK - grad(w2).Trace()) * ds, coil)))
    nxKn = math.sqrt(Integrate(InnerProduct(nxK, nxK) * ds, coil))
    nonexact = resid / (nxKn + 1e-30)
    return div_rel, nonexact


def main(no_cohomology=False):
    import demo_regcoil_fusion as D
    import calc_streamfunction as C
    from ngsolve import Mesh, H1, specialcf, TaskManager

    print("=== cohomology net-current for SF coil design (torus, b1=2) ===")
    work = os.path.join(os.environ.get("TEMP", "/tmp"), "sf_coh_netcurrent")
    os.makedirs(work, exist_ok=True)

    with TaskManager():
        coil = Mesh(D._torus_surface_vol(D.A_WIND, 0.045,
                                         os.path.join(work, "wind.vol")))
        plasma = Mesh(D._torus_surface_vol(D.A_PLASMA, 0.035,
                                           os.path.join(work, "plasma.vol")))
        fes = H1(coil, order=1, definedon=coil.Boundaries(".*"))
        n_cf = specialcf.normal(3)
        pts, nrm, theta, phi = D._plasma_points_normals(plasma, 220)
        A_n = D._A_normal(C, fes, n_cf, pts, nrm)

        # ---- (a) numeric generator validation ----
        K_list = D._secular_currents(n_cf)        # analytic h_k = n x grad(angle)
        gen = {}
        for nm, K in K_list:
            dv, nx = _div_and_nonexactness(coil, n_cf, K)
            gen[nm] = {"div_rel": dv, "nonexactness": nx}
            print(f"[gen] {nm:16s}  div_s/|K|={dv:.2e}  "
                  f"non-exactness={nx:.3f}  (1=pure cohomology)")
        b1 = None if no_cohomology else D._betti1_winding_surface(D.A_WIND, 0.06)
        print(f"[gen] Gmsh b1(winding surface) = {b1}  (expect 2)")

        # ---- TF generator field obeys Ampere B_tor*R = const ----
        radii = np.linspace(0.20, 0.40, 9)
        by_in = D._tf_field_1overR(coil, n_cf, radii)
        br = by_in * radii
        br_spread = float((br.max() - br.min()) / abs(br.mean()))
        by_out = D._tf_field_1overR(coil, n_cf, np.array([0.10, 0.55]))
        out_in = float(np.max(np.abs(by_out)) / np.max(np.abs(by_in)))
        print(f"[TF ] B_tor*R const to {br_spread*100:.2f}% inside, "
              f"outside/inside={out_in:.1e}  (Ampere)")

        # ---- (b) JOINT design: prescribe net current, fit psi to residual ----
        sec_n = D._assemble_secular_normal(coil, pts, nrm, K_list)  # (M, 2) B.n
        Bphi0 = 1.0                                # desired on-axis B_tor [T]
        Ipol = 2 * math.pi * D.R_MAJOR * Bphi0 / MU0   # net poloidal current
        c_vec = np.array([1.0, 0.0])               # prescribe TF (unit gen), no ring
        Bn_secular = sec_n @ c_vec
        Bn_target = nrm[:, 2]                       # vertical-field PF shaping target
        Bn_resid = Bn_target - Bn_secular           # psi fits what's left
        psi, res_psi, peakJ = D._design(C, fes, coil, A_n, Bn_resid, "h1", 1e-7)
        Bn_total = A_n @ psi + Bn_secular
        res_total = float(np.linalg.norm(Bn_total - Bn_target)
                          / np.linalg.norm(Bn_target))
        print(f"[joint] psi fit to (target - secular) = {res_psi:.2e}")
        print(f"[joint] total (psi + secular) vs target = {res_total:.2e}")
        print(f"[joint] prescribe B_tor={Bphi0:.0f}T -> net poloidal current "
              f"{Ipol/1e3:.0f} kA")

    result = {
        "b1_euler": b1,
        "n_secular_dofs": len(K_list),
        "gen_div_rel": {k: gen[k]["div_rel"] for k in gen},
        "gen_nonexactness": {k: gen[k]["nonexactness"] for k in gen},
        "tf_BR_spread_rel": br_spread,
        "tf_outside_over_inside": out_in,
        "joint_psi_residual": res_psi,
        "joint_total_residual": res_total,
        "net_poloidal_current_A": Ipol,
        "winding_ndof": int(fes.ndof),
    }
    ok = (gen["net_poloidal_TF"]["div_rel"] < 5e-2
          and gen["net_poloidal_TF"]["nonexactness"] > 0.8
          and gen["net_toroidal"]["nonexactness"] > 0.8
          and (no_cohomology or b1 == 2)
          and br_spread < 0.05 and out_in < 0.05
          and res_total < 1e-6)
    print(f"\n=== PASS: {ok} ===")
    return result, ok


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-cohomology", action="store_true",
                    help="skip the Gmsh b1 confirmation")
    a = ap.parse_args()
    _, ok = main(no_cohomology=a.no_cohomology)
    raise SystemExit(0 if ok else 1)
