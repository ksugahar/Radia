"""hdiv_cyoke_headtohead.py -- wall-clock head-to-head: the SCALABLE HDiv-VIM (now dense-N^2-FREE) on the
C-yoke vs the yano-type distortion-element reference (saved on this machine; numbers embedded below for
reproducibility, since the yano JSON is not committed to this repo).

The scalable HDiv-VIM path is dense-N^2-free as of 2026-06-09:
  - BUILD: O(N log N) analytic charge-Gram H-matrix with the near/far split (analytic near, monopole far),
    + sparse FE assembly (M_mass / B are scipy CSR -- no dense N^2 object anywhere).
  - SOLVE: tensor-tangent damped Newton.  The HDiv loops are ker(B) -- FIELD-NULL BY CONSTRUCTION (de Rham;
    no loop-star, no cohomology), so the tangent is well-conditioned and mu_r-independent -> 5-6 Newton iters.

yano-type reference (distortion elements + loop-star; LAB machine, saved analysis times):
  18900  DOF: 174 nonlinear iters, 99 s total
  165600 DOF: 214 nonlinear iters, 2607 s total = 582 s H-matrix build + 1953 s linear solve (2686 lin iters)

This sweeps the SAME C-yoke geometry as hdiv_cyoke_nonlinear.py over mesh density and reports BUILD + SOLVE
wall-clock, so the HDiv-VIM curve can be read directly against the two yano points.  (A C-yoke is FLAT, so
the geometry/charge setup is identical to the validated hdiv_cyoke_nonlinear.py cross-check vs Radia MMM;
this script measures TIME, not accuracy -- accuracy is locked by tests/feec/test_hdiv_vim_*.)

HONEST STATUS (2026-06-09, see docs/hdiv_vim/PRODUCTIONIZATION.md): the first runs of THIS script exposed
that the scalable nonlinear SOLVE is NOT yet mesh-robust.  BUILD is scalable + a clear win, and the SOLVE
converges to the CORRECT Mz (the method + tangent are correct -- quadratic convergence verified), BUT the
Newton iteration COUNT grows with refinement (6 -> 27 -> 37 over h 0.008 -> 0.005) because the +N MINRES
warmstart is ill-conditioned (the field-null loops ker B are near-null modes, eig ~1/chi).  The "5-6 iters"
holds only at COARSE mesh.  Fix A (committed) made the convergence test sound + fail-loud (the earlier
"Mz drift to 509k" was a FALSE-CONVERGENCE bug, now impossible -- it raises instead of returning a wrong M).
Fix B (the -N mu_r-independent material formulation / loop-space preconditioner for a mesh-independent
iteration count) is NOT done.  So this script's SOLVE timings are real but reflect the pre-Fix-B (slow,
warmstart-limited) solve; do NOT read them as the final head-to-head until Fix B lands.
"""
import json
import os
import sys

import numpy as np  # noqa: F401  (kept for interactive inspection)
import ngsolve as ng

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hdiv_cyoke_nonlinear import cyoke_mesh                      # noqa: E402  same C-yoke geometry
from radia.hdiv_vim import _nonlinear as nl                       # noqa: E402

ng.SetNumThreads(4)
HERE = os.path.dirname(os.path.abspath(__file__))

# yano-type distortion-element solver reference (saved on the LAB machine; embedded for reproducibility).
YANO_REF = {
    "method": "yano-type distortion elements + loop-star (LAB machine, saved)",
    "points": [
        {"ndof": 18900, "nonlinear_iters": 174, "t_total_s": 99.0},
        {"ndof": 165600, "nonlinear_iters": 214, "t_total_s": 2607.0,
         "t_build_s": 582.0, "t_solve_s": 1953.0, "linear_iters": 2686},
    ],
}


def run(hs, chi0=1000.0, Msat=1.0e6, H0=2.0e5, near_factor=2.0):
    rows = []
    for h in hs:
        with ng.TaskManager():                                   # caller-wraps (CLAUDE.md TaskManager policy)
            mesh = cyoke_mesh(h)
            ne = int(mesh.ne)
            Mz, nit, _D, tm = nl.solve_nonlinear_newton_scalable(
                mesh, chi0, Msat, H0, near_factor=near_factor, return_timing=True)
        rows.append(dict(h=h, ne=ne, n_charge=tm["n_charge"], ndof=tm["ndof"], Mz=Mz, iters=nit,
                         t_build_s=tm["t_build_s"], t_solve_s=tm["t_solve_s"], t_total_s=tm["t_total_s"]))
        print(f"  {h:>6.4f} {ne:>7d} {tm['n_charge']:>8d} {tm['ndof']:>7d} {Mz:>10.1f} {nit:>5d} "
              f"{tm['t_build_s']:>8.2f} {tm['t_solve_s']:>8.2f} {tm['t_total_s']:>8.2f}", flush=True)
    return rows


if __name__ == "__main__":
    hs = (0.020, 0.012, 0.008, 0.006, 0.005)
    print("Scalable HDiv-VIM C-yoke (dense-N^2-free, near/far split) -- BUILD + SOLVE wall-clock:")
    print(f"  {'maxh':>6} {'ne':>7} {'n_charge':>8} {'ndof':>7} {'Mz':>10} {'iters':>5} "
          f"{'build_s':>8} {'solve_s':>8} {'total_s':>8}")
    rows = run(hs)
    out = {"geometry": "C-yoke (Kelvin-less, iron-only volume integral -- no air, no Kelvin)",
           "chi0": 1000.0, "Msat": 1.0e6, "H0": 2.0e5, "near_factor": 2.0,
           "hdiv_vim": rows, "yano_reference": YANO_REF,
           "status_2026_06_09": ("BUILD scalable + clear win; SOLVE converges to the CORRECT Mz (Fix A: "
                                 "sound relF break + fail-loud) but iter COUNT grows with refinement "
                                 "(ill-conditioned +N warmstart) -- Fix B (-N material formulation / loop "
                                 "preconditioner) pending. NOT the final head-to-head until Fix B.")}
    with open(os.path.join(HERE, "hdiv_cyoke_headtohead.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("saved", os.path.join(HERE, "hdiv_cyoke_headtohead.json"))
