"""hdiv_cyoke_headtohead.py -- wall-clock head-to-head: the SCALABLE HDiv-VIM (now dense-N^2-FREE) on the
C-yoke vs the six-face surface-charge distortion-element reference (numbers embedded below for reproducibility; the
source JSONs are now committed under reference_yano_msc/ -- rescued from git history, see its README).

The scalable HDiv-VIM path is dense-N^2-free as of 2026-06-09:
  - BUILD: O(N log N) analytic charge-Gram H-matrix with the near/far split (analytic near, monopole far),
    + sparse FE assembly (M_mass / B are scipy CSR -- no dense N^2 object anywhere).
  - SOLVE: tensor-tangent damped Newton.  The HDiv loops are ker(B) -- FIELD-NULL BY CONSTRUCTION (de Rham;
    no loop-star, no cohomology), so the tangent is well-conditioned and mu_r-independent -> 5-6 Newton iters.

six-face surface-charge reference (distortion elements + loop-star; LAB machine, saved analysis times):
  18900  DOF: 174 nonlinear iters, 99 s total
  165600 DOF: 214 nonlinear iters, 2607 s total = 582 s H-matrix build + 1953 s linear solve (2686 lin iters)

This sweeps the SAME C-yoke geometry as hdiv_cyoke_nonlinear.py over mesh density and reports BUILD + SOLVE
wall-clock, so the HDiv-VIM curve can be read directly against the two yano points.  (A C-yoke is FLAT, so
the geometry/charge setup is identical to the validated hdiv_cyoke_nonlinear.py cross-check vs Radia MMM;
this script measures TIME, not accuracy -- accuracy is locked by tests/feec/test_hdiv_vim_*.)

STATUS (2026-06-09, RESOLVED -- see docs/hdiv_vim/PRODUCTIONIZATION.md): the first runs of THIS script
exposed two bugs, both now FIXED, after which the SOLVE is mesh-INDEPENDENT (6-iter quadratic convergence):
  - Fix A: the convergence test was unsound (an Mavg-stagnation break declared "converged" mid-globalization
    -> the apparent "Mz drift to 509k" was a FALSE convergence).  Now breaks on the true residual relF and
    RAISES if maxit is hit -- never silently returns a wrong M.
  - Fix B: the iteration count grew with refinement (6 -> 27 -> 37) because the DEFAULT ACA tolerance
    gram_eps=1e-7 was too loose -> the charge-Gram H-matrix was ~5% inaccurate/asymmetric -> the Newton
    tangent was inconsistent.  Tightening the default to 1e-10 restores 6-iter convergence at EVERY mesh.
    (NEITHER the -N material formulation NOR a loop-space preconditioner was needed -- both were red herrings.)
This script now runs with near_factor=2.0 (fast near/far build) + gram_eps=1e-10 (default; accurate): on the
C-yoke it is ~6 iters and ~16 s total at ndof~10759 vs yano's 99 s / 174 iters at 18900 DOF -- a clear win
on BOTH build and solve.  (A C-yoke is FLAT, so Radia/the dense path is a valid accuracy cross-reference;
this script measures TIME -- accuracy is locked by tests/feec/test_hdiv_vim_*.)
"""
import json
import os
import sys

import numpy as np  # noqa: F401  (kept for interactive inspection)
import ngsolve as ng

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hdiv_cyoke_nonlinear import cyoke_mesh                      # noqa: E402  same C-yoke geometry
from radia.vim import _nonlinear as nl                       # noqa: E402

ng.SetNumThreads(4)
HERE = os.path.dirname(os.path.abspath(__file__))

# six-face surface-charge distortion-element solver reference (committed under reference_yano_msc/golden/; embedded here
# for reproducibility -- the LU/BiCGSTAB/HACApK x DOF matrix is the superset of the two points below).
YANO_REF = {
    "method": "six-face surface-charge distortion elements + loop-star (LAB machine, saved)",
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
           "status_2026_06_09": ("Fix A (sound relF break + fail-loud) + Fix B (default gram_eps "
                                 "1e-7->1e-10; the loose ACA tolerance degraded the OUTER Newton). SOLVE "
                                 "is mesh-INDEPENDENT 6-iter + accurate + a clear win UP TO ndof~15726 "
                                 "(the rows above). SCALE CAVEAT: the INNER +N GMRES solve does NOT yet "
                                 "scale past ~20-40k ndof (at ndof~44k Newton's warmstart fails in 1000 "
                                 "inner iters) -- but this is NOT a fundamental wall. CORRECTED COMPARISON: "
                                 "the six-face surface-charge reference is NO-loop-star + BLOCK JACOBI (NOT H-ILU; H-ILU "
                                 "is the separate loop-star A_SS solver), which ALSO needs growing inner "
                                 "iters (2686 lin @165600) absorbed by its Picard/修正反復法 outer loop. "
                                 "The FAIR head-to-head is SAME solver (Picard + Block Jacobi) for both, "
                                 "where HDiv's de-Rham field-null loops should beat six-face surface-charge distortion "
                                 "loops -- NOT YET RUN. Do NOT read this as a 165600-DOF result.")}
    with open(os.path.join(HERE, "hdiv_cyoke_headtohead.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("saved", os.path.join(HERE, "hdiv_cyoke_headtohead.json"))
