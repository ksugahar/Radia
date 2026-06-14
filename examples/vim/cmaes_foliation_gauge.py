"""F1 remedy: CMA-ES over the foliation gauge (global, escapes the wall).

foliation_choice_wall.py established F1 honestly: choosing the foliation mu for a
multi-axis target is a BILEVEL problem -- the inner field-fit is CONVEX (exact at
every gauge), but the outer coil-complexity objective g = ||J||_L2 over foliations
is NON-CONVEX (multi-basin).  This is the documented remedy: a GLOBAL optimizer
(CMA-ES, via Optuna -- the lab does NOT re-implement CMA-ES) over the gauge.

Gauge = a 2-parameter foliation normal n(theta,phi) = (sin th cos ph,
sin th sin ph, cos th), mu = n.(x,y,z).  Objective:

    g(theta,phi) = || J(lambda*; mu) ||_L2 ,   lambda* = exact ACA+TSVD min-norm fit

A LOCAL optimizer lands in whichever basin it starts in (the F1 wall); CMA-ES
samples globally and finds the deepest basin.  We demonstrate both: map the
landscape (multi-basin), run CMA-ES (-> ~global), and run a LOCAL optimizer from
a shallow-basin start (-> stuck, higher complexity).

Reuses foliation_choice_wall.py (build_mesh, loop_field) + radia.stream_function
(aca_tsvd) UNCHANGED.  Caller wraps NGSolve work in TaskManager.

Run:  python examples/feec_vim/cmaes_foliation_gauge.py
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, "..", "..", "src", "radia"))

from ngsolve import (H1, GridFunction, CF, grad, InnerProduct, dx, x, y, z,
                     LinearForm, TaskManager, Cross, Integrate)
from stream_function import aca_tsvd, pseudo_inverse_solve

MU0 = 4e-7 * math.pi


def _n_dir(theta, phi):
    return CF((math.sin(theta) * math.cos(phi),
               math.sin(theta) * math.sin(phi),
               math.cos(theta)))


def _response(fes, targets, theta, phi):
    _, v = fes.TnT()
    gmu = _n_dir(theta, phi)
    M, N = 3 * len(targets), fes.ndof
    A = np.zeros((M, N))
    for i, (px, py, pz) in enumerate(targets):
        d = CF((px, py, pz)) - CF((x, y, z))
        K = d / (InnerProduct(d, d) ** 1.5)
        cross = Cross(Cross(grad(v), gmu), K)
        for c in range(3):
            lf = LinearForm(fes)
            lf += (MU0 / (4 * math.pi)) * cross[c] * dx
            lf.Assemble()
            A[3 * i + c, :] = lf.vec.FV().NumPy()
    return A


def _gnorm(fes, mesh, lam, theta, phi):
    gf = GridFunction(fes); gf.vec.FV().NumPy()[:] = lam
    J = Cross(grad(gf), _n_dir(theta, phi))
    return math.sqrt(Integrate(InnerProduct(J, J) * dx, mesh))


def main(n_trials=40):
    import foliation_choice_wall as F

    print("=== F1 remedy: CMA-ES over the foliation gauge ===")
    with TaskManager():
        mesh = F.build_mesh(maxh=0.06)
        fes = H1(mesh, order=2)
        gp = np.linspace(-0.5 * F.HALF, 0.5 * F.HALF, 2)       # M = 3*8 = 24
        targets = [(a, b, c) for a in gp for b in gp for c in gp]
        # two-axis target (z-loop + x-loop): no single foliation aligns with both
        B_two = F.loop_field("z", targets) + F.loop_field("x", targets)

        def cost(theta, phi):
            A = _response(fes, targets, theta, phi)
            M, N = A.shape
            res = aca_tsvd(M, N, lambda i, j: A[i, j], modes=min(M, N),
                           kmax=min(M, N), aca_eps=1e-12, method=3)
            lam = pseudo_inverse_solve(res, B_two)
            fit = float(np.linalg.norm(A @ lam - B_two)
                        / (np.linalg.norm(B_two) + 1e-30))
            return _gnorm(fes, mesh, lam, theta, phi), fit

        # ---- map the landscape (reference global + basin count) ----
        nth, nph = 7, 12
        ths = np.linspace(0.05, math.pi - 0.05, nth)
        phs = np.linspace(0, 2 * math.pi, nph, endpoint=False)
        G = np.zeros((nth, nph)); fit_max = 0.0
        for i, th in enumerate(ths):
            for j, ph in enumerate(phs):
                G[i, j], ft = cost(th, ph)
                fit_max = max(fit_max, ft)
        g_grid_min = float(G.min())
        n_minima = 0
        shallow = None
        for i in range(1, nth - 1):
            for j in range(nph):
                vv = G[i, j]
                nb = [G[i-1, j], G[i+1, j], G[i, (j-1) % nph], G[i, (j+1) % nph]]
                if all(vv < b for b in nb):
                    n_minima += 1
                    if vv > 1.05 * g_grid_min:          # a NON-global basin
                        if shallow is None or vv < shallow[0]:
                            shallow = (vv, ths[i], phs[j])
        print(f"landscape: g in [{G.min()/1e3:.0f},{G.max()/1e3:.0f}]e3, "
              f"local minima={n_minima}, fit_max={fit_max:.1e}")

        # ---- CMA-ES (global) ----
        import optuna
        from optuna.samplers import CmaEsSampler
        optuna.logging.set_verbosity(optuna.logging.WARNING)

        def objective(trial):
            th = trial.suggest_float("theta", 0.05, math.pi - 0.05)
            ph = trial.suggest_float("phi", 0.0, 2 * math.pi)
            return cost(th, ph)[0]
        study = optuna.create_study(direction="minimize",
                                    sampler=CmaEsSampler(seed=1))
        study.optimize(objective, n_trials=int(n_trials))
        g_cma = float(study.best_value)

        # ---- LOCAL optimizer from a shallow-basin start (the F1 wall) ----
        from scipy.optimize import minimize
        if shallow is not None:
            x0 = [shallow[1], shallow[2]]
        else:                                            # fallback: a fixed start
            x0 = [0.4, 0.4]
        loc = minimize(lambda p: cost(p[0], p[1])[0], x0, method="Nelder-Mead",
                       options={"xatol": 1e-3, "fatol": 1e-3, "maxiter": 60})
        g_local = float(loc.fun)

    print(f"CMA-ES (global) best g = {g_cma/1e3:.1f}e3  "
          f"(vs grid-global {g_grid_min/1e3:.1f}e3 -> ratio {g_cma/g_grid_min:.3f})")
    print(f"LOCAL  from shallow basin g = {g_local/1e3:.1f}e3  "
          f"(stuck: {g_local/g_cma:.2f}x the CMA-ES optimum)")
    print(f"-> CMA-ES escapes the F1 wall; the local optimizer does not.")

    result = {
        "ne": int(mesh.ne), "ndof": int(fes.ndof), "M": int(3 * len(targets)),
        "fit_max": fit_max, "n_local_minima": int(n_minima),
        "g_grid_min": g_grid_min, "g_cmaes": g_cma, "g_local": g_local,
        "cmaes_vs_grid": g_cma / g_grid_min, "local_vs_cmaes": g_local / g_cma,
        "n_trials": int(n_trials),
    }
    ok = (fit_max < 1e-6                      # inner field-fit exact (convex)
          and n_minima >= 2                   # outer is multi-basin (the F1 wall)
          and g_cma <= 1.05 * g_grid_min      # CMA-ES ~ global
          and g_local > 1.03 * g_cma)         # local optimizer stuck higher
    print(f"\n=== PASS: {ok} ===")
    return result, ok


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-trials", type=int, default=40)
    a = ap.parse_args()
    _, ok = main(n_trials=a.n_trials)
    raise SystemExit(0 if ok else 1)
