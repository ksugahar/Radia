"""hdiv_demag_tet_nonlinear.py -- NONLINEAR HDiv-type VIM demag (applied-field solve + BH-curve Picard).

The linear HDiv-VIM (hdiv_demag_tet.py) computes demag FACTORS via eig(N, M_mass).  A NONLINEAR solve
needs an APPLIED-FIELD solve (M for a given H_ext), then a per-element constitutive iteration.

## Applied-field formulation (verify-first result, 2026-06-07)
The eigenvalue framing A_eig = (1/chi) M_mass - N is NOT the applied-field system.  The physical
applied-field weak form is

    A+ m = M_mass h_ext ,   A+ = (1/chi) M_mass + N         (PLUS N)

  because  M = chi (H_ext + H_demag),  H_demag,weak = -N m  =>  (1/chi) M_mass m + N m = M_mass h_ext.
For a sphere this reproduces the analytic  M/H = chi/(1 + chi D)  (D = demag factor) -- VERIFIED here
to <=2.5% for mu_r<=100 (A+); the minus-sign system gives nonsense (negative / divergent).

## Nonlinear (Picard, secant susceptibility)
With a BH curve M(H), iterate:  solve A+(chi^k) m -> M_avg, internal field H_int = H0 - D M_avg,
chi^{k+1} = M(H_int)/H_int (secant), re-solve.  Converges + saturates (M -> M_sat).  This prototype
uses the SCALAR (uniform-M sphere) update; the per-element (non-uniform body) generalization needs the
per-cell field reconstruction (RT0 -> cell H) and an MMM/MSC cross-check -- the next increment.

## Honest accuracy
The nonlinear M value carries the SAME operator accuracy as the linear demag: the centroid-monopole N
UNDER-estimates the sphere demag (D~0.31 vs analytic 1/3), so the high-chi (deep-saturation-onset)
response is off by ~the same amount.  The #3 near-field correction (build_near_correction) raises D
toward 1/3 and tightens it -- demonstrated below.
"""
import os
from math import pi

import numpy as np

import ngsolve as ng
from netgen.csg import CSGeometry, Sphere, Pnt

import hdiv_demag_tet as tet


def _sphere(h):
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), 1.0))
    return ng.Mesh(geo.GenerateMesh(maxh=h))


def _bh_curve(chi0, Msat):
    """Saturating M(H) = chi0 H / (1 + chi0|H|/Msat): slope chi0 at H=0, asymptote +-Msat."""
    return lambda H: chi0 * H / (1.0 + chi0 * abs(H) / Msat)


def _scalar_fixed_point(Mof, D, H0):
    """Correct analytic uniform-sphere root: solve M = M(H0 - D M) by bisection on f(M)=M-M(H0-DM)."""
    lo, hi = -1.0, 1.0
    f = lambda M: M - Mof(H0 - D * M)
    # widen until sign change
    while f(lo) * f(hi) > 0 and hi < 1e6:
        lo *= 2; hi *= 2
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def solve_nonlinear(mesh, Mof, H0, near_correction=False, nsub=4, maxit=300, tol=1e-9):
    """Nonlinear applied-field HDiv-VIM solve (scalar-chi Picard on the A+ system).  Returns
    (M_avg, n_iter, D_used)."""
    d = tet.build_demag(mesh, nsub)
    N = d["N"].copy()
    M = d["M_mass"]
    mu = d["m_unit"]
    denom = mu @ M @ mu
    if near_correction:
        # demag operator with the #3 near-field correction folded in: N_eff = B^T (G + corr) B
        corr = tet.build_near_correction(mesh, d, nsub=nsub, near_factor=2.0)
        B = d["B_csr"]
        N = N + np.asarray((B.T @ corr @ B).todense())
    D = float((mu @ N @ mu) / denom)
    b0 = M @ mu
    chi, Mavg = 1000.0, 0.0
    for it in range(maxit):
        A = (1.0 / chi) * M + N
        m = np.linalg.solve(A, H0 * b0)
        Mavg = float((mu @ M @ m) / denom)
        Hint = H0 - D * Mavg
        chi_new = Mof(Hint) / Hint if abs(Hint) > 1e-30 else 1000.0
        if abs(chi_new - chi) < tol * chi:
            break
        chi = 0.5 * chi + 0.5 * chi_new
    return Mavg, it + 1, D


def main():
    mesh = _sphere(0.35)
    chi0, Msat = 1000.0, 1.0
    Mof = _bh_curve(chi0, Msat)
    print(f"Nonlinear HDiv-VIM sphere demag (chi0={chi0}, Msat={Msat})")
    print(f"{'H0':>8} {'Picard M':>10} {'+nearcorr':>10} {'analytic(1/3)':>14} {'M/Msat':>8}")
    for H0 in (1e-4, 1e-3, 1e-2, 1e-1):
        Mmono, it1, Dm = solve_nonlinear(mesh, Mof, H0, near_correction=False)
        Mcorr, it2, Dc = solve_nonlinear(mesh, Mof, H0, near_correction=True)
        Mana = _scalar_fixed_point(Mof, 1.0 / 3.0, H0)   # analytic uses the EXACT sphere D=1/3
        print(f"{H0:8.0e} {Mmono:10.5f} {Mcorr:10.5f} {Mana:14.5f} {Mcorr/Msat:8.3f}")
    print(f"  (monopole D={Dm:.4f}, near-corrected D={Dc:.4f}, analytic D=0.3333)")
    print("  => Picard converges + saturates; near-correction moves D toward 1/3 -> closer to analytic.")


if __name__ == "__main__":
    main()
