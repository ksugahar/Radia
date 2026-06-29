"""interior_PEC_eig_diagnostic.py — Direct eigenvalue computation for
the interior PEC eddy current problem.

DIAGNOSTIC for 2026-05-09: Tanimoto+H-H gives 225.77 us for cylinder
R=10 t=2, while BEM/v5/axifem 4-way consensus gives 218.69 us. The
hypothesis is that Tanimoto+H-H is solving the **interior PEC**
problem (n×A=0 on conductor surface) rather than vacuum-coupled.

Direct test: solve generalized eigenvalue problem
  K v = λ M v
on conductor-only mesh with HCurl(dirichlet="conductorBND"), where
  K = ∫(1/μ) curl·curl,   M = ∫σ uv
The leading eigenvalue gives τ_lead = 1/λ_min.

If τ_lead ≈ 225 us → confirms Tanimoto+H-H is interior PEC.
If τ_lead ≈ 219 us → Tanimoto+H-H matches vacuum-coupled (mystery deepens).
If τ_lead is different → something else is going on.
"""
from __future__ import annotations

import argparse
import sys
import time
from math import pi

sys.stdout.reconfigure(encoding="utf-8", errors="replace", line_buffering=True)

from netgen.occ import Cylinder, Box, Pnt, Z, OCCGeometry
from ngsolve import (
    Mesh, HCurl, BilinearForm, GridFunction, MultiVector,
    CoefficientFunction, curl, dx, TaskManager, ngsglobals,
    InnerProduct, IfPos, Conj,
)
from ngsolve.solvers import LOBPCG, PINVIT
from ngsolve.eigenvalues import LOBPCG as eig_LOBPCG, PINVIT as eig_PINVIT


mu0 = 4 * pi * 1e-7
sigma_Cu = 5.8e7
B0 = 1.0

GEOM_TABLE = {
    "cylinder":   {"kind": "cyl",  "R": 10.0,  "t":  2.0},
    "square":     {"kind": "rect", "Lx": 17.72, "Ly": 17.72, "Lz": 2.0},
    "cuboid_521": {"kind": "rect", "Lx":  5.0,  "Ly":  2.0,  "Lz": 1.0},
}


def build_geo(geom_name, h_cond):
    p = GEOM_TABLE[geom_name]
    if p["kind"] == "cyl":
        R = p["R"] * 1e-3
        t = p["t"] * 1e-3
        cond = Cylinder(Pnt(0, 0, -t/2), Z, r=R, h=t)
    else:
        Lx = p["Lx"] * 1e-3
        Ly = p["Ly"] * 1e-3
        Lz = p["Lz"] * 1e-3
        cond = Box(Pnt(-Lx/2, -Ly/2, -Lz/2),
                   Pnt( Lx/2,  Ly/2,  Lz/2))
    cond.faces.name = "conductorBND"
    cond.mat("conductor")
    cond.maxh = h_cond
    return OCCGeometry(cond)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geom", choices=list(GEOM_TABLE.keys()),
                         default="cylinder")
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--h-cond-mm", type=float, default=1.0)
    parser.add_argument("--n-eig", type=int, default=8,
                         help="Number of leading eigenvalues to compute")
    args = parser.parse_args()

    ORDER = args.order
    H_COND = args.h_cond_mm * 1e-3
    geom_name = args.geom

    ngsglobals.msg_level = 0
    print("=" * 78, flush=True)
    print(f" Interior PEC eigenvalue diagnostic (geom={geom_name})", flush=True)
    print("=" * 78, flush=True)
    print(f"  ORDER={ORDER}, H_COND={H_COND*1000}mm, "
          f"n_eig={args.n_eig}", flush=True)
    print()

    geo = build_geo(geom_name, H_COND)
    t0 = time.time()
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=H_COND, grading=0.5)
    mesh = Mesh(ngmesh)
    mesh.Curve(ORDER)
    print(f"  Mesh: ne={mesh.ne}, nv={mesh.nv} ({time.time()-t0:.1f} s)",
          flush=True)

    # Same fesA setup as Tanimoto+H-H: HCurl with Dirichlet on conductorBND
    fes = HCurl(mesh, order=ORDER, nograds=True,
                 dirichlet="conductorBND")
    print(f"  HCurl ndof = {fes.ndof}, active = {sum(fes.FreeDofs())}",
          flush=True)

    u, v = fes.TnT()
    BONUS_INT = 8

    # K = ∫ (1/μ) curl u · curl v dx
    K = BilinearForm(fes, check_unused=False)
    K += (1.0/mu0) * curl(u) * curl(v) * dx(bonus_intorder=BONUS_INT)

    # M = ∫ σ u · v dx
    M = BilinearForm(fes, check_unused=False)
    M += sigma_Cu * u * v * dx(bonus_intorder=BONUS_INT)

    print(" Assembling K, M ...", flush=True)
    t0 = time.time()
    with TaskManager():
        K.Assemble()
        M.Assemble()
    print(f"  assembly: {time.time()-t0:.1f} s", flush=True)

    # Solve generalized eigenvalue problem K v = λ M v
    # We want SMALLEST eigenvalues (largest τ = 1/λ)
    print(f" Solving for {args.n_eig} smallest eigenvalues (PINVIT) ...",
          flush=True)
    t0 = time.time()
    with TaskManager():
        # PINVIT: invert K, find smallest λ s.t. K v = λ M v
        try:
            inv_K = K.mat.Inverse(fes.FreeDofs(), inverse="pardiso")
        except Exception:
            inv_K = K.mat.Inverse(fes.FreeDofs(),
                                   inverse="sparsecholesky")

        # PINVIT iteration
        n_eig = args.n_eig
        eigs, evecs = eig_PINVIT(K.mat, M.mat, pre=inv_K,
                                   num=n_eig, maxit=200,
                                   printrates=False)
    print(f"  PINVIT: {time.time()-t0:.1f} s", flush=True)

    # τ_n = 1/λ_n
    print()
    print(" Eigenvalues and time constants:", flush=True)
    print(f"  {'n':>3} {'λ_n [1/s]':>15} {'τ_n [μs]':>15}", flush=True)
    tau_us_list = []
    for i, lam in enumerate(eigs):
        if lam > 0:
            tau_us = 1.0 / float(lam.real if hasattr(lam, "real") else lam)
            tau_us_list.append(tau_us * 1e6)
            print(f"  {i:>3} {float(lam):>15.6e} {tau_us*1e6:>15.4f}",
                  flush=True)
        else:
            print(f"  {i:>3} {float(lam):>15.6e} {'(neg)':>15}", flush=True)

    print()
    print(" Comparison with reference values:", flush=True)
    print(f"  Tanimoto+H-H stage 2 = 225.77 μs", flush=True)
    print(f"  BEM-Foster / v5 / axifem 4-way = 218.69 ± 0.5 μs", flush=True)
    if tau_us_list:
        print(f"  PEC interior eig leading τ = {tau_us_list[0]:.4f} μs",
              flush=True)


if __name__ == "__main__":
    main()
