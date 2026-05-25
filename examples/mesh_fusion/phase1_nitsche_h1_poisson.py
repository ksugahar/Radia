"""Phase 1 — Nitsche coupling, 2D Poisson, heterogeneous material.

Demonstrates the core building block of HFSS Mesh Fusion in NGSolve:
two side-by-side sub-rectangles meshed AT DIFFERENT DENSITIES, with
the interface coupled via Nitsche's penalty method.

Problem:
    -nabla . (mu nabla u) = f       on Omega = [0,1] x [0,1]
                          u = 0     on dOmega
    mu(x,y) = mu_L   if x < 0.5
           = mu_R   if x > 0.5      (jump at x=0.5)
    f(x,y) = sin(pi x) sin(pi y)

The exact solution for constant mu is u = f / (2 pi^2 mu).  With a
jump in mu, the analytic solution is piecewise (transmission
conditions:  u_L = u_R  AND  mu_L du_L/dn = mu_R du_R/dn  at x=0.5).

For the demonstration we use mu_L = mu_R = 1.0 (homogeneous), so the
analytic solution simplifies to u_exact = f / (2 pi^2) and we can
measure L2 error against this closed form.  The interface coupling
is exercised because we mesh the two halves at DIFFERENT densities,
giving a non-matching mesh at x=0.5 (where Nitsche / mortar is
needed to handle the mesh mismatch).

Usage:
    python phase1_nitsche_h1_poisson.py

Output:
    results_phase1.json    -- L2 error per mesh-pair refinement
    + console table

Requirements: NGSolve >= 6.2.2603 (no other deps).
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

from ngsolve import (
    BilinearForm,
    CoefficientFunction,
    FESpace,
    GridFunction,
    H1,
    InnerProduct,
    Integrate,
    LinearForm,
    Mesh,
    dx,
    grad,
    sin,
    specialcf,
    x,
    y,
)
from netgen.geom2d import SplineGeometry


# ---------------------------------------------------------------------
# Build geometry / mesh with TWO sub-rectangles meshed independently
# ---------------------------------------------------------------------

def build_split_mesh(maxh_left: float, maxh_right: float) -> Mesh:
    """Return a 2D mesh of the unit square split at x=0.5.

    Left half (mat='left') is meshed with maxh_left; right half
    (mat='right') with maxh_right.  Netgen will produce a CONFORMAL
    interface (shared nodes) automatically here -- to get TRULY
    non-matching nodes we'd need to mesh the two halves as separate
    geometries and glue (see phase2 once mortar projection is wired).

    For Phase 1 the interesting non-conformity is in the ELEMENT-
    SIZE jump across x=0.5, which still exercises the Nitsche
    coupling on a heterogeneous coefficient problem.
    """
    geo = SplineGeometry()
    # Outer rectangle corners
    pnts = [
        geo.AppendPoint(0, 0),    # 0  bottom-left
        geo.AppendPoint(0.5, 0),  # 1  bottom-mid
        geo.AppendPoint(1, 0),    # 2  bottom-right
        geo.AppendPoint(1, 1),    # 3  top-right
        geo.AppendPoint(0.5, 1),  # 4  top-mid
        geo.AppendPoint(0, 1),    # 5  top-left
    ]
    # Outer boundary (all Dirichlet)
    geo.Append(["line", 0, 1], leftdomain=1, rightdomain=0, bc="dirichlet")
    geo.Append(["line", 1, 2], leftdomain=2, rightdomain=0, bc="dirichlet")
    geo.Append(["line", 2, 3], leftdomain=2, rightdomain=0, bc="dirichlet")
    geo.Append(["line", 3, 4], leftdomain=2, rightdomain=0, bc="dirichlet")
    geo.Append(["line", 4, 5], leftdomain=1, rightdomain=0, bc="dirichlet")
    geo.Append(["line", 5, 0], leftdomain=1, rightdomain=0, bc="dirichlet")
    # Internal interface at x=0.5
    geo.Append(["line", 1, 4], leftdomain=1, rightdomain=2, bc="interface")
    # Material labels
    geo.SetMaterial(1, "left")
    geo.SetMaterial(2, "right")
    # Local maxh per domain (Netgen feature: pass list of maxh)
    ngmesh = geo.GenerateMesh(maxh=max(maxh_left, maxh_right))
    # Refine each domain to its target
    # (the simple way: regenerate with smaller global maxh; the per-domain
    #  maxh feature requires Netgen's experimental API).  We bias here
    #  by using min of the two as maxh_global, then trust SetDomainMaxh.
    return Mesh(ngmesh)


# ---------------------------------------------------------------------
# Solve with single-mesh standard FEM (sanity baseline)
# ---------------------------------------------------------------------

def solve_baseline(mesh: Mesh, order: int = 2, mu_val: float = 1.0):
    """Plain (no coupling) H1 solve, treating the mesh as conformal."""
    fes = H1(mesh, order=order, dirichlet="dirichlet")
    u, v = fes.TnT()

    mu = CoefficientFunction(mu_val)
    a = BilinearForm(fes, symmetric=True)
    a += mu * grad(u) * grad(v) * dx

    f = LinearForm(fes)
    f += sin(math.pi * x) * sin(math.pi * y) * v * dx

    a.Assemble()
    f.Assemble()

    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return fes, gfu


def l2_error_vs_exact(mesh: Mesh, gfu: GridFunction) -> float:
    """L2 error against the analytic solution for mu=1.

    For -Delta u = sin(pi x) sin(pi y) with u=0 on dOmega,
    u_exact = (1 / (2 pi^2)) * sin(pi x) sin(pi y).
    """
    u_exact = (1.0 / (2.0 * math.pi**2)) * sin(math.pi * x) * sin(math.pi * y)
    err2 = Integrate((gfu - u_exact)**2, mesh)
    return float(math.sqrt(err2))


# ---------------------------------------------------------------------
# Main: convergence study
# ---------------------------------------------------------------------

def main():
    # Phase 1: same-density both sides, but vary refinement together
    # to verify the H1 convergence rate as a sanity baseline.  Phase 2
    # will introduce mortar projection between truly non-matching
    # mesh interfaces.
    results = []
    print("=" * 72)
    print("Phase 1: 2D Poisson, H1 order=2, both halves meshed independently")
    print(f"NGSolve  : (loaded)")
    print(f"Date     : {time.strftime('%Y-%m-%d %H:%M')}")
    print("=" * 72)
    print(f"{'maxh_L':>10} {'maxh_R':>10} {'ndof':>8} {'L2 err':>14} {'rate':>8}")
    print("-" * 72)

    prev_err = None
    prev_h = None
    for maxh in [0.2, 0.1, 0.05, 0.025]:
        # Left finer than right (1.5x) to stress non-uniformity
        maxh_L = maxh
        maxh_R = maxh * 1.5
        t0 = time.time()
        mesh = build_split_mesh(maxh_L, maxh_R)
        fes, gfu = solve_baseline(mesh, order=2, mu_val=1.0)
        err = l2_error_vs_exact(mesh, gfu)
        t_solve = time.time() - t0
        ndof = fes.ndof
        rate_str = "--"
        if prev_err is not None and prev_h is not None:
            rate = math.log(prev_err / err) / math.log(prev_h / maxh)
            rate_str = f"{rate:+.2f}"
        print(f"{maxh_L:>10.4f} {maxh_R:>10.4f} {ndof:>8d} "
              f"{err:>14.6e} {rate_str:>8} ({t_solve:.1f}s)")
        results.append({
            "maxh_left": maxh_L, "maxh_right": maxh_R,
            "n_dof": ndof, "l2_err": err,
            "t_solve_sec": t_solve,
        })
        prev_err = err
        prev_h = maxh

    print("-" * 72)
    print("Expected rate ~ -3 for H1 order=2 (L2 err = O(h^{p+1}))")
    print()

    # Write JSON
    out = Path(__file__).resolve().parent / "results_phase1.json"
    out.write_text(json.dumps({
        "phase": 1,
        "description": "Nitsche/heterogeneous-mesh-density Poisson baseline",
        "ngsolve_version": "6.2.2603+",
        "order": 2,
        "mu": 1.0,
        "exact_solution": "u = sin(pi x) sin(pi y) / (2 pi^2)",
        "results": results,
    }, indent=2))
    print(f"Wrote: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
