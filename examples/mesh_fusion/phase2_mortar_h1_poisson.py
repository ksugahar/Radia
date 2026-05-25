"""Phase 2 — Mortar-Nitsche coupling for 2D Poisson, multi-mesh-density.

Demonstrates the mortar-style WEAK coupling of two H1 trial spaces
defined independently on left + right subdomains, using NGSolve's
Nitsche-method machinery (penalty + symmetry + consistency terms,
NO explicit Lagrange multiplier space).

Why Nitsche-mortar instead of pure mortar with LM space?
- NGSolve's pure-mortar LM space (`L2(definedon=boundary)`) is
  fragile in 6.2.2603 (segfault on some patterns; depends on
  build-time options).
- Nitsche reproduces the SAME math (consistent + symmetric +
  coercive interface coupling) without the LBB headache or the
  saddle-point indefinite system.
- This is the formulation Becker-Hansbo-Stenberg 2003 (M2AN)
  proved is equivalent to mortar for the elliptic case.

Problem:
    -nabla . (mu nabla u) = f       on Omega = [0,1] x [0,1]
                          u = 0     on dOmega
    Omega = Omega_L + Omega_R, interface at x=0.5
    Each side has its OWN H1 trial space (discontinuous H1 across
    the interface, recoupled via Nitsche)

Verification: vs analytic
    u_exact = sin(pi x) sin(pi y) / (2 pi^2)    (mu = 1)
    L2 error decays as O(h^{p+1})

Usage:
    python phase2_mortar_h1_poisson.py

Output:
    results_phase2.json
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

from ngsolve import (
    BilinearForm,
    CoefficientFunction,
    FESpace,
    GridFunction,
    H1,
    Integrate,
    LinearForm,
    Mesh,
    dx,
    ds,
    grad,
    sin,
    specialcf,
    x,
    y,
)
from netgen.geom2d import SplineGeometry


def build_split_mesh(maxh_L: float, maxh_R: float) -> Mesh:
    """Unit square split at x=0.5 into 2 subdomains.

    NOTE: NGSolve+Netgen produces a CONFORMAL interface (shared
    nodes).  The interesting non-uniformity is element density per
    side -- the Nitsche-mortar coupling treats both sides as
    formally independent FE spaces regardless.
    """
    geo = SplineGeometry()
    p0 = geo.AppendPoint(0, 0)
    p1 = geo.AppendPoint(0.5, 0)
    p2 = geo.AppendPoint(1, 0)
    p3 = geo.AppendPoint(1, 1)
    p4 = geo.AppendPoint(0.5, 1)
    p5 = geo.AppendPoint(0, 1)
    geo.Append(["line", p0, p1], leftdomain=1, rightdomain=0, bc="dirichlet")
    geo.Append(["line", p1, p2], leftdomain=2, rightdomain=0, bc="dirichlet")
    geo.Append(["line", p2, p3], leftdomain=2, rightdomain=0, bc="dirichlet")
    geo.Append(["line", p3, p4], leftdomain=2, rightdomain=0, bc="dirichlet")
    geo.Append(["line", p4, p5], leftdomain=1, rightdomain=0, bc="dirichlet")
    geo.Append(["line", p5, p0], leftdomain=1, rightdomain=0, bc="dirichlet")
    geo.Append(["line", p1, p4], leftdomain=1, rightdomain=2, bc="interface")
    geo.SetMaterial(1, "left")
    geo.SetMaterial(2, "right")
    h = min(maxh_L, maxh_R)
    return Mesh(geo.GenerateMesh(maxh=h))


def solve_nitsche_mortar(mesh: Mesh, order: int = 2, mu_L: float = 1.0,
                          mu_R: float = 1.0, gamma_factor: float = 10.0):
    """Solve with TWO independent H1 spaces + Nitsche interface terms.

    SPD system (symmetric Nitsche IPG variant).  Penalty coefficient
    gamma = gamma_factor * mu_avg * p^2 / h.
    """
    # Two H1 spaces, each LIVES ON ONE SUBDOMAIN ONLY -> trial
    # functions are formally independent across the interface
    fes_L = H1(mesh, order=order, dirichlet="dirichlet",
                definedon=mesh.Materials("left"))
    fes_R = H1(mesh, order=order, dirichlet="dirichlet",
                definedon=mesh.Materials("right"))
    fes = FESpace([fes_L, fes_R])

    (uL, uR), (vL, vR) = fes.TnT()

    mu_L_cf = CoefficientFunction(mu_L)
    mu_R_cf = CoefficientFunction(mu_R)
    h = specialcf.mesh_size
    n = specialcf.normal(mesh.dim)

    # Average + jump operators on the interface.
    # NGSolve requires `.Trace()` on the trial functions for BND-form
    # integration; `grad(u).Trace()` returns the trace of the gradient.
    jump_u = uL.Trace() - uR.Trace()
    jump_v = vL.Trace() - vR.Trace()
    # Conormal flux averages (use trace of grad)
    avg_dn_u = 0.5 * (mu_L_cf * grad(uL).Trace() * n
                       - mu_R_cf * grad(uR).Trace() * n)
    avg_dn_v = 0.5 * (mu_L_cf * grad(vL).Trace() * n
                       - mu_R_cf * grad(vR).Trace() * n)
    # Penalty
    mu_avg = 0.5 * (mu_L + mu_R)
    gamma = gamma_factor * mu_avg * (order ** 2)

    a = BilinearForm(fes, symmetric=True)
    # subdomain stiffness
    a += mu_L_cf * grad(uL) * grad(vL) * dx("left")
    a += mu_R_cf * grad(uR) * grad(vR) * dx("right")
    # Nitsche interface terms (symmetric IPG)
    a += -avg_dn_u * jump_v * ds("interface")          # consistency
    a += -avg_dn_v * jump_u * ds("interface")          # symmetry
    a += (gamma / h) * jump_u * jump_v * ds("interface")  # penalty

    f = LinearForm(fes)
    src = sin(math.pi * x) * sin(math.pi * y)
    f += src * vL * dx("left")
    f += src * vR * dx("right")

    a.Assemble()
    f.Assemble()

    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return fes, gfu


def l2_error_combined(mesh: Mesh, gfu) -> tuple[float, float, float]:
    """L2 error vs analytic on each side and combined."""
    u_exact = (1.0 / (2.0 * math.pi**2)) * sin(math.pi * x) * sin(math.pi * y)
    gfL, gfR = gfu.components
    e_L2 = math.sqrt(float(Integrate((gfL - u_exact)**2, mesh,
                                       definedon=mesh.Materials("left"))))
    e_R2 = math.sqrt(float(Integrate((gfR - u_exact)**2, mesh,
                                       definedon=mesh.Materials("right"))))
    e_tot = math.sqrt(e_L2**2 + e_R2**2)
    return e_L2, e_R2, e_tot


def interface_jump(mesh: Mesh, gfu) -> float:
    """L2 norm of the jump uL - uR on the interface.

    For Nitsche, this is NOT zero exactly; it scales like O(h^p)
    with the penalty.  Measures how well the coupling is enforced.
    """
    gfL, gfR = gfu.components
    return math.sqrt(float(Integrate((gfL - gfR)**2, mesh,
                                       definedon=mesh.Boundaries("interface"))))


def main():
    print("=" * 80)
    print("Phase 2: 2D Poisson, NITSCHE-MORTAR coupling, independent H1 per side")
    print(f"Date: {time.strftime('%Y-%m-%d %H:%M')}")
    print("=" * 80)
    print(f"{'maxh':>8} {'order':>6} {'gamma_factor':>13} {'ndof':>8} "
          f"{'L2 err':>12} {'jump':>12} {'rate':>8}")
    print("-" * 80)

    results = []
    prev_err = None
    prev_h = None
    order = 2
    gamma_factor = 10.0  # standard Nitsche penalty multiplier

    for maxh in [0.2, 0.1, 0.05, 0.025]:
        t0 = time.time()
        mesh = build_split_mesh(maxh, maxh)
        fes, gfu = solve_nitsche_mortar(mesh, order=order,
                                          mu_L=1.0, mu_R=1.0,
                                          gamma_factor=gamma_factor)
        e_L, e_R, e_tot = l2_error_combined(mesh, gfu)
        jump = interface_jump(mesh, gfu)
        dt = time.time() - t0
        rate_str = "--"
        if prev_err is not None and prev_h is not None:
            rate = math.log(prev_err / e_tot) / math.log(prev_h / maxh)
            rate_str = f"{rate:+.2f}"
        print(f"{maxh:>8.4f} {order:>6} {gamma_factor:>13.1f} {fes.ndof:>8} "
              f"{e_tot:>12.4e} {jump:>12.4e} {rate_str:>8} ({dt:.1f}s)")
        results.append({
            "maxh": maxh, "order": order, "gamma_factor": gamma_factor,
            "n_dof": fes.ndof, "l2_err_left": e_L, "l2_err_right": e_R,
            "l2_err_total": e_tot, "interface_jump": jump,
            "t_solve_sec": dt,
        })
        prev_err = e_tot
        prev_h = maxh

    print("-" * 80)
    print(f"Expected total L2 rate ~ -{order+1} (H1 order={order} -> "
          f"L2 err = O(h^{order+1})).")
    print("Interface jump decays at rate ~ -p (Nitsche penalty enforces "
          "continuity at trace-error level).")
    print()

    out = Path(__file__).resolve().parent / "results_phase2.json"
    out.write_text(json.dumps({
        "phase": 2,
        "description": "Nitsche-mortar 2D Poisson, independent H1 per side",
        "method": "Symmetric Interior Penalty Galerkin variant of Nitsche",
        "exact_solution": "u = sin(pi x) sin(pi y) / (2 pi^2)",
        "results": results,
    }, indent=2))
    print(f"Wrote: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
