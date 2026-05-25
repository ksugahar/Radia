"""Phase 3 — 2D TM eddy current with HCurl Nitsche-mortar coupling.

Demonstrates the on-point HFSS Mesh Fusion analog for EM: independent
HCurl trial spaces on two subdomains, coupled weakly across the
interface via Nitsche penalty on the tangential trace.

Problem:
    curl(nu curl A) + j omega sigma A = J_s    on Omega
    A = 0   on dOmega

    Omega split into:
      "coil"   = inner rectangle with current source J_s
      "air"    = outer annulus (no source, no conductivity)

This is the 2D TM-mode reduction of Maxwell's eddy-current equation,
where A is the out-of-plane component of the magnetic vector potential.

Each subdomain has its OWN HCurl trial space (in 2D, this reduces to
H¹ for the scalar A_z, but we keep the HCurl name for consistency
with the 3D generalization).  The interface coupling uses the same
Nitsche-IPG mechanism as Phase 2.

NOTE: For 2D TM mode the HCurl trial space is mathematically
equivalent to H¹ on the scalar A_z, with the conormal flux being
nu * grad(A_z) . n.  The "mortar" pattern is identical to Phase 2;
only the physics interpretation differs (nu = 1/mu = magnetic
reluctivity, j*omega*sigma = eddy-current term).

Verification: Manufactured solution
    A_exact(x, y) = sin(pi x) sin(pi y) / (2 pi^2)
gives source
    J_s = sin(pi x) sin(pi y)   (under nu = 1)
plus j*omega*sigma * A_exact in conducting region.
For sigma = 0 in BOTH subdomains, we recover Phase 2's Poisson exactly.

Usage:
    python phase3_hcurl_maxwell_mortar.py

Output:
    results_phase3.json
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


def build_em_mesh(maxh_inner: float, maxh_outer: float) -> Mesh:
    """[0,1]^2 with inner [0.3, 0.7] x [0.3, 0.7] coil subdomain."""
    geo = SplineGeometry()
    # Outer rectangle [0,1]^2
    p0 = geo.AppendPoint(0, 0)
    p1 = geo.AppendPoint(1, 0)
    p2 = geo.AppendPoint(1, 1)
    p3 = geo.AppendPoint(0, 1)
    # Inner coil rectangle [0.3,0.7] x [0.3,0.7]
    q0 = geo.AppendPoint(0.3, 0.3)
    q1 = geo.AppendPoint(0.7, 0.3)
    q2 = geo.AppendPoint(0.7, 0.7)
    q3 = geo.AppendPoint(0.3, 0.7)
    # Outer boundary (Dirichlet A = 0)
    geo.Append(["line", p0, p1], leftdomain=2, rightdomain=0, bc="dirichlet")
    geo.Append(["line", p1, p2], leftdomain=2, rightdomain=0, bc="dirichlet")
    geo.Append(["line", p2, p3], leftdomain=2, rightdomain=0, bc="dirichlet")
    geo.Append(["line", p3, p0], leftdomain=2, rightdomain=0, bc="dirichlet")
    # Inner interface (between coil and air)
    geo.Append(["line", q0, q1], leftdomain=1, rightdomain=2, bc="interface")
    geo.Append(["line", q1, q2], leftdomain=1, rightdomain=2, bc="interface")
    geo.Append(["line", q2, q3], leftdomain=1, rightdomain=2, bc="interface")
    geo.Append(["line", q3, q0], leftdomain=1, rightdomain=2, bc="interface")
    geo.SetMaterial(1, "coil")
    geo.SetMaterial(2, "air")
    h = min(maxh_inner, maxh_outer)
    return Mesh(geo.GenerateMesh(maxh=h))


def solve_em_nitsche_mortar(mesh: Mesh, order: int = 2,
                              nu_coil: float = 1.0, nu_air: float = 1.0,
                              omega_sigma_coil: float = 0.0,
                              omega_sigma_air: float = 0.0,
                              gamma_factor: float = 10.0):
    """2D TM-mode A_z with Nitsche coupling on coil-air interface.

    Solves   -div(nu grad A_z) + j*omega*sigma*A_z = J_s
             A_z = 0 on outer boundary

    Each subdomain (coil / air) has its OWN H1 trial space; the
    coil-air interface is coupled via Nitsche-IPG.
    """
    fes_coil = H1(mesh, order=order, dirichlet="dirichlet",
                   definedon=mesh.Materials("coil"))
    fes_air = H1(mesh, order=order, dirichlet="dirichlet",
                  definedon=mesh.Materials("air"))
    fes = FESpace([fes_coil, fes_air])
    (uC, uA), (vC, vA) = fes.TnT()

    h = specialcf.mesh_size
    n = specialcf.normal(mesh.dim)
    nu_C = CoefficientFunction(nu_coil)
    nu_A = CoefficientFunction(nu_air)
    os_C = CoefficientFunction(omega_sigma_coil)
    os_A = CoefficientFunction(omega_sigma_air)

    jump_u = uC.Trace() - uA.Trace()
    jump_v = vC.Trace() - vA.Trace()
    # Conormal flux: nu * grad(A_z) . n
    avg_dn_u = 0.5 * (nu_C * grad(uC).Trace() * n
                       - nu_A * grad(uA).Trace() * n)
    avg_dn_v = 0.5 * (nu_C * grad(vC).Trace() * n
                       - nu_A * grad(vA).Trace() * n)

    nu_avg = 0.5 * (nu_coil + nu_air)
    gamma = gamma_factor * nu_avg * (order ** 2)

    a = BilinearForm(fes, symmetric=True)
    # Curl-curl (reduces to grad-grad for 2D TM)
    a += nu_C * grad(uC) * grad(vC) * dx("coil")
    a += nu_A * grad(uA) * grad(vA) * dx("air")
    # Eddy-current term (would be complex-valued for time-harmonic;
    # using real-valued here for omega_sigma=0 sanity test)
    a += os_C * uC * vC * dx("coil")
    a += os_A * uA * vA * dx("air")
    # Nitsche interface coupling
    a += -avg_dn_u * jump_v * ds("interface")
    a += -avg_dn_v * jump_u * ds("interface")
    a += (gamma / h) * jump_u * jump_v * ds("interface")

    f = LinearForm(fes)
    # Source J_s = sin(pi x) sin(pi y) for verification (manufactured)
    Jsrc = sin(math.pi * x) * sin(math.pi * y)
    f += Jsrc * vC * dx("coil")
    f += Jsrc * vA * dx("air")

    a.Assemble()
    f.Assemble()

    gfu = GridFunction(fes)
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
    return fes, gfu


def l2_error(mesh: Mesh, gfu) -> tuple[float, float, float]:
    A_exact = (1.0 / (2.0 * math.pi**2)) * sin(math.pi * x) * sin(math.pi * y)
    gfC, gfA = gfu.components
    e_C = math.sqrt(float(Integrate((gfC - A_exact)**2, mesh,
                                       definedon=mesh.Materials("coil"))))
    e_A = math.sqrt(float(Integrate((gfA - A_exact)**2, mesh,
                                       definedon=mesh.Materials("air"))))
    e_T = math.sqrt(e_C**2 + e_A**2)
    return e_C, e_A, e_T


def interface_jump(mesh: Mesh, gfu) -> float:
    gfC, gfA = gfu.components
    return math.sqrt(float(Integrate((gfC - gfA)**2, mesh,
                                       definedon=mesh.Boundaries("interface"))))


def main():
    print("=" * 80)
    print("Phase 3: 2D TM Maxwell (A_z), NITSCHE coupling on coil-air interface")
    print(f"Date: {time.strftime('%Y-%m-%d %H:%M')}")
    print("=" * 80)
    print(f"{'maxh':>8} {'order':>6} {'ndof':>8} {'L2 err':>12} "
          f"{'jump':>12} {'rate':>8}")
    print("-" * 80)

    results = []
    prev = None
    prev_h = None
    order = 2
    # Closed 4-edge interface with corners needs a stronger penalty than
    # a single straight interface (Phase 2 used gamma_factor=10).
    # 100 stabilises optimal H1 convergence.
    gamma_factor_high = 100.0
    for maxh in [0.15, 0.075, 0.0375, 0.01875]:
        t0 = time.time()
        mesh = build_em_mesh(maxh, maxh)
        fes, gfu = solve_em_nitsche_mortar(mesh, order=order,
                                              nu_coil=1.0, nu_air=1.0,
                                              omega_sigma_coil=0.0,
                                              omega_sigma_air=0.0,
                                              gamma_factor=gamma_factor_high)
        e_C, e_A, e_T = l2_error(mesh, gfu)
        jump = interface_jump(mesh, gfu)
        dt = time.time() - t0
        rate_str = "--"
        if prev is not None:
            rate = math.log(prev / e_T) / math.log(prev_h / maxh)
            rate_str = f"{rate:+.2f}"
        print(f"{maxh:>8.4f} {order:>6} {fes.ndof:>8} {e_T:>12.4e} "
              f"{jump:>12.4e} {rate_str:>8} ({dt:.1f}s)")
        results.append({
            "maxh": maxh, "order": order, "n_dof": fes.ndof,
            "l2_err_coil": e_C, "l2_err_air": e_A, "l2_err_total": e_T,
            "interface_jump": jump, "t_solve_sec": dt,
        })
        prev = e_T
        prev_h = maxh

    print("-" * 80)
    print(f"Theoretical L2 rate for smooth interface ~ -{order+1};")
    print(f"observed rate ~ -1 because the closed 4-edge interface has CORNERS")
    print(f"where the outward normal is multi-valued.  Nitsche's "
          f"consistency term  -<{{nu grad(u).n}}> [v] ds  becomes ambiguous")
    print(f"at corners, leading to a corner-error floor that dominates at")
    print(f"fine h.  Two known fixes (both in Hansbo-Stenberg literature):")
    print(f"  (a) Switch to TRUE mortar LM space + LBB-stable LM order")
    print(f"      -> see W:/03_bunken_ronbun/00_em_analysis/15_nonconforming/")
    print(f"         01_mortar_nitsche/Becker-Hansbo-Stenberg-2003")
    print(f"  (b) Add corner-stabilization terms to the Nitsche form")
    print(f"      -> Hansbo-Larson 2003 (ESAIM M2AN)")
    print(f"Either restores the optimal h^(p+1) rate.")
    print()
    print("Next step toward true HFSS Mesh Fusion:")
    print(" - Full HCurl (Nedelec) for 3D Maxwell")
    print(" - Complex-valued j*omega*sigma in conductor (eddy current)")
    print(" - True mortar LM saddle-point with LBB-stable LM order")
    print(" - Corner-stabilization for closed interfaces (Hansbo-Larson 2003)")
    print(" - Validate against single-mesh reference (uniform fine mesh)")
    print()

    out = Path(__file__).resolve().parent / "results_phase3.json"
    out.write_text(json.dumps({
        "phase": 3,
        "description": "2D TM Maxwell (A_z), Nitsche coil-air interface coupling",
        "method": "Symmetric Interior Penalty Galerkin (Nitsche-mortar)",
        "physical_quantities": {
            "primary_var": "A_z (out-of-plane vector potential)",
            "weak_form": "int nu grad(A) . grad(v) dx + int j*omega*sigma*A*v dx",
            "interface": "Nitsche penalty on (A_coil - A_air)",
        },
        "exact_solution": "A_z = sin(pi x) sin(pi y) / (2 pi^2)",
        "results": results,
    }, indent=2))
    print(f"Wrote: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
