"""calc_axisym_volumetric.py -- 2D axisymmetric A_phi volumetric eddy current solver.

The "truth" reference solver for ESIM validation: instead of imposing a
surface-impedance Robin BC (SIBC / ESIM), this script resolves the
volumetric eddy current inside the workpiece directly via the
axisymmetric magnetic vector potential `A_phi` on a `(r, z)` mesh.

Formulation
-----------

Time-harmonic Maxwell in axisymmetric `(r, z)`, single component
`A_phi`.  Following Henrotte / Meeker FEMM convention, the discrete
unknown is the flux function `psi = 2 pi r A_phi` interpolated linearly
in `(s = r^2, z)`; element matrices follow from
``radia.radia_axifemm.{AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI}``.

Strong form:

    1/mu * (-d_r^2 A_phi + 1/r d_r A_phi - 1/r^2 A_phi - d_z^2 A_phi)
        + j w sigma A_phi = J_phi_source

where ``J_phi_source`` is the impressed coil current density (real, A/m^2)
and ``sigma > 0`` only inside the workpiece (eddy current induced).

Weak form (in flux variable psi, after Henrotte change of variable):

    int (1/mu) grad_s_z psi . grad_s_z psi' / (2 pi r)
      + j w sigma psi psi' / (2 pi r)
      = int J_phi_source psi' . (r-weighted)

(Closed-form per-element by axifemm BFIs.)

Boundary conditions
-------------------

- Axis (r = 0): A_phi = 0  (Dirichlet)
- Far field: A_phi = 0 (truncation; replace with Kelvin transform when needed)

Workflow
--------

1. Build the (r, z) mesh: workpiece + air + (optionally) Kelvin annulus.
2. Define material CFs: 1/mu(r,z), sigma(r,z), J_coil(r,z).
3. AC complex assembly: K + j w M (with K from stiffness BFI, M from
   sigma mass BFI).
4. Linear/nonlinear solve.  For nonlinear mu(|B|), an outer Picard
   iteration on the per-element |B| (max norm over quad-points).
5. Post: P_wp = (1/2) Re( int_workpiece sigma |j w A_phi|^2 dV_axisym ).

The volumetric P_wp is the "truth" against which scalar and
per-element ESIM (calc_inductance.py output) are compared.

Status (2026-05-22)
-------------------

Phase 1 (linear mu): WORK IN PROGRESS.
Phase 2 (nonlinear BH): not started.

This is task #36 in the project task list -- the validation reference
needed to determine which of {scalar ESIM, per-element ESIM} is closer
to truth on the IH benchmark.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

MU0 = 4.0e-7 * math.pi


def build_mesh(R_wp_m: float, H_wp_m: float,
               R_coil_m: float, R_outer_m: float,
               maxh_wp_m: float, maxh_air_m: float):
    """Axis-aligned 2D mesh: workpiece (cylinder cross-section in r-z) +
    surrounding air box.  Materials are tagged so a downstream CF can
    distinguish wp vs air vs coil_ring.
    """
    from ngsolve import Mesh
    from netgen.occ import OCCGeometry, MoveTo, Glue, X, Y

    # Workpiece: r in [0, R_wp], z in [-H_wp/2, H_wp/2]
    wp = MoveTo(0, -H_wp_m / 2).Rectangle(R_wp_m, H_wp_m).Face()
    wp.faces.name = "workpiece"
    wp.maxh = maxh_wp_m

    # Coil ring: thin rectangle at (R_coil, 0).  Cross-section ~ 2x2 mm.
    coil_h = 2e-3
    coil_w = 2e-3
    coil = MoveTo(R_coil_m - coil_w / 2, -coil_h / 2) \
        .Rectangle(coil_w, coil_h).Face()
    coil.faces.name = "coil"
    coil.maxh = maxh_wp_m

    # Air box: full surrounding region.
    air = MoveTo(0, -R_outer_m).Rectangle(R_outer_m, 2 * R_outer_m).Face()
    air.faces.name = "air"
    air.edges.Max(X).name = "outer"
    air.edges.Min(X).name = "axis"
    air.edges.Max(Y).name = "top"
    air.edges.Min(Y).name = "bot"

    shape = Glue([air, wp, coil])
    return Mesh(OCCGeometry(shape, dim=2).GenerateMesh(maxh=maxh_air_m))


def run_axisym_linear(args):
    """Phase 1: linear mu_r in workpiece.  Validates against Bessel
    cylinder when the field is roughly uniform along z (long-cylinder
    limit).
    """
    from ngsolve import (
        CoefficientFunction, BilinearForm, LinearForm, GridFunction,
        Integrate, dx, ds, grad, InnerProduct, Periodic, x as r_cf,
        FESpace,
    )
    import radia.radia_axifemm   # registers axihenrotte FESpace

    mesh = build_mesh(
        R_wp_m=args.R_wp, H_wp_m=args.H_wp,
        R_coil_m=args.R_coil, R_outer_m=args.R_outer,
        maxh_wp_m=args.maxh_wp, maxh_air_m=args.maxh_air,
    )
    print(f"mesh: ne={mesh.ne}, nv={mesh.nv}, mats={mesh.GetMaterials()}")

    # Material CFs.  Order matches mesh.GetMaterials().
    mu_r_d = {"workpiece": args.mu_r,
              "coil": 1.0,
              "air": 1.0}
    sigma_d = {"workpiece": args.sigma,
               "coil": 0.0,
               "air": 0.0}
    mu_cf = mesh.MaterialCF({k: v * MU0 for k, v in mu_r_d.items()},
                             default=MU0)
    sigma_cf = mesh.MaterialCF(sigma_d, default=0.0)

    # Coil source: J_phi = I / A_coil in the coil region only.
    A_coil = 2e-3 * 2e-3      # 2x2 mm cross-section
    J_phi = mesh.MaterialCF({"coil": args.current / A_coil,
                              "workpiece": 0.0,
                              "air": 0.0}, default=0.0)

    # AC complex FE problem: solve (K(1/mu) + j w sigma M) A = J.
    # Use TWO REAL FE solves coupled via:
    #   [ K   -wM ] [Re A]   [Re J]
    #   [ wM   K  ] [Im A] = [Im J]
    # because axifemm BFIs accept only real CFs.  Standard NGSolve
    # H1 with 2 pi r weighting works equivalently for axisym; the
    # Henrotte BFI is preferred for accuracy near r=0 but for our
    # workpiece interior (r > 0 mostly) standard H1 is fine.
    #
    # First pass: use standard H1 + 2 pi r weighting (simpler).  If
    # accuracy at axis is needed, switch to Henrotte later.
    p = args.order
    omega = 2 * math.pi * args.frequency

    # Complex FE problem.  NGSolve H1 with complex=True handles it.
    from ngsolve import H1
    fes = H1(mesh, order=p, complex=True,
              dirichlet="axis|outer|top|bot")
    u, v = fes.TnT()
    print(f"FES: order={p}, ndof={fes.ndof}")

    # Axisymmetric weak form for A_phi:
    #
    #   int (1/mu) [ d_r(r u) d_r(r v) / r + r d_z u d_z v ] dr dz
    #     + j w int sigma r u v dr dz
    #   = int r u_source v dr dz
    #
    # The `r d_r(A_phi) + A_phi` combination handles the r=0 boundary
    # gracefully when paired with Dirichlet u=0 on axis.  We use the
    # equivalent form  curl(A_phi e_phi) = (1/r)(d_r(r A_phi)) e_z - d_z A_phi e_r
    # and its energy ||curl||^2 = ... = (1/r)(d_r(r A))^2 + (d_z A)^2.
    # Times mu_0 (1/mu) and integrated with 2 pi r dr dz Jacobian:
    #
    #   E_field = pi int [(d_r(rA))^2 / r + r (d_z A)^2] / mu  dr dz
    #
    # In NGSolve weak-form notation:
    a = BilinearForm(fes, symmetric=False)
    a += (1.0 / mu_cf) * (1.0 / r_cf) * \
        (r_cf * grad(u)[0] + u) * (r_cf * grad(v)[0] + v) * dx
    a += (1.0 / mu_cf) * r_cf * grad(u)[1] * grad(v)[1] * dx
    a += 1j * omega * sigma_cf * r_cf * u * v * dx

    f = LinearForm(fes)
    f += J_phi * r_cf * v * dx

    t0 = time.perf_counter()
    a.Assemble()
    f.Assemble()
    print(f"assembled in {time.perf_counter()-t0:.2f}s")

    gfu = GridFunction(fes)
    t0 = time.perf_counter()
    gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f.vec
    print(f"solved in {time.perf_counter()-t0:.2f}s")

    # P_wp = (1/2) Re( int_workpiece sigma |jw A_phi|^2 2 pi r dr dz )
    A_norm_sq = InnerProduct(gfu, gfu)   # |A|^2 = A . conj(A)
    P_wp_density = 0.5 * sigma_cf * (omega ** 2) * A_norm_sq * 2 * math.pi * r_cf
    P_wp = float(Integrate(P_wp_density.real, mesh,
                            definedon=mesh.Materials("workpiece")))
    print(f"P_wp_volumetric = {P_wp:.6e} W")

    out = {
        "method": "axisym_volumetric_A_phi",
        "frequency_Hz": args.frequency,
        "current_A": args.current,
        "R_wp_m": args.R_wp,
        "H_wp_m": args.H_wp,
        "R_coil_m": args.R_coil,
        "mu_r": args.mu_r,
        "sigma_S_per_m": args.sigma,
        "fes_order": p,
        "ndof": fes.ndof,
        "P_wp_W": P_wp,
    }
    if args.output:
        with open(args.output, "w") as fh:
            json.dump(out, fh, indent=2)
        print(f"wrote {args.output}")
    return out


def main():
    parser = argparse.ArgumentParser(
        description="2D axisymmetric A_phi volumetric eddy current solver "
                    "(ESIM truth reference).")
    parser.add_argument("--R-wp", type=float, default=5e-3,
                        help="Workpiece outer radius [m]")
    parser.add_argument("--H-wp", type=float, default=10e-3,
                        help="Workpiece height (z-extent) [m]")
    parser.add_argument("--R-coil", type=float, default=20e-3,
                        help="Coil ring radius [m]")
    parser.add_argument("--R-outer", type=float, default=200e-3,
                        help="Outer air-box radius (truncation) [m]")
    parser.add_argument("--maxh-wp", type=float, default=0.2e-3,
                        help="Mesh size inside workpiece [m]")
    parser.add_argument("--maxh-air", type=float, default=5e-3,
                        help="Mesh size in air [m]")
    parser.add_argument("--mu-r", type=float, default=100.0)
    parser.add_argument("--sigma", type=float, default=2e6)
    parser.add_argument("--frequency", type=float, default=50000.0)
    parser.add_argument("--current", type=float, default=100.0)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()
    run_axisym_linear(args)


if __name__ == "__main__":
    main()
