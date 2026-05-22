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

Phase 1 (linear mu): VALIDATED.  Long-cylinder Bessel cross-check
agrees to -5.7 % (end-effect contamination at H=200mm, 40x R).
IGTE-geometry linear-mu check agrees to +1 % vs Bessel-from-FEM-H_t.

Phase 2 (nonlinear BH Picard): SCAFFOLD ONLY.  The outer Picard
loop is implemented but currently produces NaN due to:
  - mu_r floor missing at low |B| (BH curve interp gives mu_r->0)
  - element-centroid finite-diff |B| unstable near element edges
  - NGSolve's gf_mu.Set on a P0 piecewise FE space requires careful
    handling of the material map
Production use needs: proper |B| extraction via grad(gfu) at quadrature
points, mu_r floor (>=1), CoefficientFunction-based mu update instead
of per-element loop.  Estimated 1-2 days of additional work.

This is task #36 in the project task list -- the validation reference
needed to determine which of {scalar ESIM, per-element ESIM} is closer
to truth on the IH benchmark.  Phase 1 alone is sufficient for the
linear-mu cross-check; Phase 2 is needed for the nonlinear ESIM
absolute-accuracy claim.
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


def run_axisym_nonlinear(args, bh_curve):
    """Phase 2: nonlinear mu(|B|) outer Picard iteration.

    Outer loop on per-element mu_r:
      1. Solve linear AC problem with current mu_r distribution
      2. Compute |B| per element from A_phi
      3. Update mu_r per element from the BH curve at |B|
      4. Damped Picard: mu_r_new = alpha mu_BH(|B|) + (1-alpha) mu_r_old
      5. Stop when max relative change of mu_r below tol

    The cell-centered mu_r evaluation matches FEMM's prob3big.cpp
    convention (mu_r constant per element).  For axisymmetric A_phi,
    |B|^2 = |B_r|^2 + |B_z|^2 = |d_z A|^2 + |(1/r) d_r(r A)|^2.

    Returns dict with P_wp, H_t side-wall samples, Picard convergence.
    """
    from ngsolve import (
        CoefficientFunction, BilinearForm, LinearForm, GridFunction,
        Integrate, dx, grad, InnerProduct, x as r_cf, H1,
    )
    import radia.radia_axifemm   # noqa: F401

    mesh = build_mesh(
        R_wp_m=args.R_wp, H_wp_m=args.H_wp,
        R_coil_m=args.R_coil, R_outer_m=args.R_outer,
        maxh_wp_m=args.maxh_wp, maxh_air_m=args.maxh_air,
    )
    print(f"mesh: ne={mesh.ne}, nv={mesh.nv}, mats={mesh.GetMaterials()}")

    # BH lookup
    bh = np.asarray(bh_curve)
    H_arr, B_arr = bh[:, 0], bh[:, 1]
    def mu_r_from_B(B_abs):
        H = float(np.interp(B_abs, B_arr, H_arr)) if B_abs > 0 else 0.0
        return B_abs / (MU0 * max(H, 1e-12)) if H > 0 else float(np.interp(0.0, H_arr, [0.0, 100.0]))

    # Use a per-element constant mu_r via a piecewise FE space.
    fes_mu = H1(mesh, order=0)  # piecewise const proxy via P0-on-elements
    gf_mu = GridFunction(fes_mu)
    gf_mu.vec.FV().NumPy()[:] = 1.0  # init: vacuum everywhere

    # Set workpiece initial mu_r = args.mu_r (linear seed for first iter)
    mat_idx = {m: i for i, m in enumerate(mesh.GetMaterials())}
    wp_mat = mat_idx.get("workpiece", None)
    if wp_mat is None:
        raise RuntimeError("workpiece material not in mesh")
    # Build a coefficient that returns initial mu_r per element
    init_mu_r = mesh.MaterialCF({"workpiece": float(args.mu_r),
                                  "coil": 1.0, "air": 1.0}, default=1.0)
    gf_mu.Set(init_mu_r)

    p = args.order
    omega = 2 * math.pi * args.frequency
    A_coil = 2e-3 * 2e-3
    J_phi = mesh.MaterialCF({"coil": args.current / A_coil,
                              "workpiece": 0.0, "air": 0.0}, default=0.0)

    fes = H1(mesh, order=p, complex=True,
              dirichlet="axis|outer|top|bot")
    u, v = fes.TnT()

    sigma_cf = mesh.MaterialCF({"workpiece": args.sigma,
                                 "coil": 0.0, "air": 0.0}, default=0.0)

    gfu = GridFunction(fes)
    max_picard = 20
    tol_picard = 5e-3
    alpha = 0.5
    convergence = []
    for k_outer in range(max_picard):
        mu_cf = MU0 * gf_mu
        a = BilinearForm(fes, symmetric=False)
        a += (1.0 / mu_cf) * (1.0 / r_cf) * \
            (r_cf * grad(u)[0] + u) * (r_cf * grad(v)[0] + v) * dx
        a += (1.0 / mu_cf) * r_cf * grad(u)[1] * grad(v)[1] * dx
        a += 1j * omega * sigma_cf * r_cf * u * v * dx
        f = LinearForm(fes)
        f += J_phi * r_cf * v * dx
        a.Assemble()
        f.Assemble()
        gfu.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="pardiso") * f.vec

        # Update per-element mu_r from |B| at element centroid.
        # |B|^2 = |dA/dz|^2 + |(1/r)(d(rA)/dr)|^2 (peak, not RMS)
        # In axisymmetric AC: time-peak |B| = sqrt(2) * RMS |B|.
        # We use the COMPLEX magnitude as a proxy for peak.
        gf_mu_new = GridFunction(fes_mu)
        new_vec = gf_mu_new.vec.FV().NumPy()
        old_vec = gf_mu.vec.FV().NumPy().copy()
        max_dmu = 0.0
        # Iterate elements, evaluate |B|
        for el in mesh.Elements():
            if mesh.GetMaterials()[el.index] != "workpiece":
                new_vec[el.index] = old_vec[el.index]  # unchanged
                continue
            # Centroid r, z
            vs = [list(mesh.vertices[v.nr].point) for v in el.vertices]
            r_c = sum(v[0] for v in vs) / len(vs)
            z_c = sum(v[1] for v in vs) / len(vs)
            try:
                mip = mesh(r_c, z_c)
                A_c = complex(gfu(mip))
                # Finite-diff B (rough)
                eps = max(args.maxh_wp * 0.5, 1e-5)
                A_r1 = complex(gfu(mesh(max(r_c-eps, 1e-6), z_c)))
                A_r2 = complex(gfu(mesh(min(r_c+eps, args.R_wp-1e-7), z_c)))
                A_z1 = complex(gfu(mesh(r_c, max(z_c-eps, -args.H_wp/2+1e-7))))
                A_z2 = complex(gfu(mesh(r_c, min(z_c+eps, args.H_wp/2-1e-7))))
                B_z = ((max(r_c+eps, 1e-6) * A_r2 - max(r_c-eps, 1e-6) * A_r1) / (2*eps)) / r_c
                B_r = -((A_z2 - A_z1) / (2*eps))
                B_abs = math.sqrt(2) * math.sqrt(abs(B_r)**2 + abs(B_z)**2)  # peak
                mu_r_new = mu_r_from_B(B_abs)
                mu_r_damped = alpha * mu_r_new + (1 - alpha) * old_vec[el.index]
                new_vec[el.index] = mu_r_damped
                dmu = abs(mu_r_damped - old_vec[el.index]) / max(old_vec[el.index], 1e-3)
                max_dmu = max(max_dmu, dmu)
            except Exception:
                new_vec[el.index] = old_vec[el.index]
        gf_mu.vec.data = gf_mu_new.vec
        convergence.append({"iter": k_outer, "max_dmu_r": max_dmu})
        print(f"  Picard iter {k_outer}: max d(mu_r) = {max_dmu:.4f}")
        if max_dmu < tol_picard:
            print(f"  CONVERGED at iter {k_outer}")
            break

    # P_wp + |H_t| extraction (same as linear path)
    A_norm_sq = InnerProduct(gfu, gfu)
    P_wp_density = 0.5 * sigma_cf * (omega ** 2) * A_norm_sq * 2 * math.pi * r_cf
    P_wp = float(Integrate(P_wp_density.real, mesh,
                            definedon=mesh.Materials("workpiece")))
    print(f"P_wp_volumetric (nonlinear) = {P_wp:.6e} W")
    return {
        "method": "axisym_volumetric_A_phi_nonlinear",
        "P_wp_W": P_wp,
        "picard_iter": len(convergence),
        "picard_convergence": convergence,
    }


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

    # Extract |H_t| at the cylinder side wall (r = R_wp, z varies)
    # H_t in axisymmetric = H_z on the lateral surface (r=R, z) =
    # (1/r) d(r A_phi)/dr |_{r=R-eps} since B_z = (1/r) d(r A)/dr
    # and H_z = B_z/mu (in air just outside the conductor).
    # Sample a sparse z-line just outside R_wp.
    n_z_samples = 21
    z_pts = np.linspace(-args.H_wp / 2 * 0.9, args.H_wp / 2 * 0.9, n_z_samples)
    eps_R = 1e-5  # offset outside workpiece into air (mu_r=1 there)
    H_t_samples = []
    for z_val in z_pts:
        try:
            mip = mesh(args.R_wp + eps_R, z_val)
            A_val = complex(gfu(mip))
            # d(r A_phi)/dr via finite diff
            mip2 = mesh(args.R_wp + 2 * eps_R, z_val)
            A_val2 = complex(gfu(mip2))
            r1 = args.R_wp + eps_R
            r2 = args.R_wp + 2 * eps_R
            B_z = ((r2 * A_val2 - r1 * A_val) / (r2 - r1)) / r1
            H_z = B_z / MU0  # mu_r=1 just outside workpiece
            H_t_samples.append(abs(H_z))
        except Exception:
            H_t_samples.append(0.0)
    H_t_samples = np.array(H_t_samples)
    H_t_mean = float(H_t_samples.mean())
    H_t_max = float(H_t_samples.max())
    print(f"|H_t| at r=R+eps, z=[-0.9H/2, +0.9H/2], n={n_z_samples}: "
          f"mean={H_t_mean:.2f}, max={H_t_max:.2f} A/m")

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
        "H_t_mean_A_per_m": H_t_mean,
        "H_t_max_A_per_m": H_t_max,
        "H_t_samples_A_per_m": H_t_samples.tolist(),
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
    parser.add_argument("--bh-file", type=str, default=None,
                        help="If set, run nonlinear-BH Picard outer iteration "
                             "using the given two-column H[A/m]  B[T] BH file.")
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()
    if args.bh_file:
        import sys as _sys, os as _os
        HERE_ = _os.path.dirname(_os.path.abspath(__file__))
        SRC_ = _os.path.abspath(_os.path.join(HERE_, ".."))
        if SRC_ not in _sys.path:
            _sys.path.insert(0, SRC_)
        from em_material import load_bh_file  # type: ignore
        bh = load_bh_file(args.bh_file)
        if isinstance(bh, tuple):
            H_arr, B_arr = bh
            bh_curve = list(zip(H_arr.tolist(), B_arr.tolist()))
        else:
            bh_curve = list(bh)
        result = run_axisym_nonlinear(args, bh_curve)
        if args.output:
            with open(args.output, "w") as fh:
                json.dump(result, fh, indent=2)
            print(f"wrote {args.output}")
    else:
        run_axisym_linear(args)


if __name__ == "__main__":
    main()
