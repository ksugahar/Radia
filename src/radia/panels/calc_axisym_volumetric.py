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
``radia.axifem.{AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI}``.

Strong form:

    1/mu * (-d_r^2 A_phi + 1/r d_r A_phi - 1/r^2 A_phi - d_z^2 A_phi)
        + j w sigma A_phi = J_phi_source

where ``J_phi_source`` is the impressed coil current density (real, A/m^2)
and ``sigma > 0`` only inside the workpiece (eddy current induced).

Weak form (in flux variable psi, after Henrotte change of variable):

    int (1/mu) grad_s_z psi . grad_s_z psi' / (2 pi r)
      + j w sigma psi psi' / (2 pi r)
      = int J_phi_source psi' . (r-weighted)

(Closed-form per-element by axifem BFIs.)

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
    from ngsolve import Mesh, TaskManager
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
      1. Solve AC linear problem with current mu_r distribution.
      2. Build a |B| CoefficientFunction from grad(gfu) (no finite diff).
      3. Project |B| onto L2(order=0) -> per-element constant.
      4. Update mu_r per element via BH inverse: mu_r = B / (mu_0 H(B))
         with mu_r floor at 1.0 (vacuum minimum).
      5. Damped Picard: mu_r_new = alpha mu_BH + (1-alpha) mu_r_old.
      6. Stop when max relative change of mu_r drops below tol.

    For axisymmetric A_phi:
        B_z = (1/r) d(r A)/dr = A/r + dA/dr
        B_r = -dA/dz
        |B|^2 = |B_r|^2 + |B_z|^2  (complex; peak = sqrt(2) * RMS)
    """
    from ngsolve import (
        BilinearForm, LinearForm, GridFunction, Integrate, dx, grad,
        Conj, sqrt as ng_sqrt, x as r_cf, H1, L2,
    )
    import radia.axifem   # noqa: F401

    mesh = build_mesh(
        R_wp_m=args.R_wp, H_wp_m=args.H_wp,
        R_coil_m=args.R_coil, R_outer_m=args.R_outer,
        maxh_wp_m=args.maxh_wp, maxh_air_m=args.maxh_air,
    )
    print(f"mesh: ne={mesh.ne}, nv={mesh.nv}, mats={mesh.GetMaterials()}")

    # BH lookup -- given |B|, return mu_r = B / (mu_0 H).
    bh = np.asarray(bh_curve)
    H_arr_bh, B_arr_bh = bh[:, 0], bh[:, 1]
    # Initial-slope mu_r for low |B| (B ~ 0 floor).
    mu_r_init = float(args.mu_r)

    def mu_r_from_B_array(B_arr_input):
        """Vectorised: given per-element |B|, return per-element mu_r."""
        out = np.full_like(B_arr_input, mu_r_init, dtype=float)
        # For B above the first BH-curve point: invert BH numerically.
        mask = B_arr_input > B_arr_bh[1] if len(B_arr_bh) > 1 else B_arr_input > 1e-6
        if mask.any():
            H_vals = np.interp(B_arr_input[mask], B_arr_bh, H_arr_bh)
            mu_r = B_arr_input[mask] / (MU0 * np.maximum(H_vals, 1e-9))
            out[mask] = mu_r
        # Floor at 1 (vacuum).
        return np.maximum(out, 1.0)

    # Per-element mu_r via L2(order=0) (true piecewise-constant per element).
    fes_mu = L2(mesh, order=0)
    gf_mu = GridFunction(fes_mu)
    init_mu_r_cf = mesh.MaterialCF({"workpiece": float(args.mu_r),
                                     "coil": 1.0, "air": 1.0}, default=1.0)
    with TaskManager():
        gf_mu.Set(init_mu_r_cf)

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
        max_picard = 25
        tol_picard = 1e-2
        alpha = 0.5
        convergence = []
        old_vec = gf_mu.vec.FV().NumPy().copy()
        # Map workpiece elements: index list
        wp_elem_idx = np.array(
            [i for i, el in enumerate(mesh.Elements())
             if mesh.GetMaterials()[el.index] == "workpiece"],
            dtype=int,
        )
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

            # |B| as a CoefficientFunction from grad(gfu).
            # B_z = grad[0] + A/r,  B_r = -grad[1]
            # |B|^2 = Re(B_z conj(B_z)) + Re(B_r conj(B_r)) is the COMPLEX
            # magnitude squared; peak |B| = sqrt(2) * sqrt(|B|^2).
            B_z_cf = grad(gfu)[0] + gfu / r_cf
            B_r_cf = -grad(gfu)[1]
            B_abs_sq_cf = (B_z_cf * Conj(B_z_cf) + B_r_cf * Conj(B_r_cf)).real
            # peak |B| = sqrt(2) * sqrt(|B|_complex^2)
            B_peak_cf = ng_sqrt(2.0 * B_abs_sq_cf)

            # Project |B| onto L2(order=0) -> per-element constant.
            gf_B = GridFunction(fes_mu)
            gf_B.Set(B_peak_cf)
            B_per_elem = gf_B.vec.FV().NumPy()

            # Update mu_r only on workpiece elements.
            new_vec = old_vec.copy()
            if len(wp_elem_idx) > 0:
                B_wp = B_per_elem[wp_elem_idx]
                mu_r_bh = mu_r_from_B_array(B_wp)
                mu_r_damped = alpha * mu_r_bh + (1.0 - alpha) * old_vec[wp_elem_idx]
                mu_r_damped = np.maximum(mu_r_damped, 1.0)
                new_vec[wp_elem_idx] = mu_r_damped

            max_dmu = float(np.max(np.abs(new_vec - old_vec)
                                    / np.maximum(old_vec, 1.0)))
            gf_mu.vec.FV().NumPy()[:] = new_vec
            old_vec = new_vec.copy()
            convergence.append({"iter": k_outer, "max_dmu_r": max_dmu,
                                 "mu_r_wp_mean": float(new_vec[wp_elem_idx].mean()),
                                 "B_wp_max": float(B_per_elem[wp_elem_idx].max()),
                                 "B_wp_mean": float(B_per_elem[wp_elem_idx].mean())})
            print(f"  Picard iter {k_outer}: max d(mu_r)={max_dmu:.4f}, "
                  f"<mu_r_wp>={float(new_vec[wp_elem_idx].mean()):.1f}, "
                  f"|B|_max={float(B_per_elem[wp_elem_idx].max()):.3f} T")
            if max_dmu < tol_picard and k_outer > 0:
                print(f"  CONVERGED at iter {k_outer}")
                break

        # P_wp via volumetric integration.
        from ngsolve import InnerProduct
        A_norm_sq = InnerProduct(gfu, gfu)
        P_wp_density = 0.5 * sigma_cf * (omega ** 2) * A_norm_sq * 2 * math.pi * r_cf
        P_wp = float(Integrate(P_wp_density.real, mesh,
                                definedon=mesh.Materials("workpiece")))
        print(f"P_wp_volumetric (nonlinear) = {P_wp:.6e} W")

        # |H_t| at side wall (same as linear path).
        n_z = 21
        z_pts = np.linspace(-args.H_wp / 2 * 0.9, args.H_wp / 2 * 0.9, n_z)
        eps_R = 1e-5
        H_t_samples = []
        for z_val in z_pts:
            try:
                r1 = args.R_wp + eps_R
                r2 = args.R_wp + 2 * eps_R
                A1 = complex(gfu(mesh(r1, z_val)))
                A2 = complex(gfu(mesh(r2, z_val)))
                B_z = ((r2 * A2 - r1 * A1) / (r2 - r1)) / r1
                H_z = B_z / MU0  # vacuum just outside workpiece
                H_t_samples.append(abs(H_z))
            except Exception:
                H_t_samples.append(0.0)
        H_t_samples = np.array(H_t_samples)

        return {
            "method": "axisym_volumetric_A_phi_nonlinear",
            "frequency_Hz": args.frequency,
            "current_A": args.current,
            "R_wp_m": args.R_wp,
            "H_wp_m": args.H_wp,
            "R_coil_m": args.R_coil,
            "sigma_S_per_m": args.sigma,
            "mu_r_init": float(args.mu_r),
            "fes_order": p,
            "ndof": fes.ndof,
            "P_wp_W": P_wp,
            "H_t_mean_A_per_m": float(H_t_samples.mean()),
            "H_t_max_A_per_m": float(H_t_samples.max()),
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
        FESpace, TaskManager,
    )
    import radia.axifem   # registers axihenrotte FESpace

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
    # because axifem BFIs accept only real CFs.  Standard NGSolve
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

    # TaskManager-Only policy: wrap the heavy FE work (assembly + pardiso
    # solve + power integral) so it runs PARALLEL on the NGSolve threadpool,
    # like run_axisym_nonlinear (L178).  Without this the LINEAR solver ran
    # SERIALLY -- calc_main does not open a TaskManager context, and this
    # function's wrap was missing (run_axisym_nonlinear had it, the linear
    # twin did not).
    with TaskManager():
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
