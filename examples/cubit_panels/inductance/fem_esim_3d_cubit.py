"""
3D FEM-ESIM with Cubit high-order mesh: induction heating with Karl iteration.

Cubit -> export_NGSolveCurvedMesh -> HCurl FEM + SIBC Robin BC.

Pipeline:
  1. Load Cubit model (.cub5) or auto-create Go-Tech toymodel
  2. export_NGSolveCurvedMesh(order=2): high-order curved tet mesh
  3. HCurl FEM: curl(nu*curl(A)) = J0 in coil, nu_K in Kelvin exterior domain
  4. SIBC Robin BC: -(jw/Z_s) * A_t . v_t  on wp_surface
  5. Karl iteration: solve -> sample H_t -> update Z_s -> repeat
  6. Power: P = sum(P'(H_t_i) * dA_i) via ESIM per panel

Blocks (from Cubit):
  coil        - volume, J = J0 * e_theta
  air         - volume, nu = nu0
  kelvin      - volume, nu = nu0 * (r / R_air)^4
  wp_surface  - surface, Robin BC
  outer       - surface, Dirichlet A = 0

Usage:
    python fem_esim_3d_cubit.py                              # auto toymodel
    python fem_esim_3d_cubit.py --cub5 induction_model.cub5  # existing model
    python fem_esim_3d_cubit.py --freq 7000 --material steel
    python fem_esim_3d_cubit.py --order 3                    # cubic elements
"""

import argparse
import math
import os
import sys
import time
import numpy as np

# Pre-import scipy BEFORE Cubit (Cubit bundles a broken scipy that shadows system)
import scipy  # noqa: F401
import scipy.interpolate  # noqa: F401

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4e-7 * np.pi
NU_0 = 1.0 / MU_0

STEEL_BH = [
    [0, 0], [50, 0.1], [100, 0.25], [200, 0.55],
    [500, 0.95], [1000, 1.2], [2000, 1.4], [5000, 1.55],
    [10000, 1.65], [20000, 1.75], [50000, 1.9], [100000, 2.0],
]


# ============================================================
# Cubit model loading + mesh export
# ============================================================
def _init_cubit():
    """Import and initialize Cubit (batch mode).

    Cubit's bin/ directory is appended (not prepended) to sys.path to
    avoid its bundled scipy/numpy from shadowing the system packages.
    """
    cubit_path = os.environ.get("CUBIT_PATH")
    if cubit_path and cubit_path not in sys.path:
        sys.path.append(cubit_path)
    # Also append Cubit's internal Python site-packages at the END
    # so that system scipy/numpy take priority
    cubit_pylib = os.path.join(cubit_path, "python3", "lib", "site-packages") \
        if cubit_path else None
    if cubit_pylib and cubit_pylib in sys.path:
        sys.path.remove(cubit_pylib)
        sys.path.append(cubit_pylib)
    import cubit
    try:
        cubit.init(['cubit', '-nojournal', '-batch'])
    except:
        pass  # Already initialized
    return cubit


def load_cubit_mesh(cub5_path, order=2):
    """Load Cubit .cub5 and export NGSolve curved mesh.

    Returns:
        ngsolve.Mesh (with CallbackGeometry curving)
    """
    # NGSolve MUST be imported BEFORE Cubit (DLL conflict on Windows)
    import ngsolve  # noqa: F401 - ensure DLL loaded first
    cubit = _init_cubit()
    cubit.cmd(f'open "{cub5_path}"')

    from cubit_mesh_export import extract_curved_mesh
    from ngsolve import Mesh
    return Mesh(extract_curved_mesh(cubit, order=order))


def create_and_load(order=2, mesh_size=0.008):
    """Auto-create toymodel in Cubit and export mesh.

    Returns:
        (ngsolve.Mesh, dict with geometry params)
    """
    import ngsolve  # noqa: F401
    cubit = _init_cubit()

    from create_induction_model import create_induction_model
    info = create_induction_model(cubit, mesh_size=mesh_size)

    from cubit_mesh_export import extract_curved_mesh
    from ngsolve import Mesh
    mesh = Mesh(extract_curved_mesh(cubit, order=order))
    return mesh, info


# ============================================================
# Solver
# ============================================================
def solve(mesh, R_coil=0.030, a_coil=0.005,
          R_wp=0.015, H_wp=0.040, R_air=0.150,
          sigma=2e6, frequency=7000, material="steel",
          I_total=1.0, order=2,
          max_iter=15, tol=1e-3, relax=0.5,
          solver="auto"):
    """3D FEM-ESIM with Karl iteration.

    Args:
        mesh: NGSolve Mesh with materials {coil, air, kelvin} and
              boundaries {wp_surface, outer}
        R_coil, a_coil: Coil geometry [m]
        R_wp, H_wp: Workpiece geometry [m]
        R_air: Kelvin boundary radius [m]
        sigma: Workpiece conductivity [S/m]
        frequency: Operating frequency [Hz]
        material: "steel", "copper", or "aluminum"
        I_total: Coil current [A]
        order: FE polynomial order
        max_iter: Max Karl iterations
        tol: Z_s convergence tolerance
        relax: Under-relaxation factor (0.5 = 50% new + 50% old)
        solver: "pardiso" (direct), "bddc" (iterative), or "auto"

    Returns:
        dict with P_total, Q_total, L, Karl convergence, panel data
    """
    from ngsolve import (HCurl, BilinearForm, LinearForm, GridFunction,
                         Integrate, curl, dx, ds, CF, TaskManager,
                         Conj, VOL, BND, Preconditioner)
    from ngsolve.krylovspace import CGSolver, GMResSolver
    from ngsolve import x, y, z, sqrt, IfPos
    from esim_cell_problem import ESIMFiniteSlabSolver

    omega = 2 * np.pi * frequency
    J0 = I_total / (np.pi * a_coil**2)

    bh_curve = STEEL_BH if material == "steel" else None
    mu_r_lin = 1.0 if material in ("copper", "aluminum") else None

    rho = 1.0 / sigma
    mu_r_eff = 1000 if material == "steel" else 1.0
    delta = np.sqrt(2 * rho / (omega * MU_0 * mu_r_eff))

    print("=" * 65)
    print("3D FEM-ESIM (Cubit mesh) + Karl Iteration")
    print("=" * 65)
    print(f"Coil:      R={R_coil*1e3:.0f}mm, a={a_coil*1e3:.1f}mm")
    print(f"Workpiece: R={R_wp*1e3:.0f}mm, H={H_wp*1e3:.0f}mm (SIBC hole)")
    print(f"Material:  {material}, sigma={sigma:.0e} S/m, f={frequency} Hz")
    print(f"Skin depth: delta={delta*1e3:.3f}mm, xi=R/delta={R_wp/delta:.1f}")
    print()

    materials = mesh.GetMaterials()
    boundaries = mesh.GetBoundaries()
    ne = mesh.GetNE(VOL)
    print(f"Mesh: {ne} elements, order={order}")
    print(f"  Materials:  {materials}")
    print(f"  Boundaries: {boundaries}")

    # --- Coefficient functions ---
    # Kelvin weight: ALL material constants scale by (a/r')^2
    # (conformal symmetry of Maxwell's equations, Sugahara IEEE Trans. Magn. 2022)
    # nu_kelvin = nu0 * (R_air/r)^2 where r = distance from origin
    r_sq = x*x + y*y + z*z
    kelvin_fac = R_air**2 / (r_sq + 1e-20)  # (a/r)^2

    nu_dict = {}
    for m in materials:
        if "kelvin" in m.lower():
            nu_dict[m] = NU_0 * kelvin_fac
        else:
            nu_dict[m] = NU_0
    nu_cf = mesh.MaterialCF(nu_dict, default=NU_0)

    # Source current: J = J0 * e_theta (tangential to torus)
    r_xy = sqrt(x*x + y*y)
    r_xy_safe = IfPos(r_xy - 1e-10, r_xy, 1e-10)
    J_source = J0 * CF((-y / r_xy_safe, x / r_xy_safe, 0))

    # --- ESIM solver ---
    # When bh_curve is provided (steel), mu_r must be None so ESIM uses
    # the BH curve for nonlinear iteration. Passing mu_r=1.0 with bh_curve
    # causes ESIM to ignore the BH curve and use linear mu_r=1.
    esim_solver = ESIMFiniteSlabSolver(
        half_thickness=R_wp, bh_curve=bh_curve, sigma=sigma,
        frequency=frequency,
        mu_r=mu_r_lin if bh_curve is None else None, n_nodes=200,
        geometry='cylinder')

    # Initial Z_s
    sol0 = esim_solver.solve(5.0)
    Z_s = sol0['Z']
    print(f"\nInitial Z_s = {Z_s:.4e}")
    print(f"  delta (ESIM) = {sol0['delta']*1e3:.3f} mm")

    # --- Sample points on workpiece surface ---
    n_phi, n_z = 24, 12
    n_r_cap = 6
    eps = 0.001  # 1mm offset from surface

    def build_sample_panels():
        """Build sample panels on workpiece surface."""
        panels = []
        # Lateral surface (cylinder at r = R_wp)
        for i in range(n_phi):
            phi = (i + 0.5) * 2 * np.pi / n_phi
            cos_p, sin_p = np.cos(phi), np.sin(phi)
            r_eval = R_wp + eps
            for j in range(n_z):
                z_j = -H_wp / 2 + (j + 0.5) * H_wp / n_z
                px, py, pz = r_eval * cos_p, r_eval * sin_p, z_j
                dA = R_wp * (2 * np.pi / n_phi) * (H_wp / n_z)
                nr = np.array([cos_p, sin_p, 0.0])
                panels.append({
                    'type': 'lateral', 'pt': (px, py, pz),
                    'normal': nr, 'dA': dA, 'z': z_j, 'phi': phi,
                })
        # Top/bottom caps
        for sign in [+1, -1]:
            for i in range(n_phi):
                phi = (i + 0.5) * 2 * np.pi / n_phi
                cos_p, sin_p = np.cos(phi), np.sin(phi)
                for k in range(n_r_cap):
                    r_k = (k + 0.5) * R_wp / n_r_cap
                    if r_k < 1e-6:
                        continue
                    px = r_k * cos_p
                    py = r_k * sin_p
                    pz = sign * (H_wp / 2 + eps)
                    dA = r_k * (2 * np.pi / n_phi) * (R_wp / n_r_cap)
                    nr = np.array([0.0, 0.0, float(sign)])
                    panels.append({
                        'type': 'top' if sign > 0 else 'bottom',
                        'pt': (px, py, pz), 'normal': nr,
                        'dA': dA, 'r': r_k, 'phi': phi,
                    })
        return panels

    sample_panels = build_sample_panels()
    n_panels = len(sample_panels)
    print(f"Sample panels: {n_panels} ({n_phi}x{n_z} lateral + "
          f"2x{n_phi}x{n_r_cap} caps)")

    # --- FE space ---
    fes = HCurl(mesh, order=order, complex=True, dirichlet="outer")
    u, v = fes.TnT()
    print(f"DOFs: {fes.ndof}")

    # Select solver: direct (pardiso) for small, iterative (bddc) for large
    if solver == "auto":
        solver = "pardiso" if fes.ndof < 200000 else "bddc"
    print(f"Solver: {solver}")

    # --- Karl iteration ---
    print()
    print(f"{'Iter':>4s} {'|Z_s|':>12s} {'Re(Z_s)':>12s} "
          f"{'H_t_rms':>10s} {'dZ/Z':>10s} {'t_solve':>8s}")
    print("-" * 70)

    t0_total = time.perf_counter()
    gfu = None
    history = []

    for iteration in range(max_iter):
        t0_iter = time.perf_counter()

        # Robin BC coefficient
        robin = -1j * omega / Z_s

        # Assemble (gauge regularization for curl-curl kernel stability)
        # bonus_intorder=4 required for Kelvin domain (spatially varying nu)
        a_bf = BilinearForm(fes)
        a_bf += nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)
        a_bf += 1e-6 * NU_0 * u * v * dx  # gauge regularization
        a_bf += robin * u.Trace() * v.Trace() * ds("wp_surface")

        if solver != "pardiso":
            pre = Preconditioner(a_bf, "bddc")  # must register BEFORE Assemble

        a_bf.Assemble()

        f_lf = LinearForm(fes)
        f_lf += J_source * v * dx("coil")
        f_lf.Assemble()

        # Solve
        gfu = GridFunction(fes)
        with TaskManager():
            if solver == "pardiso":
                gfu.vec.data = a_bf.mat.Inverse(
                    fes.FreeDofs(), inverse="pardiso") * f_lf.vec
            else:
                inv = GMResSolver(a_bf.mat, pre.mat,
                                  maxiter=2000, tol=1e-8,
                                  printrates=False)
                gfu.vec.data = inv * f_lf.vec
        t_solve = time.perf_counter() - t0_iter

        # H_t from SIBC relation on boundary (NOT offset sampling)
        # SIBC: |H_t| = w*|A_t|/|Z_s|, so H_t_rms from int|A_t|^2 dS:
        #   H_rms^2 = w^2/|Z_s|^2 * int|A_t|^2 dS / A_wp
        int_A2 = Integrate(
            gfu * Conj(gfu), mesh, BND,
            definedon=mesh.Boundaries("wp_surface")).real
        wp_area = 2 * np.pi * R_wp * H_wp + 2 * np.pi * R_wp**2  # lat + caps
        H_t_rms = omega / abs(Z_s) * np.sqrt(
            max(int_A2 / wp_area, 0)) if int_A2 > 0 else 1e-3

        # Update Z_s (under-relaxation)
        Z_s_old = Z_s
        sol_new = esim_solver.solve(max(float(H_t_rms), 1e-3))
        Z_s = relax * sol_new['Z'] + (1 - relax) * Z_s_old

        dZ = abs(Z_s - Z_s_old) / max(abs(Z_s_old), 1e-30)
        history.append({
            'iteration': iteration,
            'Z_s': complex(Z_s), 'H_t_rms': float(H_t_rms),
            'dZ': float(dZ), 't_solve': t_solve,
        })
        print(f"{iteration:4d} {abs(Z_s):12.4e} {Z_s.real:12.4e} "
              f"{H_t_rms:10.2f} {dZ:10.4e} {t_solve:8.1f}s")

        if dZ < tol and iteration > 0:
            print(f"  Converged at iteration {iteration}")
            break

    t_total = time.perf_counter() - t0_total

    # --- Post-process: P from Poynting flux (self-consistent) ---
    # P = 0.5 * w^2 * Re(Z_s)/|Z_s|^2 * int_S |A_t|^2 dS
    from ngsolve import BND
    int_A2_final = Integrate(
        gfu * Conj(gfu), mesh, BND,
        definedon=mesh.Boundaries("wp_surface")).real
    P_total = 0.5 * omega**2 * Z_s.real / abs(Z_s)**2 * int_A2_final
    # Q from imaginary part of Z_s
    Q_total = 0.5 * omega**2 * Z_s.imag / abs(Z_s)**2 * int_A2_final
    panel_results = []  # per-panel not available from boundary integral

    # Energy method for L (high integration order for Kelvin nu)
    W_mag = Integrate(
        0.5 * nu_cf * curl(gfu) * Conj(curl(gfu)), mesh,
        order=10).real
    L = 2 * W_mag / I_total**2

    # --- Summary ---
    print()
    print("=" * 65)
    print("Results: FEM-ESIM with Robin BC (Poynting flux)")
    print("=" * 65)
    print(f"  P (Poynting)     = {P_total:.6e} W")
    print(f"  Q                = {Q_total:.6e} var")
    print(f"  L (energy)       = {L*1e9:.2f} nH")
    print(f"  Z_s (converged)  = {Z_s:.4e}")
    print(f"  H_t rms (SIBC)   = {H_t_rms:.2f} A/m")
    print(f"  Karl iterations  = {len(history)}")
    print(f"  DOFs             = {fes.ndof}")
    print(f"  Total time       = {t_total:.1f} s")
    print("=" * 65)

    return {
        'P_total': float(P_total),
        'Q_total': float(Q_total),
        'L': float(L),
        'Z_s': complex(Z_s),
        'H_t_rms': float(H_t_rms),
        'ndof': fes.ndof,
        'ne': ne,
        'iterations': len(history),
        'history': history,
        'panels': panel_results,
        't_total': t_total,
        'delta': float(delta),
        'frequency': frequency,
        'sigma': sigma,
        'material': material,
    }


# ============================================================
# Main
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description="3D FEM-ESIM with Cubit high-order mesh")
    parser.add_argument("--cub5", default="", help="Cubit .cub5 file")
    parser.add_argument("--toymodel", action="store_true",
                        help="Auto-create Go-Tech toymodel")
    parser.add_argument("--freq", type=float, default=7000,
                        help="Frequency [Hz] (default: 7000)")
    parser.add_argument("--sigma", type=float, default=2e6,
                        help="Conductivity [S/m] (default: 2e6)")
    parser.add_argument("--material", default="steel",
                        choices=["steel", "copper", "aluminum"])
    parser.add_argument("--order", type=int, default=1,
                        help="FE polynomial order (default: 1)")
    parser.add_argument("--mesh-size", type=float, default=0.008,
                        help="Mesh size [m] (toymodel only)")
    parser.add_argument("--max-iter", type=int, default=15,
                        help="Max Karl iterations")
    parser.add_argument("--current", type=float, default=1.0,
                        help="Coil current [A] (default: 1.0)")
    parser.add_argument("--solver", default="auto",
                        choices=["pardiso", "bddc", "auto"],
                        help="Linear solver (default: auto)")
    args = parser.parse_args()

    # Geometry parameters (Go-Tech toymodel)
    R_coil, a_coil = 0.030, 0.005
    R_wp, H_wp = 0.015, 0.040
    R_air = 0.150

    if args.cub5:
        print(f"Loading Cubit model: {args.cub5}")
        mesh = load_cubit_mesh(args.cub5, order=args.order)
    else:
        print("Creating toymodel in Cubit...")
        mesh, info = create_and_load(order=args.order,
                                     mesh_size=args.mesh_size)
        R_coil = info.get('R_coil', R_coil)
        a_coil = info.get('a_coil', a_coil)
        R_wp = info.get('R_wp', R_wp)
        H_wp = info.get('H_wp', H_wp)
        R_air = info.get('R_air', R_air)

    result = solve(
        mesh, R_coil=R_coil, a_coil=a_coil,
        R_wp=R_wp, H_wp=H_wp, R_air=R_air,
        sigma=args.sigma, frequency=args.freq,
        material=args.material, I_total=args.current,
        order=args.order, max_iter=args.max_iter,
        solver=args.solver)

    # Comparison with existing BEM-ESIM (if available)
    try:
        print("\nRunning BEM-ESIM for comparison...")
        from impedance_esim import run as bem_run
        r_bem = bem_run(
            material=args.material, frequency=args.freq,
            R=R_coil, a=a_coil,
            wp_radius=R_wp, wp_height=H_wp,
            n_phi=24, n_z=12, verbose=False)
        print()
        print("Comparison: FEM-ESIM (Cubit 3D) vs BEM-ESIM")
        print("-" * 55)
        P_f, P_b = result['P_total'], r_bem['P_total']
        L_f, L_b = result['L'], r_bem['L_coil']
        print(f"{'':>15s} {'FEM-ESIM':>14s} {'BEM-ESIM':>14s} {'diff':>8s}")
        print(f"{'P [W]':>15s} {P_f:14.4e} {P_b:14.4e} "
              f"{(P_f-P_b)/max(abs(P_b),1e-30)*100:+8.1f}%")
        print(f"{'L [nH]':>15s} {L_f*1e9:14.2f} {L_b*1e9:14.2f} "
              f"{(L_f-L_b)/max(abs(L_b),1e-30)*100:+8.1f}%")
    except Exception as e:
        print(f"  (BEM-ESIM comparison skipped: {e})")


if __name__ == "__main__":
    main()
