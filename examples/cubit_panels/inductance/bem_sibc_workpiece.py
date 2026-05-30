"""
BEM-SIBC workpiece eddy current: coil (Biot-Savart) + cylindrical workpiece (scalar BIE).

Full pipeline:
  1. Coil: filamentary current loop -> Biot-Savart -> H_inc
  2. phi_inc at workpiece surface via path integration
  3. Scalar BIE + SIBC on workpiece surface -> J_s, P_loss, H_t_rms
  4. Compare with FEM-SIBC (fem_esim_3d.py)

No C++ code needed - uses only existing ngsolve.bem operators.

Usage:
    python bem_sibc_workpiece.py
    python bem_sibc_workpiece.py --freq 7000 --material steel
    python bem_sibc_workpiece.py --freq 1000 --material copper
"""

import math
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4e-7 * np.pi


def run(R_coil=0.030, a_coil=0.003, gap_deg=5,
        R_wp=0.010, H_wp=0.020,
        sigma=2e6, frequency=7000, material="steel",
        maxh_factor=4, h1_order=1):
    """BEM-SIBC workpiece solver.

    Args:
        R_coil, a_coil: Coil major/minor radius [m]
        gap_deg: Gap angle [degrees]
        R_wp, H_wp: Workpiece cylinder radius/height [m]
        sigma: Conductivity [S/m]
        frequency: [Hz]
        material: steel/copper/aluminum
        maxh_factor: mesh density (maxh = R_wp / factor)
        h1_order: H1 polynomial order
    """
    from ngsolve import (Mesh, Integrate, CF, ds, BND, TaskManager)
    from netgen.occ import (Cylinder, Pnt, Dir, OCCGeometry, Glue)
    from radia.bem_sibc_solver import (ScalarBIESIBCSolver, compute_phi_inc_from_loop)
    from esim_cell_problem import ESIMFiniteSlabSolver

    omega = 2 * np.pi * frequency
    I_total = 1.0  # 1 A total current (same as fem_esim_3d)

    # Material properties
    STEEL_BH = [
        [0, 0], [50, 0.1], [100, 0.25], [200, 0.55],
        [500, 0.95], [1000, 1.2], [2000, 1.4], [5000, 1.55],
        [10000, 1.65], [20000, 1.75], [50000, 1.9], [100000, 2.0],
    ]
    bh_curve = STEEL_BH if material == "steel" else None
    mu_r = 1.0 if material in ("copper", "aluminum") else None

    # Skin depth estimate
    mu_eff = MU_0 * (mu_r if mu_r else 100.0)
    delta = math.sqrt(2.0 / (omega * mu_eff * sigma))

    print("=" * 65)
    print("BEM-SIBC: Coil (Biot-Savart) + Workpiece (Scalar BIE)")
    print("=" * 65)
    print(f"Coil:      R={R_coil*1e3:.0f}mm, a={a_coil*1e3:.1f}mm, gap={gap_deg}deg")
    print(f"Workpiece: R={R_wp*1e3:.0f}mm, H={H_wp*1e3:.0f}mm")
    print(f"Material:  {material}, sigma={sigma:.0e}, f={frequency}Hz")
    print(f"Skin depth: delta={delta*1e3:.2f}mm")
    print()

    # === 1. Workpiece surface mesh ===
    print("[1/4] Meshing workpiece surface...")
    t0 = time.perf_counter()

    wp_cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)
    for f in wp_cyl.faces:
        f.name = "wp_surface"

    geo = OCCGeometry(Glue(wp_cyl.faces))
    maxh = R_wp / maxh_factor
    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)
    with TaskManager():
        mesh.Curve(3)

        ne = mesh.GetNE(BND)
        nv = mesh.nv
        area = Integrate(CF(1), mesh, BND)
        t_mesh = time.perf_counter() - t0
        print(f"  {ne} elements, {nv} vertices, area = {abs(area):.6e} m^2 ({t_mesh:.1f}s)")

        # === 2. Assemble BEM operators ===
        print("[2/4] Assembling BEM operators...")
        solver = ScalarBIESIBCSolver(mesh, order=h1_order)
        print(f"  {solver.ndof} DOFs, assembly: {solver.t_assembly:.1f}s")

        # === 3. Compute phi_inc from coil ===
        print("[3/4] Computing phi_inc (Biot-Savart + path integration)...")
        t0 = time.perf_counter()

        # Get mesh vertex coordinates
        node_coords = np.array([[mesh.vertices[i].point[j] for j in range(3)]
                                for i in range(mesh.nv)])

        phi_inc_nodes = compute_phi_inc_from_loop(
            node_coords, loop_center=[0, 0, 0], loop_radius=R_coil,
            current=I_total, n_quad=30, gap_deg=gap_deg)

        phi_inc_vec = phi_inc_nodes  # H1 order=1: vertex DOFs match

        t_phi = time.perf_counter() - t0
        print(f"  phi_inc range: [{phi_inc_nodes.min():.4f}, {phi_inc_nodes.max():.4f}] ({t_phi:.1f}s)")

        # === 4. ESIM surface impedance + Karl iteration ===
        print("[4/4] Solving scalar BIE + SIBC...")

        esim_solver = ESIMFiniteSlabSolver(
            half_thickness=R_wp, bh_curve=bh_curve, sigma=sigma,
            frequency=frequency,
            mu_r=mu_r if bh_curve is None else None,
            n_nodes=200, geometry='cylinder')

        # Initial Z_s
        H_t_init = 5.0
        Z_s = esim_solver.solve(H_t_init)['Z']

        max_iter = 15
        tol = 1e-3
        relax = 0.5

        print(f"  |Z_s| init = {abs(Z_s):.4e}")
        print(f"  {'Iter':>4s} {'|Z_s|':>12s} {'H_t_rms':>10s} "
              f"{'P [W]':>12s} {'dZ/Z':>10s}")

        for iteration in range(max_iter):
            result = solver.solve(phi_inc_vec, Z_s=Z_s, omega=omega)
            H_t_rms = result['H_t_rms']
            P_density = result['P_density']
            P_total = P_density * result['area']

            # ESIM update
            Z_s_old = Z_s
            sol = esim_solver.solve(max(H_t_rms, 1e-3))
            Z_s_new = sol['Z']
            Z_s = relax * Z_s_new + (1 - relax) * Z_s_old

            dZ = abs(Z_s - Z_s_old) / max(abs(Z_s_old), 1e-30)
            print(f"  {iteration:4d} {abs(Z_s):12.4e} {H_t_rms:10.2f} "
                  f"{P_total:12.4e} {dZ:10.4e}")

            if dZ < tol and iteration > 0:
                print(f"  Converged at iteration {iteration}")
                break

        # Final result with converged Z_s
        result = solver.solve(phi_inc_vec, Z_s=Z_s, omega=omega)
        H_t_rms = result['H_t_rms']
        P_density = result['P_density']
        P_total = P_density * result['area']
        Q_density = 0.5 * Z_s.imag * H_t_rms**2 if Z_s != 0 else 0
        Q_total = Q_density * result['area']

        print()
        print("=" * 65)
        print("Results: BEM-SIBC Workpiece")
        print("-" * 65)
        print(f"  P (BEM-SIBC)       = {P_total:.6e} W")
        print(f"  Q                  = {Q_total:.6e} var")
        print(f"  H_t_rms            = {H_t_rms:.2f} A/m")
        print(f"  |Z_s|              = {abs(Z_s):.4e} Ohm")
        print(f"  BEM DOFs: {solver.ndof}")
        print(f"  Time: mesh={t_mesh:.1f}s, assembly={solver.t_assembly:.1f}s, "
              f"phi_inc={t_phi:.1f}s, solve={result['t_solve']:.3f}s")
        print("-" * 65)

        return {
            'P_total': float(P_total),
            'Q_total': float(Q_total),
            'H_t_rms': float(H_t_rms),
            'Z_s': complex(Z_s),
            'ndof': solver.ndof,
            'area': result['area'],
        }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="BEM-SIBC workpiece solver")
    parser.add_argument("--freq", type=float, default=7000)
    parser.add_argument("--material", default="steel")
    parser.add_argument("--maxh", type=int, default=4,
                        help="Mesh density: maxh = R_wp / this value")
    args = parser.parse_args()

    sigma_map = {'steel': 2e6, 'copper': 5.8e7, 'aluminum': 3.5e7}
    sigma = sigma_map.get(args.material, 2e6)

    r_bem = run(frequency=args.freq, material=args.material, sigma=sigma,
                maxh_factor=args.maxh)

    # Compare with FEM-SIBC
    print()
    print("Running FEM-SIBC for comparison...")
    from fem_esim_3d import run as fem_run
    r_fem = fem_run(frequency=args.freq, material=args.material, sigma=sigma,
                    approach='hole')

    if r_fem:
        print()
        print("=" * 65)
        print(f"Comparison: BEM-SIBC vs FEM-SIBC ({args.material}, {args.freq}Hz)")
        print("=" * 65)
        print(f"{'':>20s} {'BEM-SIBC':>14s} {'FEM-SIBC':>14s} {'diff':>8s}")
        P_b, P_f = r_bem['P_total'], r_fem['P_total']
        if abs(P_f) > 1e-30:
            print(f"{'P [W]':>20s} {P_b:14.4e} {P_f:14.4e} "
                  f"{(P_b-P_f)/P_f*100:+8.1f}%")
        print(f"{'H_t_rms [A/m]':>20s} {r_bem['H_t_rms']:14.2f} "
              f"{r_fem['H_t_rms']:14.2f}")
        print(f"{'|Z_s| [Ohm]':>20s} {abs(r_bem['Z_s']):14.4e} "
              f"{abs(r_fem['Z_s']):14.4e}")
        print(f"{'DOFs':>20s} {r_bem['ndof']:14d} {r_fem['ndof']:14d}")
