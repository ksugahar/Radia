"""Test two-body coupled BEM: coil (EFIE) + workpiece (Scalar BIE + SIBC).

Compares:
  - Uncoupled BEM: EFIE coil alone + one-way BEM-SIBC workpiece
  - Coupled BEM: iterative coil-workpiece coupling
  - FEM-ESIM: 3D reference (coil + workpiece together)
"""
import math
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))
sys.path.insert(0, os.path.dirname(__file__))

MU_0 = 4e-7 * np.pi


def run():
    from ngsolve import Mesh, Integrate, CF, BND, TaskManager
    from netgen.occ import Cylinder, Pnt, Dir, OCCGeometry, Glue
    from impedance_esim import make_gapped_torus_mesh
    from bem_coupled_solver import CoupledBEMSolver

    R_coil, a_coil, gap_deg = 0.030, 0.003, 5
    R_wp, H_wp = 0.010, 0.020
    sigma = 5.8e7  # copper
    freq = 1000
    omega = 2 * math.pi * freq
    delta = math.sqrt(2.0 / (omega * MU_0 * sigma))
    Z_s = complex(1, 1) / (sigma * delta)

    print("=" * 65)
    print("Coupled BEM: Coil (EFIE) + Workpiece (Scalar BIE + SIBC)")
    print("=" * 65)
    print(f"Coil: R={R_coil*1e3:.0f}mm, a={a_coil*1e3:.1f}mm")
    print(f"Workpiece: R={R_wp*1e3:.0f}mm, H={H_wp*1e3:.0f}mm, copper")
    print(f"f={freq}Hz, delta={delta*1e3:.2f}mm, |Z_s|={abs(Z_s):.4e}")
    print()

    # --- Meshes ---
    print("[1/3] Building meshes...")
    mesh_coil = make_gapped_torus_mesh(R_coil, a_coil, gap_deg=gap_deg,
                                        maxh=a_coil)

    wp_cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)
    for f in wp_cyl.faces:
        f.name = "wp"
    mesh_wp = Mesh(OCCGeometry(Glue(wp_cyl.faces)).GenerateMesh(maxh=R_wp / 3))
    with TaskManager():
        mesh_wp.Curve(3)

        print(f"  Coil: {mesh_coil.GetNE(BND)} elems")
        print(f"  WP:   {mesh_wp.GetNE(BND)} elems")

        # --- Coupled solver ---
        print("[2/3] Assembling coupled solver...")
        t0 = time.perf_counter()
        solver = CoupledBEMSolver(mesh_coil, mesh_wp)
        t_asm = time.perf_counter() - t0
        print(f"  Coil: {solver.n_J} J DOFs, assembly: {solver.t_coil_assembly:.1f}s")
        print(f"  WP:   {solver.wp_solver.ndof} phi DOFs, assembly: {solver.wp_solver.t_assembly:.1f}s")
        print(f"  Total assembly: {t_asm:.1f}s")

        # --- Solve ---
        print("[3/3] Coupled iteration...")
        result = solver.solve(Z_s=Z_s, omega=omega, max_iter=10, tol=1e-3)

        print()
        print("=" * 65)
        print("Coupled BEM Results")
        print("-" * 65)
        print(f"  L_self     = {result['L_self']*1e9:.2f} nH")
        print(f"  Delta_L    = {result['Delta_L']*1e9:.2f} nH")
        print(f"  L_total    = {result['L_total']*1e9:.2f} nH")
        print(f"  P_total    = {result['P_total']:.4e} W")
        print(f"  H_t_rms    = {result['H_t_rms']:.2f} A/m")
        print(f"  Iterations = {result['iterations']}")
        print("-" * 65)

        # --- FEM reference ---
        print()
        print("Running FEM-ESIM for comparison...")
        from fem_esim_3d import run as fem_run
        r_fem = fem_run(frequency=freq, material="copper", sigma=sigma, approach='hole')

        if r_fem:
            print()
            print("=" * 65)
            print("Comparison: Coupled BEM vs FEM-ESIM")
            print("=" * 65)
            print(f"{'':>15s} {'Coupled BEM':>14s} {'FEM-ESIM':>14s} {'diff':>8s}")
            print("-" * 55)
            L_c = result['L_total']
            L_f = r_fem.get('L', 0)
            print(f"{'L [nH]':>15s} {L_c*1e9:14.2f} {L_f*1e9:14.2f} "
                  f"{(L_c/L_f-1)*100:+8.1f}%")
            P_c = result['P_total']
            P_f = r_fem['P_total']
            print(f"{'P [W]':>15s} {P_c:14.4e} {P_f:14.4e} "
                  f"{(P_c/P_f-1)*100:+8.1f}%")
            print(f"{'H_t_rms':>15s} {result['H_t_rms']:14.2f} "
                  f"{r_fem['H_t_rms']:14.2f}")


if __name__ == "__main__":
    run()
