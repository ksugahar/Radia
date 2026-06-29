"""Coupled BEM test with steel workpiece (nonlinear mu, Karl iteration).

Steel has mu_r >> 1 at low H, so Delta_L should be significant.
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
    from ngsolve import Mesh, Integrate, CF, BND, TaskManager, GridFunction, grad
    from netgen.occ import Cylinder, Pnt, Dir, OCCGeometry, Glue
    from impedance_esim import make_gapped_torus_mesh
    from bem_coupled_solver import (CoupledBEMSolver, extract_element_J,
                                     vector_potential_from_J)
    from radia.bem_sibc_solver import compute_phi_inc_from_surface_J
    from esim_cell_problem import ESIMFiniteSlabSolver
    from scipy.linalg import lu_solve

    R_coil, a_coil, gap_deg = 0.030, 0.003, 5
    R_wp, H_wp = 0.010, 0.020
    sigma = 2e6  # steel
    freq = 7000
    omega = 2 * math.pi * freq

    STEEL_BH = [
        [0, 0], [50, 0.1], [100, 0.25], [200, 0.55],
        [500, 0.95], [1000, 1.2], [2000, 1.4], [5000, 1.55],
        [10000, 1.65], [20000, 1.75], [50000, 1.9], [100000, 2.0],
    ]

    # Estimate delta at mu_r ~ 100 (moderate H)
    mu_eff = MU_0 * 100
    delta = math.sqrt(2.0 / (omega * mu_eff * sigma))

    print("=" * 65)
    print("Coupled BEM: Steel Workpiece (nonlinear, Karl iteration)")
    print("=" * 65)
    print(f"Coil: R={R_coil*1e3:.0f}mm, a={a_coil*1e3:.1f}mm")
    print(f"Workpiece: R={R_wp*1e3:.0f}mm, H={H_wp*1e3:.0f}mm, steel")
    print(f"f={freq}Hz, sigma={sigma:.0e}, delta~{delta*1e3:.2f}mm (mu_r~100)")
    print()

    # === Meshes ===
    print("[1/3] Building meshes...")
    mesh_coil = make_gapped_torus_mesh(R_coil, a_coil, gap_deg=gap_deg,
                                        maxh=a_coil)
    wp_cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)
    for f in wp_cyl.faces:
        f.name = "wp"
    mesh_wp = Mesh(OCCGeometry(Glue(wp_cyl.faces)).GenerateMesh(maxh=R_wp / 4))
    with TaskManager():
        mesh_wp.Curve(3)
        print(f"  Coil: {mesh_coil.GetNE(BND)} elems, WP: {mesh_wp.GetNE(BND)} elems")

        # === Assemble ===
        print("[2/3] Assembling coupled solver...")
        solver = CoupledBEMSolver(mesh_coil, mesh_wp)
        print(f"  Coil: {solver.n_J} DOFs ({solver.t_coil_assembly:.1f}s)")
        print(f"  WP:   {solver.wp_solver.ndof} DOFs ({solver.wp_solver.t_assembly:.1f}s)")

        # === ESIM for Z_s (nonlinear) ===
        esim = ESIMFiniteSlabSolver(
            half_thickness=R_wp, bh_curve=STEEL_BH, sigma=sigma,
            frequency=freq, n_nodes=200, geometry='cylinder')

        # === Manual coupled iteration with ESIM Karl loop ===
        print("[3/3] Coupled iteration (ESIM Karl + BEM coupling)...")

        # Initial EFIE (no back-reaction)
        n_J = solver.n_J
        n_c = solver.n_constraint
        rhs0 = np.zeros(n_J + n_c)
        rhs0[n_J:] = solver.g_red
        x = lu_solve(solver.K_lu, rhs0)
        J_coil = x[:n_J]
        L_self = MU_0 * J_coil @ solver.SL_coil @ J_coil

        gf_J = GridFunction(solver.fes_J)
        gf_J.vec.FV().NumPy()[:] = J_coil
        coil_c, coil_a, coil_J = extract_element_J(mesh_coil, gf_J)

        wp_nodes = solver.wp_nodes

        # Initial Z_s
        Z_s = esim.solve(5.0)['Z']
        max_outer = 5   # coupling iterations
        max_karl = 10   # Karl iterations per coupling step
        tol_karl = 1e-3
        relax_karl = 0.5
        tol_outer = 1e-3
        L_prev = L_self

        print(f"  L_self = {L_self*1e9:.2f} nH, |Z_s| init = {abs(Z_s):.4e}")
        print()
        print(f"  {'outer':>5s} {'karl':>5s} {'L_total':>10s} {'DeltaL':>10s} "
              f"{'H_t':>8s} {'P [W]':>12s} {'|Z_s|':>12s} {'dL':>10s}")
        print("-" * 85)

        for outer in range(max_outer):
            # Forward: coil -> phi_inc at workpiece
            phi_inc = compute_phi_inc_from_surface_J(
                wp_nodes, coil_c, coil_a, coil_J, n_quad=20)

            # Karl iteration for Z_s
            for karl in range(max_karl):
                result = solver.wp_solver.solve(phi_inc, Z_s=Z_s, omega=omega)
                H_t_rms = result['H_t_rms']

                Z_s_old = Z_s
                sol = esim.solve(max(H_t_rms, 1e-3))
                Z_s = relax_karl * sol['Z'] + (1 - relax_karl) * Z_s_old
                dZ = abs(Z_s - Z_s_old) / max(abs(Z_s_old), 1e-30)
                if dZ < tol_karl and karl > 0:
                    break

            P_total = result['P_density'] * result['area']

            # Back-reaction: SCATTERED surface current from phi
            # phi_scat = phi_total - phi_inc (large, represents workpiece reflection)
            # J_s_scat = n x (-grad_s(phi_scat))
            phi_vec = result['phi_vec']
            gf_scat_re = GridFunction(solver.wp_solver.fes)
            gf_scat_im = GridFunction(solver.wp_solver.fes)
            gf_scat_re.vec.FV().NumPy()[:] = phi_vec.real - phi_inc
            gf_scat_im.vec.FV().NumPy()[:] = phi_vec.imag

            elem_A_wp = Integrate(CF(1), mesh_wp, BND, element_wise=True)
            grad_re = [Integrate(grad(gf_scat_re)[k], mesh_wp, BND, element_wise=True)
                       for k in range(3)]
            grad_im = [Integrate(grad(gf_scat_im)[k], mesh_wp, BND, element_wise=True)
                       for k in range(3)]

            wp_c, wp_a, wp_Js_re, wp_Js_im = [], [], [], []
            for el in mesh_wp.Elements(BND):
                area = abs(elem_A_wp[el.nr])
                if area < 1e-30:
                    continue
                ht_re = np.array([-grad_re[k][el.nr] / area for k in range(3)])
                ht_im = np.array([-grad_im[k][el.nr] / area for k in range(3)])
                # Compute outward normal from vertex cross product
                verts = [np.array(mesh_wp.vertices[v.nr].point) for v in el.vertices]
                e1 = verts[1] - verts[0]
                e2 = verts[2] - verts[0]
                n_vec = np.cross(e1, e2)
                n_mag = np.linalg.norm(n_vec)
                if n_mag > 1e-30:
                    n_vec /= n_mag
                # Ensure outward: normal should point away from origin for convex shape
                centroid = np.mean(verts, axis=0)
                if np.dot(n_vec, centroid) < 0:
                    n_vec = -n_vec
                # J_s = n x H_t
                js_re = np.cross(n_vec, ht_re)
                js_im = np.cross(n_vec, ht_im)
                verts = [mesh_wp.vertices[v.nr].point for v in el.vertices]
                c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
                wp_c.append(c); wp_a.append(area)
                wp_Js_re.append(js_re); wp_Js_im.append(js_im)
            wp_c = np.array(wp_c); wp_a = np.array(wp_a)
            wp_Js_re = np.array(wp_Js_re); wp_Js_im = np.array(wp_Js_im)

            # Delta_L from complex coupling:
            # Delta_Z = -jw * mu0 * J_coil . A(J_wp_complex)
            # Delta_L = -mu0 * J_coil . A(J_wp_re)
            A_back_re = vector_potential_from_J(coil_c, wp_c, wp_a, wp_Js_re)
            coupling_re = MU_0 * np.sum(np.sum(coil_J * A_back_re, axis=1) * coil_a)
            Delta_L = -coupling_re
            L_total = L_self + Delta_L

            dL = abs(L_total - L_prev) / max(abs(L_prev), 1e-30)
            print(f"  {outer:5d} {karl:5d} {L_total*1e9:10.2f} {Delta_L*1e9:10.4f} "
                  f"{H_t_rms:8.2f} {P_total:12.4e} {abs(Z_s):12.4e} {dL:10.4e}")

            L_prev = L_total

            # Re-solve EFIE with back-reaction (approximate)
            alpha = Delta_L / max(abs(L_self), 1e-30)
            f_back = alpha * (solver.SL_coil @ J_coil)
            rhs = np.zeros(n_J + n_c)
            rhs[:n_J] = -f_back
            rhs[n_J:] = solver.g_red
            x = lu_solve(solver.K_lu, rhs)
            J_coil = x[:n_J]
            gf_J.vec.FV().NumPy()[:] = J_coil
            coil_c, coil_a, coil_J = extract_element_J(mesh_coil, gf_J)
            L_self = MU_0 * J_coil @ solver.SL_coil @ J_coil

            if dL < tol_outer and outer > 0:
                print(f"  Converged at outer iteration {outer}")
                break

        print()
        print("=" * 65)
        print("Coupled BEM Results (Steel)")
        print("-" * 65)
        print(f"  L_self     = {L_self*1e9:.2f} nH")
        print(f"  Delta_L    = {Delta_L*1e9:.4f} nH")
        print(f"  L_total    = {L_total*1e9:.2f} nH")
        print(f"  P_total    = {P_total:.4e} W")
        print(f"  H_t_rms    = {H_t_rms:.2f} A/m")
        print(f"  |Z_s|      = {abs(Z_s):.4e}")
        print("-" * 65)

        # FEM reference
        print()
        print("Running FEM-ESIM for comparison...")
        from fem_esim_3d import run as fem_run
        r_fem = fem_run(frequency=freq, material="steel", sigma=sigma, approach='hole')

        if r_fem:
            print()
            print("=" * 65)
            print(f"Comparison: Coupled BEM vs FEM-ESIM (steel, {freq}Hz)")
            print("=" * 65)
            print(f"{'':>15s} {'Coupled BEM':>14s} {'FEM-ESIM':>14s} {'diff':>8s}")
            print("-" * 55)
            L_f = r_fem.get('L', 0)
            print(f"{'L [nH]':>15s} {L_total*1e9:14.2f} {L_f*1e9:14.2f} "
                  f"{(L_total/L_f-1)*100:+8.1f}%")
            P_f = r_fem['P_total']
            print(f"{'P [W]':>15s} {P_total:14.4e} {P_f:14.4e} "
                  f"{(P_total/P_f-1)*100:+8.1f}%")
            print(f"{'H_t_rms':>15s} {H_t_rms:14.2f} {r_fem['H_t_rms']:14.2f}")
            print(f"{'|Z_s|':>15s} {abs(Z_s):14.4e} {abs(r_fem['Z_s']):14.4e}")


if __name__ == "__main__":
    run()
