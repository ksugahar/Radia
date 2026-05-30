"""
Precision test: filament approximation vs actual surface current Biot-Savart
for BEM-SIBC workpiece solver.

Compares phi_inc computed two ways:
  (a) compute_phi_inc_from_loop   -- filament (line current) approximation
  (b) compute_phi_inc_from_surface_J -- actual surface current from EFIE solve

Then solves BEM-SIBC with each phi_inc (copper workpiece, 1 kHz) and compares
H_t_rms and P_total.  Also runs FEM-ESIM as a reference.

Usage:
    python -u test_surface_J_accuracy.py
"""

import math
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))
sys.path.insert(0, os.path.dirname(__file__))

MU_0 = 4e-7 * np.pi


def main():
    from ngsolve import Mesh, Integrate, CF, ds, BND, TaskManager
    from netgen.occ import Cylinder, Pnt, Dir, OCCGeometry, Glue
    from radia.bem_sibc_solver import (ScalarBIESIBCSolver,
                                  compute_phi_inc_from_loop,
                                  compute_phi_inc_from_surface_J)
    from bem_inductance import compute_inductance_source_sink
    from esim_cell_problem import ESIMFiniteSlabSolver
    from impedance_esim import make_gapped_torus_mesh

    # ---- Parameters ----
    R_coil = 0.030      # 30 mm
    a_coil = 0.003      # 3 mm
    gap_deg = 5
    R_wp = 0.010        # 10 mm
    H_wp = 0.020        # 20 mm
    sigma = 5.8e7       # copper
    frequency = 1000    # 1 kHz
    material = "copper"
    I_total = 1.0       # 1 A

    omega = 2 * np.pi * frequency
    delta = math.sqrt(2.0 / (omega * MU_0 * sigma))

    print("=" * 72)
    print("Precision Test: Filament vs Surface-J Biot-Savart for BEM-SIBC")
    print("=" * 72)
    print(f"Coil:      R={R_coil*1e3:.0f}mm, a={a_coil*1e3:.1f}mm, gap={gap_deg}deg")
    print(f"Workpiece: R={R_wp*1e3:.0f}mm, H={H_wp*1e3:.0f}mm")
    print(f"Material:  {material}, sigma={sigma:.1e} S/m, f={frequency} Hz")
    print(f"Skin depth: delta={delta*1e3:.3f}mm, delta/R={delta/R_wp:.3f}")
    print()

    # ==================================================================
    # Step 1: Coil EFIE solve -> surface current J
    # ==================================================================
    print("[1/6] Coil BEM solve (EFIE source/sink)...")
    t0 = time.perf_counter()
    mesh_coil = make_gapped_torus_mesh(R_coil, a_coil, gap_deg=gap_deg,
                                        maxh=a_coil)
    sol_coil = compute_inductance_source_sink(mesh_coil)
    if 'error' in sol_coil:
        print(f"  ERROR: {sol_coil['error']}")
        return
    L_coil = sol_coil['L']
    gf_J_coil = sol_coil['gf_J']
    t_coil = time.perf_counter() - t0
    print(f"  L_coil = {L_coil*1e9:.2f} nH, "
          f"n_J = {sol_coil['n_J']}, n_f = {sol_coil['n_f']} ({t_coil:.1f}s)")

    # Extract per-element J data from coil mesh for surface-J Biot-Savart
    print("  Extracting per-element coil surface current...")
    elem_A = Integrate(CF(1), mesh_coil, VOL_or_BND=BND, element_wise=True)
    elem_Jx = Integrate(gf_J_coil[0], mesh_coil, VOL_or_BND=BND,
                         element_wise=True)
    elem_Jy = Integrate(gf_J_coil[1], mesh_coil, VOL_or_BND=BND,
                         element_wise=True)
    elem_Jz = Integrate(gf_J_coil[2], mesh_coil, VOL_or_BND=BND,
                         element_wise=True)

    src_centroids, src_areas, src_J_vecs = [], [], []
    for el in mesh_coil.Elements(BND):
        area = abs(elem_A[el.nr])
        if area < 1e-30:
            continue
        jvec = np.array([elem_Jx[el.nr], elem_Jy[el.nr],
                         elem_Jz[el.nr]]) / area
        verts = [mesh_coil.vertices[v.nr].point for v in el.vertices]
        c = np.mean([(v[0], v[1], v[2]) for v in verts], axis=0)
        src_centroids.append(c)
        src_areas.append(area)
        src_J_vecs.append(jvec)

    src_centroids = np.array(src_centroids)
    src_areas = np.array(src_areas)
    src_J_vecs = np.array(src_J_vecs)
    print(f"  {len(src_areas)} coil elements, "
          f"total area = {np.sum(src_areas):.6e} m^2")

    # ==================================================================
    # Step 2: Workpiece surface mesh
    # ==================================================================
    print()
    print("[2/6] Meshing workpiece surface...")
    t0 = time.perf_counter()

    wp_cyl = Cylinder(Pnt(0, 0, -H_wp / 2), Dir(0, 0, 1), R_wp, H_wp)
    for f in wp_cyl.faces:
        f.name = "wp_surface"
    geo = OCCGeometry(Glue(wp_cyl.faces))
    maxh_wp = R_wp / 4
    ngmesh = geo.GenerateMesh(maxh=maxh_wp)
    mesh_wp = Mesh(ngmesh)
    with TaskManager():
        mesh_wp.Curve(3)

    ne_wp = mesh_wp.GetNE(BND)
    nv_wp = mesh_wp.nv
    area_wp = abs(Integrate(CF(1), mesh_wp, BND))
    t_mesh = time.perf_counter() - t0
    print(f"  {ne_wp} elements, {nv_wp} vertices, "
          f"area = {area_wp:.6e} m^2 ({t_mesh:.1f}s)")

    # ==================================================================
    # Step 3: Assemble BEM operators on workpiece
    # ==================================================================
    print()
    print("[3/6] Assembling BEM operators (ScalarBIESIBCSolver)...")
    solver = ScalarBIESIBCSolver(mesh_wp, order=1)
    print(f"  {solver.ndof} DOFs, assembly: {solver.t_assembly:.1f}s")

    # Get vertex coordinates for phi_inc evaluation
    node_coords = np.array([[mesh_wp.vertices[i].point[j] for j in range(3)]
                            for i in range(mesh_wp.nv)])

    # ==================================================================
    # Step 4a: phi_inc from filament approximation
    # ==================================================================
    print()
    print("[4/6] Computing phi_inc...")
    print("  (a) Filament approximation (compute_phi_inc_from_loop)...")
    t0 = time.perf_counter()
    phi_inc_filament = compute_phi_inc_from_loop(
        node_coords, loop_center=[0, 0, 0], loop_radius=R_coil,
        current=I_total, n_quad=30, gap_deg=gap_deg)
    t_fil = time.perf_counter() - t0
    print(f"      range: [{phi_inc_filament.min():.6f}, "
          f"{phi_inc_filament.max():.6f}] ({t_fil:.1f}s)")

    # ==================================================================
    # Step 4b: phi_inc from actual surface current
    # ==================================================================
    print("  (b) Surface current (compute_phi_inc_from_surface_J)...")
    t0 = time.perf_counter()
    phi_inc_surfJ = compute_phi_inc_from_surface_J(
        node_coords, src_centroids, src_areas, src_J_vecs, n_quad=20)
    t_surf = time.perf_counter() - t0
    print(f"      range: [{phi_inc_surfJ.min():.6f}, "
          f"{phi_inc_surfJ.max():.6f}] ({t_surf:.1f}s)")

    # Compare phi_inc directly
    diff_phi = phi_inc_surfJ - phi_inc_filament
    rel_diff_phi = np.abs(diff_phi) / np.maximum(np.abs(phi_inc_filament), 1e-30)
    print()
    print("  phi_inc comparison (surface-J vs filament):")
    print(f"    max |diff|          = {np.max(np.abs(diff_phi)):.6e}")
    print(f"    mean |diff|         = {np.mean(np.abs(diff_phi)):.6e}")
    print(f"    max |rel diff|      = {np.max(rel_diff_phi)*100:.2f}%")
    print(f"    mean |rel diff|     = {np.mean(rel_diff_phi)*100:.2f}%")
    print(f"    RMS(filament)       = "
          f"{np.sqrt(np.mean(phi_inc_filament**2)):.6e}")
    print(f"    RMS(surface-J)      = "
          f"{np.sqrt(np.mean(phi_inc_surfJ**2)):.6e}")

    # ==================================================================
    # Step 5: Solve BEM-SIBC with each phi_inc (Karl iteration)
    # ==================================================================
    print()
    print("[5/6] Solving BEM-SIBC with ESIM (Karl iteration)...")

    esim_solver = ESIMFiniteSlabSolver(
        half_thickness=R_wp, bh_curve=None, sigma=sigma,
        frequency=frequency, mu_r=1.0, n_nodes=200,
        geometry='cylinder')

    results = {}
    for label, phi_inc_vec in [("filament", phi_inc_filament),
                                ("surface-J", phi_inc_surfJ)]:
        print(f"\n  --- {label} ---")
        Z_s = esim_solver.solve(5.0)['Z']
        max_iter = 15
        tol = 1e-3
        relax_k = 0.5

        print(f"  |Z_s| init = {abs(Z_s):.4e}")
        print(f"  {'Iter':>4s} {'|Z_s|':>12s} {'H_t_rms':>10s} "
              f"{'P [W]':>12s} {'dZ/Z':>10s}")

        for iteration in range(max_iter):
            result = solver.solve(phi_inc_vec, Z_s=Z_s, omega=omega)
            H_t_rms = result['H_t_rms']
            P_density = result['P_density']
            P_total = P_density * result['area']

            Z_s_old = Z_s
            sol = esim_solver.solve(max(H_t_rms, 1e-3))
            Z_s_new = sol['Z']
            Z_s = relax_k * Z_s_new + (1 - relax_k) * Z_s_old

            dZ = abs(Z_s - Z_s_old) / max(abs(Z_s_old), 1e-30)
            print(f"  {iteration:4d} {abs(Z_s):12.4e} {H_t_rms:10.2f} "
                  f"{P_total:12.4e} {dZ:10.4e}")

            if dZ < tol and iteration > 0:
                print(f"  Converged at iteration {iteration}")
                break

        # Final solve with converged Z_s
        result = solver.solve(phi_inc_vec, Z_s=Z_s, omega=omega)
        H_t_rms = result['H_t_rms']
        P_total = result['P_density'] * result['area']
        Q_total = 0.5 * Z_s.imag * H_t_rms**2 * result['area']

        results[label] = {
            'P_total': P_total,
            'Q_total': Q_total,
            'H_t_rms': H_t_rms,
            'Z_s': Z_s,
            'area': result['area'],
        }

    # ==================================================================
    # Step 6: FEM-ESIM reference
    # ==================================================================
    print()
    print("[6/6] FEM-ESIM reference (may take 100-200s)...")
    from fem_esim_3d import run as fem_run
    r_fem = fem_run(R_coil=R_coil, a_coil=a_coil, gap_deg=gap_deg,
                     R_wp=R_wp, H_wp=H_wp,
                     sigma=sigma, frequency=frequency, material=material,
                     R_air=0.12, maxh_air=0.015, maxh_coil=0.005,
                     order=1, esim_geometry='cylinder',
                     approach='hole')

    # ==================================================================
    # Comparison table
    # ==================================================================
    print()
    print("=" * 72)
    print("COMPARISON TABLE: Filament vs Surface-J vs FEM-ESIM")
    print("=" * 72)
    print(f"  Material: {material}, f={frequency} Hz, "
          f"delta={delta*1e3:.3f} mm, delta/R={delta/R_wp:.3f}")
    print()

    r_fil = results["filament"]
    r_sj = results["surface-J"]

    headers = ["", "Filament", "Surface-J", "FEM-ESIM",
               "fil-sJ diff", "fil-FEM diff"]
    row_fmt = "{:>20s} {:>14s} {:>14s} {:>14s} {:>12s} {:>12s}"
    print(row_fmt.format(*headers))
    print("-" * 88)

    def fmt_row(label, v_fil, v_sj, v_fem):
        s_fil = f"{v_fil:.4e}"
        s_sj = f"{v_sj:.4e}"
        s_fem = f"{v_fem:.4e}" if v_fem is not None else "N/A"
        if abs(v_sj) > 1e-30:
            d_fs = f"{(v_fil - v_sj) / v_sj * 100:+.2f}%"
        else:
            d_fs = "N/A"
        if v_fem is not None and abs(v_fem) > 1e-30:
            d_ff = f"{(v_fil - v_fem) / v_fem * 100:+.2f}%"
        else:
            d_ff = "N/A"
        print(row_fmt.format(label, s_fil, s_sj, s_fem, d_fs, d_ff))

    P_fem = r_fem['P_total'] if r_fem else None
    H_fem = r_fem['H_t_rms'] if r_fem else None
    Zs_fem = abs(r_fem['Z_s']) if r_fem else None

    fmt_row("P_total [W]",
            r_fil['P_total'], r_sj['P_total'], P_fem)
    fmt_row("H_t_rms [A/m]",
            r_fil['H_t_rms'], r_sj['H_t_rms'], H_fem)
    fmt_row("|Z_s| [Ohm]",
            abs(r_fil['Z_s']), abs(r_sj['Z_s']), Zs_fem)

    print()
    print("phi_inc statistics:")
    print(f"  max |phi_fil - phi_sJ| / |phi_fil| = "
          f"{np.max(rel_diff_phi)*100:.2f}%")
    print(f"  mean |phi_fil - phi_sJ| / |phi_fil| = "
          f"{np.mean(rel_diff_phi)*100:.2f}%")

    # Physics interpretation
    print()
    print("Interpretation:")
    dp_fs = ((r_fil['P_total'] - r_sj['P_total'])
             / r_sj['P_total'] * 100)
    print(f"  Filament vs Surface-J power difference: {dp_fs:+.2f}%")
    if abs(dp_fs) < 1:
        print("  -> Filament approximation is adequate for this geometry.")
        print("     (Workpiece is far enough from coil cross-section.)")
    elif abs(dp_fs) < 5:
        print("  -> Small difference; filament is reasonable but not exact.")
    else:
        print("  -> Significant difference; surface current distribution matters.")
        print("     (Workpiece is close to coil cross-section.)")

    if r_fem:
        dp_ff = (r_fil['P_total'] - P_fem) / P_fem * 100
        print(f"  Filament BEM vs FEM-ESIM power difference: {dp_ff:+.2f}%")
        dp_sf = (r_sj['P_total'] - P_fem) / P_fem * 100
        print(f"  Surface-J BEM vs FEM-ESIM power difference: {dp_sf:+.2f}%")

    print()
    print("Done.")


if __name__ == "__main__":
    main()
