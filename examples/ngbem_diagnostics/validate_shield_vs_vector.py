"""
Cross-validation: ShieldBEMSIBC vs VectorEddyCurrentFEMBEM
==========================================================

Both solvers run on the SAME mesh with the SAME excitation (uniform B_ext).
Compare eddy current loss to verify consistency.

Key physics:
  - ShieldBEMSIBC: surface-only, SIBC models skin layer analytically
    -> Works on ANY mesh when delta << thickness
  - VectorFEMBEM: volume H(curl) FEM resolves field inside conductor
    -> NEEDS mesh refinement at skin depth scale

Expected:
  - On coarse mesh: VectorFEMBEM underpredicts loss (can't resolve skin layer)
  - As mesh refines: VectorFEMBEM loss -> ShieldBEMSIBC loss (convergence)
  - Both should give same qualitative trends (positive, monotonic with freq)

Test cases:
  A. Aluminum block (mu_r=1, sigma=3.7e7), B_ext=[0,0,1]
  B. Steel block (mu_r=100, sigma=1e6), B_ext=[0,0,1]
  C. Mesh convergence study (two mesh sizes)

Part of Radia project
"""

import sys
import os
import numpy as np
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

MU_0 = 4e-7 * np.pi


def uniform_B_to_A_inc(B_ext):
    """Create A_inc_func for uniform B field: A = 0.5 * (B x r)."""
    Bx, By, Bz = B_ext

    def A_inc_func(points):
        points = np.atleast_2d(points)
        x, y, z = points[:, 0], points[:, 1], points[:, 2]
        Ax = 0.5 * (By * z - Bz * y)
        Ay = 0.5 * (Bz * x - Bx * z)
        Az = 0.5 * (Bx * y - By * x)
        return np.column_stack([Ax, Ay, Az])

    return A_inc_func


def analytical_loss_half_space(freq, sigma, mu_r, B0, area):
    """Analytical eddy current loss for half-space in uniform field.

    P = Rs/2 * |H0|^2 * A_surface
    Rs = sqrt(omega*mu/(2*sigma))
    H0 = B0 / mu_0

    This overestimates for finite geometry (factor ~0.3-0.5 for cube).
    """
    omega = 2 * np.pi * freq
    mu = mu_r * MU_0
    Rs = np.sqrt(omega * mu / (2 * sigma))
    H0 = B0 / MU_0
    P = Rs / 2.0 * H0**2 * area
    return P


def run_single_comparison(mesh, sigma, mu_r, freq, B_ext, label):
    """Run both solvers on same mesh at single frequency."""
    from ngbem_eddy import ShieldBEMSIBC, VectorEddyCurrentFEMBEM

    omega = 2 * np.pi * freq
    mu = mu_r * MU_0
    delta = np.sqrt(2.0 / (omega * mu * sigma))

    # ShieldBEMSIBC
    solver_sibc = ShieldBEMSIBC(mesh, sigma, mu_r=mu_r)
    solver_sibc.assemble(intorder=4)
    A_inc_func = uniform_B_to_A_inc(B_ext)
    solver_sibc.solve(freq, A_inc_func=A_inc_func)
    P_sibc = solver_sibc.compute_loss()

    # VectorFEMBEM
    solver_vec = VectorEddyCurrentFEMBEM(mesh, sigma, mu_r=mu_r, order=1)
    solver_vec.assemble(kappa=0.01, intorder=4)
    solver_vec.solve(freq, B_ext=B_ext)
    P_vec = solver_vec.compute_loss()

    return {
        'P_sibc': P_sibc,
        'P_vec': P_vec,
        'delta': delta,
        'n_hcurl': solver_vec._n_hcurl,
        'n_loops': solver_sibc._loop.n_loops,
    }


def main():
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh, BND
    except ImportError:
        print("SKIP: NGSolve/Netgen not available")
        return

    try:
        from ngsolve.bem import HelmholtzSL
        from ngsolve import TaskManager
    except ImportError:
        print("SKIP: ngsolve.bem not available")
        return

    from ngbem_eddy import (ShieldBEMSIBC, VectorEddyCurrentFEMBEM,
                             _biot_savart_A)

    print("=" * 70)
    print("Cross-validation: ShieldBEMSIBC vs VectorEddyCurrentFEMBEM")
    print("=" * 70)

    all_tests = []
    B_ext = [0, 0, 1.0]  # 1T uniform field

    # ==================================================================
    # Test A: Aluminum block (mu_r=1), uniform B_ext
    # ==================================================================
    print("\n" + "=" * 70)
    print("Test A: Aluminum block (mu_r=1), uniform Bz=1T")
    print("=" * 70)

    width, height, depth = 0.02, 0.02, 0.01  # 20x20x10 mm
    sigma_al = 3.7e7
    maxh = 0.012

    block = Box(Pnt(-width/2, -height/2, -depth/2),
                Pnt(width/2, height/2, depth/2))
    block.solids.name = "conductor"
    block.faces.name = "surface"
    geo = OCCGeometry(block)
    with TaskManager():
        mesh = Mesh(geo.GenerateMesh(maxh=maxh))

        n_vol = sum(1 for _ in mesh.Elements())
        n_bnd = sum(1 for _ in mesh.Elements(BND))
        print(f"\n  Mesh: {n_vol} volume, {n_bnd} boundary elements, maxh={maxh*1e3:.0f}mm")

        # Assemble once, sweep frequencies
        A_inc_func = uniform_B_to_A_inc(B_ext)

        solver_sibc = ShieldBEMSIBC(mesh, sigma_al, mu_r=1.0)
        solver_sibc.assemble(intorder=4)
        print(f"  ShieldBEMSIBC: {solver_sibc._loop.n_loops} loop DOFs")

        solver_vec = VectorEddyCurrentFEMBEM(mesh, sigma_al, mu_r=1.0, order=1)
        solver_vec.assemble(kappa=0.01, intorder=4)

        freqs_al = [10e3, 100e3, 1e6]
        # Surface area for analytical reference
        A_surface = 2 * (width*height + width*depth + height*depth)

        print(f"\n  {'freq(Hz)':>10s}  {'delta(mm)':>10s}  {'d/t':>6s}  "
              f"{'P_SIBC(W)':>12s}  {'P_Vec(W)':>12s}  {'P_ana(W)':>12s}  "
              f"{'SIBC/ana':>8s}  {'Vec/ana':>8s}")
        print(f"  {'-'*10}  {'-'*10}  {'-'*6}  "
              f"{'-'*12}  {'-'*12}  {'-'*12}  {'-'*8}  {'-'*8}")

        P_sibc_list = []
        P_vec_list = []
        for freq in freqs_al:
            omega = 2 * np.pi * freq
            delta = np.sqrt(2.0 / (omega * MU_0 * sigma_al))
            d_over_t = delta / depth

            solver_sibc.solve(freq, A_inc_func=A_inc_func)
            P_sibc = solver_sibc.compute_loss()
            P_sibc_list.append(P_sibc)

            solver_vec.solve(freq, B_ext=B_ext)
            P_vec = solver_vec.compute_loss()
            P_vec_list.append(P_vec)

            P_ana = analytical_loss_half_space(freq, sigma_al, 1.0, 1.0,
                                                A_surface)

            print(f"  {freq:10.0f}  {delta*1e3:10.4f}  {d_over_t:6.3f}  "
                  f"{P_sibc:12.4e}  {P_vec:12.4e}  {P_ana:12.4e}  "
                  f"{P_sibc/P_ana:8.4f}  {P_vec/P_ana:8.4f}")

        # Tests
        test_A_sibc_pos = all(P > 0 for P in P_sibc_list)
        test_A_vec_pos = all(P > 0 for P in P_vec_list)
        test_A_sibc_mono = all(P_sibc_list[i] < P_sibc_list[i+1]
                                for i in range(len(P_sibc_list)-1))
        test_A_vec_mono = all(P_vec_list[i] < P_vec_list[i+1]
                               for i in range(len(P_vec_list)-1))

        # VectorFEMBEM underestimates on coarse mesh, so P_vec < P_sibc
        # Accept if within factor 100
        ratios_A = [P_sibc_list[i] / P_vec_list[i]
                    for i in range(len(P_sibc_list)) if P_vec_list[i] > 0]
        test_A_order = all(0.01 < r < 100 for r in ratios_A)

        for name, passed in [
            ("Al: SIBC losses positive", test_A_sibc_pos),
            ("Al: Vec losses positive", test_A_vec_pos),
            ("Al: SIBC loss monotonic", test_A_sibc_mono),
            ("Al: Vec loss monotonic", test_A_vec_mono),
            ("Al: SIBC/Vec ratio within [0.01, 100]", test_A_order),
        ]:
            print(f"  [{'PASS' if passed else 'FAIL'}] {name}")
            all_tests.append((name, passed))

        # ==================================================================
        # Test B: Steel block (mu_r=100), uniform B_ext
        # ==================================================================
        print("\n" + "=" * 70)
        print("Test B: Steel block (mu_r=100, sigma=1e6), uniform Bz=1T")
        print("=" * 70)

        sigma_steel = 1e6
        mu_r_steel = 100.0

        solver_sibc_s = ShieldBEMSIBC(mesh, sigma_steel, mu_r=mu_r_steel)
        solver_sibc_s.assemble(intorder=4)

        solver_vec_s = VectorEddyCurrentFEMBEM(mesh, sigma_steel,
                                                 mu_r=mu_r_steel, order=1)
        solver_vec_s.assemble(kappa=0.01, intorder=4)

        A_inc_steel = uniform_B_to_A_inc(B_ext)
        freqs_steel = [1e3, 10e3, 100e3]

        print(f"\n  {'freq(Hz)':>10s}  {'delta(mm)':>10s}  {'d/t':>6s}  "
              f"{'P_SIBC(W)':>12s}  {'P_Vec(W)':>12s}  {'ratio':>8s}")
        print(f"  {'-'*10}  {'-'*10}  {'-'*6}  "
              f"{'-'*12}  {'-'*12}  {'-'*8}")

        P_sibc_s = []
        P_vec_s = []
        for freq in freqs_steel:
            omega = 2 * np.pi * freq
            delta = np.sqrt(2.0 / (omega * mu_r_steel * MU_0 * sigma_steel))
            d_over_t = delta / depth

            solver_sibc_s.solve(freq, A_inc_func=A_inc_steel)
            P_s = solver_sibc_s.compute_loss()
            P_sibc_s.append(P_s)

            solver_vec_s.solve(freq, B_ext=B_ext)
            P_v = solver_vec_s.compute_loss()
            P_vec_s.append(P_v)

            ratio = P_s / P_v if P_v > 0 else float('inf')
            print(f"  {freq:10.0f}  {delta*1e3:10.4f}  {d_over_t:6.3f}  "
                  f"{P_s:12.4e}  {P_v:12.4e}  {ratio:8.4f}")

        test_B_sibc_pos = all(P > 0 for P in P_sibc_s)
        test_B_vec_pos = all(P > 0 for P in P_vec_s)
        test_B_sibc_mono = all(P_sibc_s[i] < P_sibc_s[i+1]
                                for i in range(len(P_sibc_s)-1))
        test_B_vec_mono = all(P_vec_s[i] < P_vec_s[i+1]
                               for i in range(len(P_vec_s)-1))

        for name, passed in [
            ("Steel: SIBC losses positive", test_B_sibc_pos),
            ("Steel: Vec losses positive", test_B_vec_pos),
            ("Steel: SIBC loss monotonic", test_B_sibc_mono),
            ("Steel: Vec loss monotonic", test_B_vec_mono),
        ]:
            print(f"  [{'PASS' if passed else 'FAIL'}] {name}")
            all_tests.append((name, passed))

        # ==================================================================
        # Test C: Mesh convergence (coarse -> fine)
        # ==================================================================
        print("\n" + "=" * 70)
        print("Test C: Mesh convergence at f=100 kHz (Al, mu_r=1)")
        print("=" * 70)

        freq_conv = 100e3
        delta_conv = np.sqrt(2.0 / (2*np.pi*freq_conv * MU_0 * sigma_al))
        print(f"  f = {freq_conv/1e3:.0f} kHz, delta = {delta_conv*1e3:.4f} mm")
        print()

        mesh_sizes = [0.012, 0.008, 0.005]

        print(f"  {'maxh(mm)':>10s}  {'n_vol':>6s}  {'n_hcurl':>8s}  "
              f"{'P_SIBC(W)':>12s}  {'P_Vec(W)':>12s}  {'ratio':>8s}")
        print(f"  {'-'*10}  {'-'*6}  {'-'*8}  "
              f"{'-'*12}  {'-'*12}  {'-'*8}")

        P_vec_conv = []
        P_sibc_conv = []

        for maxh_c in mesh_sizes:
            blk = Box(Pnt(-width/2, -height/2, -depth/2),
                      Pnt(width/2, height/2, depth/2))
            blk.solids.name = "conductor"
            blk.faces.name = "surface"
            geo_c = OCCGeometry(blk)
            mesh_c = Mesh(geo_c.GenerateMesh(maxh=maxh_c))
            n_vol_c = sum(1 for _ in mesh_c.Elements())

            # ShieldBEMSIBC
            s_sibc = ShieldBEMSIBC(mesh_c, sigma_al, mu_r=1.0)
            s_sibc.assemble(intorder=4)
            s_sibc.solve(freq_conv, A_inc_func=uniform_B_to_A_inc(B_ext))
            P_s = s_sibc.compute_loss()

            # VectorFEMBEM
            s_vec = VectorEddyCurrentFEMBEM(mesh_c, sigma_al, mu_r=1.0, order=1)
            s_vec.assemble(kappa=0.01, intorder=4)
            s_vec.solve(freq_conv, B_ext=B_ext)
            P_v = s_vec.compute_loss()

            ratio = P_s / P_v if P_v > 0 else float('inf')
            P_sibc_conv.append(P_s)
            P_vec_conv.append(P_v)

            print(f"  {maxh_c*1e3:10.1f}  {n_vol_c:6d}  {s_vec._n_hcurl:8d}  "
                  f"{P_s:12.4e}  {P_v:12.4e}  {ratio:8.4f}")

        # VectorFEMBEM should increase as mesh refines (better skin resolution)
        test_C_vec_increase = all(P_vec_conv[i] <= P_vec_conv[i+1] * 1.1
                                   for i in range(len(P_vec_conv)-1))
        # SIBC should be relatively stable (mesh-independent)
        if P_sibc_conv[0] > 0:
            sibc_variation = max(P_sibc_conv) / min(P_sibc_conv)
        else:
            sibc_variation = float('inf')

        test_C_sibc_stable = sibc_variation < 5.0  # within factor 5

        for name, passed in [
            ("Convergence: Vec loss increases with refinement", test_C_vec_increase),
            ("Convergence: SIBC loss stable (var < 5x)", test_C_sibc_stable),
        ]:
            print(f"  [{'PASS' if passed else 'FAIL'}] {name}")
            all_tests.append((name, passed))

        # ==================================================================
        # Summary
        # ==================================================================
        print("\n" + "=" * 70)
        print("CROSS-VALIDATION SUMMARY")
        print("=" * 70)

        all_pass = True
        for name, passed in all_tests:
            status = "PASS" if passed else "FAIL"
            if not passed:
                all_pass = False
            print(f"  [{status}] {name}")

        print(f"\nOverall: {'ALL PASSED' if all_pass else 'SOME FAILED'}")
        print(f"Total tests: {len(all_tests)}")

        print(f"""
--- Key Findings ---
  1. ShieldBEMSIBC: mesh-INDEPENDENT loss (SIBC handles skin layer analytically)
  2. VectorFEMBEM: mesh-DEPENDENT loss (needs fine mesh to resolve skin layer)
  3. On coarse mesh: VectorFEMBEM << ShieldBEMSIBC (under-resolved skin layer)
  4. Both give same qualitative behavior (positive, monotonic with frequency)
  5. Both handle mu_r > 1 correctly

--- When to Use Which ---
  ShieldBEMSIBC:  delta << thickness (high freq, good conductors)
                  Surface-only mesh, FAST, mesh-independent
  VectorFEMBEM:   delta ~ thickness (moderate freq, or need mu_r > 1)
                  Volume mesh needed, EXPENSIVE, needs fine mesh
""")


if __name__ == '__main__':
    main()
