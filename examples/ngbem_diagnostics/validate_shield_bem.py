"""
Validation: ShieldBEMSIBC eddy current solver
==============================================

Test case: 1A current wire above a 100x100x2mm aluminum plate at 100 kHz.

Analytical reference (infinite plate, image method):
  P_eddy = (mu_0^2 * omega^2 * I^2) / (32 * pi * sigma * delta^2)
         * integral over wire length

For a finite wire and finite plate, the analytical value serves as
order-of-magnitude check. We verify:
  1. Loop basis construction (sign sums = 0, correct loop count)
  2. V_LL eigenvalues positive (no V_div leakage)
  3. Loss magnitude in physically reasonable range
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

MU_0 = 4e-7 * np.pi


def main():
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh, HDivSurface, BND, TaskManager
    except ImportError:
        print("SKIP: NGSolve/Netgen not available")
        return

    try:
        from ngsolve.bem import MaxwellSingleLayerPotentialOperator
    except ImportError:
        print("SKIP: ngsolve.bem not available")
        return

    from ngbem_eddy import ShieldBEMSIBC, _biot_savart_A

    print("=" * 70)
    print("Validation: ShieldBEMSIBC with corrected divergence matrix")
    print("=" * 70)

    # --- Problem setup ---
    width, height, depth = 0.1, 0.1, 0.002  # 100x100x2 mm aluminum plate
    sigma_al = 3.7e7  # S/m (aluminum)
    freq = 100e3      # 100 kHz
    I_wire = 1.0      # 1 A
    wire_height = 0.01  # 10 mm above plate

    omega = 2 * np.pi * freq
    delta = np.sqrt(2 / (omega * MU_0 * sigma_al))
    print(f"\nProblem parameters:")
    print(f"  Plate: {width*1e3:.0f} x {height*1e3:.0f} x {depth*1e3:.1f} mm")
    print(f"  sigma = {sigma_al:.1e} S/m (aluminum)")
    print(f"  f = {freq/1e3:.0f} kHz")
    print(f"  Skin depth = {delta*1e3:.4f} mm")
    print(f"  Wire: I={I_wire} A at z={wire_height*1e3:.0f} mm")

    # --- Step 1: Create mesh ---
    maxh = 0.025  # 25mm element size
    block = Box(Pnt(-width/2, -height/2, -depth/2),
                Pnt(width/2, height/2, depth/2))
    block.solids.name = "conductor"
    block.faces.name = "surface"
    geo = OCCGeometry(block)
    with TaskManager():
        ngmesh = geo.GenerateMesh(maxh=maxh)
        mesh = Mesh(ngmesh)

    n_bnd = sum(1 for _ in mesh.Elements(BND))
    n_vert = mesh.nv
    print(f"\n  Mesh: {n_bnd} boundary faces, {n_vert} vertices")

    # --- Step 2: Build loop basis and check ---
    print("\n--- Loop basis diagnostics ---")
    fes = HDivSurface(mesh, order=0)
    from ngbem_eddy import LoopBasisBuilder
    loop = LoopBasisBuilder(mesh, fes)

    print(f"  Active DOFs: {loop.n_active}")
    print(f"  Faces: {loop.n_faces}")
    print(f"  Loops: {loop.n_loops}")

    # Check sign sums (should be 0 for closed surface)
    D = loop._D
    sign_sums = []
    for idx in range(loop.n_active):
        col = D[:, idx]
        nz = col[col != 0]
        sign_sums.append(np.sum(nz))

    min_ss = min(sign_sums)
    max_ss = max(sign_sums)
    n_zero = sum(1 for s in sign_sums if abs(s) < 0.01)
    print(f"  D sign sums: min={min_ss:.0f}, max={max_ss:.0f}, "
          f"#zero={n_zero}/{loop.n_active}")

    # Expected: V - 1 loops (genus-0 closed surface)
    expected_loops = n_vert - 1
    print(f"  Expected loops (V-1): {expected_loops}")

    test_sign_sums = abs(min_ss) < 0.01 and abs(max_ss) < 0.01
    test_loop_count = loop.n_loops == expected_loops
    print(f"  Sign sum check: {'PASS' if test_sign_sums else 'FAIL'}")
    print(f"  Loop count check: {'PASS' if test_loop_count else 'FAIL'}")

    # --- Step 3: Assemble and check V_LL eigenvalues ---
    print("\n--- Assembling BEM operators ---")
    solver = ShieldBEMSIBC(mesh, sigma_al, mu_r=1.0)
    solver.assemble(intorder=4)
    print(f"  Assembly time: {solver.t_assemble:.2f} s")

    V_LL = solver._V_LL
    # In loop subspace: V_LL = V_vec_LL (positive definite)
    # The V_div term vanishes for divergence-free currents
    V_LL_real_sym = 0.5 * (V_LL.real + V_LL.real.T)
    eigvals = np.linalg.eigvalsh(V_LL_real_sym)
    print(f"  V_LL eigenvalues: min={eigvals.min():.6e}, "
          f"max={eigvals.max():.6e}")
    print(f"  mu_0 * V_LL -> L range: "
          f"[{MU_0*eigvals.min():.4e}, {MU_0*eigvals.max():.4e}] H")

    test_V_positive = eigvals.min() > 0
    print(f"  V_LL positive definite: "
          f"{'PASS' if test_V_positive else 'FAIL'}")

    M_LL = solver._M_LL
    M_eigvals = np.linalg.eigvalsh(0.5 * (M_LL + M_LL.T))
    print(f"  M_LL eigenvalues: min={M_eigvals.min():.6e}, "
          f"max={M_eigvals.max():.6e}")
    test_M_positive = M_eigvals.min() > 0
    print(f"  M_LL positive definite: "
          f"{'PASS' if test_M_positive else 'FAIL'}")

    # --- Step 4: Solve with wire excitation ---
    print("\n--- Solving at f={:.0f} kHz ---".format(freq/1e3))

    wire_center = np.array([0, 0, wire_height + depth/2])
    wire_dir = np.array([1.0, 0, 0])  # x-directed wire
    wire_length = 0.2  # 200mm wire (longer than plate)

    def A_inc(points):
        return I_wire * _biot_savart_A(points, wire_center, wire_dir,
                                       wire_length)

    result = solver.solve(freq, A_inc_func=A_inc)
    print(f"  Solve time: {result['t_solve']:.4f} s")
    print(f"  |J_loop|_max = {result['J_loop_max']:.6e}")

    # --- Step 5: Compute loss ---
    P = solver.compute_loss()
    print(f"  Eddy current loss P = {P:.6e} W")

    # --- Analytical reference (order of magnitude) ---
    # For infinite plate with infinite wire at height h:
    # P/L = mu_0^2 * omega^2 * I^2 / (32*pi*sigma*delta^2) * h_factor
    # This is just for order-of-magnitude comparison.
    # Expected: ~10^-6 to 10^-4 W for this geometry.
    print(f"\n--- Results summary ---")
    print(f"  Loss: {P:.4e} W")
    print(f"  |Loss|: {abs(P):.4e} W")

    test_loss_positive = P > 0
    test_loss_order = 1e-8 < abs(P) < 1e-1
    print(f"  Loss positive: {'PASS' if test_loss_positive else 'FAIL'}")
    print(f"  Loss in [1e-8, 1e-1] W: "
          f"{'PASS' if test_loss_order else 'FAIL'}")

    # --- Step 6: Frequency sweep ---
    print("\n--- Frequency sweep ---")
    freqs = [1e3, 10e3, 100e3, 1e6]
    losses = []
    for f in freqs:
        solver.solve(f, A_inc_func=A_inc)
        P_f = solver.compute_loss()
        losses.append(P_f)
        d_f = np.sqrt(2 / (2*np.pi*f * MU_0 * sigma_al))
        print(f"  f={f/1e3:8.1f} kHz, delta={d_f*1e3:.4f} mm, "
              f"P={P_f:.4e} W")

    # Loss should increase with frequency (more eddy currents)
    test_monotonic = all(losses[i] < losses[i+1]
                         for i in range(len(losses)-1))
    print(f"  Loss increases with frequency: "
          f"{'PASS' if test_monotonic else 'FAIL'}")

    # --- Summary ---
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)
    tests = [
        ("D sign sums all zero", test_sign_sums),
        ("Loop count = V-1", test_loop_count),
        ("V_LL positive definite", test_V_positive),
        ("M_LL positive definite", test_M_positive),
        ("Loss > 0", test_loss_positive),
        ("Loss in reasonable range", test_loss_order),
        ("Loss monotonic with freq", test_monotonic),
    ]
    all_pass = True
    for name, passed in tests:
        status = "PASS" if passed else "FAIL"
        if not passed:
            all_pass = False
        print(f"  [{status}] {name}")

    print(f"\nOverall: {'ALL PASSED' if all_pass else 'SOME FAILED'}")


if __name__ == '__main__':
    main()
