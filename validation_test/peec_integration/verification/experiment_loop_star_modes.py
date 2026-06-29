"""
test_loop_star_modes.py

Compare all three solver modes of NGBEMPEECSolver:
  1. 'mqs':        Loop-only (R + jwL)
  2. 'full':       Reformulated Schur complement (fixed P/(jw) issue)
  3. 'stabilized': Weggler's stabilized BEM (Q_0 from LaplaceSL)

Key findings:
  - 'full' and 'stabilized' include displacement current (capacitive) effects
  - For good conductors (Cu), displacement current is negligible -> MQS is correct
  - full/stabilized show higher Z because single-plate capacitance is very small
    (1/(jwC) >> jwL for single conductor without return path)
  - For multi-conductor systems, capacitive effects are more reasonable
  - full and stabilized should agree (both correct, different numerics)

Part of Radia project
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

from ngbem_peec import NGBEMPEECSolver, create_plate_mesh

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7
EPS_0 = 8.854187817e-12
C_0 = 1.0 / np.sqrt(MU_0 * EPS_0)


def main():
    print("=" * 70)
    print("  Loop-Star Mode Comparison Test")
    print("  (mqs vs full vs stabilized)")
    print("=" * 70)

    # --- Geometry: 100mm x 10mm plate conductor ---
    width = 0.1      # 100mm
    height = 0.01    # 10mm
    maxh = 0.008     # 8mm mesh
    sigma = 5.8e7    # Copper
    thickness = 1e-3  # 1mm

    print(f"\n  Geometry: {width*1e3:.0f}mm x {height*1e3:.0f}mm plate")
    print(f"  Thickness: {thickness*1e3:.1f}mm, Sigma: {sigma:.1e} S/m")
    print(f"  Mesh: maxh={maxh*1e3:.1f}mm")

    # --- Create mesh and solver ---
    mesh = create_plate_mesh(width, height, maxh, label="conductor")
    solver = NGBEMPEECSolver(mesh, conductor_label="conductor",
                              sigma=sigma, thickness=thickness,
                              order=0, intorder=5)

    print("\n  Assembling matrices...")
    solver.assemble()
    solver.print_summary()

    # --- Diagnostic: check matrix scales ---
    print("  Matrix Scale Analysis:")
    P_inv_M = solver._P_inv_M_LS
    if P_inv_M is not None:
        cap_matrix = solver.M_LS.T @ P_inv_M
        cap_eigs = np.linalg.eigvalsh(cap_matrix.real)
        L_eigs = np.linalg.eigvalsh(solver.L.real)
        print(f"    L eigenvalues: [{L_eigs[0]:.4e}, {L_eigs[-1]:.4e}]")
        print(f"    M^T*P^-1*M eigenvalues: [{cap_eigs[0]:.4e}, {cap_eigs[-1]:.4e}]")
        ratio = cap_eigs[-1] / L_eigs[-1] if L_eigs[-1] > 0 else float('inf')
        print(f"    Ratio (cap/L max eig): {ratio:.1f}")
        if ratio > 100:
            print(f"    -> Cap correction >> L: Schur complement is ill-scaled")
            print(f"    -> Expected: MQS != full/stabilized (single conductor)")
    print()

    if solver.Q_0 is not None:
        print(f"    Q_0 vs M_LS (operator structure comparison):")
        print(f"      M_LS: sparse, entries in {{-1,0,+1}}, "
              f"nnz={np.count_nonzero(np.abs(solver.M_LS) > 1e-15)}")
        print(f"      Q_0:  dense, entries ~ {np.mean(np.abs(solver.Q_0)):.4e}, "
              f"nnz={np.count_nonzero(np.abs(solver.Q_0) > 1e-15)}")
        print()

    # --- Frequency sweep ---
    freqs = np.array([100, 1e3, 1e4, 1e5, 1e6, 1e7])

    print("  Solving for each mode...")
    Z_mqs = solver.solve_frequency(freqs, mode='mqs')
    Z_full = solver.solve_frequency(freqs, mode='full')
    Z_stab = solver.solve_frequency(freqs, mode='stabilized')

    R_mqs = np.real(Z_mqs)
    R_full = np.real(Z_full)
    R_stab = np.real(Z_stab)
    L_mqs = np.imag(Z_mqs) / (2 * np.pi * freqs)
    L_full = np.imag(Z_full) / (2 * np.pi * freqs)
    L_stab = np.imag(Z_stab) / (2 * np.pi * freqs)

    # --- Results table ---
    print("\n" + "=" * 90)
    print("  Inductance [nH]")
    print(f"  {'Freq':>10s} | {'kappa':>10s} | "
          f"{'L_mqs':>10s} | {'L_full':>10s} | {'L_stab':>10s} | "
          f"{'full/stab':>10s}")
    print("  " + "-" * 70)

    for i, f in enumerate(freqs):
        kappa = 2 * np.pi * f / C_0
        ratio_fs = L_full[i] / L_stab[i] if abs(L_stab[i]) > 0 else float('inf')
        print(f"  {f:10.0e} | {kappa:10.3e} | "
              f"{L_mqs[i]*1e9:10.2f} | {L_full[i]*1e9:10.2f} | "
              f"{L_stab[i]*1e9:10.2f} | "
              f"{ratio_fs:10.4f}")

    # --- Resistance comparison ---
    print(f"\n  Resistance [mOhm]")
    print(f"  {'Freq':>10s} | "
          f"{'R_mqs':>12s} | {'R_full':>12s} | {'R_stab':>12s}")
    print("  " + "-" * 52)
    for i, f in enumerate(freqs):
        print(f"  {f:10.0e} | "
              f"{R_mqs[i]*1e3:12.4f} | {R_full[i]*1e3:12.4f} | "
              f"{R_stab[i]*1e3:12.4f}")

    # --- Validation ---
    print("\n" + "=" * 70)
    print("  Validation")
    print("=" * 70)
    n_pass = 0
    n_total = 0

    # Test 1: full vs stabilized should agree (both include cap effects)
    max_diff_fs = np.max(
        np.abs(L_full - L_stab) / (np.abs(L_stab) + 1e-30)) * 100
    n_total += 1
    if max_diff_fs < 5:
        print(f"  PASS: full vs stabilized L: max {max_diff_fs:.2f}% (< 5%)")
        n_pass += 1
    else:
        print(f"  FAIL: full vs stabilized L: max {max_diff_fs:.2f}% (>= 5%)")

    # Test 2: R should agree between full and stabilized
    max_R_diff = np.max(
        np.abs(R_full - R_stab) / (np.abs(R_stab) + 1e-30)) * 100
    n_total += 1
    if max_R_diff < 5:
        print(f"  PASS: full vs stabilized R: max {max_R_diff:.2f}% (< 5%)")
        n_pass += 1
    else:
        print(f"  FAIL: full vs stabilized R: max {max_R_diff:.2f}% (>= 5%)")

    # Test 3: All L values positive
    n_total += 1
    if np.all(L_mqs > 0) and np.all(L_full > 0) and np.all(L_stab > 0):
        print(f"  PASS: All L values positive")
        n_pass += 1
    else:
        print(f"  FAIL: Negative L values detected")

    # Test 4: L_full >= L_mqs (cap adds impedance for single conductor)
    n_total += 1
    if np.all(L_full >= L_mqs * 0.95):
        print(f"  PASS: L_full >= L_mqs (cap correction adds impedance)")
        n_pass += 1
    else:
        print(f"  FAIL: L_full < L_mqs at some frequencies")

    # Test 5: MQS inductance is physically reasonable
    # For 100mm strip: L ~ 50-150 nH expected
    n_total += 1
    L_mqs_low = L_mqs[0] * 1e9
    if 10 < L_mqs_low < 500:
        print(f"  PASS: L_mqs(100 Hz) = {L_mqs_low:.1f} nH (reasonable)")
        n_pass += 1
    else:
        print(f"  FAIL: L_mqs(100 Hz) = {L_mqs_low:.1f} nH (unexpected)")

    print(f"\n  Result: {n_pass}/{n_total} tests passed")

    # --- Physics interpretation ---
    print("\n" + "=" * 70)
    print("  Physics Interpretation")
    print("=" * 70)
    print("""
  MQS vs Full/Stabilized difference is EXPECTED for single conductor:
  - Single plate has very low self-capacitance C ~ eps_0 * Area / d_eff
  - Capacitive impedance 1/(jwC) >> jwL at moderate frequencies
  - This ADDS to total impedance, making L_full > L_mqs
  - For multi-conductor systems (with return path), C is larger
    and the capacitive effect is more moderate

  Recommendation:
  - Use mode='mqs' for conductor impedance extraction (standard PEEC)
  - Use mode='full' or 'stabilized' when capacitive coupling matters
    (e.g., interconnect analysis with parasitic capacitance)
  - Full and stabilized should agree; prefer 'stabilized' for low freq
""")

    return n_pass == n_total


if __name__ == '__main__':
    main()
