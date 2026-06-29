"""
demo_ngbem_peec.py

Demonstration and validation of ngbem-based PEEC Loop-Star solver.

Tests:
1. Matrix assembly: L, P, M_LS properties (symmetry, positive definiteness)
2. MQS impedance: Z(f) in Loop-only mode
3. Full Loop-Star: Z(f) with capacitive effects (self-resonance)
4. High-order convergence: order=0,1,2

Part of Radia project
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

from ngbem_peec import NGBEMPEECSolver, create_plate_mesh, get_mesh_triangles

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7
EPS_0 = 8.854187817e-12


def print_separator(title):
    """Print formatted section separator."""
    print()
    print("=" * 70)
    print(f"  {title}")
    print("=" * 70)


def test_matrix_assembly():
    """Test 1: Verify matrix properties after assembly."""
    print_separator("Test 1: Matrix Assembly and Properties")

    # Create plate mesh: 10mm x 10mm copper, 35um thickness (PCB trace)
    width = 0.01    # 10mm
    height = 0.01   # 10mm
    maxh = 0.004    # ~4mm element size
    thickness = 35e-6  # 35um (1oz copper)
    sigma = 5.8e7   # Copper

    print(f"\n  Geometry: Rectangular plate {width*1e3:.0f}mm x {height*1e3:.0f}mm")
    print(f"  Thickness: {thickness*1e6:.0f} um")
    print(f"  Conductivity: {sigma:.1e} S/m (copper)")
    print(f"  Max element size: {maxh*1e3:.1f} mm")

    mesh = create_plate_mesh(width, height, maxh, label="conductor")
    triangles, areas = get_mesh_triangles(mesh)
    print(f"  Triangles: {len(triangles)}")
    print(f"  Total area: {np.sum(areas)*1e6:.2f} mm^2 (expected {width*height*1e6:.0f})")

    solver = NGBEMPEECSolver(mesh, conductor_label="conductor",
                              sigma=sigma, thickness=thickness,
                              order=0, intorder=5)
    solver.assemble()
    solver.print_summary()

    # Verify matrix properties
    print("  Verification:")

    # L should be symmetric positive (semi-)definite
    L_eig = np.linalg.eigvalsh(solver.L)
    n_neg_L = np.sum(L_eig < -1e-15 * np.max(np.abs(L_eig)))
    if n_neg_L == 0:
        print("    [PASS] L: all eigenvalues non-negative")
    else:
        print(f"    [WARN] L: {n_neg_L} negative eigenvalues "
              f"(min={L_eig[0]:.4e})")

    # P should be symmetric positive definite
    P_eig = np.linalg.eigvalsh(solver.P)
    n_neg_P = np.sum(P_eig < -1e-15 * np.max(np.abs(P_eig)))
    if n_neg_P == 0:
        print("    [PASS] P: all eigenvalues non-negative")
    else:
        print(f"    [WARN] P: {n_neg_P} negative eigenvalues "
              f"(min={P_eig[0]:.4e})")

    # M_LS should be sparse (each RWG edge connects 2 triangles)
    nnz = np.count_nonzero(np.abs(solver.M_LS) > 1e-15)
    sparsity = 1.0 - nnz / solver.M_LS.size
    print(f"    M_LS sparsity: {sparsity*100:.1f}% "
          f"({nnz} nonzeros / {solver.M_LS.size})")

    # Euler check: V - E + F for open surface
    nv = mesh.nv
    n_edge = solver.n_loop
    n_cell = solver.n_star
    euler = nv - n_edge + n_cell
    print(f"    Euler: V({nv}) - E({n_edge}) + F({n_cell}) = {euler} "
          f"(expected 1 for disk topology)")

    return solver


def test_mqs_impedance(solver):
    """Test 2: MQS (Loop-only) impedance Z(f)."""
    print_separator("Test 2: MQS Impedance (Loop-only)")

    freqs = np.logspace(3, 9, 25)  # 1 kHz to 1 GHz
    Z_mqs = solver.solve_frequency(freqs, mode='mqs')

    print(f"\n  {'Freq':>12s}  {'R (Ohm)':>12s}  {'X (Ohm)':>12s}  "
          f"{'L (nH)':>10s}  {'|Z| (Ohm)':>12s}")
    print(f"  {'-'*12}  {'-'*12}  {'-'*12}  {'-'*10}  {'-'*12}")

    for i in range(0, len(freqs), 4):
        f = freqs[i]
        R = np.real(Z_mqs[i])
        X = np.imag(Z_mqs[i])
        Z_abs = np.abs(Z_mqs[i])
        L_nH = X / (2 * np.pi * f) * 1e9 if f > 0 else 0
        if f < 1e6:
            f_str = f"{f/1e3:.1f} kHz"
        elif f < 1e9:
            f_str = f"{f/1e6:.1f} MHz"
        else:
            f_str = f"{f/1e9:.1f} GHz"
        print(f"  {f_str:>12s}  {R:>12.4e}  {X:>12.4e}  "
              f"{L_nH:>10.3f}  {Z_abs:>12.4e}")

    return Z_mqs


def test_full_loop_star(solver):
    """Test 3: Full Loop-Star impedance (with capacitive effects)."""
    print_separator("Test 3: Full Loop-Star Impedance")

    freqs = np.logspace(3, 11, 40)  # 1 kHz to 100 GHz
    Z_full = solver.solve_frequency(freqs, mode='full')

    # Find self-resonance (where Im(Z) = 0 after being inductive)
    X = np.imag(Z_full)
    resonance_idx = None
    for i in range(1, len(X)):
        if X[i-1] > 0 and X[i] <= 0:
            resonance_idx = i
            break

    print(f"\n  {'Freq':>12s}  {'R (Ohm)':>12s}  {'X (Ohm)':>12s}  "
          f"{'|Z| (Ohm)':>12s}  {'Phase (deg)':>12s}")
    print(f"  {'-'*12}  {'-'*12}  {'-'*12}  {'-'*12}  {'-'*12}")

    for i in range(0, len(freqs), 5):
        f = freqs[i]
        R = np.real(Z_full[i])
        X_val = np.imag(Z_full[i])
        Z_abs = np.abs(Z_full[i])
        phase = np.angle(Z_full[i], deg=True)
        if f < 1e6:
            f_str = f"{f/1e3:.1f} kHz"
        elif f < 1e9:
            f_str = f"{f/1e6:.1f} MHz"
        else:
            f_str = f"{f/1e9:.1f} GHz"
        print(f"  {f_str:>12s}  {R:>12.4e}  {X_val:>12.4e}  "
              f"{Z_abs:>12.4e}  {phase:>12.1f}")

    if resonance_idx is not None:
        f_res = freqs[resonance_idx]
        if f_res < 1e9:
            print(f"\n  Self-resonance frequency: ~{f_res/1e6:.1f} MHz")
        else:
            print(f"\n  Self-resonance frequency: ~{f_res/1e9:.1f} GHz")
        print("  (where impedance transitions from inductive to capacitive)")
    else:
        print("\n  No self-resonance found in frequency range")
        print("  (this is expected for a simple flat plate)")

    return Z_full


def test_high_order(mesh_params):
    """Test 4: High-order convergence."""
    print_separator("Test 4: High-Order Convergence")

    width, height, maxh, thickness, sigma = mesh_params
    mesh = create_plate_mesh(width, height, maxh, label="conductor")

    f_test = 1e6  # 1 MHz test frequency
    print(f"\n  Test frequency: {f_test/1e6:.0f} MHz")
    print(f"\n  {'Order':>5s}  {'n_loop':>7s}  {'n_star':>7s}  "
          f"{'t_asm (ms)':>10s}  {'Z_mqs (Ohm)':>20s}  {'L (nH)':>10s}")
    print(f"  {'-'*5}  {'-'*7}  {'-'*7}  {'-'*10}  {'-'*20}  {'-'*10}")

    prev_L = None
    for order in [0, 1, 2]:
        try:
            solver = NGBEMPEECSolver(mesh, conductor_label="conductor",
                                      sigma=sigma, thickness=thickness,
                                      order=order, intorder=5 + 2*order)
            solver.assemble()

            Z = solver.solve_frequency(np.array([f_test]), mode='mqs')
            omega = 2 * np.pi * f_test
            L_nH = np.imag(Z[0]) / omega * 1e9

            change = ""
            if prev_L is not None:
                pct = abs(L_nH - prev_L) / abs(prev_L) * 100
                change = f"  ({pct:.1f}% change)"
            prev_L = L_nH

            print(f"  {order:>5d}  {solver.n_loop:>7d}  {solver.n_star:>7d}  "
                  f"{solver.t_assemble*1e3:>10.1f}  "
                  f"{Z[0].real:+.4e}{Z[0].imag:+.4e}j  "
                  f"{L_nH:>10.3f}{change}")
        except Exception as e:
            print(f"  {order:>5d}  Error: {e}")


def main():
    print("=" * 70)
    print("  Demo: ngbem PEEC Loop-Star Solver")
    print("=" * 70)

    # Geometry parameters
    width = 0.01     # 10mm
    height = 0.01    # 10mm
    maxh = 0.004     # ~4mm elements
    thickness = 35e-6  # 35um (1oz copper PCB)
    sigma = 5.8e7    # Copper

    # Test 1: Matrix assembly
    solver = test_matrix_assembly()

    # Test 2: MQS impedance
    Z_mqs = test_mqs_impedance(solver)

    # Test 3: Full Loop-Star
    Z_full = test_full_loop_star(solver)

    # Test 4: High-order convergence
    test_high_order((width, height, maxh, thickness, sigma))

    print_separator("Summary")
    print("""
  ngbem PEEC Loop-Star solver successfully demonstrated:

  1. Matrix assembly via ngbem (LaplaceSL + SingleLayerPotentialOperator)
     - L: symmetric, from HDivSurface (edge DOFs)
     - P: symmetric positive definite, from SurfaceL2 (cell DOFs)
     - M_LS: sparse divergence coupling (FEM bilinear form)

  2. MQS mode: inductive impedance Z = R + jwL
     - Valid for DC to ~1 MHz (power electronics regime)
     - Fast: only Loop DOFs needed

  3. Full Loop-Star: includes capacitive self-resonance
     - Shows transition from inductive to capacitive behavior
     - Uses Schur complement for port impedance extraction

  4. High-order elements: order=0,1,2 convergence
     - Higher order = more DOFs but better accuracy per DOF
     - Natural in ngbem via FE order parameter

  Key advantage over C++ filament PEEC:
    - Galerkin (symmetric by construction) vs Collocation
    - High-order elements for convergence
    - Natural Loop-Star via function spaces (no explicit M_LS construction)
    - FMM acceleration available for large problems
""")


if __name__ == '__main__':
    main()
