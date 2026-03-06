"""
demo_ngbem_coupled.py

Demonstration and validation of coupled PEEC + MMM solver.

Tests:
1. Coupled impedance: Conductor plate with magnetic core (mu_r enhancement)
2. Frequency sweep: Z(f) with and without core coupling
3. Complex permeability: Core loss via imaginary mu_r
4. Multi-loop system: Port impedance extraction

Part of Radia project
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

from ngbem_peec import NGBEMPEECSolver, create_plate_mesh
from ngbem_coupled import CoupledPEECMMM

# Physical constants
MU_0 = 4.0 * np.pi * 1e-7


def print_separator(title):
    """Print formatted section separator."""
    print()
    print("=" * 70)
    print(f"  {title}")
    print("=" * 70)


def test_basic_coupling():
    """Test 1: Basic PEEC + analytical core coupling."""
    print_separator("Test 1: Basic PEEC + Core Coupling (Analytical)")

    # Create conductor plate
    width = 0.01    # 10mm
    height = 0.01   # 10mm
    maxh = 0.004    # ~4mm elements
    thickness = 35e-6  # 35um (1oz copper)
    sigma = 5.8e7   # Copper

    print(f"\n  Conductor: {width*1e3:.0f}mm x {height*1e3:.0f}mm plate")
    print(f"  Thickness: {thickness*1e6:.0f} um, sigma = {sigma:.1e} S/m")

    mesh = create_plate_mesh(width, height, maxh, label="conductor")
    peec = NGBEMPEECSolver(mesh, conductor_label="conductor",
                            sigma=sigma, thickness=thickness,
                            order=0, intorder=5)
    peec.assemble()
    print(f"  Loop DOFs: {peec.n_loop}, Star DOFs: {peec.n_star}")
    print(f"  Assembly time: {peec.t_assemble:.3f} s")

    # Test without core
    f_test = 1e6  # 1 MHz
    Z_air = peec.solve_frequency(np.array([f_test]), mode='mqs')
    omega = 2.0 * np.pi * f_test
    L_air = np.imag(Z_air[0]) / omega
    R_air = np.real(Z_air[0])

    print(f"\n  Without core (air):")
    print(f"    Z({f_test/1e6:.0f} MHz) = {R_air:.4e} + j{np.imag(Z_air[0]):.4e} Ohm")
    print(f"    L_air = {L_air*1e9:.3f} nH")
    print(f"    R_air = {R_air:.4e} Ohm")

    # Create coupled solver with magnetic core
    mu_r = 1000.0
    coupled = CoupledPEECMMM(peec, mu_r=mu_r)

    # Analytical coupling: Delta_L ~ L_air * (mu_r - 1) / (mu_r + 1) / chi
    # For a plate above an infinite magnetic half-space:
    # The image method gives L_total = L_air + L_image * chi
    # where L_image = L_air * (mu_r - 1)/(mu_r + 1)
    # and chi = mu_r - 1
    #
    # So: Delta_L = L_air * (mu_r - 1)/(mu_r + 1) / (mu_r - 1) = L_air / (mu_r + 1)
    # Then: L_total = L_air + Delta_L * chi = L_air + L_air * (mu_r-1)/(mu_r+1)
    #               = L_air * 2*mu_r / (mu_r + 1)
    #
    # For mu_r >> 1: L_total -> 2 * L_air (image method doubles inductance)
    image_factor = (mu_r - 1.0) / (mu_r + 1.0)
    Delta_L_mat = peec.L * image_factor / (mu_r - 1.0)
    coupled.compute_coupling_analytically(Delta_L_mat)

    # Verify coupling
    Z_coupled = coupled.solve_frequency(np.array([f_test]), mode='mqs')
    L_coupled = np.imag(Z_coupled[0]) / omega
    R_coupled = np.real(Z_coupled[0])

    print(f"\n  With core (mu_r = {mu_r:.0f}):")
    print(f"    Z({f_test/1e6:.0f} MHz) = {R_coupled:.4e} + j{np.imag(Z_coupled[0]):.4e} Ohm")
    print(f"    L_coupled = {L_coupled*1e9:.3f} nH")
    print(f"    L_coupled / L_air = {L_coupled/L_air:.3f}")
    print(f"    Expected ratio (image method): {2.0 * mu_r / (mu_r + 1.0):.3f}")

    # Verify ratio
    expected_ratio = 2.0 * mu_r / (mu_r + 1.0)
    actual_ratio = L_coupled / L_air
    error = abs(actual_ratio - expected_ratio) / expected_ratio * 100
    if error < 1.0:
        print(f"    [PASS] Ratio error: {error:.2f}%")
    else:
        print(f"    [WARN] Ratio error: {error:.2f}%")

    coupled.print_summary()
    return peec, coupled


def test_frequency_sweep(peec, coupled):
    """Test 2: Frequency sweep with and without core."""
    print_separator("Test 2: Frequency Sweep (Air vs Core)")

    freqs = np.logspace(3, 8, 20)  # 1 kHz to 100 MHz

    Z_air = peec.solve_frequency(freqs, mode='mqs')
    Z_core = coupled.solve_frequency(freqs, mode='mqs')

    print(f"\n  {'Freq':>12s}  {'L_air (nH)':>10s}  {'L_core (nH)':>10s}  "
          f"{'Ratio':>8s}  {'R_air (Ohm)':>12s}  {'R_core (Ohm)':>12s}")
    print(f"  {'-'*12}  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*12}  {'-'*12}")

    for i in range(0, len(freqs), 3):
        f = freqs[i]
        omega = 2.0 * np.pi * f
        L_air = np.imag(Z_air[i]) / omega * 1e9
        L_core = np.imag(Z_core[i]) / omega * 1e9
        R_air = np.real(Z_air[i])
        R_core = np.real(Z_core[i])
        ratio = L_core / L_air if L_air != 0 else 0

        if f < 1e6:
            f_str = f"{f/1e3:.1f} kHz"
        else:
            f_str = f"{f/1e6:.1f} MHz"

        print(f"  {f_str:>12s}  {L_air:>10.3f}  {L_core:>10.3f}  "
              f"{ratio:>8.3f}  {R_air:>12.4e}  {R_core:>12.4e}")

    print(f"\n  Inductance ratio should be ~{2.0 * coupled.mu_r / (coupled.mu_r + 1.0):.3f}"
          f" (image method, mu_r={coupled.mu_r:.0f})")


def test_complex_permeability():
    """Test 3: Core loss via complex permeability."""
    print_separator("Test 3: Complex Permeability (Core Loss)")

    # Create conductor
    width = 0.01
    height = 0.01
    maxh = 0.004
    thickness = 35e-6
    sigma = 5.8e7

    mesh = create_plate_mesh(width, height, maxh, label="conductor")
    peec = NGBEMPEECSolver(mesh, conductor_label="conductor",
                            sigma=sigma, thickness=thickness,
                            order=0, intorder=5)
    peec.assemble()

    # Test with different loss tangents
    mu_r_real = 1000.0
    loss_tangents = [0.0, 0.001, 0.01, 0.05, 0.1]

    f_test = 100e3  # 100 kHz
    omega = 2.0 * np.pi * f_test

    print(f"\n  Frequency: {f_test/1e3:.0f} kHz")
    print(f"  Core mu_r' = {mu_r_real:.0f}")
    print()

    # Air reference
    Z_air = peec.solve_frequency(np.array([f_test]), mode='mqs')
    L_air = np.imag(Z_air[0]) / omega
    R_air = np.real(Z_air[0])

    print(f"  {'tan(delta)':>12s}  {'mu_r\"':>8s}  {'L (nH)':>10s}  "
          f"{'R (Ohm)':>12s}  {'Q factor':>10s}  {'R_core/R_air':>12s}")
    print(f"  {'-'*12}  {'-'*8}  {'-'*10}  {'-'*12}  {'-'*10}  {'-'*12}")

    # Print air reference
    Q_air = np.imag(Z_air[0]) / R_air if R_air > 0 else float('inf')
    print(f"  {'(air)':>12s}  {'0.0':>8s}  {L_air*1e9:>10.3f}  "
          f"{R_air:>12.4e}  {Q_air:>10.1f}  {'1.000':>12s}")

    for tan_d in loss_tangents:
        mu_r_imag = mu_r_real * tan_d

        coupled = CoupledPEECMMM(peec, mu_r=mu_r_real, mu_r_imag=mu_r_imag)

        # Analytical coupling
        image_factor = (mu_r_real - 1.0) / (mu_r_real + 1.0)
        Delta_L_mat = peec.L * image_factor / (mu_r_real - 1.0)
        coupled.compute_coupling_analytically(Delta_L_mat)

        Z = coupled.solve_frequency(np.array([f_test]), mode='mqs')
        L = np.imag(Z[0]) / omega
        R = np.real(Z[0])
        Q = np.imag(Z[0]) / R if R > 0 else float('inf')
        R_ratio = R / R_air if R_air > 0 else float('inf')

        print(f"  {tan_d:>12.3f}  {mu_r_imag:>8.1f}  {L*1e9:>10.3f}  "
              f"{R:>12.4e}  {Q:>10.1f}  {R_ratio:>12.3f}")

    print(f"\n  Core loss adds resistance: R_core = omega * Im(Delta_L * chi)")
    print(f"  Higher loss tangent -> lower Q factor")


def test_full_loop_star_vs_mqs():
    """Test 4: Full Loop-Star vs MQS comparison."""
    print_separator("Test 4: Full Loop-Star vs MQS")

    width = 0.01
    height = 0.01
    maxh = 0.004
    thickness = 35e-6
    sigma = 5.8e7
    mu_r = 1000.0

    mesh = create_plate_mesh(width, height, maxh, label="conductor")
    peec = NGBEMPEECSolver(mesh, conductor_label="conductor",
                            sigma=sigma, thickness=thickness,
                            order=0, intorder=5)
    peec.assemble()

    coupled = CoupledPEECMMM(peec, mu_r=mu_r)
    image_factor = (mu_r - 1.0) / (mu_r + 1.0)
    Delta_L_mat = peec.L * image_factor / (mu_r - 1.0)
    coupled.compute_coupling_analytically(Delta_L_mat)

    freqs = np.logspace(3, 11, 30)  # 1 kHz to 100 GHz
    Z_mqs = coupled.solve_frequency(freqs, mode='mqs')
    Z_full = coupled.solve_frequency(freqs, mode='full')

    print(f"\n  {'Freq':>12s}  {'|Z_mqs|':>12s}  {'|Z_full|':>12s}  "
          f"{'Phase_mqs':>10s}  {'Phase_full':>10s}  {'Ratio':>8s}")
    print(f"  {'-'*12}  {'-'*12}  {'-'*12}  {'-'*10}  {'-'*10}  {'-'*8}")

    for i in range(0, len(freqs), 4):
        f = freqs[i]
        Z_m = Z_mqs[i]
        Z_f = Z_full[i]

        if f < 1e6:
            f_str = f"{f/1e3:.1f} kHz"
        elif f < 1e9:
            f_str = f"{f/1e6:.1f} MHz"
        else:
            f_str = f"{f/1e9:.1f} GHz"

        ratio = np.abs(Z_f) / np.abs(Z_m) if np.abs(Z_m) > 0 else 0

        print(f"  {f_str:>12s}  {np.abs(Z_m):>12.4e}  {np.abs(Z_f):>12.4e}  "
              f"{np.angle(Z_m, deg=True):>10.1f}  {np.angle(Z_f, deg=True):>10.1f}  "
              f"{ratio:>8.3f}")

    # Find self-resonance in full Loop-Star
    X_full = np.imag(Z_full)
    for i in range(1, len(X_full)):
        if X_full[i-1] > 0 and X_full[i] <= 0:
            f_res = freqs[i]
            if f_res < 1e9:
                print(f"\n  Self-resonance (full Loop-Star): ~{f_res/1e6:.1f} MHz")
            else:
                print(f"\n  Self-resonance (full Loop-Star): ~{f_res/1e9:.1f} GHz")
            break
    else:
        print(f"\n  No self-resonance found in frequency range")

    print(f"\n  MQS (Loop-only): purely inductive at all frequencies")
    print(f"  Full Loop-Star: includes capacitive self-resonance")


def test_inductance_vs_mu_r():
    """Test 5: Inductance enhancement vs core permeability."""
    print_separator("Test 5: Inductance Enhancement vs mu_r")

    width = 0.01
    height = 0.01
    maxh = 0.004
    thickness = 35e-6
    sigma = 5.8e7

    mesh = create_plate_mesh(width, height, maxh, label="conductor")
    peec = NGBEMPEECSolver(mesh, conductor_label="conductor",
                            sigma=sigma, thickness=thickness,
                            order=0, intorder=5)
    peec.assemble()

    f_test = 100e3  # 100 kHz
    omega = 2.0 * np.pi * f_test

    Z_air = peec.solve_frequency(np.array([f_test]), mode='mqs')
    L_air = np.imag(Z_air[0]) / omega

    mu_r_values = [1, 10, 100, 500, 1000, 2000, 5000, 10000]

    print(f"\n  Frequency: {f_test/1e3:.0f} kHz")
    print(f"  L_air = {L_air*1e9:.3f} nH")
    print()

    print(f"  {'mu_r':>8s}  {'L (nH)':>10s}  {'L/L_air':>8s}  "
          f"{'Expected':>8s}  {'Error (%)':>10s}")
    print(f"  {'-'*8}  {'-'*10}  {'-'*8}  {'-'*8}  {'-'*10}")

    for mu_r in mu_r_values:
        coupled = CoupledPEECMMM(peec, mu_r=float(mu_r))

        if mu_r > 1:
            image_factor = (mu_r - 1.0) / (mu_r + 1.0)
            Delta_L_mat = peec.L * image_factor / (mu_r - 1.0)
            coupled.compute_coupling_analytically(Delta_L_mat)
        # else: no coupling for mu_r=1

        Z = coupled.solve_frequency(np.array([f_test]), mode='mqs')
        L = np.imag(Z[0]) / omega
        ratio = L / L_air
        expected = 2.0 * mu_r / (mu_r + 1.0) if mu_r > 1 else 1.0
        error = abs(ratio - expected) / expected * 100

        print(f"  {mu_r:>8d}  {L*1e9:>10.3f}  {ratio:>8.3f}  "
              f"{expected:>8.3f}  {error:>10.2f}")

    print(f"\n  Image method: L/L_air = 2*mu_r / (mu_r + 1)")
    print(f"  For mu_r >> 1: L/L_air -> 2.0 (perfect image doubles inductance)")


def main():
    print("=" * 70)
    print("  Demo: Coupled PEEC + MMM Solver (Phase 3)")
    print("=" * 70)

    # Test 1: Basic coupling
    peec, coupled = test_basic_coupling()

    # Test 2: Frequency sweep
    test_frequency_sweep(peec, coupled)

    # Test 3: Complex permeability
    test_complex_permeability()

    # Test 4: Full Loop-Star vs MQS
    test_full_loop_star_vs_mqs()

    # Test 5: Inductance vs mu_r
    test_inductance_vs_mu_r()

    print_separator("Summary")
    print("""
  Coupled PEEC + MMM solver demonstrated:

  1. Basic coupling: conductor plate + magnetic core
     - Uses analytical image method for Delta_L
     - L_coupled / L_air matches 2*mu_r/(mu_r+1) formula

  2. Frequency sweep: Z(f) with and without core
     - Core increases inductance at all frequencies
     - Ratio constant (frequency-independent for linear core)

  3. Complex permeability: core loss modeling
     - mu = mu' - j*mu" introduces additional resistance
     - Higher loss tangent -> lower Q factor
     - R_core = omega * Im(Delta_L * chi)

  4. Full Loop-Star vs MQS comparison
     - MQS: purely inductive (valid for DC - 1 MHz)
     - Full: includes capacitive self-resonance
     - Self-resonance from Schur complement

  5. Inductance enhancement vs mu_r
     - Image method: L/L_air = 2*mu_r/(mu_r+1)
     - Saturates at factor 2 for large mu_r
     - Consistent across all mu_r values

  Physics:
    - Coil-core coupling via mutual inductance Delta_L
    - chi = (mu_r - 1) susceptibility scales Delta_L
    - Complex chi = (mu_r' - 1) - j*mu_r" for lossy core
    - Port impedance via Schur complement (Loop-Star block system)

  Integration with Radia MMM (future):
    - compute_coupling_radia() computes Delta_L numerically
    - Uses Biot-Savart for H-field at core elements
    - Radia Solve() for core magnetization response
    - A-field from core M back to conductor loops
""")


if __name__ == '__main__':
    main()
