"""
validate_shield_delta_r.py

Validation of PEEC shield conductor modeling.

Test: Square loop coil above an aluminum shield plate.
Expected behavior:
  - L decreases with shield (field screening)
  - R increases with shield (eddy current losses)
  - Effect increases with frequency

Physics:
  Shield eddy currents oppose the coil field -> reduced flux linkage -> lower L.
  Eddy current losses dissipate power -> reflected resistance -> higher R.
  At low frequency, skin depth >> shield -> weak shielding.
  At high frequency, skin depth << shield -> strong shielding.

Author: Sugahara
Date: 2026-02-17
"""

import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

from peec_matrices import PEECBuilder
from peec_topology import PEECCircuitSolver
from peec_shield import (add_shield_mesh, add_coil_from_segments,
                         compute_shield_effect)


MU_0 = 4 * np.pi * 1e-7
SIGMA_AL = 3.5e7  # S/m aluminum
SIGMA_CU = 5.8e7  # S/m copper


def create_square_loop(side_m, z_m=0.0, n_seg_per_side=10):
    """Create a square loop coil as segment arrays.

    Args:
        side_m: Side length [m]
        z_m: Z position [m]
        n_seg_per_side: Segments per side

    Returns:
        p1, p2: (N, 3) arrays of segment endpoints [m]
    """
    half = side_m / 2
    # Four sides: bottom, right, top, left
    corners = [
        [-half, -half, z_m],
        [+half, -half, z_m],
        [+half, +half, z_m],
        [-half, +half, z_m],
    ]

    p1_list = []
    p2_list = []
    for ic in range(4):
        c1 = np.array(corners[ic])
        c2 = np.array(corners[(ic + 1) % 4])
        for i in range(n_seg_per_side):
            t1 = i / n_seg_per_side
            t2 = (i + 1) / n_seg_per_side
            p1_list.append(c1 + t1 * (c2 - c1))
            p2_list.append(c1 + t2 * (c2 - c1))

    return np.array(p1_list), np.array(p2_list)


def test_basic_shield():
    """Test 1: Basic shield effect - L decreases, R increases."""
    print("=" * 60)
    print("Test 1: Basic shield effect (square loop + Al plate)")
    print("=" * 60)

    # Square loop coil: 50mm side, at z=0
    side = 0.05  # 50 mm
    wire_w = 1e-3
    wire_h = 1e-3
    p1, p2 = create_square_loop(side, z_m=0.0, n_seg_per_side=8)
    n_coil = len(p1)
    print(f"  Coil: {side*1e3:.0f}mm square loop, {n_coil} segments")

    # Shield: 80mm x 80mm Al plate at z=-5mm, 1mm thick
    shield_z = -0.005
    shield_w = 0.08
    shield_t = 1e-3
    nx, ny = 6, 6
    print(f"  Shield: {shield_w*1e3:.0f}x{shield_w*1e3:.0f}mm Al, "
          f"t={shield_t*1e3:.1f}mm, z={shield_z*1e3:.1f}mm")
    print(f"  Shield mesh: {nx}x{ny}")

    # Frequency sweep
    freqs = np.logspace(2, 6, 20)  # 100 Hz to 1 MHz

    t0 = time.time()
    result = compute_shield_effect(
        p1, p2, wire_w, wire_h,
        shield_center=[0, 0, shield_z],
        shield_size_x=shield_w, shield_size_y=shield_w,
        shield_thickness=shield_t, shield_sigma=SIGMA_AL,
        shield_nx=nx, shield_ny=ny,
        freqs=freqs, coil_sigma=SIGMA_CU)
    t_elapsed = time.time() - t0

    print(f"  Segments: coil={result['n_coil_segments']}, "
          f"shield={result['n_shield_segments']}, "
          f"total={result['n_coil_segments'] + result['n_shield_segments']}")
    print(f"  Computation time: {t_elapsed:.1f} s")

    # Check key properties
    L0 = result['L_without'][0]  # Low-frequency L
    L_hi = result['L_with'][-1]  # High-frequency L with shield
    L_no_hi = result['L_without'][-1]  # High-frequency L without shield

    Delta_R_hi = result['Delta_R'][-1]
    Delta_L_hi = result['Delta_L'][-1]

    print(f"\n  Results:")
    print(f"  {'Freq [Hz]':>12s} {'L_no [nH]':>10s} {'L_with [nH]':>11s} "
          f"{'DeltaL [nH]':>11s} {'DeltaR [mOhm]':>13s}")
    print(f"  {'-'*60}")

    for i in [0, 5, 10, 15, -1]:
        f = result['freqs'][i]
        Ln = result['L_without'][i] * 1e9
        Lw = result['L_with'][i] * 1e9
        dL = result['Delta_L'][i] * 1e9
        dR = result['Delta_R'][i] * 1e3
        print(f"  {f:12.0f} {Ln:10.2f} {Lw:11.2f} {dL:11.2f} {dR:13.3f}")

    # Validation checks
    n_pass = 0
    n_total = 3

    # Check 1: L decreases with shield at high frequency
    passed = Delta_L_hi < 0
    status = "PASS" if passed else "FAIL"
    n_pass += passed
    print(f"\n  Check 1: L decreases with shield at high freq ... {status}")
    print(f"    Delta_L = {Delta_L_hi*1e9:.2f} nH (should be < 0)")

    # Check 2: R increases with shield (positive Delta_R)
    passed = Delta_R_hi > 0
    status = "PASS" if passed else "FAIL"
    n_pass += passed
    print(f"  Check 2: R increases with shield (Delta_R > 0)  ... {status}")
    print(f"    Delta_R = {Delta_R_hi*1e3:.3f} mOhm at {freqs[-1]/1e3:.0f} kHz")

    # Check 3: Shield effect increases with frequency
    Delta_R_lo = result['Delta_R'][0]
    passed = abs(Delta_R_hi) > abs(Delta_R_lo)
    status = "PASS" if passed else "FAIL"
    n_pass += passed
    print(f"  Check 3: Shield effect increases with frequency  ... {status}")
    print(f"    Delta_R(low)={Delta_R_lo*1e3:.4f} mOhm, "
          f"Delta_R(high)={Delta_R_hi*1e3:.3f} mOhm")

    print(f"\n  Result: {n_pass}/{n_total} PASSED")
    return n_pass == n_total


def test_shield_distance():
    """Test 2: Shield effect decreases with distance."""
    print("\n" + "=" * 60)
    print("Test 2: Shield effect vs distance")
    print("=" * 60)

    side = 0.05
    wire_w = 1e-3
    wire_h = 1e-3
    p1, p2 = create_square_loop(side, z_m=0.0, n_seg_per_side=8)

    shield_w = 0.08
    shield_t = 1e-3
    nx, ny = 5, 5

    freq = 100e3  # 100 kHz

    distances = [2e-3, 5e-3, 10e-3, 20e-3]  # 2, 5, 10, 20 mm
    delta_Ls = []
    delta_Rs = []

    for d in distances:
        result = compute_shield_effect(
            p1, p2, wire_w, wire_h,
            shield_center=[0, 0, -d],
            shield_size_x=shield_w, shield_size_y=shield_w,
            shield_thickness=shield_t, shield_sigma=SIGMA_AL,
            shield_nx=nx, shield_ny=ny,
            freqs=np.array([freq]), coil_sigma=SIGMA_CU)

        delta_Ls.append(result['Delta_L'][0])
        delta_Rs.append(result['Delta_R'][0])

    print(f"  Freq = {freq/1e3:.0f} kHz")
    print(f"  {'Distance [mm]':>14s} {'Delta_L [nH]':>12s} {'Delta_R [mOhm]':>14s}")
    print(f"  {'-'*44}")
    for d, dL, dR in zip(distances, delta_Ls, delta_Rs):
        print(f"  {d*1e3:14.1f} {dL*1e9:12.2f} {dR*1e3:14.3f}")

    # Check: closer shield -> larger effect
    passed = (abs(delta_Ls[0]) > abs(delta_Ls[-1]) and
              abs(delta_Rs[0]) > abs(delta_Rs[-1]))
    status = "PASS" if passed else "FAIL"
    print(f"\n  Check: closer shield has larger effect ... {status}")
    print(f"    |Delta_L(2mm)| = {abs(delta_Ls[0])*1e9:.2f} nH "
          f"> |Delta_L(20mm)| = {abs(delta_Ls[-1])*1e9:.2f} nH")

    return passed


def test_no_shield_baseline():
    """Test 3: Without shield, Delta_Z should be zero."""
    print("\n" + "=" * 60)
    print("Test 3: No-shield baseline (Delta_Z = 0)")
    print("=" * 60)

    side = 0.05
    wire_w = 1e-3
    wire_h = 1e-3
    p1, p2 = create_square_loop(side, z_m=0.0, n_seg_per_side=8)

    # Build model without shield
    builder = PEECBuilder()
    ports = add_coil_from_segments(builder, p1, p2, wire_w, wire_h, SIGMA_CU)
    builder.add_port(ports[0], ports[1])
    topo = builder.build_topology()
    solver = PEECCircuitSolver(topo)

    freq = 100e3
    Z = solver.compute_port_impedance(freq)
    omega = 2 * np.pi * freq
    L = Z.imag / omega
    R = Z.real

    print(f"  Freq = {freq/1e3:.0f} kHz")
    print(f"  Z = {R*1e3:.3f} + j{Z.imag*1e3:.3f} mOhm")
    print(f"  L = {L*1e9:.2f} nH")
    print(f"  R = {R*1e3:.3f} mOhm")

    # Analytical self-inductance of square loop (approximate)
    a = side
    r_w = wire_w / 2
    L_analytical = 2 * MU_0 * a / np.pi * (np.log(a / r_w) - 0.774)
    print(f"\n  Analytical L (square loop): {L_analytical*1e9:.2f} nH")
    err = abs(L - L_analytical) / L_analytical * 100
    print(f"  Error: {err:.1f}%")

    passed = err < 30  # 30% tolerance (analytical formula is itself approximate)
    status = "PASS" if passed else "FAIL"
    print(f"  Check: L matches analytical within 30% ... {status}")

    return passed


def main():
    print("PEEC Shield Conductor Validation")
    print("=" * 60)

    results = []
    results.append(test_basic_shield())
    results.append(test_shield_distance())
    results.append(test_no_shield_baseline())

    n_pass = sum(results)
    n_total = len(results)

    print("\n" + "=" * 60)
    print(f"SUMMARY: {n_pass}/{n_total} tests PASSED")
    print("=" * 60)

    return 0 if n_pass == n_total else 1


if __name__ == '__main__':
    sys.exit(main())
