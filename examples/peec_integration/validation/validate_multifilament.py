"""
validate_multifilament.py

Validation of PEEC Multi-Filament (nwinc/nhinc) implementation.

Tests:
  0. API: nwinc/nhinc parameters accepted, correct number of filaments created
  1. DC resistance: nwinc*nhinc sub-filaments in parallel -> R = R_single / (nwinc*nhinc)
  2. Inductance: Multi-filament L < single-filament L (mutual coupling reduces effective L)
  3. Skin effect: At high frequency, R(f) increases due to non-uniform current distribution
  4. Convergence: Increasing nwinc/nhinc should converge

Part of Radia project
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))
from peec_matrices import PEECBuilder
from peec_topology import PEECCircuitSolver

MU_0 = 4e-7 * np.pi


def test_0_api():
    """Test: nwinc/nhinc parameters work, correct number of filaments."""
    print("=" * 60)
    print("Test 0: Multi-filament API")
    print("=" * 60)

    # Single filament (nwinc=1, nhinc=1)
    builder1 = PEECBuilder()
    n1 = builder1.add_node_at(0, 0, 0)
    n2 = builder1.add_node_at(0.1, 0, 0)
    builder1.add_connected_segment(n1, n2, 1e-3, 1e-3)
    builder1.add_port(n1, n2)
    topo1 = builder1.build_topology()

    n_fil_1x1 = topo1['n_loop']
    print(f"  1x1 filaments: {n_fil_1x1}")
    assert n_fil_1x1 == 1, f"Expected 1 filament, got {n_fil_1x1}"

    # 3x3 = 9 sub-filaments
    builder9 = PEECBuilder()
    n1 = builder9.add_node_at(0, 0, 0)
    n2 = builder9.add_node_at(0.1, 0, 0)
    builder9.add_connected_segment(n1, n2, 3e-3, 3e-3, nwinc=3, nhinc=3)
    builder9.add_port(n1, n2)
    topo9 = builder9.build_topology()

    n_fil_3x3 = topo9['n_loop']
    print(f"  3x3 filaments: {n_fil_3x3}")
    assert n_fil_3x3 == 9, f"Expected 9 filaments, got {n_fil_3x3}"

    # 2x4 = 8 sub-filaments
    builder8 = PEECBuilder()
    n1 = builder8.add_node_at(0, 0, 0)
    n2 = builder8.add_node_at(0.1, 0, 0)
    builder8.add_connected_segment(n1, n2, 2e-3, 4e-3, nwinc=2, nhinc=4)
    builder8.add_port(n1, n2)
    topo8 = builder8.build_topology()

    n_fil_2x4 = topo8['n_loop']
    print(f"  2x4 filaments: {n_fil_2x4}")
    assert n_fil_2x4 == 8, f"Expected 8 filaments, got {n_fil_2x4}"

    # Segment connectivity: all sub-filaments share same node_from/node_to
    seg_nodes = topo9['segment_nodes']
    for i in range(n_fil_3x3):
        assert seg_nodes[i, 0] == 0, f"Sub-filament {i} node_from should be 0"
        assert seg_nodes[i, 1] == 1, f"Sub-filament {i} node_to should be 1"

    print("  All sub-filaments share same node_from/node_to: OK")
    print("  PASSED\n")
    return True


def test_1_dc_resistance():
    """Test: DC resistance of multi-filament = R_single / (nwinc*nhinc)."""
    print("=" * 60)
    print("Test 1: DC Resistance (parallel sub-filaments)")
    print("=" * 60)

    length = 0.1  # m
    width = 3e-3  # m
    height = 3e-3  # m
    sigma = 5.8e7  # S/m (copper)
    area = width * height  # 9e-6 m^2

    # Analytical DC resistance
    R_analytical = length / (sigma * area)
    print(f"  Analytical R_dc = {R_analytical*1e6:.4f} uOhm")

    # Single filament
    builder1 = PEECBuilder()
    n1 = builder1.add_node_at(0, 0, 0)
    n2 = builder1.add_node_at(length, 0, 0)
    builder1.add_connected_segment(n1, n2, width, height)
    builder1.add_port(n1, n2)
    topo1 = builder1.build_topology()
    solver1 = PEECCircuitSolver(topo1)
    Z_dc_1 = solver1.compute_port_impedance(0.001)  # near-DC
    R_dc_1 = Z_dc_1.real
    print(f"  1x1 R_dc = {R_dc_1*1e6:.4f} uOhm")

    # 3x3 = 9 sub-filaments (same total cross-section)
    builder9 = PEECBuilder()
    n1 = builder9.add_node_at(0, 0, 0)
    n2 = builder9.add_node_at(length, 0, 0)
    builder9.add_connected_segment(n1, n2, width, height, nwinc=3, nhinc=3)
    builder9.add_port(n1, n2)
    topo9 = builder9.build_topology()
    solver9 = PEECCircuitSolver(topo9)
    Z_dc_9 = solver9.compute_port_impedance(0.001)
    R_dc_9 = Z_dc_9.real
    print(f"  3x3 R_dc = {R_dc_9*1e6:.4f} uOhm")

    # Error vs analytical
    err_1 = abs(R_dc_1 - R_analytical) / R_analytical * 100
    err_9 = abs(R_dc_9 - R_analytical) / R_analytical * 100
    print(f"  1x1 error vs analytical: {err_1:.2f}%")
    print(f"  3x3 error vs analytical: {err_9:.2f}%")

    # Both should match analytical (parallel sub-filaments have same total area)
    if err_1 < 1.0 and err_9 < 1.0:
        print("  PASSED\n")
        return True
    else:
        print("  FAILED\n")
        return False


def test_2_inductance_reduction():
    """Test: Multi-filament inductance < single-filament inductance."""
    print("=" * 60)
    print("Test 2: Inductance Reduction (multi-filament)")
    print("=" * 60)

    length = 0.1  # m
    width = 3e-3  # m
    height = 3e-3  # m

    # Single filament
    builder1 = PEECBuilder()
    n1 = builder1.add_node_at(0, 0, 0)
    n2 = builder1.add_node_at(length, 0, 0)
    builder1.add_connected_segment(n1, n2, width, height)
    builder1.add_port(n1, n2)
    topo1 = builder1.build_topology()
    solver1 = PEECCircuitSolver(topo1)

    freq = 1000  # Low frequency (no skin effect)
    Z_1 = solver1.compute_port_impedance(freq)
    L_1 = Z_1.imag / (2 * np.pi * freq)
    print(f"  1x1 L = {L_1*1e9:.4f} nH")

    # 3x3 = 9 sub-filaments
    builder9 = PEECBuilder()
    n1 = builder9.add_node_at(0, 0, 0)
    n2 = builder9.add_node_at(length, 0, 0)
    builder9.add_connected_segment(n1, n2, width, height, nwinc=3, nhinc=3)
    builder9.add_port(n1, n2)
    topo9 = builder9.build_topology()
    solver9 = PEECCircuitSolver(topo9)

    Z_9 = solver9.compute_port_impedance(freq)
    L_9 = Z_9.imag / (2 * np.pi * freq)
    print(f"  3x3 L = {L_9*1e9:.4f} nH")

    # Multi-filament should have LOWER inductance due to mutual coupling
    # (current distributes across sub-filaments, partial flux cancellation)
    reduction = (L_1 - L_9) / L_1 * 100
    print(f"  L reduction: {reduction:.2f}%")

    if L_9 < L_1:
        print("  L_3x3 < L_1x1: OK (mutual coupling reduces effective L)")
        print("  PASSED\n")
        return True
    else:
        print("  L_3x3 >= L_1x1: UNEXPECTED")
        print("  FAILED\n")
        return False


def test_3_ac_resistance():
    """Test: At high frequency, R increases due to skin effect (non-uniform current)."""
    print("=" * 60)
    print("Test 3: AC Resistance Increase (skin effect)")
    print("=" * 60)

    length = 0.1  # m
    width = 3e-3  # m
    height = 3e-3  # m
    sigma = 5.8e7  # S/m

    # Skin depth
    freq_high = 1e6  # 1 MHz
    delta = np.sqrt(2 / (2 * np.pi * freq_high * MU_0 * sigma))
    print(f"  Skin depth at {freq_high/1e6:.0f} MHz: {delta*1e6:.1f} um")
    print(f"  Conductor size: {width*1e3:.1f} mm x {height*1e3:.1f} mm")
    print(f"  Ratio d/delta: {min(width, height)/delta:.1f}")

    # 5x5 = 25 sub-filaments (reasonable resolution for skin effect)
    builder = PEECBuilder()
    n1 = builder.add_node_at(0, 0, 0)
    n2 = builder.add_node_at(length, 0, 0)
    builder.add_connected_segment(n1, n2, width, height, nwinc=5, nhinc=5)
    builder.add_port(n1, n2)
    topo = builder.build_topology()
    solver = PEECCircuitSolver(topo)

    # DC resistance
    Z_dc = solver.compute_port_impedance(0.001)
    R_dc = Z_dc.real

    # AC resistance at high frequency
    Z_ac = solver.compute_port_impedance(freq_high)
    R_ac = Z_ac.real

    ratio = R_ac / R_dc
    print(f"  R_dc = {R_dc*1e6:.4f} uOhm")
    print(f"  R_ac({freq_high/1e6:.0f} MHz) = {R_ac*1e6:.4f} uOhm")
    print(f"  R_ac/R_dc = {ratio:.2f}")

    # At high frequency, R should increase (skin effect pushes current to surface)
    # With only 5x5 filaments, the effect may be modest
    if ratio > 1.0:
        print("  R_ac > R_dc: OK (skin effect present)")
        print("  PASSED\n")
        return True
    else:
        print("  R_ac <= R_dc: UNEXPECTED (no skin effect detected)")
        print("  NOTE: With only 5x5 sub-filaments, skin effect may be weak")
        print("  FAILED\n")
        return False


def test_4_convergence():
    """Test: Increasing nwinc/nhinc should converge."""
    print("=" * 60)
    print("Test 4: Convergence with increasing subdivision")
    print("=" * 60)

    length = 0.1  # m
    width = 3e-3  # m
    height = 3e-3  # m
    freq = 10000  # 10 kHz

    subdivisions = [1, 2, 3, 4, 5]
    Z_values = []

    for n in subdivisions:
        builder = PEECBuilder()
        n1 = builder.add_node_at(0, 0, 0)
        n2 = builder.add_node_at(length, 0, 0)
        builder.add_connected_segment(n1, n2, width, height, nwinc=n, nhinc=n)
        builder.add_port(n1, n2)
        topo = builder.build_topology()
        solver = PEECCircuitSolver(topo)
        Z = solver.compute_port_impedance(freq)
        Z_values.append(Z)

        n_fil = topo['n_loop']
        print(f"  {n}x{n} ({n_fil:3d} fils): Z = {Z.real*1e6:.4f} + j{Z.imag*1e6:.4f} uOhm")

    # Check convergence: difference between last two should be smaller than earlier
    diffs = []
    for i in range(1, len(Z_values)):
        diff = abs(Z_values[i] - Z_values[i-1])
        diffs.append(diff)

    if len(diffs) >= 2:
        converging = diffs[-1] < diffs[0]
        print(f"  First diff: {diffs[0]*1e6:.6f} uOhm")
        print(f"  Last diff:  {diffs[-1]*1e6:.6f} uOhm")
        if converging:
            print("  Convergence: OK (differences decreasing)")
            print("  PASSED\n")
            return True
        else:
            print("  Convergence: differences NOT decreasing")
            print("  FAILED\n")
            return False
    else:
        print("  Not enough data for convergence check")
        print("  FAILED\n")
        return False


def test_5_series_multifilament():
    """Test: Series connection with multi-filament segments."""
    print("=" * 60)
    print("Test 5: Series + Multi-filament")
    print("=" * 60)

    length = 0.05  # m per segment
    width = 2e-3   # m
    height = 2e-3  # m

    # Two series segments, each with 2x2 sub-filaments
    builder = PEECBuilder()
    n1 = builder.add_node_at(0, 0, 0)
    n2 = builder.add_node_at(length, 0, 0)
    n3 = builder.add_node_at(2*length, 0, 0)
    builder.add_connected_segment(n1, n2, width, height, nwinc=2, nhinc=2)
    builder.add_connected_segment(n2, n3, width, height, nwinc=2, nhinc=2)
    builder.add_port(n1, n3)
    topo = builder.build_topology()

    n_fil = topo['n_loop']
    print(f"  2 series segments x 2x2 = {n_fil} filaments")
    assert n_fil == 8, f"Expected 8 filaments (2 segments x 4 sub-fils), got {n_fil}"

    # Should have 1 junction node (n2 connects 4+4 filaments)
    n_junction = topo['n_junction']
    print(f"  Junction nodes: {n_junction}")

    # Solve
    solver = PEECCircuitSolver(topo)
    Z = solver.compute_port_impedance(1000)  # 1 kHz

    # Compare with single-filament series
    builder_1 = PEECBuilder()
    n1 = builder_1.add_node_at(0, 0, 0)
    n2 = builder_1.add_node_at(length, 0, 0)
    n3 = builder_1.add_node_at(2*length, 0, 0)
    builder_1.add_connected_segment(n1, n2, width, height)
    builder_1.add_connected_segment(n2, n3, width, height)
    builder_1.add_port(n1, n3)
    topo_1 = builder_1.build_topology()
    solver_1 = PEECCircuitSolver(topo_1)
    Z_1 = solver_1.compute_port_impedance(1000)

    print(f"  1x1 Z = {Z_1.real*1e6:.4f} + j{Z_1.imag*1e6:.4f} uOhm")
    print(f"  2x2 Z = {Z.real*1e6:.4f} + j{Z.imag*1e6:.4f} uOhm")

    # DC resistance should match (same total area)
    R_ratio = Z.real / Z_1.real
    print(f"  R ratio (2x2/1x1): {R_ratio:.4f}")

    # Inductance: multi-filament should be lower
    L_1 = Z_1.imag / (2 * np.pi * 1000)
    L_mf = Z.imag / (2 * np.pi * 1000)
    print(f"  L_1x1 = {L_1*1e9:.4f} nH, L_2x2 = {L_mf*1e9:.4f} nH")

    if abs(R_ratio - 1.0) < 0.05 and L_mf < L_1:
        print("  PASSED\n")
        return True
    else:
        print("  FAILED\n")
        return False


if __name__ == '__main__':
    results = []
    results.append(("Test 0: API", test_0_api()))
    results.append(("Test 1: DC Resistance", test_1_dc_resistance()))
    results.append(("Test 2: Inductance Reduction", test_2_inductance_reduction()))
    results.append(("Test 3: AC Resistance", test_3_ac_resistance()))
    results.append(("Test 4: Convergence", test_4_convergence()))
    results.append(("Test 5: Series + Multi-filament", test_5_series_multifilament()))

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    all_passed = True
    for name, passed in results:
        status = "PASSED" if passed else "FAILED"
        print(f"  {name}: {status}")
        if not passed:
            all_passed = False

    print()
    if all_passed:
        print("All tests PASSED!")
    else:
        print("Some tests FAILED!")
    print()
