"""
Test N-Port Block Lanczos SPICE Generation

Tests the NPortBlockLanczosSPICE class for multi-port model order reduction
and SPICE netlist generation.

Tests:
1. 2-port WPT (Wireless Power Transfer) system
2. 3-port system
3. Frequency response accuracy
4. SPICE netlist validity
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import numpy as np
from lanczos_reduction import NPortBlockLanczosSPICE


def test_2port_wpt():
    """Test 2-port WPT system."""
    print("=" * 60)
    print("Test 1: 2-Port WPT System")
    print("=" * 60)

    # Create 2-port WPT system
    n = 20  # Total DOFs
    p = 2   # Ports

    # Inductance matrix
    L = np.eye(n) * 10e-6  # 10 uH self inductance

    # TX region coupling (DOFs 0-9)
    for i in range(10):
        for j in range(10):
            if i != j:
                L[i, j] = 1e-6 / (1 + abs(i - j))

    # RX region coupling (DOFs 10-19)
    for i in range(10, 20):
        for j in range(10, 20):
            if i != j:
                L[i, j] = 1e-6 / (1 + abs(i - j))

    # TX-RX mutual inductance
    k_coupling = 0.3
    M = k_coupling * 10e-6
    for i in range(10):
        for j in range(10, 20):
            dist = abs(i - (j - 10))
            L[i, j] = M / (1 + 0.5 * dist)
            L[j, i] = L[i, j]

    # Resistance
    R = np.ones(n) * 0.1

    port_indices = [0, 10]

    # Reduce
    reducer = NPortBlockLanczosSPICE(n_ports=p, n_stages=5)
    result = reducer.reduce(L, R, port_indices)

    # Generate SPICE
    netlist = reducer.to_spice(result, "WPT_2PORT", ["tx_p", "tx_n", "rx_p", "rx_n"])

    # Verify frequency response
    frequencies = np.logspace(3, 6, 20)
    Z_full = reducer.compute_impedance_full(L, R, port_indices, frequencies)
    Z_red = reducer.compute_impedance(result, frequencies)

    err_Z11 = np.max(np.abs(Z_full[:, 0, 0] - Z_red[:, 0, 0]) / np.abs(Z_full[:, 0, 0])) * 100
    err_Z12 = np.max(np.abs(Z_full[:, 0, 1] - Z_red[:, 0, 1]) / (np.abs(Z_full[:, 0, 1]) + 1e-15)) * 100

    print(f"  Reduction: {n} DOF -> {result.n_stages * p} DOF ({100*(1 - result.n_stages*p/n):.0f}% compression)")
    print(f"  Z_11 max error: {err_Z11:.2f}%")
    print(f"  Z_12 max error: {err_Z12:.2f}%")
    print(f"  SPICE netlist: {len(netlist.split(chr(10)))} lines")

    passed = err_Z11 < 1.0 and err_Z12 < 1.0
    print(f"  Result: {'PASSED' if passed else 'FAILED'}")
    return passed


def test_3port_system():
    """Test 3-port system."""
    print("\n" + "=" * 60)
    print("Test 2: 3-Port System")
    print("=" * 60)

    # Create 3-port system
    n = 30
    p = 3

    # Inductance matrix with block structure
    L = np.eye(n) * 5e-6

    # Intra-block coupling
    for block in range(3):
        start = block * 10
        for i in range(start, start + 10):
            for j in range(start, start + 10):
                if i != j:
                    L[i, j] = 0.5e-6 / (1 + abs(i - j))

    # Inter-block coupling
    for b1 in range(3):
        for b2 in range(b1 + 1, 3):
            k = 0.2 / (1 + abs(b1 - b2))  # Coupling decreases with distance
            M = k * 5e-6
            for i in range(10):
                for j in range(10):
                    L[b1*10 + i, b2*10 + j] = M / (1 + 0.3 * (abs(i - j)))
                    L[b2*10 + j, b1*10 + i] = L[b1*10 + i, b2*10 + j]

    R = np.ones(n) * 0.05

    port_indices = [0, 10, 20]

    # Reduce
    reducer = NPortBlockLanczosSPICE(n_ports=p, n_stages=5)
    result = reducer.reduce(L, R, port_indices)

    # Generate SPICE
    netlist = reducer.to_spice(result, "THREE_PORT")

    # Verify frequency response
    frequencies = np.logspace(3, 6, 20)
    Z_full = reducer.compute_impedance_full(L, R, port_indices, frequencies)
    Z_red = reducer.compute_impedance(result, frequencies)

    max_err = 0.0
    for i in range(p):
        for j in range(p):
            err = np.max(np.abs(Z_full[:, i, j] - Z_red[:, i, j]) / (np.abs(Z_full[:, i, j]) + 1e-15)) * 100
            max_err = max(max_err, err)
            print(f"  Z_{i+1}{j+1} max error: {err:.2f}%")

    print(f"  Reduction: {n} DOF -> {result.n_stages * p} DOF ({100*(1 - result.n_stages*p/n):.0f}% compression)")
    print(f"  SPICE netlist: {len(netlist.split(chr(10)))} lines")

    passed = max_err < 5.0
    print(f"  Result: {'PASSED' if passed else 'FAILED'}")
    return passed


def test_spice_netlist_format():
    """Test SPICE netlist format validity."""
    print("\n" + "=" * 60)
    print("Test 3: SPICE Netlist Format")
    print("=" * 60)

    n = 10
    p = 2

    L = np.eye(n) * 1e-6
    L[0, 5] = L[5, 0] = 0.3e-6  # Mutual between ports
    R = np.ones(n) * 0.01

    port_indices = [0, 5]

    reducer = NPortBlockLanczosSPICE(n_ports=p, n_stages=3)
    result = reducer.reduce(L, R, port_indices)
    netlist = reducer.to_spice(result, "TEST_FORMAT", ["in1", "out1", "in2", "out2"])

    lines = netlist.strip().split('\n')

    # Check required elements
    has_subckt = any('.SUBCKT' in line.upper() for line in lines)
    has_ends = any('.ENDS' in line.upper() for line in lines)
    has_inductors = any(line.strip().startswith('L') for line in lines)
    has_k_statements = any(line.strip().startswith('K') for line in lines)
    has_resistors = any(line.strip().startswith('R') for line in lines)
    has_ports = 'in1' in netlist and 'out1' in netlist and 'in2' in netlist and 'out2' in netlist

    print(f"  .SUBCKT present: {has_subckt}")
    print(f"  .ENDS present: {has_ends}")
    print(f"  Inductors (L): {has_inductors}")
    print(f"  Coupling (K): {has_k_statements}")
    print(f"  Resistors (R): {has_resistors}")
    print(f"  Port names: {has_ports}")

    passed = all([has_subckt, has_ends, has_inductors, has_resistors, has_ports])
    print(f"  Result: {'PASSED' if passed else 'FAILED'}")

    # Print sample netlist
    print("\n  Sample netlist (first 20 lines):")
    for line in lines[:20]:
        print(f"    {line}")

    return passed


def test_convergence_with_stages():
    """Test convergence with increasing Lanczos stages."""
    print("\n" + "=" * 60)
    print("Test 4: Convergence with Stages")
    print("=" * 60)

    n = 40
    p = 2

    # Dense coupling matrix
    L = np.eye(n) * 10e-6
    for i in range(n):
        for j in range(i + 1, n):
            coupling = 3e-6 * np.exp(-0.1 * abs(i - j))
            L[i, j] = L[j, i] = coupling

    R = np.ones(n) * 0.1
    port_indices = [0, 20]

    frequencies = np.array([100e3])  # Single frequency for simplicity

    Z_full = None
    errors = []

    for stages in [2, 4, 6, 8, 10]:
        reducer = NPortBlockLanczosSPICE(n_ports=p, n_stages=stages)
        result = reducer.reduce(L, R, port_indices)

        if Z_full is None:
            Z_full = reducer.compute_impedance_full(L, R, port_indices, frequencies)

        Z_red = reducer.compute_impedance(result, frequencies)

        err = np.abs(Z_full[0] - Z_red[0]) / np.abs(Z_full[0]) * 100
        max_err = np.max(err)
        errors.append(max_err)

        print(f"  Stages={stages:2d}: reduced DOF={stages*p:3d}, max error={max_err:.4f}%")

    # Check monotonic decrease (with some tolerance)
    monotonic = all(errors[i] >= errors[i+1] * 0.5 for i in range(len(errors) - 1))
    converged = errors[-1] < 1.0

    passed = converged
    print(f"  Convergence achieved: {converged}")
    print(f"  Result: {'PASSED' if passed else 'FAILED'}")
    return passed


def main():
    """Run all tests."""
    print("N-Port Block Lanczos SPICE Generation Test Suite")
    print("=" * 60)

    results = {
        'test_2port_wpt': test_2port_wpt(),
        'test_3port_system': test_3port_system(),
        'test_spice_netlist_format': test_spice_netlist_format(),
        'test_convergence_with_stages': test_convergence_with_stages(),
    }

    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)

    all_passed = True
    for name, passed in results.items():
        status = "PASSED" if passed else "FAILED"
        print(f"  {name}: {status}")
        if not passed:
            all_passed = False

    print("-" * 60)
    overall = "ALL TESTS PASSED" if all_passed else "SOME TESTS FAILED"
    print(f"  Overall: {overall}")

    return 0 if all_passed else 1


if __name__ == '__main__':
    sys.exit(main())
