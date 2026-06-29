"""
validate_fasthenry.py

Validation script for FastHenry .inp parser.

Tests:
  0. Parser: .Units, N, E, .external, .freq, .default, .equiv directives
  1. Single wire: Compare parsed model vs manual PEECBuilder
  2. Parallel wires: Two-conductor bus bar from .inp
  3. Multi-filament: nwinc/nhinc from .inp matches manual build
  4. Series chain: N1-N2-N3-N4 with .external N1 N4
  5. Frequency sweep: Parse .freq, compute Z(f), verify R and L

Part of Radia project
"""

import sys
import os
import numpy as np

# Add src/radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

from fasthenry_parser import FastHenryParser

MU_0 = 4e-7 * np.pi

def test_parser():
    """Test 0: Parser directives."""
    print("Test 0: Parser directives")
    print("-" * 40)

    inp_text = """
* Test FastHenry input file
.Units mm
.default sigma=5.8e7 w=1 h=0.5 nwinc=1 nhinc=1

N1 x=0 y=0 z=0
N2 x=10 y=0 z=0
N3 x=20 y=0 z=0

E1 N1 N2 w=2 h=1
E2 N2 N3

.external N1 N3

.freq fmin=1e3 fmax=1e6 ndec=5

.end
"""
    parser = FastHenryParser()
    parser.parse_string(inp_text)

    summary = parser.get_summary()

    # Check units
    assert summary['units'] == 'mm', f"Expected mm, got {summary['units']}"
    print(f"  Units: {summary['units']} -> scale={parser.unit_scale}")

    # Check nodes
    assert summary['n_nodes'] == 3, f"Expected 3 nodes, got {summary['n_nodes']}"
    # N1 at (0, 0, 0) in meters
    assert abs(parser.nodes['N1']['x'] - 0.0) < 1e-10
    # N2 at (0.01, 0, 0) in meters (10mm)
    assert abs(parser.nodes['N2']['x'] - 0.01) < 1e-10
    # N3 at (0.02, 0, 0) in meters (20mm)
    assert abs(parser.nodes['N3']['x'] - 0.02) < 1e-10
    print(f"  Nodes: {summary['n_nodes']} (positions in meters: OK)")

    # Check segments
    assert summary['n_segments'] == 2, f"Expected 2 segments, got {summary['n_segments']}"
    # E1: w=2mm, h=1mm
    assert abs(parser.segments['E1']['w'] - 0.002) < 1e-10
    assert abs(parser.segments['E1']['h'] - 0.001) < 1e-10
    # E2: uses defaults w=1mm, h=0.5mm
    assert abs(parser.segments['E2']['w'] - 0.001) < 1e-10
    assert abs(parser.segments['E2']['h'] - 0.0005) < 1e-10
    print(f"  Segments: {summary['n_segments']} (dimensions in meters: OK)")

    # Check ports
    assert summary['n_ports'] == 1
    assert parser.ports[0] == ('N1', 'N3')
    print(f"  Ports: {summary['n_ports']} ({parser.ports[0]})")

    # Check frequency
    freqs = parser.get_frequencies()
    assert freqs is not None
    assert abs(freqs[0] - 1e3) < 1e-6
    assert abs(freqs[-1] - 1e6) < 1e-3
    n_decades = 3  # 1e3 to 1e6
    assert len(freqs) == int(n_decades * 5) + 1, f"Expected {int(n_decades*5)+1} freq points, got {len(freqs)}"
    print(f"  Frequencies: {len(freqs)} points ({freqs[0]:.0f} Hz to {freqs[-1]:.0f} Hz)")

    print("  PASSED\n")
    return True


def test_equiv():
    """Test .equiv directive."""
    print("Test 0b: .equiv directive")
    print("-" * 40)

    inp_text = """
.Units m
N1 x=0 y=0 z=0
N2 x=0.1 y=0 z=0
N3 x=0 y=0 z=0
N4 x=0.2 y=0 z=0

E1 N1 N2 w=1e-3 h=1e-3 sigma=5.8e7
E2 N3 N4 w=1e-3 h=1e-3 sigma=5.8e7

* Merge N1 and N3 (same position)
.equiv N1 N3

.external N1 N4
.end
"""
    parser = FastHenryParser()
    parser.parse_string(inp_text)

    # After equiv, E2 should reference N1 instead of N3
    assert parser.segments['E2']['n1'] == 'N1', f"Expected N1, got {parser.segments['E2']['n1']}"
    print(f"  .equiv N1 N3: E2.n1 = {parser.segments['E2']['n1']} (OK)")
    print("  PASSED\n")
    return True


def test_single_wire():
    """Test 1: Single wire - compare parsed vs manual builder."""
    print("Test 1: Single wire (parsed vs manual)")
    print("-" * 40)

    from peec_matrices import PEECBuilder as PyPEECBuilder
    from peec_topology import PEECCircuitSolver

    # FastHenry input
    inp_text = """
.Units m
N1 x=0 y=0 z=0
N2 x=0.1 y=0 z=0

E1 N1 N2 w=1e-3 h=1e-3 sigma=5.8e7

.external N1 N2
.end
"""

    # Parse and build
    parser = FastHenryParser()
    parser.parse_string(inp_text)
    builder_parsed = parser.to_peec_builder()
    topo_parsed = builder_parsed.build_topology()
    solver_parsed = PEECCircuitSolver(topo_parsed)

    # Manual build
    builder_manual = PyPEECBuilder()
    n1 = builder_manual.add_node_at(0, 0, 0)
    n2 = builder_manual.add_node_at(0.1, 0, 0)
    builder_manual.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
    builder_manual.add_port(n1, n2)
    topo_manual = builder_manual.build_topology()
    solver_manual = PEECCircuitSolver(topo_manual)

    # Compare at 1 MHz
    freq = 1e6
    Z_parsed = solver_parsed.compute_port_impedance(freq)
    Z_manual = solver_manual.compute_port_impedance(freq)

    R_parsed = np.real(Z_parsed)
    R_manual = np.real(Z_manual)
    L_parsed = np.imag(Z_parsed) / (2 * np.pi * freq)
    L_manual = np.imag(Z_manual) / (2 * np.pi * freq)

    R_err = abs(R_parsed - R_manual) / abs(R_manual) * 100 if abs(R_manual) > 0 else 0
    L_err = abs(L_parsed - L_manual) / abs(L_manual) * 100 if abs(L_manual) > 0 else 0

    print(f"  R_parsed = {R_parsed*1e3:.4f} mOhm, R_manual = {R_manual*1e3:.4f} mOhm (err={R_err:.4f}%)")
    print(f"  L_parsed = {L_parsed*1e9:.2f} nH, L_manual = {L_manual*1e9:.2f} nH (err={L_err:.4f}%)")

    passed = R_err < 0.01 and L_err < 0.01
    print(f"  {'PASSED' if passed else 'FAILED'}\n")
    return passed


def test_parallel_wires():
    """Test 2: Parallel wires (two conductors in parallel)."""
    print("Test 2: Parallel wires")
    print("-" * 40)

    from peec_topology import PEECCircuitSolver

    inp_text = """
.Units m
N1 x=0 y=0 z=0
N2 x=0.1 y=0 z=0

* Two parallel wires (same nodes = parallel)
E1 N1 N2 w=1e-3 h=1e-3 sigma=5.8e7
E2 N1 N2 w=1e-3 h=1e-3 sigma=5.8e7

.external N1 N2
.end
"""

    parser = FastHenryParser()
    parser.parse_string(inp_text)
    builder = parser.to_peec_builder()
    topo = builder.build_topology()
    solver = PEECCircuitSolver(topo)

    # At DC, two identical parallel resistors -> R_total = R_single / 2
    Z_dc = solver.compute_port_impedance(1.0)  # Near-DC
    R_dc = np.real(Z_dc)

    # Single wire resistance
    l = 0.1
    w = 1e-3
    h = 1e-3
    sigma = 5.8e7
    R_single = l / (sigma * w * h)
    R_expected = R_single / 2

    R_err = abs(R_dc - R_expected) / R_expected * 100

    print(f"  R_parallel = {R_dc*1e3:.4f} mOhm, Expected = {R_expected*1e3:.4f} mOhm (err={R_err:.2f}%)")

    # At 1 MHz, check inductance
    freq = 1e6
    Z_ac = solver.compute_port_impedance(freq)
    L_parallel = np.imag(Z_ac) / (2 * np.pi * freq)

    print(f"  L_parallel = {L_parallel*1e9:.2f} nH at {freq/1e6:.0f} MHz")

    passed = R_err < 0.1
    print(f"  {'PASSED' if passed else 'FAILED'}\n")
    return passed


def test_multifilament():
    """Test 3: Multi-filament from .inp matches manual build."""
    print("Test 3: Multi-filament (nwinc/nhinc from .inp)")
    print("-" * 40)

    from peec_matrices import PEECBuilder as PyPEECBuilder
    from peec_topology import PEECCircuitSolver

    inp_text = """
.Units m
N1 x=0 y=0 z=0
N2 x=0.1 y=0 z=0

E1 N1 N2 w=3e-3 h=3e-3 sigma=5.8e7 nwinc=3 nhinc=3

.external N1 N2
.end
"""

    # Parse and build
    parser = FastHenryParser()
    parser.parse_string(inp_text)
    builder_parsed = parser.to_peec_builder()
    topo_parsed = builder_parsed.build_topology()
    solver_parsed = PEECCircuitSolver(topo_parsed)

    # Manual build with same nwinc/nhinc
    builder_manual = PyPEECBuilder()
    n1 = builder_manual.add_node_at(0, 0, 0)
    n2 = builder_manual.add_node_at(0.1, 0, 0)
    builder_manual.add_connected_segment(n1, n2, 3e-3, 3e-3, sigma=5.8e7, nwinc=3, nhinc=3)
    builder_manual.add_port(n1, n2)
    topo_manual = builder_manual.build_topology()
    solver_manual = PEECCircuitSolver(topo_manual)

    # Compare at 1 MHz
    freq = 1e6
    Z_parsed = solver_parsed.compute_port_impedance(freq)
    Z_manual = solver_manual.compute_port_impedance(freq)

    Z_err = abs(Z_parsed - Z_manual) / abs(Z_manual) * 100

    print(f"  Z_parsed = {Z_parsed:.6f} Ohm")
    print(f"  Z_manual = {Z_manual:.6f} Ohm")
    print(f"  Vector error: {Z_err:.4f}%")

    # Also check filament count
    n_fil_parsed = topo_parsed['n_loop']
    n_fil_manual = topo_manual['n_loop']
    assert n_fil_parsed == n_fil_manual == 9, \
        f"Expected 9 filaments, got parsed={n_fil_parsed}, manual={n_fil_manual}"
    print(f"  Filaments: parsed={n_fil_parsed}, manual={n_fil_manual} (both 9)")

    passed = Z_err < 0.01
    print(f"  {'PASSED' if passed else 'FAILED'}\n")
    return passed


def test_series_chain():
    """Test 4: Series chain N1-N2-N3-N4."""
    print("Test 4: Series chain (4 segments)")
    print("-" * 40)

    from peec_topology import PEECCircuitSolver

    inp_text = """
.Units m
.default w=1e-3 h=1e-3 sigma=5.8e7

N1 x=0   y=0 z=0
N2 x=0.025 y=0 z=0
N3 x=0.05  y=0 z=0
N4 x=0.075 y=0 z=0
N5 x=0.1   y=0 z=0

E1 N1 N2
E2 N2 N3
E3 N3 N4
E4 N4 N5

.external N1 N5
.end
"""

    parser = FastHenryParser()
    parser.parse_string(inp_text)
    builder = parser.to_peec_builder()
    topo = builder.build_topology()
    solver = PEECCircuitSolver(topo)

    # DC resistance: 4 series segments, each 25mm
    Z_dc = solver.compute_port_impedance(1.0)
    R_dc = np.real(Z_dc)

    l_total = 0.1
    sigma = 5.8e7
    w = 1e-3
    h = 1e-3
    R_expected = l_total / (sigma * w * h)

    R_err = abs(R_dc - R_expected) / R_expected * 100

    print(f"  R_series = {R_dc*1e3:.4f} mOhm, Expected = {R_expected*1e3:.4f} mOhm (err={R_err:.2f}%)")

    passed = R_err < 0.1
    print(f"  {'PASSED' if passed else 'FAILED'}\n")
    return passed


def test_frequency_sweep():
    """Test 5: Frequency sweep with .freq directive."""
    print("Test 5: Frequency sweep (.freq)")
    print("-" * 40)

    inp_text = """
.Units m
N1 x=0 y=0 z=0
N2 x=0.1 y=0 z=0

E1 N1 N2 w=1e-3 h=1e-3 sigma=5.8e7

.external N1 N2
.freq fmin=100 fmax=1e6 ndec=3
.end
"""

    parser = FastHenryParser()
    parser.parse_string(inp_text)

    # Use solve() convenience method
    result = parser.solve()

    freqs = result['freqs']
    Z_port = result['Z_port']
    R = result['R']
    L = result['L']

    print(f"  Frequency range: {freqs[0]:.0f} Hz to {freqs[-1]:.0f} Hz ({len(freqs)} points)")
    print(f"  DC: R = {R[0]*1e3:.4f} mOhm, L = {L[0]*1e9:.2f} nH")
    print(f"  1 MHz: R = {R[-1]*1e3:.4f} mOhm, L = {L[-1]*1e9:.2f} nH")

    # R should be positive
    assert np.all(R > 0), f"R has negative values!"

    # L should be positive
    valid_L = L[freqs > 0]
    assert np.all(valid_L > 0), f"L has negative values!"

    # R should be approximately constant (no skin effect with 1x1 filament)
    R_variation = (np.max(R) - np.min(R)) / np.mean(R) * 100
    print(f"  R variation across freq: {R_variation:.2f}%")

    passed = R_variation < 1.0  # R should be nearly constant for 1x1
    print(f"  {'PASSED' if passed else 'FAILED'}\n")
    return passed


def test_default_params():
    """Test 6: .default parameters inheritance."""
    print("Test 6: .default parameters")
    print("-" * 40)

    inp_text = """
.Units mm
.default w=2 h=1 sigma=1e7 nwinc=2 nhinc=2

N1 x=0 y=0 z=0
N2 x=100 y=0 z=0

* E1 inherits all defaults
E1 N1 N2

.external N1 N2
.end
"""

    parser = FastHenryParser()
    parser.parse_string(inp_text)

    seg = parser.segments['E1']
    # w=2mm = 0.002m, h=1mm = 0.001m
    assert abs(seg['w'] - 0.002) < 1e-10, f"w={seg['w']}, expected 0.002"
    assert abs(seg['h'] - 0.001) < 1e-10, f"h={seg['h']}, expected 0.001"
    assert abs(seg['sigma'] - 1e7) < 1, f"sigma={seg['sigma']}, expected 1e7"
    assert seg['nwinc'] == 2, f"nwinc={seg['nwinc']}, expected 2"
    assert seg['nhinc'] == 2, f"nhinc={seg['nhinc']}, expected 2"

    print(f"  w = {seg['w']*1e3:.1f} mm (expected 2.0 mm)")
    print(f"  h = {seg['h']*1e3:.1f} mm (expected 1.0 mm)")
    print(f"  sigma = {seg['sigma']:.0e} (expected 1e7)")
    print(f"  nwinc = {seg['nwinc']} (expected 2)")
    print(f"  nhinc = {seg['nhinc']} (expected 2)")
    print("  PASSED\n")
    return True


def test_continuation_lines():
    """Test 7: Line continuation with '+'."""
    print("Test 7: Continuation lines")
    print("-" * 40)

    inp_text = """
.Units m
N1 x=0 y=0 z=0
N2 x=0.1 y=0 +
z=0

E1 N1 N2 w=1e-3 +
h=1e-3 sigma=5.8e7

.external N1 N2
.end
"""

    parser = FastHenryParser()
    parser.parse_string(inp_text)

    # Check that continuation worked
    assert abs(parser.nodes['N2']['x'] - 0.1) < 1e-10
    assert abs(parser.nodes['N2']['z'] - 0.0) < 1e-10
    assert abs(parser.segments['E1']['w'] - 1e-3) < 1e-10
    assert abs(parser.segments['E1']['h'] - 1e-3) < 1e-10

    print(f"  N2: x={parser.nodes['N2']['x']}, z={parser.nodes['N2']['z']} (OK)")
    print(f"  E1: w={parser.segments['E1']['w']}, h={parser.segments['E1']['h']} (OK)")
    print("  PASSED\n")
    return True


if __name__ == '__main__':
    print("=" * 60)
    print("FastHenry Parser Validation")
    print("=" * 60)
    print()

    results = []
    results.append(("Test 0: Parser directives", test_parser()))
    results.append(("Test 0b: .equiv", test_equiv()))
    results.append(("Test 1: Single wire", test_single_wire()))
    results.append(("Test 2: Parallel wires", test_parallel_wires()))
    results.append(("Test 3: Multi-filament", test_multifilament()))
    results.append(("Test 4: Series chain", test_series_chain()))
    results.append(("Test 5: Frequency sweep", test_frequency_sweep()))
    results.append(("Test 6: .default params", test_default_params()))
    results.append(("Test 7: Continuation lines", test_continuation_lines()))

    print("=" * 60)
    print("Summary")
    print("=" * 60)
    n_pass = sum(1 for _, p in results if p)
    n_total = len(results)
    for name, passed in results:
        print(f"  {name}: {'PASSED' if passed else 'FAILED'}")
    print(f"\n  {n_pass}/{n_total} tests passed")
    print("=" * 60)
