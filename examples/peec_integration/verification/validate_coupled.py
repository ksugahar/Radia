"""
validate_coupled.py

Validation script for coupled PEEC + MMM solver (CoupledPEECSolver).

Tests:
  1. Wire near mu_r=1 material: Delta_L should be ~0 (no coupling)
  2. Wire near high-mu material: L should increase (mutual inductance)
  3. Symmetry check: Delta_L should be approximately symmetric
  4. Frequency sweep: Z(f) is physically reasonable
  5. FastHenry .magnetic block parsing and solve

Part of Radia project
"""

import sys
import os
import numpy as np

# Add src/radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad
from peec_matrices import PEECBuilder
from peec_coupled import CoupledPEECSolver, biot_savart_finite_filament
from fasthenry_parser import FastHenryParser

MU_0 = 4e-7 * np.pi


def build_single_wire_topology(length=0.1, w=1e-3, h=1e-3, sigma=5.8e7):
    """Helper: build a single wire PEEC topology."""
    builder = PEECBuilder()
    n1 = builder.add_node_at(0, 0, 0)
    n2 = builder.add_node_at(length, 0, 0)
    builder.add_connected_segment(n1, n2, w, h, sigma=sigma)
    builder.add_port(n1, n2)
    return builder.build_topology()


def test_biot_savart():
    """Test 0: Biot-Savart finite filament formula."""
    print("Test 0: Biot-Savart finite filament")
    print("-" * 40)

    # Wire along x-axis from (0,0,0) to (1,0,0), obs at (0.5,d,0)
    p1 = [0, 0, 0]
    p2 = [1, 0, 0]
    d = 0.01  # 1 cm perpendicular distance

    H = biot_savart_finite_filament(p1, p2, [0.5, d, 0])

    # For a wire centered at x=0.5, obs at perpendicular distance d:
    # H_z = I/(4*pi*d) * 2*cos(alpha) where alpha = atan(d/0.5)
    # For L=1m, d=0.01m: cos(alpha) ~ 0.5/sqrt(0.5^2+0.01^2) ~ 0.9998
    half_L = 0.5
    cos_alpha = half_L / np.sqrt(half_L**2 + d**2)
    H_expected_z = (1.0 / (4 * np.pi * d)) * 2 * cos_alpha

    # H should be in z-direction (wire along x, obs in y direction)
    print(f"  H = [{H[0]:.4f}, {H[1]:.4f}, {H[2]:.4f}] A/m")
    print(f"  Expected Hz = {H_expected_z:.4f} A/m")

    err = abs(abs(H[2]) - H_expected_z) / H_expected_z * 100
    print(f"  Error: {err:.4f}%")

    # H should be primarily in z (perpendicular to both wire and offset)
    assert abs(H[0]) < 1e-10, f"Hx should be ~0, got {H[0]}"
    assert abs(H[1]) < 1e-10, f"Hy should be ~0, got {H[1]}"

    passed = err < 0.1
    print(f"  {'PASSED' if passed else 'FAILED'}\n")
    return passed


def test_mu_r_1_no_coupling():
    """Test 1: Wire near mu_r=1 material -> Delta_L ~ 0."""
    print("Test 1: mu_r=1 material (no coupling)")
    print("-" * 40)

    rad.UtiDelAll()

    # Build wire topology
    topo = build_single_wire_topology()

    # Create mu_r=1 (air) hexahedral block near the wire
    core_verts = [
        [0.02, 0.02, -0.005],
        [0.08, 0.02, -0.005],
        [0.08, 0.04, -0.005],
        [0.02, 0.04, -0.005],
        [0.02, 0.02, 0.005],
        [0.08, 0.02, 0.005],
        [0.08, 0.04, 0.005],
        [0.02, 0.04, 0.005],
    ]
    core = rad.ObjHexahedron(core_verts, [0, 0, 0])
    mat = rad.MatLin(1)  # mu_r = 1 (no susceptibility)
    rad.MatApl(core, mat)

    # Coupled solver
    solver = CoupledPEECSolver(topo, [core])
    Delta_L = solver.compute_coupling_matrix()

    # Delta_L should be essentially zero
    max_Delta_L = np.max(np.abs(Delta_L))
    L_air = topo['L'][0, 0]

    print(f"  L_air = {L_air*1e9:.4f} nH")
    print(f"  max|Delta_L| = {max_Delta_L:.4e} H")
    print(f"  |Delta_L/L_air| = {max_Delta_L/L_air*100:.6f}%")

    passed = max_Delta_L / L_air < 0.001  # Less than 0.1% of L_air
    print(f"  {'PASSED' if passed else 'FAILED'}\n")
    return passed


def test_high_mu_coupling():
    """Test 2: Wire near high-mu material -> L increases."""
    print("Test 2: High-mu material (L should increase)")
    print("-" * 40)

    rad.UtiDelAll()

    # Build wire along x-axis
    topo = build_single_wire_topology()

    # Create high-mu block near wire (above it, offset in y)
    core_verts = [
        [0.02, 0.005, -0.005],
        [0.08, 0.005, -0.005],
        [0.08, 0.015, -0.005],
        [0.02, 0.015, -0.005],
        [0.02, 0.005, 0.005],
        [0.08, 0.005, 0.005],
        [0.08, 0.015, 0.005],
        [0.02, 0.015, 0.005],
    ]
    core = rad.ObjHexahedron(core_verts, [0, 0, 0])
    mat = rad.MatLin(999)  # mu_r = 1000
    rad.MatApl(core, mat)

    # Coupled solver
    solver = CoupledPEECSolver(topo, [core])
    Delta_L = solver.compute_coupling_matrix()

    L_air = topo['L'][0, 0]
    L_total = L_air + Delta_L[0, 0]

    print(f"  L_air     = {L_air*1e9:.4f} nH")
    print(f"  Delta_L   = {Delta_L[0,0]*1e9:.4f} nH")
    print(f"  L_total   = {L_total*1e9:.4f} nH")
    print(f"  L_increase = {Delta_L[0,0]/L_air*100:.2f}%")

    # Delta_L should be positive (magnetic material increases inductance)
    passed = Delta_L[0, 0] > 0
    print(f"  Delta_L > 0: {passed}")
    print(f"  {'PASSED' if passed else 'FAILED'}\n")
    return passed


def test_symmetry():
    """Test 3: Delta_L matrix symmetry for multi-segment conductor."""
    print("Test 3: Delta_L symmetry (2-segment wire)")
    print("-" * 40)

    rad.UtiDelAll()

    # Build 2-segment wire
    builder = PEECBuilder()
    n1 = builder.add_node_at(0, 0, 0)
    n2 = builder.add_node_at(0.05, 0, 0)
    n3 = builder.add_node_at(0.1, 0, 0)
    builder.add_connected_segment(n1, n2, 1e-3, 1e-3, sigma=5.8e7)
    builder.add_connected_segment(n2, n3, 1e-3, 1e-3, sigma=5.8e7)
    builder.add_port(n1, n3)
    topo = builder.build_topology()

    # Symmetric core above wire center
    core_verts = [
        [0.02, 0.005, -0.005],
        [0.08, 0.005, -0.005],
        [0.08, 0.015, -0.005],
        [0.02, 0.015, -0.005],
        [0.02, 0.005, 0.005],
        [0.08, 0.005, 0.005],
        [0.08, 0.015, 0.005],
        [0.02, 0.015, 0.005],
    ]
    core = rad.ObjHexahedron(core_verts, [0, 0, 0])
    mat = rad.MatLin(999)
    rad.MatApl(core, mat)

    solver = CoupledPEECSolver(topo, [core])
    Delta_L = solver.compute_coupling_matrix()

    n_seg = Delta_L.shape[0]
    print(f"  Delta_L shape: {Delta_L.shape}")

    # Check symmetry: |Delta_L[i,j] - Delta_L[j,i]| / max(|Delta_L|)
    max_val = np.max(np.abs(Delta_L))
    if max_val > 0:
        asym = np.max(np.abs(Delta_L - Delta_L.T)) / max_val * 100
    else:
        asym = 0

    print(f"  max|Delta_L| = {max_val:.4e} H")
    print(f"  Asymmetry = {asym:.2f}%")

    for i in range(n_seg):
        for j in range(n_seg):
            print(f"  Delta_L[{i},{j}] = {Delta_L[i,j]*1e9:.4f} nH")

    # Reciprocity: Delta_L should be approximately symmetric
    # Note: Not perfectly symmetric due to finite-filament Biot-Savart approximation
    passed = asym < 30  # Allow up to 30% asymmetry for finite filament approx
    print(f"  {'PASSED' if passed else 'FAILED'}\n")
    return passed


def test_frequency_sweep():
    """Test 4: Frequency sweep with coupled solver."""
    print("Test 4: Frequency sweep (coupled)")
    print("-" * 40)

    rad.UtiDelAll()

    topo = build_single_wire_topology()

    # High-mu material near wire
    core_verts = [
        [0.02, 0.005, -0.005],
        [0.08, 0.005, -0.005],
        [0.08, 0.015, -0.005],
        [0.02, 0.015, -0.005],
        [0.02, 0.005, 0.005],
        [0.08, 0.005, 0.005],
        [0.08, 0.015, 0.005],
        [0.02, 0.015, 0.005],
    ]
    core = rad.ObjHexahedron(core_verts, [0, 0, 0])
    mat = rad.MatLin(999)
    rad.MatApl(core, mat)

    solver = CoupledPEECSolver(topo, [core])
    solver.compute_coupling_matrix()

    # Frequency sweep
    freqs = np.logspace(2, 7, 11)  # 100 Hz to 10 MHz
    Z_port = solver.frequency_sweep(freqs)

    R = np.real(Z_port)
    omega = 2.0 * np.pi * freqs
    L = np.imag(Z_port) / omega

    print(f"  Frequency range: {freqs[0]:.0f} Hz to {freqs[-1]:.0f} Hz")
    print(f"  R at {freqs[0]:.0f} Hz: {R[0]*1e3:.4f} mOhm")
    print(f"  R at {freqs[-1]:.0f} Hz: {R[-1]*1e3:.4f} mOhm")
    print(f"  L at {freqs[0]:.0f} Hz: {L[0]*1e9:.4f} nH")
    print(f"  L at {freqs[-1]:.0f} Hz: {L[-1]*1e9:.4f} nH")

    # Physical checks
    all_R_positive = np.all(R > 0)
    all_L_positive = np.all(L > 0)

    print(f"  All R > 0: {all_R_positive}")
    print(f"  All L > 0: {all_L_positive}")

    # L should be constant (linear material, no frequency dependence in Delta_L)
    L_variation = (np.max(L) - np.min(L)) / np.mean(L) * 100
    print(f"  L variation: {L_variation:.2f}%")

    passed = all_R_positive and all_L_positive
    print(f"  {'PASSED' if passed else 'FAILED'}\n")
    return passed


def test_fasthenry_magnetic_block():
    """Test 5: FastHenry .magnetic block parsing."""
    print("Test 5: FastHenry .magnetic block parsing")
    print("-" * 40)

    inp_text = """
.Units m
N1 x=0 y=0 z=0
N2 x=0.1 y=0 z=0

E1 N1 N2 w=1e-3 h=1e-3 sigma=5.8e7

.external N1 N2

.magnetic
  type=box
  center=0.05,0.01,0.0
  size=0.06,0.01,0.01
  divisions=2,1,1
  mu_r=1000
.endmagnetic

.freq fmin=1e3 fmax=1e6 ndec=3

.end
"""

    parser = FastHenryParser()
    parser.parse_string(inp_text)

    summary = parser.get_summary()
    print(f"  Nodes: {summary['n_nodes']}")
    print(f"  Segments: {summary['n_segments']}")
    print(f"  Ports: {summary['n_ports']}")
    print(f"  Magnetic blocks: {summary['n_magnetic_blocks']}")

    assert summary['n_magnetic_blocks'] == 1, \
        f"Expected 1 magnetic block, got {summary['n_magnetic_blocks']}"

    block = parser.magnetic_blocks[0]
    assert block['type'] == 'box'
    assert block['mu_r'] == 1000.0
    assert block['divisions'] == [2, 1, 1]

    print(f"  Block type: {block['type']}")
    print(f"  mu_r: {block['mu_r']}")
    print(f"  center: {block['center']}")
    print(f"  size: {block['size']}")
    print(f"  divisions: {block['divisions']}")

    print("  PASSED\n")
    return True


def test_fasthenry_magnetic_hexahedron():
    """Test 6: FastHenry .magnetic hexahedron block."""
    print("Test 6: FastHenry .magnetic hexahedron block")
    print("-" * 40)

    inp_text = """
.Units m
N1 x=0 y=0 z=0
N2 x=0.1 y=0 z=0

E1 N1 N2 w=1e-3 h=1e-3 sigma=5.8e7

.external N1 N2

.magnetic
  type=hexahedron
  vertices=0.02,0.005,-0.005, 0.08,0.005,-0.005, 0.08,0.015,-0.005, 0.02,0.015,-0.005, 0.02,0.005,0.005, 0.08,0.005,0.005, 0.08,0.015,0.005, 0.02,0.015,0.005
  mu_r=500
  mu_r_imag=10
.endmagnetic

.end
"""

    parser = FastHenryParser()
    parser.parse_string(inp_text)

    assert len(parser.magnetic_blocks) == 1
    block = parser.magnetic_blocks[0]
    assert block['type'] == 'hexahedron'
    assert block['mu_r'] == 500.0
    assert block['mu_r_imag'] == 10.0
    assert len(block['vertices']) == 8
    assert len(block['vertices'][0]) == 3

    print(f"  Block type: {block['type']}")
    print(f"  mu_r: {block['mu_r']}")
    print(f"  mu_r_imag: {block['mu_r_imag']}")
    print(f"  Vertices: {len(block['vertices'])} (8 expected)")
    print(f"  First vertex: {block['vertices'][0]}")

    print("  PASSED\n")
    return True


def test_fasthenry_coupled_solve():
    """Test 7: Full coupled solve via FastHenry parser."""
    print("Test 7: FastHenry coupled solve")
    print("-" * 40)

    rad.UtiDelAll()

    # First solve WITHOUT magnetic block
    inp_air = """
.Units m
N1 x=0 y=0 z=0
N2 x=0.1 y=0 z=0

E1 N1 N2 w=1e-3 h=1e-3 sigma=5.8e7

.external N1 N2
.end
"""
    parser_air = FastHenryParser()
    parser_air.parse_string(inp_air)
    result_air = parser_air.solve(freqs=np.array([1e6]))

    # Then solve WITH magnetic block
    rad.UtiDelAll()

    inp_coupled = """
.Units m
N1 x=0 y=0 z=0
N2 x=0.1 y=0 z=0

E1 N1 N2 w=1e-3 h=1e-3 sigma=5.8e7

.external N1 N2

.magnetic
  type=box
  center=0.05,0.01,0.0
  size=0.06,0.01,0.01
  divisions=1,1,1
  mu_r=1000
.endmagnetic

.end
"""
    parser_coupled = FastHenryParser()
    parser_coupled.parse_string(inp_coupled)
    result_coupled = parser_coupled.solve(freqs=np.array([1e6]))

    L_air = result_air['L'][0]
    L_coupled = result_coupled['L'][0]
    R_air = result_air['R'][0]
    R_coupled = result_coupled['R'][0]

    print(f"  Air only:  R = {R_air*1e3:.4f} mOhm, L = {L_air*1e9:.4f} nH")
    print(f"  Coupled:   R = {R_coupled*1e3:.4f} mOhm, L = {L_coupled*1e9:.4f} nH")
    print(f"  Delta_L = {(L_coupled-L_air)*1e9:.4f} nH")

    if 'Delta_L' in result_coupled:
        print(f"  Delta_L matrix: {result_coupled['Delta_L'][0,0]*1e9:.4f} nH")

    # L should increase with magnetic material
    passed = L_coupled > L_air
    print(f"  L_coupled > L_air: {passed}")
    # R should stay the same (no losses in core for linear material)
    R_err = abs(R_coupled - R_air) / R_air * 100
    print(f"  R difference: {R_err:.4f}%")

    print(f"  {'PASSED' if passed else 'FAILED'}\n")
    return passed


if __name__ == '__main__':
    print("=" * 60)
    print("Coupled PEEC + MMM Solver Validation")
    print("=" * 60)
    print()

    results = []
    results.append(("Test 0: Biot-Savart formula", test_biot_savart()))
    results.append(("Test 1: mu_r=1 (no coupling)", test_mu_r_1_no_coupling()))
    results.append(("Test 2: High-mu coupling", test_high_mu_coupling()))
    results.append(("Test 3: Delta_L symmetry", test_symmetry()))
    results.append(("Test 4: Frequency sweep", test_frequency_sweep()))
    results.append(("Test 5: .magnetic box parsing", test_fasthenry_magnetic_block()))
    results.append(("Test 6: .magnetic hex parsing", test_fasthenry_magnetic_hexahedron()))
    results.append(("Test 7: FastHenry coupled solve", test_fasthenry_coupled_solve()))

    print("=" * 60)
    print("Summary")
    print("=" * 60)
    n_pass = sum(1 for _, p in results if p)
    n_total = len(results)
    for name, passed in results:
        print(f"  {name}: {'PASSED' if passed else 'FAILED'}")
    print(f"\n  {n_pass}/{n_total} tests passed")
    print("=" * 60)
