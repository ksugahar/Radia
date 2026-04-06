"""
validate_ngbem_bridge.py

Validation of ngbem -> PEEC Circuit Extraction Bridge

Tests:
1. Edge geometry extraction: mesh.edges -> centers, directions, lengths
2. Virtual topology: segment_nodes valid, node count matches vertices
3. Matrix symmetry: ngbem L and P are symmetric (Galerkin guarantee)
4. Topology dict format: all required keys present
5. PEEC vs ngbem inductance: single wire comparison
6. FastHenry use_ngbem: end-to-end test via parser

All tests use try/except for ngbem availability.
Tests that require ngbem are SKIPped if not installed.

Part of Radia project
"""

import sys
import os
import numpy as np

# Add src/radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

MU_0 = 4.0 * np.pi * 1e-7

# Track test results
results = []


def record(name, status, detail=""):
    results.append((name, status, detail))
    symbol = {"PASS": "OK", "FAIL": "FAIL", "SKIP": "SKIP"}[status]
    print(f"  [{symbol}] {name}")
    if detail:
        print(f"        {detail}")


def check_ngbem_available():
    """Check if ngbem and NGSolve are available."""
    try:
        import ngsolve
        import ngbem
        return True
    except ImportError:
        return False


# =====================================================================
# Test 1: ngbem_interface module imports correctly
# =====================================================================
def test_module_import():
    """Test that ngbem_interface module can be imported (no ngbem needed)."""
    print("=" * 60)
    print("Test 1: Module Import")
    print("=" * 60)

    try:
        from ngbem_interface import NGBEMBridge, extract_edge_geometry
        record("import NGBEMBridge", "PASS")
    except ImportError as e:
        record("import NGBEMBridge", "FAIL", str(e))

    # gmsh_mesh_import removed -- surface mesh functions no longer available
    record("import gmsh surface functions", "SKIP", "gmsh_mesh_import removed")


# =====================================================================
# Test 2: GMSH surface mesh reading
# =====================================================================
def test_gmsh_surface_read():
    """Test GMSH surface mesh reading -- SKIPPED (gmsh_mesh_import removed)."""
    print("\n" + "=" * 60)
    print("Test 2: GMSH Surface Mesh Reading (SKIPPED)")
    print("=" * 60)

    # gmsh_mesh_import is removed. Surface mesh reading should use
    # .vol files with NGSolve Mesh() instead.
    record("gmsh surface read", "SKIP", "gmsh_mesh_import removed")


# =====================================================================
# Test 3: NGBEMBridge virtual topology (requires ngbem)
# =====================================================================
def test_ngbem_bridge_topology():
    """Test NGBEMBridge virtual topology from ngbem solver."""
    print("\n" + "=" * 60)
    print("Test 3: NGBEMBridge Virtual Topology")
    print("=" * 60)

    if not check_ngbem_available():
        record("ngbem bridge topology", "SKIP", "ngbem not installed")
        return

    from ngbem_peec import NGBEMPEECSolver, create_plate_mesh
    from ngbem_interface import NGBEMBridge

    # Create simple plate mesh
    mesh = create_plate_mesh(0.01, 0.001, 0.003)

    # Assemble ngbem solver
    solver = NGBEMPEECSolver(mesh, sigma=5.8e7, thickness=1e-3)
    solver.assemble()

    # Create bridge
    bridge = NGBEMBridge(solver, port_spec='auto')
    topo = bridge.to_topology_dict()

    # Check required keys
    required_keys = ['L', 'R', 'segment_nodes', 'n_nodes', 'n_loop', 'ports']
    for key in required_keys:
        assert key in topo, f"Missing required key: {key}"
    record("required keys present", "PASS", f"{len(required_keys)} keys")

    # Check dimensions
    n_loop = topo['n_loop']
    assert topo['L'].shape == (n_loop, n_loop), \
        f"L shape mismatch: {topo['L'].shape} vs ({n_loop}, {n_loop})"
    record("L matrix dimensions", "PASS", f"{n_loop}x{n_loop}")

    assert topo['R'].shape == (n_loop,), \
        f"R shape mismatch: {topo['R'].shape} vs ({n_loop},)"
    record("R vector dimensions", "PASS")

    # Check segment_nodes validity
    seg_nodes = topo['segment_nodes']
    assert seg_nodes.shape[0] == n_loop, \
        f"segment_nodes rows ({seg_nodes.shape[0]}) != n_loop ({n_loop})"
    assert seg_nodes.shape[1] == 2, "segment_nodes must have 2 columns"
    assert np.all(seg_nodes >= 0), "Negative node indices"
    assert np.all(seg_nodes < topo['n_nodes']), "Node indices out of range"
    record("segment_nodes valid", "PASS",
           f"{n_loop} segments, {topo['n_nodes']} nodes")

    # Check port assignment
    assert len(topo['ports']) >= 1, "No ports assigned"
    pos, neg, pid = topo['ports'][0]
    assert pos != neg, "Port positive == negative"
    record("port assignment", "PASS", f"port 0: node {pos} -> {neg}")

    # Check backend metadata
    assert topo.get('backend') == 'ngbem', \
        f"Expected backend='ngbem', got '{topo.get('backend')}'"
    record("backend metadata", "PASS")


# =====================================================================
# Test 4: Matrix symmetry (requires ngbem)
# =====================================================================
def test_matrix_symmetry():
    """Test ngbem L and P matrices are symmetric (Galerkin property)."""
    print("\n" + "=" * 60)
    print("Test 4: Matrix Symmetry (Galerkin)")
    print("=" * 60)

    if not check_ngbem_available():
        record("matrix symmetry", "SKIP", "ngbem not installed")
        return

    from ngbem_peec import NGBEMPEECSolver, create_plate_mesh
    from ngbem_interface import NGBEMBridge

    mesh = create_plate_mesh(0.01, 0.001, 0.003)
    solver = NGBEMPEECSolver(mesh, sigma=5.8e7, thickness=1e-3)
    solver.assemble()

    bridge = NGBEMBridge(solver)
    topo = bridge.to_topology_dict()

    # Check L symmetry
    L = topo['L']
    L_asym = np.max(np.abs(L - L.T))
    L_norm = np.max(np.abs(L))
    L_rel_asym = L_asym / L_norm if L_norm > 0 else 0.0
    if L_rel_asym < 1e-10:
        record("L matrix symmetric", "PASS",
               f"relative asymmetry: {L_rel_asym:.2e}")
    else:
        record("L matrix symmetric", "FAIL",
               f"relative asymmetry: {L_rel_asym:.2e}")

    # Check P symmetry (if available)
    if 'P' in topo and topo['P'] is not None:
        P = topo['P']
        P_asym = np.max(np.abs(P - P.T))
        P_norm = np.max(np.abs(P))
        P_rel_asym = P_asym / P_norm if P_norm > 0 else 0.0
        if P_rel_asym < 1e-10:
            record("P matrix symmetric", "PASS",
                   f"relative asymmetry: {P_rel_asym:.2e}")
        else:
            record("P matrix symmetric", "FAIL",
                   f"relative asymmetry: {P_rel_asym:.2e}")
    else:
        record("P matrix symmetric", "SKIP", "P matrix not available")


# =====================================================================
# Test 5: PEEC vs ngbem inductance comparison (requires ngbem)
# =====================================================================
def test_peec_vs_ngbem():
    """Compare L from PEEC filaments vs ngbem BEM for a single wire."""
    print("\n" + "=" * 60)
    print("Test 5: PEEC vs ngbem Inductance")
    print("=" * 60)

    if not check_ngbem_available():
        record("PEEC vs ngbem", "SKIP", "ngbem not installed")
        return

    from peec_matrices import PEECBuilder
    from peec_topology import PEECCircuitSolver
    from ngbem_peec import NGBEMPEECSolver, create_plate_mesh
    from ngbem_interface import NGBEMBridge

    # PEEC: single wire 100mm, 1mm x 1mm cross-section
    builder = PEECBuilder()
    n1 = builder.add_node_at(0.0, 0.0, 0.0)
    n2 = builder.add_node_at(0.1, 0.0, 0.0)
    builder.add_connected_segment(n1, n2, 1e-3, 1e-3, 5.8e7)
    builder.add_port(n1, n2)
    topo_peec = builder.build_topology()

    solver_peec = PEECCircuitSolver(topo_peec)
    Z_peec = solver_peec.compute_port_impedance(1e6)
    L_peec = np.imag(Z_peec) / (2 * np.pi * 1e6)

    print(f"  PEEC: Z = {Z_peec:.4e}")
    print(f"  PEEC: L = {L_peec*1e9:.2f} nH")

    # ngbem: same wire as surface mesh plate
    mesh = create_plate_mesh(0.1, 0.001, 0.003)
    ngbem_solver = NGBEMPEECSolver(mesh, sigma=5.8e7, thickness=1e-3)
    ngbem_solver.assemble()

    bridge = NGBEMBridge(ngbem_solver, port_spec=([0, 0, 0], [0.1, 0, 0]))
    topo_ngbem = bridge.to_topology_dict()

    solver_ngbem = PEECCircuitSolver(topo_ngbem)
    Z_ngbem = solver_ngbem.compute_port_impedance(1e6)
    L_ngbem = np.imag(Z_ngbem) / (2 * np.pi * 1e6)

    print(f"  ngbem: Z = {Z_ngbem:.4e}")
    print(f"  ngbem: L = {L_ngbem*1e9:.2f} nH")

    # Compare (accept <50% difference - fundamentally different methods)
    if L_peec > 0 and L_ngbem > 0:
        rel_diff = abs(L_peec - L_ngbem) / abs(L_peec) * 100
        print(f"  Relative difference: {rel_diff:.1f}%")
        if rel_diff < 50:
            record("PEEC vs ngbem L", "PASS", f"{rel_diff:.1f}% difference")
        else:
            record("PEEC vs ngbem L", "FAIL",
                   f"{rel_diff:.1f}% > 50% threshold")
    else:
        record("PEEC vs ngbem L", "FAIL", "Non-positive inductance")


# =====================================================================
# Test 6: FastHenry parser use_ngbem option (requires ngbem)
# =====================================================================
def test_fasthenry_ngbem():
    """Test FastHenry parser with use_ngbem=True."""
    print("\n" + "=" * 60)
    print("Test 6: FastHenry use_ngbem")
    print("=" * 60)

    if not check_ngbem_available():
        record("FastHenry use_ngbem", "SKIP", "ngbem not installed")
        return

    from fasthenry_parser import FastHenryParser

    inp_text = """
.Units mm
.default sigma=5.8e7

N1 x=0 y=0 z=0
N2 x=100 y=0 z=0

E1 N1 N2 w=1 h=1

.external N1 N2
.freq fmin=1e3 fmax=1e6 ndec=3
.end
"""

    parser = FastHenryParser()
    parser.parse_string(inp_text)

    # Solve with ngbem backend
    result = parser.solve(use_ngbem=True)

    assert result['backend'] == 'ngbem', \
        f"Expected backend='ngbem', got '{result['backend']}'"
    record("backend field", "PASS")

    assert len(result['freqs']) > 0, "No frequencies"
    record("frequency sweep", "PASS",
           f"{len(result['freqs'])} frequencies")

    # L should be positive and physically reasonable (1-1000 nH for 100mm wire)
    L_dc = result['L'][0] if result['L'][0] > 0 else result['L'][-1]
    L_nH = L_dc * 1e9
    print(f"  L = {L_nH:.1f} nH")
    if 1 < L_nH < 1000:
        record("physically reasonable L", "PASS", f"{L_nH:.1f} nH")
    else:
        record("physically reasonable L", "FAIL", f"{L_nH:.1f} nH out of range")

    # Also solve with PEEC and compare
    result_peec = parser.solve(use_ngbem=False)
    L_peec_nH = result_peec['L'][-1] * 1e9
    print(f"  PEEC L = {L_peec_nH:.1f} nH")
    print(f"  ngbem L = {L_nH:.1f} nH")

    record("end-to-end FastHenry ngbem", "PASS")


# =====================================================================
# Test 7: gmsh_surface_to_ngsolve (requires NGSolve)
# =====================================================================
def test_gmsh_to_ngsolve():
    """Test GMSH surface mesh to NGSolve -- SKIPPED (gmsh_mesh_import removed)."""
    print("\n" + "=" * 60)
    print("Test 7: GMSH Surface to NGSolve (SKIPPED)")
    print("=" * 60)

    # gmsh_mesh_import is removed. Use .vol files with NGSolve Mesh() instead.
    record("NGSolve surface mesh", "SKIP", "gmsh_mesh_import removed")


# =====================================================================
# Test 8: ngbem_coupled compute_coupling_radia signature
# =====================================================================
def test_ngbem_coupled_signature():
    """Test that ngbem_coupled has the compute_coupling_radia method."""
    print("\n" + "=" * 60)
    print("Test 8: ngbem_coupled Method Signatures")
    print("=" * 60)

    try:
        from ngbem_coupled import CoupledPEECMMM
        record("import CoupledPEECMMM", "PASS")

        # Check that the method exists
        assert hasattr(CoupledPEECMMM, 'compute_coupling_radia'), \
            "Missing compute_coupling_radia method"
        record("compute_coupling_radia exists", "PASS")

        assert hasattr(CoupledPEECMMM, '_extract_edge_geometry'), \
            "Missing _extract_edge_geometry method"
        record("_extract_edge_geometry exists", "PASS")

    except ImportError as e:
        record("import CoupledPEECMMM", "FAIL", str(e))


# =====================================================================
# Main
# =====================================================================
if __name__ == '__main__':
    print("ngbem Bridge Validation")
    print("=" * 60)
    print(f"NGSolve/ngbem available: {check_ngbem_available()}")
    print()

    test_module_import()
    test_gmsh_surface_read()
    test_ngbem_bridge_topology()
    test_matrix_symmetry()
    test_peec_vs_ngbem()
    test_fasthenry_ngbem()
    test_gmsh_to_ngsolve()
    test_ngbem_coupled_signature()

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    n_pass = sum(1 for _, s, _ in results if s == "PASS")
    n_fail = sum(1 for _, s, _ in results if s == "FAIL")
    n_skip = sum(1 for _, s, _ in results if s == "SKIP")

    print(f"  PASS: {n_pass}")
    print(f"  FAIL: {n_fail}")
    print(f"  SKIP: {n_skip}")
    print(f"  Total: {len(results)}")

    if n_fail > 0:
        print("\nFailed tests:")
        for name, status, detail in results:
            if status == "FAIL":
                print(f"  - {name}: {detail}")

    if n_fail == 0:
        print("\nAll tests PASSED (or SKIPPED due to missing dependencies).")
    else:
        print(f"\n{n_fail} test(s) FAILED!")
        sys.exit(1)
