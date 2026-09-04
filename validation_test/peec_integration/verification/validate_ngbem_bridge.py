"""
validate_ngbem_bridge.py

Validation of the NGSolve BEM -> PEEC circuit extraction bridge.

Tests:
1. Virtual topology and port assignment
2. Matrix dimensions and metadata
3. Galerkin symmetry of the surface-current L and P matrices

All tests use try/except for NGSolve BEM availability.
Tests that require NGSolve BEM are SKIPped if not installed.

Part of Radia project
"""

import sys
import os
import numpy as np

# Add the repository package root to the import path.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src'))

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
    """Check if NGSolve and its BEM module are available."""
    try:
        import ngsolve
        from ngsolve import bem
        return True
    except ImportError:
        return False


# =====================================================================
# Test 1: ngsbem_interface module imports correctly
# =====================================================================
def test_module_import():
    """Test that ngsbem_interface can be imported without assembling BEM."""
    print("=" * 60)
    print("Test 1: Module Import")
    print("=" * 60)

    try:
        from radia.ngsbem_interface import NGBEMBridge, extract_edge_geometry
        record("import NGBEMBridge", "PASS")
    except ImportError as e:
        record("import NGBEMBridge", "FAIL", str(e))

# =====================================================================
# Test 2: NGBEMBridge virtual topology (requires NGSolve BEM)
# =====================================================================
def test_ngbem_bridge_topology():
    """Test NGBEMBridge virtual topology from the NGSolve BEM solver."""
    print("\n" + "=" * 60)
    print("Test 2: NGBEMBridge Virtual Topology")
    print("=" * 60)

    if not check_ngbem_available():
        record("NGSolve BEM bridge topology", "SKIP", "NGSolve BEM unavailable")
        return

    from radia.ngsbem_peec import NGBEMPEECSolver, create_plate_mesh
    from radia.ngsbem_interface import NGBEMBridge

    # Create simple plate mesh
    mesh = create_plate_mesh(0.01, 0.001, 0.003)

    # Assemble the NGSolve BEM solver.
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
    assert topo.get('backend') == 'ngsbem', \
        f"Expected backend='ngsbem', got '{topo.get('backend')}'"
    record("backend metadata", "PASS")


# =====================================================================
# Test 3: Matrix symmetry (requires NGSolve BEM)
# =====================================================================
def test_matrix_symmetry():
    """Test NGSolve BEM L and P matrices are symmetric (Galerkin property)."""
    print("\n" + "=" * 60)
    print("Test 3: Matrix Symmetry (Galerkin)")
    print("=" * 60)

    if not check_ngbem_available():
        record("matrix symmetry", "SKIP", "NGSolve BEM unavailable")
        return

    from radia.ngsbem_peec import NGBEMPEECSolver, create_plate_mesh
    from radia.ngsbem_interface import NGBEMBridge

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
# Main
# =====================================================================
if __name__ == '__main__':
    print("NGSolve BEM Bridge Validation")
    print("=" * 60)
    print(f"NGSolve BEM available: {check_ngbem_available()}")
    print()

    test_module_import()
    test_ngbem_bridge_topology()
    test_matrix_symmetry()

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
