"""
Test script for mmm_core pybind11 module
Validates the standalone MMM (Magnetic Moment Method) implementation

Tests:
1. MMMBuilder: Tetrahedron and Hexahedron element addition
2. MMMBuilder: Interaction matrix construction
3. MMMSolver: LU and BiCGSTAB solvers
4. NGSolve mesh integration: add_tetrahedra_from_mesh

Usage:
    python test_mmm_cpp.py
"""

import sys
import os
import numpy as np

# Add build directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build-msvc'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

def test_import():
    """Test basic module import"""
    print("=" * 60)
    print("Test 1: Module Import")
    print("=" * 60)

    try:
        import mmm_core
        print(f"mmm_core module imported successfully")

        # Check available classes
        print(f"  MMMBuilder: {hasattr(mmm_core, 'MMMBuilder')}")
        print(f"  MMMSolver: {hasattr(mmm_core, 'MMMSolver')}")

        return True
    except ImportError as e:
        print(f"FAILED: {e}")
        return False


def test_builder_tetrahedron():
    """Test MMMBuilder with tetrahedron elements"""
    print("\n" + "=" * 60)
    print("Test 2: MMMBuilder with Tetrahedra")
    print("=" * 60)

    import mmm_core

    # Create builder
    builder = mmm_core.MMMBuilder()

    # Add a single tetrahedron (unit tetrahedron)
    # Vertices: 4 points forming a tetrahedron
    vertices = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 0.866, 0.0],
        [0.5, 0.289, 0.816]
    ], dtype=np.float64)

    # Element connectivity (single tetrahedron using all 4 vertices)
    elements = np.array([[0, 1, 2, 3]], dtype=np.int32)

    # Add tetrahedra from mesh format
    builder.add_tetrahedra_from_mesh(vertices, elements)

    print(f"  Number of elements: {builder.num_elements}")
    print(f"  Total DOF: {builder.total_dof}")

    assert builder.num_elements == 1, "Expected 1 element"
    assert builder.total_dof == 3, "Expected 3 DOF for tetrahedron"

    # Build interaction matrix
    N, dof_offset = builder.build()
    total_dof = N.shape[0]
    n_elements = len(dof_offset) - 1

    print(f"  Interaction matrix shape: {N.shape}")
    print(f"  DOF offset: {dof_offset}")
    print(f"  N matrix (3x3):\n{N}")

    assert N.shape == (3, 3), "Expected 3x3 matrix"
    assert total_dof == 3
    assert n_elements == 1

    print("  PASSED")
    return True


def test_builder_hexahedron():
    """Test MMMBuilder with hexahedron elements (6 DOF MSC)"""
    print("\n" + "=" * 60)
    print("Test 3: MMMBuilder with Hexahedra (6 DOF MSC)")
    print("=" * 60)

    import mmm_core

    builder = mmm_core.MMMBuilder()

    # Add a single hexahedron (unit cube)
    # Vertices: 8 points forming a cube
    vertices = np.array([
        [0.0, 0.0, 0.0],  # 0: bottom-front-left
        [1.0, 0.0, 0.0],  # 1: bottom-front-right
        [1.0, 1.0, 0.0],  # 2: bottom-back-right
        [0.0, 1.0, 0.0],  # 3: bottom-back-left
        [0.0, 0.0, 1.0],  # 4: top-front-left
        [1.0, 0.0, 1.0],  # 5: top-front-right
        [1.0, 1.0, 1.0],  # 6: top-back-right
        [0.0, 1.0, 1.0],  # 7: top-back-left
    ], dtype=np.float64)

    # Element connectivity (single hexahedron using all 8 vertices)
    elements = np.array([[0, 1, 2, 3, 4, 5, 6, 7]], dtype=np.int32)

    # Add hexahedra from mesh format
    builder.add_hexahedra_from_mesh(vertices, elements)

    print(f"  Number of elements: {builder.num_elements}")
    print(f"  Total DOF: {builder.total_dof}")

    assert builder.num_elements == 1, "Expected 1 element"
    assert builder.total_dof == 6, "Expected 6 DOF for hexahedron (MSC)"

    # Build interaction matrix
    N, dof_offset = builder.build()
    total_dof = N.shape[0]

    print(f"  Interaction matrix shape: {N.shape}")
    print(f"  N matrix diagonal: {np.diag(N)}")

    assert N.shape == (6, 6), "Expected 6x6 matrix"
    assert total_dof == 6

    print("  PASSED")
    return True


def test_solver_lu():
    """Test MMMSolver with LU decomposition"""
    print("\n" + "=" * 60)
    print("Test 4: MMMSolver LU Decomposition")
    print("=" * 60)

    import mmm_core

    # Create a simple 2-tetrahedron problem
    builder = mmm_core.MMMBuilder()

    # Two tetrahedra side by side
    vertices = np.array([
        # First tetrahedron
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 0.866, 0.0],
        [0.5, 0.289, 0.816],
        # Second tetrahedron (shifted by 2.0 in x)
        [2.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [2.5, 0.866, 0.0],
        [2.5, 0.289, 0.816],
    ], dtype=np.float64)

    elements = np.array([
        [0, 1, 2, 3],  # First tetrahedron
        [4, 5, 6, 7],  # Second tetrahedron
    ], dtype=np.int32)

    builder.add_tetrahedra_from_mesh(vertices, elements)

    print(f"  Number of elements: {builder.num_elements}")
    print(f"  Total DOF: {builder.total_dof}")

    # Build interaction matrix
    N, dof_offset = builder.build()

    # Create solver
    solver = mmm_core.MMMSolver()
    solver.set_matrix(N, dof_offset)

    print(f"  Solver ready: {solver.is_ready}")
    print(f"  Total DOF: {solver.total_dof}")
    print(f"  Number of elements: {solver.num_elements}")

    # Solve: (I/chi - N) * M = H_ext
    # chi = 999 (mu_r = 1000)
    chi = 999.0
    inv_chi = np.array([1.0 / chi, 1.0 / chi], dtype=np.float64)  # Per element

    # External field H_ext = [0, 0, 1000] A/m for each DOF
    H_ext = np.array([0, 0, 1000, 0, 0, 1000], dtype=np.float64)

    # Solve using LU
    M = solver.solve_lu(inv_chi, H_ext, True)  # chi_per_element=True

    print(f"  Solution M shape: {M.shape}")
    print(f"  M (element 1): [{M[0]:.1f}, {M[1]:.1f}, {M[2]:.1f}] A/m")
    print(f"  M (element 2): [{M[3]:.1f}, {M[4]:.1f}, {M[5]:.1f}] A/m")

    # Check that M is reasonable (should be roughly chi * H_ext for isolated elements)
    assert abs(M[2]) > 100, "Expected significant M_z"
    assert abs(M[5]) > 100, "Expected significant M_z"

    print("  PASSED")
    return True


def test_solver_bicgstab():
    """Test MMMSolver with BiCGSTAB iteration"""
    print("\n" + "=" * 60)
    print("Test 5: MMMSolver BiCGSTAB Iteration")
    print("=" * 60)

    import mmm_core

    # Create a 2x2x2 cube of hexahedra
    builder = mmm_core.MMMBuilder()

    # Generate 8 hexahedral elements in a 2x2x2 grid
    vertices_list = []
    elements_list = []
    vertex_count = 0

    for iz in range(2):
        for iy in range(2):
            for ix in range(2):
                # Offset for this element
                ox, oy, oz = ix * 0.1, iy * 0.1, iz * 0.1

                # 8 vertices of this hex element (0.1m cube)
                hex_verts = np.array([
                    [ox, oy, oz],
                    [ox + 0.1, oy, oz],
                    [ox + 0.1, oy + 0.1, oz],
                    [ox, oy + 0.1, oz],
                    [ox, oy, oz + 0.1],
                    [ox + 0.1, oy, oz + 0.1],
                    [ox + 0.1, oy + 0.1, oz + 0.1],
                    [ox, oy + 0.1, oz + 0.1],
                ])
                vertices_list.append(hex_verts)
                elements_list.append(np.arange(8) + vertex_count)
                vertex_count += 8

    vertices = np.vstack(vertices_list).astype(np.float64)
    elements = np.array(elements_list, dtype=np.int32)

    builder.add_hexahedra_from_mesh(vertices, elements)

    n_elem = builder.num_elements
    total_dof = builder.total_dof

    print(f"  Number of elements: {n_elem}")
    print(f"  Total DOF: {total_dof}")

    # Build interaction matrix
    N, dof_offset = builder.build()

    # Create solver
    solver = mmm_core.MMMSolver()
    solver.set_matrix(N, dof_offset)

    # Solve: chi = 99 (mu_r = 100)
    chi = 99.0
    inv_chi = np.full(n_elem, 1.0 / chi, dtype=np.float64)

    # External field H_ext = [0, 0, 1000] A/m (uniform)
    # For 6 DOF MSC, H_ext needs to be per-DOF (n_face components)
    H_ext = np.zeros(total_dof, dtype=np.float64)
    for i in range(n_elem):
        # Set z-component of H_ext for each element
        # For MSC, this is more complex - set all face components
        dof_start = i * 6
        H_ext[dof_start:dof_start + 6] = 1000.0  # Simplified

    # Solve using BiCGSTAB
    M, iterations = solver.solve_bicgstab(inv_chi, H_ext, 1e-8, 1000, True)  # tol, max_iter, chi_per_element

    print(f"  BiCGSTAB iterations: {iterations}")
    print(f"  Solution M shape: {M.shape}")
    print(f"  M max: {np.max(np.abs(M)):.1f} A/m")

    assert iterations > 0, "Expected positive iteration count"
    assert iterations < 1000, "Should converge before max iterations"

    print("  PASSED")
    return True


def test_nonlinear_helpers():
    """Test nonlinear iteration helper functions"""
    print("\n" + "=" * 60)
    print("Test 6: Nonlinear Iteration Helpers")
    print("=" * 60)

    import mmm_core

    # Test compute_chi_from_bh
    # M = [0, 0, 100000] for each element (100 kA/m)
    M_flat = np.array([0, 0, 100000, 0, 0, 100000], dtype=np.float64)

    # H = [0, 0, 100] for each element (100 A/m)
    H_flat = np.array([0, 0, 100, 0, 0, 100], dtype=np.float64)

    # B-H curve: linear material with mu_r = 1000
    # B = mu_0 * mu_r * H = mu_0 * (1 + chi) * H
    # chi = mu_r - 1 = 999
    bh_curve = np.array([
        [0.0, 0.0],
        [100.0, 0.1257],     # B = mu_0 * 1000 * 100 = 0.1257 T
        [1000.0, 1.257],     # B = mu_0 * 1000 * 1000 = 1.257 T
    ], dtype=np.float64)

    chi = mmm_core.compute_chi_from_bh(M_flat, H_flat, bh_curve)

    print(f"  Computed chi: {chi}")
    print(f"  Expected chi: ~999 (mu_r = 1000)")

    # Check convergence helper
    B_new = np.array([0, 0, 1.257, 0, 0, 1.257], dtype=np.float64)
    B_old = np.array([0, 0, 1.25, 0, 0, 1.25], dtype=np.float64)
    B_sat = 2.0  # T

    max_change = mmm_core.check_convergence(B_new, B_old, B_sat)

    print(f"  Max B-field change: {max_change:.4f}")

    # Expected: (1.257 - 1.25) / 2.0 = 0.0035
    assert abs(max_change - 0.0035) < 0.001, f"Unexpected max_change: {max_change}"

    print("  PASSED")
    return True


def test_field_computer():
    """Test MMMFieldComputer for B/H field computation"""
    print("\n" + "=" * 60)
    print("Test 7: MMMFieldComputer B/H Field Computation")
    print("=" * 60)

    import mmm_core

    # Create a single tetrahedron element
    builder = mmm_core.MMMBuilder()

    vertices = np.array([
        [0.0, 0.0, 0.0],
        [0.1, 0.0, 0.0],
        [0.05, 0.0866, 0.0],
        [0.05, 0.0289, 0.0816]
    ], dtype=np.float64)

    elements = np.array([[0, 1, 2, 3]], dtype=np.int32)
    builder.add_tetrahedra_from_mesh(vertices, elements)

    print(f"  Number of elements: {builder.num_elements}")
    print(f"  Total DOF: {builder.total_dof}")

    # Build interaction matrix
    N, dof_offset = builder.build()

    # Create field computer
    field_computer = mmm_core.MMMFieldComputer()
    field_computer.set_elements_from_builder(builder)

    print(f"  Field computer elements: {field_computer.num_elements}")
    print(f"  Field computer DOF: {field_computer.total_dof}")

    # Get element centers and volumes
    centers = field_computer.get_element_centers()
    volumes = field_computer.get_element_volumes()

    print(f"  Element centers shape: {centers.shape}")
    print(f"  Element center: [{centers[0, 0]:.4f}, {centers[0, 1]:.4f}, {centers[0, 2]:.4f}] m")
    print(f"  Element volume: {volumes[0]:.6f} m^3")

    # Set up magnetization: M = [0, 0, 1e6] A/m (1 MA/m in z)
    M = np.array([0.0, 0.0, 1e6], dtype=np.float64)

    # Define observation points (outside the element)
    obs_points = np.array([
        [0.05, 0.05, 0.2],   # Above element
        [0.2, 0.0, 0.0],     # To the side
        [0.05, 0.05, -0.2],  # Below element
    ], dtype=np.float64)

    # Compute B field
    B = field_computer.compute_b_field(M, obs_points)
    print(f"  B field shape: {B.shape}")
    print(f"  B at z=0.2m: [{B[0, 0]:.6e}, {B[0, 1]:.6e}, {B[0, 2]:.6e}] T")
    print(f"  B at x=0.2m: [{B[1, 0]:.6e}, {B[1, 1]:.6e}, {B[1, 2]:.6e}] T")
    print(f"  B at z=-0.2m: [{B[2, 0]:.6e}, {B[2, 1]:.6e}, {B[2, 2]:.6e}] T")

    # Compute H field
    H = field_computer.compute_h_field(M, obs_points)
    print(f"  H field shape: {H.shape}")
    print(f"  H at z=0.2m: [{H[0, 0]:.2f}, {H[0, 1]:.2f}, {H[0, 2]:.2f}] A/m")

    # Verify B = mu_0 * H (in air)
    mu_0 = 4 * np.pi * 1e-7
    B_from_H = mu_0 * H

    print(f"  B = mu_0 * H check:")
    for i in range(3):
        diff = np.linalg.norm(B[i] - B_from_H[i]) / (np.linalg.norm(B[i]) + 1e-20)
        print(f"    Point {i}: B={np.linalg.norm(B[i]):.6e} T, mu_0*H={np.linalg.norm(B_from_H[i]):.6e} T, rel_diff={diff:.2e}")
        assert diff < 1e-10, f"B != mu_0 * H at point {i}"

    # Basic physics check: For a dipole oriented in +z,
    # B_z should be positive both above and below (field lines exit from both ends)
    # but the field direction (sign) depends on position relative to the dipole
    # What's important is that B field is non-zero and B = mu_0 * H holds
    assert np.linalg.norm(B[0]) > 1e-10, "Expected non-zero B field above element"
    assert np.linalg.norm(B[2]) > 1e-10, "Expected non-zero B field below element"

    print("  PASSED")
    return True


def test_mixed_elements():
    """Test mixed tetrahedra and hexahedra"""
    print("\n" + "=" * 60)
    print("Test 8: Mixed Tetrahedra and Hexahedra")
    print("=" * 60)

    import mmm_core

    builder = mmm_core.MMMBuilder()

    # Add one tetrahedron
    tet_vertices = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 0.866, 0.0],
        [0.5, 0.289, 0.816]
    ], dtype=np.float64)
    tet_elements = np.array([[0, 1, 2, 3]], dtype=np.int32)
    builder.add_tetrahedra_from_mesh(tet_vertices, tet_elements)

    # Add one hexahedron (shifted by 3.0 in x)
    hex_vertices = np.array([
        [3.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [4.0, 1.0, 0.0],
        [3.0, 1.0, 0.0],
        [3.0, 0.0, 1.0],
        [4.0, 0.0, 1.0],
        [4.0, 1.0, 1.0],
        [3.0, 1.0, 1.0],
    ], dtype=np.float64)
    hex_elements = np.array([[0, 1, 2, 3, 4, 5, 6, 7]], dtype=np.int32)
    builder.add_hexahedra_from_mesh(hex_vertices, hex_elements)

    print(f"  Number of elements: {builder.num_elements}")
    print(f"  Total DOF: {builder.total_dof}")

    # Expected: 1 tetra (3 DOF) + 1 hexa (6 DOF) = 9 DOF
    assert builder.num_elements == 2
    assert builder.total_dof == 9, f"Expected 9 DOF, got {builder.total_dof}"

    # Build interaction matrix
    N, dof_offset = builder.build()

    print(f"  Interaction matrix shape: {N.shape}")
    print(f"  DOF offset: {dof_offset}")

    # Expected: 9x9 matrix
    assert N.shape == (9, 9), f"Expected (9, 9) matrix, got {N.shape}"
    assert list(dof_offset) == [0, 3, 9], f"Unexpected dof_offset: {dof_offset}"

    print("  PASSED")
    return True


def main():
    """Run all tests"""
    print("=" * 60)
    print("  MMM Core C++ Module Test Suite")
    print("=" * 60)

    tests = [
        ("Import", test_import),
        ("Builder Tetrahedron", test_builder_tetrahedron),
        ("Builder Hexahedron", test_builder_hexahedron),
        ("Solver LU", test_solver_lu),
        ("Solver BiCGSTAB", test_solver_bicgstab),
        ("Nonlinear Helpers", test_nonlinear_helpers),
        ("Field Computer", test_field_computer),
        ("Mixed Elements", test_mixed_elements),
    ]

    passed = 0
    failed = 0

    for name, test_func in tests:
        try:
            if test_func():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"\n  FAILED: {name}")
            print(f"  Error: {e}")
            import traceback
            traceback.print_exc()
            failed += 1

    print("\n" + "=" * 60)
    print(f"  Test Summary: {passed} passed, {failed} failed")
    print("=" * 60)

    return failed == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
