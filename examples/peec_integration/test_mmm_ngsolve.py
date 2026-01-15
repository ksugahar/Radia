"""
Test script for mmm_ngsolve module (NGSolve integration for mmm_core)

Tests:
1. MMMSystem.from_ngsolve_mesh: Build system from NGSolve mesh
2. MMMSystem.solve: Solve with external field
3. MMMFieldCF: CoefficientFunction evaluation
4. Comparison with Radia (optional)

Usage:
    python test_mmm_ngsolve.py
"""

import sys
import os
import numpy as np

# Add paths
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build-msvc'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))


def test_import():
    """Test basic module import"""
    print("=" * 60)
    print("Test 1: Module Import")
    print("=" * 60)

    try:
        from mmm_ngsolve import MMMSystem, MMMFieldCF, create_mmm_from_mesh
        print("  mmm_ngsolve module imported successfully")
        print(f"    MMMSystem: {MMMSystem}")
        print(f"    MMMFieldCF: {MMMFieldCF}")
        print(f"    create_mmm_from_mesh: {create_mmm_from_mesh}")
        return True
    except ImportError as e:
        print(f"  FAILED: {e}")
        return False


def test_from_arrays():
    """Test MMMSystem creation from numpy arrays"""
    print("\n" + "=" * 60)
    print("Test 2: MMMSystem from Arrays")
    print("=" * 60)

    from mmm_ngsolve import MMMSystem

    # Create simple tetrahedral mesh (single element)
    vertices = np.array([
        [0.0, 0.0, 0.0],
        [0.1, 0.0, 0.0],
        [0.05, 0.0866, 0.0],
        [0.05, 0.0289, 0.0816]
    ], dtype=np.float64)

    elements = np.array([[0, 1, 2, 3]], dtype=np.int32)

    # Create system
    system = MMMSystem.from_arrays(vertices, elements, mu_r=1000)

    print(f"  Number of elements: {system.num_elements}")
    print(f"  Total DOF: {system.total_dof}")

    assert system.num_elements == 1
    assert system.total_dof == 3

    # Check element properties
    centers = system.get_element_centers()
    volumes = system.get_element_volumes()

    print(f"  Element center: [{centers[0, 0]:.4f}, {centers[0, 1]:.4f}, {centers[0, 2]:.4f}] m")
    print(f"  Element volume: {volumes[0]:.6e} m^3")

    assert volumes[0] > 0, "Volume must be positive"

    print("  PASSED")
    return True


def test_solve():
    """Test solving with external field"""
    print("\n" + "=" * 60)
    print("Test 3: MMMSystem Solve")
    print("=" * 60)

    from mmm_ngsolve import MMMSystem

    # Create 2x2x2 grid of tetrahedral elements (8 elements)
    vertices_list = []
    elements_list = []
    vertex_count = 0

    # Create grid of 8 tetrahedra
    grid_size = 0.1  # 10 cm
    for iz in range(2):
        for iy in range(2):
            for ix in range(2):
                ox = ix * grid_size
                oy = iy * grid_size
                oz = iz * grid_size

                # Simple tetrahedron at each grid position
                tet_verts = np.array([
                    [ox, oy, oz],
                    [ox + grid_size, oy, oz],
                    [ox + grid_size/2, oy + grid_size*0.866, oz],
                    [ox + grid_size/2, oy + grid_size*0.289, oz + grid_size*0.816]
                ])
                vertices_list.append(tet_verts)
                elements_list.append(np.arange(4) + vertex_count)
                vertex_count += 4

    vertices = np.vstack(vertices_list).astype(np.float64)
    elements = np.array(elements_list, dtype=np.int32)

    # Create and solve system
    system = MMMSystem.from_arrays(vertices, elements, mu_r=100)

    print(f"  Number of elements: {system.num_elements}")
    print(f"  Total DOF: {system.total_dof}")

    # Solve with uniform external field
    H_ext = np.array([0, 0, 1000], dtype=np.float64)  # 1000 A/m in z
    M = system.solve(H_ext, method='lu')

    print(f"  Solution M shape: {M.shape}")
    print(f"  M max: {np.max(np.abs(M)):.1f} A/m")

    # Check that magnetization is reasonable
    # For chi = 99, M should be roughly chi * H = 99 * 1000 = 99000 A/m
    # But with demagnetization, actual M will be lower
    assert np.max(np.abs(M)) > 100, "Magnetization too small"
    assert np.max(np.abs(M)) < 1e8, "Magnetization too large"

    print("  PASSED")
    return True


def test_field_computation():
    """Test B/H field computation"""
    print("\n" + "=" * 60)
    print("Test 4: Field Computation")
    print("=" * 60)

    from mmm_ngsolve import MMMSystem

    # Create single tetrahedron
    vertices = np.array([
        [0.0, 0.0, 0.0],
        [0.1, 0.0, 0.0],
        [0.05, 0.0866, 0.0],
        [0.05, 0.0289, 0.0816]
    ], dtype=np.float64)

    elements = np.array([[0, 1, 2, 3]], dtype=np.int32)

    system = MMMSystem.from_arrays(vertices, elements, mu_r=1000)
    system.solve([0, 0, 1e5])  # Strong field to get significant M

    # Observation points
    obs_points = np.array([
        [0.05, 0.05, 0.2],   # Above element
        [0.2, 0.0, 0.0],     # To the side
        [0.05, 0.05, -0.2],  # Below element
    ], dtype=np.float64)

    # Compute fields
    B = system.compute_b_field(obs_points)
    H = system.compute_h_field(obs_points)

    print(f"  B field shape: {B.shape}")
    print(f"  H field shape: {H.shape}")
    print(f"  B at z=0.2m: [{B[0, 0]:.6e}, {B[0, 1]:.6e}, {B[0, 2]:.6e}] T")
    print(f"  H at z=0.2m: [{H[0, 0]:.2f}, {H[0, 1]:.2f}, {H[0, 2]:.2f}] A/m")

    # Verify B = mu_0 * H in air
    mu_0 = 4 * np.pi * 1e-7
    for i in range(3):
        B_from_H = mu_0 * H[i]
        diff = np.linalg.norm(B[i] - B_from_H) / (np.linalg.norm(B[i]) + 1e-20)
        assert diff < 1e-10, f"B != mu_0 * H at point {i}"

    print("  B = mu_0 * H verified")
    print("  PASSED")
    return True


def test_coefficient_function():
    """Test MMMFieldCF coefficient function"""
    print("\n" + "=" * 60)
    print("Test 5: MMMFieldCF Coefficient Function")
    print("=" * 60)

    from mmm_ngsolve import MMMSystem, MMMFieldCF

    # Create and solve system
    vertices = np.array([
        [0.0, 0.0, 0.0],
        [0.1, 0.0, 0.0],
        [0.05, 0.0866, 0.0],
        [0.05, 0.0289, 0.0816]
    ], dtype=np.float64)

    elements = np.array([[0, 1, 2, 3]], dtype=np.int32)

    system = MMMSystem.from_arrays(vertices, elements, mu_r=1000)
    system.solve([0, 0, 1e5])

    # Get CoefficientFunctions
    B_cf = system.get_b_field_cf()
    H_cf = system.get_h_field_cf()

    print(f"  B_cf: {B_cf}")
    print(f"  H_cf: {H_cf}")

    # Test point evaluation
    test_point = [0.05, 0.05, 0.2]
    B_val = B_cf(test_point)
    H_val = H_cf(test_point)

    print(f"  B at {test_point}: {B_val}")
    print(f"  H at {test_point}: {H_val}")

    # Test batch evaluation
    obs_points = np.array([
        [0.05, 0.05, 0.2],
        [0.2, 0.0, 0.0],
    ])
    B_batch = B_cf.Evaluate(obs_points)
    H_batch = H_cf.Evaluate(obs_points)

    print(f"  Batch B shape: {B_batch.shape}")
    print(f"  Batch H shape: {H_batch.shape}")

    # Verify consistency between point and batch evaluation
    assert np.allclose(B_val, B_batch[0]), "Inconsistent point/batch evaluation"

    print("  PASSED")
    return True


def test_ngsolve_mesh():
    """Test with actual NGSolve mesh (requires ngsolve)"""
    print("\n" + "=" * 60)
    print("Test 6: NGSolve Mesh Integration")
    print("=" * 60)

    try:
        from ngsolve import Mesh, VOL
        from netgen.occ import Box, Pnt, OCCGeometry
        from netgen.meshing import MeshingParameters
    except ImportError as e:
        print(f"  SKIPPED: NGSolve not available ({e})")
        return True  # Don't fail if NGSolve not installed

    from mmm_ngsolve import MMMSystem

    # Create simple cube mesh
    print("  Creating NGSolve mesh...")
    box = Box(Pnt(0, 0, 0), Pnt(0.1, 0.1, 0.1))
    box.mat('magnetic')
    geo = OCCGeometry(box)

    mp = MeshingParameters(maxh=0.05)
    ngmesh = geo.GenerateMesh(mp)
    mesh = Mesh(ngmesh)

    print(f"  Mesh vertices: {mesh.nv}")
    print(f"  Mesh elements: {sum(1 for _ in mesh.Elements(VOL))}")

    # Create MMM system from NGSolve mesh
    system = MMMSystem.from_ngsolve_mesh(mesh, mu_r=1000, material_filter='magnetic')

    print(f"  MMM elements: {system.num_elements}")
    print(f"  MMM DOF: {system.total_dof}")

    # Solve
    H_ext = [0, 0, 1e4]  # 10 kA/m
    M = system.solve(H_ext)

    print(f"  Solution M max: {np.max(np.abs(M)):.1f} A/m")

    # Compute field at a point outside the cube
    obs_point = np.array([[0.05, 0.05, 0.15]])
    B = system.compute_b_field(obs_point)
    H = system.compute_h_field(obs_point)

    print(f"  B at (0.05, 0.05, 0.15): |B| = {np.linalg.norm(B[0]):.6e} T")
    print(f"  H at (0.05, 0.05, 0.15): |H| = {np.linalg.norm(H[0]):.2f} A/m")

    print("  PASSED")
    return True


def test_convenience_function():
    """Test create_mmm_from_mesh convenience function"""
    print("\n" + "=" * 60)
    print("Test 7: Convenience Function")
    print("=" * 60)

    try:
        from ngsolve import Mesh
        from netgen.occ import Box, Pnt, OCCGeometry
    except ImportError as e:
        print(f"  SKIPPED: NGSolve not available ({e})")
        return True

    from mmm_ngsolve import create_mmm_from_mesh

    # Create mesh
    box = Box(Pnt(0, 0, 0), Pnt(0.1, 0.1, 0.1))
    box.mat('iron')
    geo = OCCGeometry(box)
    mesh = Mesh(geo.GenerateMesh(maxh=0.05))

    # Use convenience function
    system = create_mmm_from_mesh(
        mesh,
        mu_r=500,
        H_ext=[0, 0, 5000],
        material_filter='iron'
    )

    print(f"  System elements: {system.num_elements}")
    print(f"  System solved: {system._is_solved}")

    # Get field CF
    B_cf = system.get_b_field_cf()
    B_val = B_cf([0.05, 0.05, 0.15])

    print(f"  B at (0.05, 0.05, 0.15): {B_val}")

    print("  PASSED")
    return True


def main():
    """Run all tests"""
    print("=" * 60)
    print("  MMM NGSolve Integration Test Suite")
    print("=" * 60)

    tests = [
        ("Import", test_import),
        ("From Arrays", test_from_arrays),
        ("Solve", test_solve),
        ("Field Computation", test_field_computation),
        ("Coefficient Function", test_coefficient_function),
        ("NGSolve Mesh", test_ngsolve_mesh),
        ("Convenience Function", test_convenience_function),
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
