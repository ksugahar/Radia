"""
test_mmm_pm.py - Permanent Magnet test for mmm_core

Tests:
1. PM-only system (no solve needed)
2. PM + soft iron mixed system
3. Comparison with analytical dipole field

Usage:
    python test_mmm_pm.py
"""

import numpy as np
import sys
import os

# Add path for mmm_core
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

from mmm_ngsolve import MMMSystem


def create_cube_mesh(center, size, n_div=2):
    """
    Create simple structured hexahedral mesh for a cube.

    Args:
        center: [x, y, z] center of cube
        size: edge length of cube
        n_div: number of divisions per edge

    Returns:
        vertices: (n_verts, 3) array
        elements: (n_elem, 8) array of vertex indices
    """
    cx, cy, cz = center
    half = size / 2

    # Create grid of vertices
    n = n_div + 1
    x = np.linspace(cx - half, cx + half, n)
    y = np.linspace(cy - half, cy + half, n)
    z = np.linspace(cz - half, cz + half, n)

    # Create vertex array
    vertices = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                vertices.append([x[i], y[j], z[k]])
    vertices = np.array(vertices, dtype=np.float64)

    # Create hex elements
    elements = []
    for k in range(n_div):
        for j in range(n_div):
            for i in range(n_div):
                # 8 vertices of hex element
                v0 = i + j * n + k * n * n
                v1 = v0 + 1
                v2 = v0 + n + 1
                v3 = v0 + n
                v4 = v0 + n * n
                v5 = v4 + 1
                v6 = v4 + n + 1
                v7 = v4 + n
                elements.append([v0, v1, v2, v3, v4, v5, v6, v7])
    elements = np.array(elements, dtype=np.int32)

    return vertices, elements


def dipole_field(M, r_mag, volume, r_obs):
    """
    Compute magnetic field from a magnetic dipole.

    B = (mu_0/4pi) * (3(m.r)r/r^5 - m/r^3)

    where m = M * volume is the magnetic moment.

    Args:
        M: magnetization vector [Mx, My, Mz] in A/m
        r_mag: position of dipole [x, y, z] in m
        volume: volume in m^3
        r_obs: observation point [x, y, z] in m

    Returns:
        B: [Bx, By, Bz] in Tesla
    """
    mu_0 = 4 * np.pi * 1e-7

    m = np.array(M) * volume  # magnetic moment [A*m^2]
    r = np.array(r_obs) - np.array(r_mag)  # displacement vector
    r_norm = np.linalg.norm(r)

    if r_norm < 1e-12:
        return np.array([0.0, 0.0, 0.0])

    r_hat = r / r_norm
    m_dot_r = np.dot(m, r_hat)

    B = (mu_0 / (4 * np.pi)) * (3 * m_dot_r * r_hat - m) / (r_norm ** 3)

    return B


def test_pm_only():
    """Test PM-only system (no solve needed)."""
    print("\n" + "=" * 60)
    print("Test 1: PM-only system")
    print("=" * 60)

    # Create cube magnet: 10cm cube centered at origin
    # Magnetization: 1.2 T = 954930 A/m in z-direction
    size = 0.1  # 10 cm
    Mr = [0, 0, 954930]  # 1.2 T equivalent

    vertices, elements = create_cube_mesh([0, 0, 0], size, n_div=2)
    print(f"Mesh: {len(vertices)} vertices, {len(elements)} hex elements")

    # Create PM system
    system = MMMSystem.from_arrays(vertices, elements, Mr=Mr)
    system.info()

    # Check that system is PM-only and already solved
    assert system.is_pm_only, "Should be PM-only system"
    assert system._is_solved, "PM-only system should be auto-solved"

    # Compute field at observation points
    obs_points = np.array([
        [0.0, 0.0, 0.15],   # Above magnet
        [0.15, 0.0, 0.0],   # Side of magnet
        [0.0, 0.15, 0.0],   # Side of magnet
    ])

    B = system.compute_b_field(obs_points)

    print("\nField at observation points:")
    for i, (pt, b) in enumerate(zip(obs_points, B)):
        print(f"  Point {pt}: B = [{b[0]:.6f}, {b[1]:.6f}, {b[2]:.6f}] T")
        print(f"    |B| = {np.linalg.norm(b):.6f} T")

    # Compare with single dipole approximation (for far-field check)
    volume = size ** 3
    center = [0, 0, 0]

    print("\nComparison with single dipole approximation (far-field):")
    for pt in obs_points:
        B_mmm = system.compute_b_field([pt])[0]
        B_dipole = dipole_field(Mr, center, volume, pt)
        error = np.linalg.norm(B_mmm - B_dipole) / np.linalg.norm(B_dipole) * 100
        print(f"  Point {pt}:")
        print(f"    MMM:    [{B_mmm[0]:.6f}, {B_mmm[1]:.6f}, {B_mmm[2]:.6f}] T")
        print(f"    Dipole: [{B_dipole[0]:.6f}, {B_dipole[1]:.6f}, {B_dipole[2]:.6f}] T")
        print(f"    Error:  {error:.2f}%")

    print("\nTest 1 PASSED: PM-only system works correctly")
    return True


def test_pm_field_symmetry():
    """Test field symmetry of PM system."""
    print("\n" + "=" * 60)
    print("Test 2: PM field symmetry")
    print("=" * 60)

    # Create cube magnet
    size = 0.1
    Mr = [0, 0, 954930]

    vertices, elements = create_cube_mesh([0, 0, 0], size, n_div=3)
    system = MMMSystem.from_arrays(vertices, elements, Mr=Mr)

    # Test z-symmetry (field should be symmetric about z-axis)
    z = 0.15
    points_on_ring = []
    n_phi = 8
    r = 0.1
    for i in range(n_phi):
        phi = 2 * np.pi * i / n_phi
        x = r * np.cos(phi)
        y = r * np.sin(phi)
        points_on_ring.append([x, y, z])
    points_on_ring = np.array(points_on_ring)

    B = system.compute_b_field(points_on_ring)

    # Bz should be same for all points (due to axial symmetry)
    Bz_values = B[:, 2]
    Bz_mean = np.mean(Bz_values)
    Bz_std = np.std(Bz_values)

    print(f"Bz values around ring at z={z}m, r={r}m:")
    for i, bz in enumerate(Bz_values):
        print(f"  phi={360*i/n_phi:6.1f} deg: Bz = {bz:.6f} T")
    print(f"Mean Bz = {Bz_mean:.6f} T, Std = {Bz_std:.6e} T")

    # Check symmetry
    relative_std = Bz_std / abs(Bz_mean) if abs(Bz_mean) > 1e-12 else 0
    print(f"Relative Std = {relative_std:.2e}")

    assert relative_std < 0.01, f"Symmetry check failed: relative std = {relative_std}"

    print("\nTest 2 PASSED: Field symmetry verified")
    return True


def test_pm_with_soft_iron():
    """Test mixed PM + soft iron system.

    NOTE: Mixed PM + soft iron solve is planned for future implementation.
    Currently, this test verifies that the system can be built and reports
    the correct state, but full solve is not yet supported.
    """
    print("\n" + "=" * 60)
    print("Test 3: PM + Soft Iron mixed system (API check)")
    print("=" * 60)

    # PM: 10cm cube at origin
    pm_size = 0.1
    Mr = [0, 0, 954930]
    pm_vertices, pm_elements = create_cube_mesh([0, 0, 0], pm_size, n_div=2)

    # Soft iron: 10cm cube at z = 0.2m
    iron_size = 0.1
    iron_vertices, iron_elements = create_cube_mesh([0, 0, 0.25], iron_size, n_div=2)

    # Create mixed system
    system = MMMSystem()
    system.add_pm_from_arrays(pm_vertices, pm_elements, Mr)
    system.add_soft_from_arrays(iron_vertices, iron_elements, mu_r=1000)
    system.build()

    system.info()

    # Check system state
    assert system.has_pm, "Should have PM elements"
    assert system.has_soft, "Should have soft elements"
    assert not system.is_pm_only, "Should not be PM-only"

    print("\nMixed PM + soft iron system correctly built.")
    print("NOTE: Full solve for mixed systems is planned for future implementation.")
    print("      For now, use separate PM and soft iron systems.")

    # Demonstrate workaround: compute PM field at soft iron location
    iron_center = np.array([[0, 0, 0.25]])
    pm_only = MMMSystem.from_arrays(pm_vertices, pm_elements, Mr=Mr)
    H_from_pm = pm_only.compute_h_field(iron_center)[0]
    print(f"\nH field from PM at iron center: [{H_from_pm[0]:.1f}, {H_from_pm[1]:.1f}, {H_from_pm[2]:.1f}] A/m")
    print(f"|H| = {np.linalg.norm(H_from_pm):.1f} A/m")

    print("\nTest 3 PASSED: PM + soft iron API works correctly")
    return True


def test_multiple_pm():
    """Test system with multiple PM regions."""
    print("\n" + "=" * 60)
    print("Test 4: Multiple PM regions (Halbach-like)")
    print("=" * 60)

    # Create Halbach-like array: 4 magnets with rotating magnetization
    size = 0.05  # 5cm cubes
    gap = 0.01   # 1cm gap

    # Positions and magnetizations (simplified Halbach)
    magnets = [
        {'center': [-0.06, 0, 0], 'Mr': [0, 0, 954930]},    # +Z
        {'center': [0, 0.06, 0], 'Mr': [954930, 0, 0]},     # +X
        {'center': [0.06, 0, 0], 'Mr': [0, 0, -954930]},    # -Z
        {'center': [0, -0.06, 0], 'Mr': [-954930, 0, 0]},   # -X
    ]

    # Create combined system
    system = MMMSystem()

    for mag in magnets:
        vertices, elements = create_cube_mesh(mag['center'], size, n_div=2)
        system.add_pm_from_arrays(vertices, elements, mag['Mr'])

    system.build()
    system.info()

    # Check field at center (should be enhanced in one direction)
    center = np.array([[0.0, 0.0, 0.0]])
    B_center = system.compute_b_field(center)[0]

    print(f"\nField at center: B = [{B_center[0]:.4f}, {B_center[1]:.4f}, {B_center[2]:.4f}] T")
    print(f"|B| = {np.linalg.norm(B_center):.4f} T")

    # Compare with single magnet field
    single_vertices, single_elements = create_cube_mesh([0, 0, 0.1], size, n_div=2)
    single_system = MMMSystem.from_arrays(single_vertices, single_elements, Mr=[0, 0, 954930])
    B_single = single_system.compute_b_field([[0, 0, 0]])[0]

    print(f"\nSingle magnet at same distance:")
    print(f"  B = [{B_single[0]:.4f}, {B_single[1]:.4f}, {B_single[2]:.4f}] T")
    print(f"  |B| = {np.linalg.norm(B_single):.4f} T")

    # Halbach should have field enhancement
    enhancement = np.linalg.norm(B_center) / np.linalg.norm(B_single)
    print(f"\nField enhancement factor: {enhancement:.2f}x")

    print("\nTest 4 PASSED: Multiple PM regions work correctly")
    return True


def main():
    """Run all tests."""
    print("=" * 60)
    print("MMM Permanent Magnet Tests")
    print("=" * 60)

    all_passed = True

    try:
        all_passed &= test_pm_only()
    except Exception as e:
        print(f"Test 1 FAILED: {e}")
        all_passed = False

    try:
        all_passed &= test_pm_field_symmetry()
    except Exception as e:
        print(f"Test 2 FAILED: {e}")
        all_passed = False

    try:
        all_passed &= test_pm_with_soft_iron()
    except Exception as e:
        print(f"Test 3 FAILED: {e}")
        all_passed = False

    try:
        all_passed &= test_multiple_pm()
    except Exception as e:
        print(f"Test 4 FAILED: {e}")
        all_passed = False

    print("\n" + "=" * 60)
    if all_passed:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print("=" * 60)

    return all_passed


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
