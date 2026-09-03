"""
Test analytical magnet field computations against magpylib.

This test verifies the Python implementation of analytical magnetic field
computation for various magnet shapes.
"""

import sys

import numpy as np

from radia.analytical_magnet import CuboidMagnet, CurrentLoop, SphericalMagnet

# Try to import magpylib for comparison
try:
    import magpylib as magpy
    HAS_MAGPYLIB = True
except ImportError:
    HAS_MAGPYLIB = False
    print("Warning: magpylib not installed, skipping comparison tests")

# Physical constants
PI = np.pi
MU0 = 4.0 * np.pi * 1.0e-7  # H/m


def test_spherical_magnet():
    """Test spherical magnet B-field."""
    print("=" * 60)
    print("Test: Spherical Magnet B-field")
    print("=" * 60)

    # Sphere parameters
    D = 20.0  # mm, diameter
    Mz = 1e6  # A/m, magnetization

    J = MU0 * Mz  # Polarization in Tesla
    print(f"Sphere: D={D:.0f} mm, Mz={Mz:.0e} A/m (J={J:.4f} T)")

    # Create our sphere
    sphere_ours = SphericalMagnet(
        center=[0, 0, 0],
        diameter=D,
        magnetization=[0, 0, Mz]
    )

    # Test points in mm
    test_points_mm = [
        [0, 0, 15],      # On axis, outside
        [0, 0, 5],       # On axis, inside
        [15, 0, 0],      # Off axis, outside
        [10, 10, 10],    # Diagonal, outside
        [5, 5, 0],       # Diagonal, inside
    ]

    if HAS_MAGPYLIB:
        print(f"\n{'Point [mm]':>20} {'Bz_ours':>12} {'Bz_magpy':>12} {'Error [%]':>12}")
        print("-" * 65)

        sphere_magpy = magpy.magnet.Sphere(
            polarization=(0, 0, J),
            diameter=D / 1000.0,  # meters
            position=(0, 0, 0)
        )

        max_error = 0.0
        for pt in test_points_mm:
            B_ours = sphere_ours.get_B(pt)
            pt_m = [p / 1000.0 for p in pt]
            B_magpy = sphere_magpy.getB(pt_m)

            # Compute relative error
            B_mag = np.linalg.norm(B_magpy)
            if B_mag > 1e-15:
                error = np.linalg.norm(np.array(B_ours) - np.array(B_magpy)) / B_mag * 100
                max_error = max(max_error, error)
            else:
                error = 0.0

            print(f"({pt[0]:5.1f}, {pt[1]:5.1f}, {pt[2]:5.1f}) "
                  f"{B_ours[2]:>12.6e} {B_magpy[2]:>12.6e} {error:>12.2f}")

        print("-" * 65)
        assert max_error < 1.0, f"TEST FAILED: Maximum error {max_error:.2f}%"
        print(f"TEST PASSED: Maximum error {max_error:.4f}%")
    else:
        print(f"\n{'Point [mm]':>25} {'Bx [T]':>12} {'By [T]':>12} {'Bz [T]':>12}")
        print("-" * 70)

        for pt in test_points_mm:
            B = sphere_ours.get_B(pt)
            print(f"({pt[0]:6.1f}, {pt[1]:6.1f}, {pt[2]:6.1f}) "
                  f"{B[0]:>12.6e} {B[1]:>12.6e} {B[2]:>12.6e}")

        print("-" * 70)
        print("TEST INFO: Values shown (no magpylib for comparison)")


def test_cuboid_magnet():
    """Test cuboid magnet B-field."""
    print("\n" + "=" * 60)
    print("Test: Cuboid Magnet B-field")
    print("=" * 60)

    # Cuboid parameters
    dims = [20.0, 20.0, 10.0]  # mm
    Mz = 1e6  # A/m, magnetization

    J = MU0 * Mz
    print(f"Cuboid: {dims[0]:.0f}x{dims[1]:.0f}x{dims[2]:.0f} mm, Mz={Mz:.0e} A/m (J={J:.4f} T)")

    # Create our cuboid
    cuboid_ours = CuboidMagnet(
        center=[0, 0, 0],
        dimensions=dims,
        magnetization=[0, 0, Mz]
    )

    # Test points in mm (outside the cuboid)
    test_points_mm = [
        [0, 0, 15],      # On axis, above
        [15, 0, 0],      # Side
        [15, 15, 0],     # Corner region
        [10, 10, 15],    # Diagonal
    ]

    if HAS_MAGPYLIB:
        print(f"\n{'Point [mm]':>20} {'Bz_ours':>12} {'Bz_magpy':>12} {'Error [%]':>12}")
        print("-" * 65)

        cuboid_magpy = magpy.magnet.Cuboid(
            polarization=(0, 0, J),
            dimension=[d / 1000.0 for d in dims],  # meters
            position=(0, 0, 0)
        )

        max_error = 0.0
        for pt in test_points_mm:
            B_ours = cuboid_ours.get_B(pt)
            pt_m = [p / 1000.0 for p in pt]
            B_magpy = cuboid_magpy.getB(pt_m)

            B_mag = np.linalg.norm(B_magpy)
            if B_mag > 1e-15:
                error = np.linalg.norm(np.array(B_ours) - np.array(B_magpy)) / B_mag * 100
                max_error = max(max_error, error)
            else:
                error = 0.0

            print(f"({pt[0]:5.1f}, {pt[1]:5.1f}, {pt[2]:5.1f}) "
                  f"{B_ours[2]:>12.6e} {B_magpy[2]:>12.6e} {error:>12.2f}")

        print("-" * 65)
        assert max_error < 1.0, f"TEST FAILED: Maximum error {max_error:.2f}%"
        print(f"TEST PASSED: Maximum error {max_error:.4f}%")
    else:
        print(f"\n{'Point [mm]':>25} {'Bx [T]':>12} {'By [T]':>12} {'Bz [T]':>12}")
        print("-" * 70)

        for pt in test_points_mm:
            B = cuboid_ours.get_B(pt)
            print(f"({pt[0]:6.1f}, {pt[1]:6.1f}, {pt[2]:6.1f}) "
                  f"{B[0]:>12.6e} {B[1]:>12.6e} {B[2]:>12.6e}")

        print("-" * 70)
        print("TEST INFO: Values shown (no magpylib for comparison)")


def test_current_loop():
    """Test current loop B-field."""
    print("\n" + "=" * 60)
    print("Test: Current Loop B-field")
    print("=" * 60)

    # Loop parameters
    D = 50.0  # mm, diameter
    I = 100.0  # A, current

    print(f"Current Loop: D={D:.0f} mm, I={I:.0f} A")

    # Create our loop
    loop_ours = CurrentLoop(
        center=[0, 0, 0],
        diameter=D,
        current=I,
        axis='z'
    )

    # Test points in mm
    test_points_mm = [
        [0, 0, 0],       # Center
        [0, 0, 25],      # On axis
        [0, 0, 50],      # Far on axis
        [20, 0, 25],     # Off axis
        [15, 15, 30],    # Diagonal
    ]

    if HAS_MAGPYLIB:
        print(f"\n{'Point [mm]':>20} {'Bz_ours':>12} {'Bz_magpy':>12} {'Error [%]':>12}")
        print("-" * 65)

        loop_magpy = magpy.current.Circle(
            current=I,
            diameter=D / 1000.0,  # meters
            position=(0, 0, 0)
        )

        max_error = 0.0
        for pt in test_points_mm:
            B_ours = loop_ours.get_B(pt)
            pt_m = [p / 1000.0 for p in pt]
            B_magpy = loop_magpy.getB(pt_m)

            B_mag = np.linalg.norm(B_magpy)
            if B_mag > 1e-15:
                error = np.linalg.norm(np.array(B_ours) - np.array(B_magpy)) / B_mag * 100
                max_error = max(max_error, error)
            else:
                error = 0.0

            print(f"({pt[0]:5.1f}, {pt[1]:5.1f}, {pt[2]:5.1f}) "
                  f"{B_ours[2]:>12.6e} {B_magpy[2]:>12.6e} {error:>12.2f}")

        print("-" * 65)
        assert max_error < 1.0, f"TEST FAILED: Maximum error {max_error:.2f}%"
        print(f"TEST PASSED: Maximum error {max_error:.4f}%")
    else:
        print(f"\n{'Point [mm]':>25} {'Bx [T]':>12} {'By [T]':>12} {'Bz [T]':>12}")
        print("-" * 70)

        for pt in test_points_mm:
            B = loop_ours.get_B(pt)
            print(f"({pt[0]:6.1f}, {pt[1]:6.1f}, {pt[2]:6.1f}) "
                  f"{B[0]:>12.6e} {B[1]:>12.6e} {B[2]:>12.6e}")

        print("-" * 70)
        print("TEST INFO: Values shown (no magpylib for comparison)")


def test_vector_potentials():
    """Test vector potential implementations via curl(A) = B."""
    print("\n" + "=" * 60)
    print("Test: Vector Potentials (curl(A) = B)")
    print("=" * 60)

    h = 0.1  # Step for numerical differentiation in mm
    h_m = h / 1000.0  # Step in meters (for A in T*m)

    def numerical_curl(get_A, pt, use_meters=False):
        """Compute curl(A) numerically at point pt.

        If use_meters=True, assumes A is in T*m and divides by h in meters.
        Otherwise assumes A is in T*mm and divides by h in mm.
        """
        A_px = get_A([pt[0] + h, pt[1], pt[2]])
        A_mx = get_A([pt[0] - h, pt[1], pt[2]])
        A_py = get_A([pt[0], pt[1] + h, pt[2]])
        A_my = get_A([pt[0], pt[1] - h, pt[2]])
        A_pz = get_A([pt[0], pt[1], pt[2] + h])
        A_mz = get_A([pt[0], pt[1], pt[2] - h])

        # Choose denominator based on A units
        denom = 2 * h_m if use_meters else 2 * h

        dAz_dy = (A_py[2] - A_my[2]) / denom
        dAy_dz = (A_pz[1] - A_mz[1]) / denom
        dAx_dz = (A_pz[0] - A_mz[0]) / denom
        dAz_dx = (A_px[2] - A_mx[2]) / denom
        dAy_dx = (A_px[1] - A_mx[1]) / denom
        dAx_dy = (A_py[0] - A_my[0]) / denom

        return [
            dAz_dy - dAy_dz,
            dAx_dz - dAz_dx,
            dAy_dx - dAx_dy
        ]

    all_passed = True

    # Test 1: Spherical Magnet (A in T*m, so use meters for curl)
    print("\n1. Spherical Magnet at [15, 0, 0]:")
    sphere = SphericalMagnet(center=[0, 0, 0], diameter=20.0, magnetization=[0, 0, 1e6])
    pt = [15, 0, 0]
    B = sphere.get_B(pt)
    curl_A = numerical_curl(sphere.get_A, pt, use_meters=True)
    B_mag = np.linalg.norm(B)
    error = np.linalg.norm(np.array(curl_A) - np.array(B)) / B_mag * 100 if B_mag > 1e-15 else 0
    print(f"   B       = [{B[0]:.6e}, {B[1]:.6e}, {B[2]:.6e}]")
    print(f"   curl(A) = [{curl_A[0]:.6e}, {curl_A[1]:.6e}, {curl_A[2]:.6e}]")
    print(f"   Error: {error:.2f}%")
    if error > 5.0:
        all_passed = False

    # Test 2: Current Loop (A in T*m, so use meters for curl)
    print("\n2. Current Loop at [10, 0, 25]:")
    loop = CurrentLoop(center=[0, 0, 0], diameter=50.0, current=100.0)
    pt = [10, 0, 25]
    B = loop.get_B(pt)
    curl_A = numerical_curl(loop.get_A, pt, use_meters=True)
    B_mag = np.linalg.norm(B)
    error = np.linalg.norm(np.array(curl_A) - np.array(B)) / B_mag * 100 if B_mag > 1e-15 else 0
    print(f"   B       = [{B[0]:.6e}, {B[1]:.6e}, {B[2]:.6e}]")
    print(f"   curl(A) = [{curl_A[0]:.6e}, {curl_A[1]:.6e}, {curl_A[2]:.6e}]")
    print(f"   Error: {error:.2f}%")
    if error > 5.0:
        all_passed = False

    # Test 3: Cuboid Magnet (A in T*m, exact analytical solution)
    print("\n3. Cuboid Magnet at [50, 0, 0] (exact analytical):")
    cuboid = CuboidMagnet(center=[0, 0, 0], dimensions=[20, 20, 10], magnetization=[0, 0, 1e6])
    pt = [50, 0, 0]
    B = cuboid.get_B(pt)
    curl_A = numerical_curl(cuboid.get_A, pt, use_meters=True)
    B_mag = np.linalg.norm(B)
    error = np.linalg.norm(np.array(curl_A) - np.array(B)) / B_mag * 100 if B_mag > 1e-15 else 0
    print(f"   B       = [{B[0]:.6e}, {B[1]:.6e}, {B[2]:.6e}]")
    print(f"   curl(A) = [{curl_A[0]:.6e}, {curl_A[1]:.6e}, {curl_A[2]:.6e}]")
    print(f"   Error: {error:.2f}%")
    if error > 5.0:
        all_passed = False

    print("-" * 60)
    assert all_passed, "TEST FAILED: Some vector potentials do not satisfy curl(A) = B"
    print("TEST PASSED: All vector potentials satisfy curl(A) ~= B")


if __name__ == "__main__":
    print("Analytical Magnet Field Test Suite")
    print("=" * 60)
    print(f"magpylib available: {HAS_MAGPYLIB}")
    print()

    tests = [
        ("Spherical Magnet", test_spherical_magnet),
        ("Cuboid Magnet", test_cuboid_magnet),
        ("Current Loop", test_current_loop),
        ("Vector Potentials", test_vector_potentials),
    ]

    results = []
    for name, func in tests:
        try:
            func()
            results.append((name, True))
        except AssertionError as e:
            print(str(e))
            results.append((name, False))

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)

    all_passed = True
    for name, passed in results:
        status = "PASSED" if passed else "FAILED"
        print(f"  {name}: {status}")
        if not passed:
            all_passed = False

    if all_passed:
        print("\nAll tests completed!")
        sys.exit(0)
    else:
        print("\nSome tests failed!")
        sys.exit(1)
