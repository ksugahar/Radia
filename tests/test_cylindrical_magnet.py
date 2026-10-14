"""
Test cylindrical magnet analytical field computation against magpylib.

This test verifies the Python implementation of cylindrical magnet field
computation using the Derby-Olbert (2010) and Caciagli (2018) formulations.

References:
    [1] Derby, N., Olbert, S., "Cylindrical Magnets and Ideal Solenoids",
        American Journal of Physics, Vol. 78(3), pp. 229-235, 2010.
    [2] Caciagli, A., et al., "Exact expression for the magnetic field of
        a finite cylinder with arbitrary uniform magnetization",
        Journal of Magnetism and Magnetic Materials, 456, 423-432, 2018.
"""

import sys

import numpy as np

from radia.cylindrical_magnet import CylindricalMagnet, RingMagnet

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


def test_axial_cylinder_on_axis():
    """Test axially magnetized cylinder B-field on axis."""
    print("=" * 60)
    print("Test: Axially magnetized cylinder on-axis B-field")
    print("=" * 60)

    # Cylinder parameters
    R = 10.0   # mm, radius
    L = 20.0   # mm, total height (half-height = 10mm)
    Mz = 1e6   # A/m, magnetization (Br ~ 1.26 T)

    print(f"Cylinder: R={R} mm, L={L} mm, Mz={Mz:.0e} A/m")

    # Create our cylindrical magnet
    cyl_ours = CylindricalMagnet(
        center=[0, 0, 0],
        radius=R,
        height=L,
        magnetization=[0, 0, Mz]
    )

    # On-axis test points (x=0, y=0)
    z_values = [0, 5, 10, 15, 20, 30, 50]  # mm

    print(f"\n{'z [mm]':>10} {'Bz_ours [T]':>15} {'Bz_magpylib [T]':>18} {'Error [%]':>12}")
    print("-" * 65)

    half_L = L / 2.0
    J = MU0 * Mz  # Polarization in Tesla

    max_error = 0.0
    for z in z_values:
        # Our implementation
        B_ours = cyl_ours.get_B([0, 0, z])
        Bz_ours = B_ours[2]

        # Expected from magpylib if available
        if HAS_MAGPYLIB:
            cyl = magpy.magnet.Cylinder(
                polarization=(0, 0, J),
                dimension=(R * 2 / 1000.0, L / 1000.0),  # diameter, height in meters
                position=(0, 0, 0)
            )
            B_magpy = cyl.getB((0, 0, z / 1000.0))  # magpylib uses meters
            Bz_expected = B_magpy[2]
        else:
            # Use analytical on-axis formula as reference
            z_m = z / 1000.0
            R_m = R / 1000.0
            half_L_m = half_L / 1000.0
            zp = z_m + half_L_m
            zm = z_m - half_L_m
            Bz_expected = J / 2.0 * (
                zp / np.sqrt(zp**2 + R_m**2) -
                zm / np.sqrt(zm**2 + R_m**2)
            )

        if abs(Bz_expected) > 1e-15:
            error = abs(Bz_ours - Bz_expected) / abs(Bz_expected) * 100
        else:
            error = 0.0
        max_error = max(max_error, error)

        print(f"{z:>10.1f} {Bz_ours:>15.6e} {Bz_expected:>18.6e} {error:>12.2f}")

    print("-" * 65)
    assert max_error < 1.0, f"TEST FAILED: Maximum error {max_error:.2f}%"
    print(f"TEST PASSED: Maximum error {max_error:.4f}%")


def test_axial_cylinder_off_axis():
    """Test axially magnetized cylinder B-field off axis."""
    print("\n" + "=" * 60)
    print("Test: Axially magnetized cylinder off-axis B-field")
    print("=" * 60)

    # Cylinder parameters
    R = 10.0   # mm, radius
    L = 20.0   # mm, total height
    Mz = 1e6   # A/m, magnetization

    J = MU0 * Mz  # Polarization in Tesla
    print(f"Cylinder: R={R:.0f} mm, L={L:.0f} mm, M={Mz:.0e} A/m (J={J:.4f} T)")

    # Create our cylinder
    cyl_ours = CylindricalMagnet(
        center=[0, 0, 0],
        radius=R,
        height=L,
        magnetization=[0, 0, Mz]
    )

    # Test points in mm
    test_points_mm = [
        (5, 0, 15),
        (10, 0, 20),
        (15, 0, 10),
        (0, 10, 15),
        (7, 7, 25),
    ]

    if HAS_MAGPYLIB:
        print(f"\n{'Point [mm]':>20} {'Bx_ours':>12} {'Bx_magpy':>12} {'Bz_ours':>12} {'Bz_magpy':>12}")
        print("-" * 75)

        cyl_magpy = magpy.magnet.Cylinder(
            polarization=(0, 0, J),
            dimension=(R * 2 / 1000.0, L / 1000.0),  # diameter, height in meters
            position=(0, 0, 0)
        )

        max_error = 0.0
        for pt in test_points_mm:
            B_ours = cyl_ours.get_B(list(pt))
            pt_m = tuple(p / 1000.0 for p in pt)
            B_magpy = cyl_magpy.getB(pt_m)

            # Compute relative error on Bz (usually largest component)
            if abs(B_magpy[2]) > 1e-15:
                error = abs(B_ours[2] - B_magpy[2]) / abs(B_magpy[2]) * 100
                max_error = max(max_error, error)

            print(f"({pt[0]:5.1f}, {pt[1]:5.1f}, {pt[2]:5.1f}) "
                  f"{B_ours[0]:>12.6e} {B_magpy[0]:>12.6e} "
                  f"{B_ours[2]:>12.6e} {B_magpy[2]:>12.6e}")

        print("-" * 75)
        assert max_error < 1.0, f"TEST FAILED: Maximum Bz error {max_error:.2f}%"
        print(f"TEST PASSED: Maximum Bz error {max_error:.4f}%")
    else:
        print(f"\n{'Point [mm]':>25} {'Bx [T]':>12} {'By [T]':>12} {'Bz [T]':>12}")
        print("-" * 70)

        for pt in test_points_mm:
            B = cyl_ours.get_B(list(pt))
            print(f"({pt[0]:6.1f}, {pt[1]:6.1f}, {pt[2]:6.1f}) "
                  f"{B[0]:>12.6e} {B[1]:>12.6e} {B[2]:>12.6e}")

        print("-" * 70)
        print("TEST INFO: Values shown (no magpylib for comparison)")


def test_diametric_cylinder():
    """Test diametrically magnetized cylinder B-field."""
    print("\n" + "=" * 60)
    print("Test: Diametrically magnetized cylinder B-field")
    print("=" * 60)

    # Cylinder parameters
    R = 10.0   # mm, radius
    L = 20.0   # mm, total height
    Mx = 1e6   # A/m, magnetization in x-direction

    J = MU0 * Mx  # Polarization in Tesla
    print(f"Cylinder: R={R:.0f} mm, L={L:.0f} mm, Mx={Mx:.0e} A/m (Jx={J:.4f} T)")

    # Create our cylinder with diametric magnetization
    cyl_ours = CylindricalMagnet(
        center=[0, 0, 0],
        radius=R,
        height=L,
        magnetization=[Mx, 0, 0]
    )

    # Test points in mm
    test_points_mm = [
        (20, 0, 0),      # Along magnetization direction
        (0, 20, 0),      # Perpendicular to magnetization
        (0, 0, 20),      # On axis
        (15, 0, 15),
        (10, 10, 20),
    ]

    if HAS_MAGPYLIB:
        print(f"\n{'Point [mm]':>20} {'Bx_ours':>12} {'Bx_magpy':>12} {'Bz_ours':>12} {'Bz_magpy':>12}")
        print("-" * 75)

        cyl_magpy = magpy.magnet.Cylinder(
            polarization=(J, 0, 0),  # Magnetized in x-direction
            dimension=(R * 2 / 1000.0, L / 1000.0),
            position=(0, 0, 0)
        )

        max_error = 0.0
        for pt in test_points_mm:
            B_ours = cyl_ours.get_B(list(pt))
            pt_m = tuple(p / 1000.0 for p in pt)
            B_magpy = cyl_magpy.getB(pt_m)

            # Compute relative error on Bx (expected to be largest)
            if abs(B_magpy[0]) > 1e-15:
                error = abs(B_ours[0] - B_magpy[0]) / abs(B_magpy[0]) * 100
                max_error = max(max_error, error)

            print(f"({pt[0]:5.1f}, {pt[1]:5.1f}, {pt[2]:5.1f}) "
                  f"{B_ours[0]:>12.6e} {B_magpy[0]:>12.6e} "
                  f"{B_ours[2]:>12.6e} {B_magpy[2]:>12.6e}")

        print("-" * 75)
        assert max_error < 5.0, f"TEST FAILED: Maximum Bx error {max_error:.2f}%"  # Allow 5% error for diametric case (Taylor expansion)
        print(f"TEST PASSED: Maximum Bx error {max_error:.4f}%")
    else:
        print(f"\n{'Point [mm]':>25} {'Bx [T]':>12} {'By [T]':>12} {'Bz [T]':>12}")
        print("-" * 70)

        for pt in test_points_mm:
            B = cyl_ours.get_B(list(pt))
            print(f"({pt[0]:6.1f}, {pt[1]:6.1f}, {pt[2]:6.1f}) "
                  f"{B[0]:>12.6e} {B[1]:>12.6e} {B[2]:>12.6e}")

        print("-" * 70)
        print("TEST INFO: Values shown (no magpylib for comparison)")


def test_ring_magnet():
    """Test ring magnet (hollow cylinder) B-field."""
    print("\n" + "=" * 60)
    print("Test: Ring magnet (hollow cylinder) B-field")
    print("=" * 60)

    # Ring magnet parameters
    R_inner = 5.0   # mm, inner radius
    R_outer = 15.0  # mm, outer radius
    L = 20.0        # mm, total height
    Mz = 1e6        # A/m, magnetization

    J = MU0 * Mz  # Polarization in Tesla
    print(f"Ring: R_inner={R_inner:.0f} mm, R_outer={R_outer:.0f} mm, "
          f"L={L:.0f} mm, Mz={Mz:.0e} A/m")

    # Create our ring magnet
    ring_ours = RingMagnet(
        center=[0, 0, 0],
        inner_radius=R_inner,
        outer_radius=R_outer,
        height=L,
        magnetization=[0, 0, Mz]
    )

    # Test points on axis
    z_values_mm = [0, 10, 20, 30, 50]

    if HAS_MAGPYLIB:
        print(f"\n{'z [mm]':>10} {'Bz_ours [T]':>15} {'Bz_magpy [T]':>15} {'Error [%]':>12}")
        print("-" * 55)

        ring_magpy = magpy.magnet.CylinderSegment(
            polarization=(0, 0, J),
            dimension=(R_inner / 1000.0, R_outer / 1000.0, L / 1000.0, 0, 360),
            position=(0, 0, 0)
        )

        max_error = 0.0
        for z in z_values_mm:
            B_ours = ring_ours.get_B([0, 0, z])
            B_magpy = ring_magpy.getB((0, 0, z / 1000.0))

            if abs(B_magpy[2]) > 1e-15:
                error = abs(B_ours[2] - B_magpy[2]) / abs(B_magpy[2]) * 100
                max_error = max(max_error, error)
            else:
                error = 0.0

            print(f"{z:>10.1f} {B_ours[2]:>15.6e} {B_magpy[2]:>15.6e} {error:>12.2f}")

        print("-" * 55)
        assert max_error < 1.0, f"TEST FAILED: Maximum error {max_error:.2f}%"
        print(f"TEST PASSED: Maximum error {max_error:.4f}%")
    else:
        print(f"\n{'z [mm]':>10} {'Bz [T]':>15}")
        print("-" * 30)

        for z in z_values_mm:
            B = ring_ours.get_B([0, 0, z])
            print(f"{z:>10.1f} {B[2]:>15.6e}")

        print("-" * 30)
        print("TEST INFO: Values shown (no magpylib for comparison)")


def test_vector_potential():
    """Test vector potential A for axially magnetized cylinder."""
    print("\n" + "=" * 60)
    print("Test: Vector Potential (A) for Axially Magnetized Cylinder")
    print("=" * 60)

    # Cylinder parameters
    R = 10.0   # mm, radius
    L = 20.0   # mm, total height
    Mz = 1e6   # A/m, magnetization

    print(f"Cylinder: R={R:.0f} mm, L={L:.0f} mm, Mz={Mz:.0e} A/m")

    # Create our cylinder
    cyl = CylindricalMagnet(
        center=[0, 0, 0],
        radius=R,
        height=L,
        magnetization=[0, 0, Mz]
    )

    # Test points
    test_points = [
        [0, 0, 20],      # On axis (should be zero)
        [5, 0, 20],      # Off-axis
        [10, 0, 15],     # At cylinder surface level
        [15, 0, 0],      # In the midplane
        [7, 7, 25],      # Off-axis, diagonal
    ]

    print(f"\n{'Point [mm]':>25} {'Ax [T*mm]':>15} {'Ay [T*mm]':>15} {'Az [T*mm]':>15}")
    print("-" * 75)

    all_reasonable = True
    for pt in test_points:
        A = cyl.get_A(pt)

        # Check reasonableness
        # For axial magnetization, Az should be ~0
        # A should be purely azimuthal (Ax = -A_phi*sin(phi), Ay = A_phi*cos(phi))
        if abs(A[2]) > 1e-15:  # Az should be zero for axial magnetization
            all_reasonable = False

        # On-axis (rho=0), A should be zero
        rho = np.sqrt(pt[0]**2 + pt[1]**2)
        if rho < 1e-10 and (abs(A[0]) > 1e-15 or abs(A[1]) > 1e-15):
            all_reasonable = False

        print(f"({pt[0]:6.1f}, {pt[1]:6.1f}, {pt[2]:6.1f}) "
              f"{A[0]:>15.6e} {A[1]:>15.6e} {A[2]:>15.6e}")

    print("-" * 75)

    # Verify curl(A) ~ B at a test point (numerical curl)
    print("\nVerifying curl(A) ~= B at test point [10, 0, 20]:")
    pt = [10, 0, 20]
    h = 0.1  # Small step for numerical differentiation

    # Get A at neighboring points
    A_px = cyl.get_A([pt[0] + h, pt[1], pt[2]])
    A_mx = cyl.get_A([pt[0] - h, pt[1], pt[2]])
    A_py = cyl.get_A([pt[0], pt[1] + h, pt[2]])
    A_my = cyl.get_A([pt[0], pt[1] - h, pt[2]])
    A_pz = cyl.get_A([pt[0], pt[1], pt[2] + h])
    A_mz = cyl.get_A([pt[0], pt[1], pt[2] - h])

    # curl(A) = (dAz/dy - dAy/dz, dAx/dz - dAz/dx, dAy/dx - dAx/dy)
    dAz_dy = (A_py[2] - A_my[2]) / (2 * h)
    dAy_dz = (A_pz[1] - A_mz[1]) / (2 * h)
    dAx_dz = (A_pz[0] - A_mz[0]) / (2 * h)
    dAz_dx = (A_px[2] - A_mx[2]) / (2 * h)
    dAy_dx = (A_px[1] - A_mx[1]) / (2 * h)
    dAx_dy = (A_py[0] - A_my[0]) / (2 * h)

    curl_A = [
        dAz_dy - dAy_dz,
        dAx_dz - dAz_dx,
        dAy_dx - dAx_dy
    ]

    # Get B at same point
    B = cyl.get_B(pt)

    print(f"  curl(A) = [{curl_A[0]:.6e}, {curl_A[1]:.6e}, {curl_A[2]:.6e}] T")
    print(f"  B       = [{B[0]:.6e}, {B[1]:.6e}, {B[2]:.6e}] T")

    # Check relative error
    B_mag = np.sqrt(B[0]**2 + B[1]**2 + B[2]**2)
    if B_mag > 1e-15:
        error_x = abs(curl_A[0] - B[0]) / B_mag * 100
        error_y = abs(curl_A[1] - B[1]) / B_mag * 100
        error_z = abs(curl_A[2] - B[2]) / B_mag * 100
        max_error = max(error_x, error_y, error_z)
        print(f"  Max relative error: {max_error:.2f}%")

        # Allow 5% error due to numerical differentiation
        assert max_error < 5.0 and all_reasonable, "TEST FAILED: curl(A) != B or unreasonable A values"
        print("TEST PASSED: curl(A) ~= B within 5%")
    else:
        print("TEST INFO: B magnitude too small for meaningful comparison")
        assert all_reasonable, "TEST FAILED: unreasonable A values"


def test_vector_potential_diametric_raises():
    """No-Fallback: get_A() must RAISE on diametric (transverse) magnetization
    rather than SILENTLY returning an axial-only (wrong) / zero A -- the
    diametric contribution is not implemented.  (get_B handles all directions;
    see test_diametric_cylinder.)"""
    import pytest
    # pure diametric (Mx) -> raise instead of silently returning zero A
    cyl = CylindricalMagnet(center=[0, 0, 0], radius=10.0, height=20.0,
                            magnetization=[1e6, 0, 0])
    with pytest.raises(NotImplementedError, match="diametric"):
        cyl.get_A([5.0, 0.0, 5.0])
    # mixed axial + diametric -> still raise (the Mx part would be dropped)
    cyl_mix = CylindricalMagnet(center=[0, 0, 0], radius=10.0, height=20.0,
                                magnetization=[1e6, 0, 1e6])
    with pytest.raises(NotImplementedError, match="diametric"):
        cyl_mix.get_A([5.0, 0.0, 5.0])
    # pure axial still works (regression guard)
    cyl_ax = CylindricalMagnet(center=[0, 0, 0], radius=10.0, height=20.0,
                               magnetization=[0, 0, 1e6])
    A = cyl_ax.get_A([5.0, 0.0, 5.0])
    assert abs(A[2]) < 1e-12          # Az = 0 for axial magnetization


if __name__ == "__main__":
    print("Cylindrical Magnet Field Test Suite")
    print("=" * 60)
    print(f"magpylib available: {HAS_MAGPYLIB}")
    print()

    tests = [
        ("On-axis B-field", test_axial_cylinder_on_axis),
        ("Off-axis B-field", test_axial_cylinder_off_axis),
        ("Diametric B-field", test_diametric_cylinder),
        ("Ring magnet", test_ring_magnet),
        ("Vector potential A", test_vector_potential),
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
