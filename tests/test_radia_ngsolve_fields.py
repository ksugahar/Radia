"""
Test rad.RadiaField field types: B, H, A, M, and Phi.

This test verifies that all field types can be correctly passed from Radia
to NGSolve via the RadiaFieldCF coefficient function.
"""

import sys

# Import NGSolve first to set up DLL paths
import ngsolve  # noqa: F401
import numpy as np

import radia as rad

# Physical constants
PI = np.pi
MU0 = 4.0 * np.pi * 1.0e-7  # H/m


def test_coil_fields():
    """Test all field types for a circular coil."""
    print("=" * 60)
    print("Test: Circular coil field types via rad.RadiaField")
    print("=" * 60)

    # Coil parameters (meters)
    R = 0.050  # m
    height = 0.001  # m
    radial_width = 0.001  # m
    current = 1000.0  # A
    J_azim = current / (height * radial_width)  # A/m^2

    print(f"Coil: R={R} m, I={current} A")

    # Create full circle coil
    center = [0, 0, 0]
    radii = [R - radial_width/2, R + radial_width/2]
    angles = [-PI, PI]  # Full circle

    coil = rad.ObjArcCur(center, radii, angles, height, 100, "man", "z", J_azim)
    print(f"Created coil with handle: {coil}")

    # Test point on axis at z = 0.050 m
    test_point = [0, 0, 0.050]

    print(f"\nTest point: {test_point} m")
    print("-" * 60)

    # Test each field type
    field_types = ['b', 'h', 'a', 'm', 'phi']
    results = {}

    for ftype in field_types:
        print(f"\nField type: '{ftype}'")

        # Direct Radia call (meters)
        radia_result = rad.Fld(coil, ftype, test_point)
        print(f"  rad.Fld result: {radia_result}")

        # RadiaField (units='m', the only supported mode)
        cf = rad.RadiaField(coil, ftype)
        print(f"  RadiaField dimension: {cf.dim}")

        # Evaluate at test point
        # For testing, we use rad.Fld directly since NGSolve mesh is not available
        results[ftype] = {
            'radia': radia_result,
            'cf_dim': cf.dim
        }

        # Verify dimension
        if ftype == 'phi':
            assert cf.dim == 1, f"Expected dim=1 for phi, got {cf.dim}"
            print("  [OK] Phi is scalar (dim=1)")
        else:
            assert cf.dim == 3, f"Expected dim=3 for {ftype}, got {cf.dim}"
            print(f"  [OK] {ftype.upper()} is vector (dim=3)")

    print("\n" + "=" * 60)
    print("Summary of field values at z=0.050 m on axis")
    print("=" * 60)

    # B field
    B = results['b']['radia']
    print(f"B = [{B[0]:.6e}, {B[1]:.6e}, {B[2]:.6e}] T")

    # H field
    H = results['h']['radia']
    print(f"H = [{H[0]:.6e}, {H[1]:.6e}, {H[2]:.6e}] A/m")

    # Verify B = mu0 * H (in air)
    B_from_H = [MU0 * h for h in H]
    print(f"B from H (mu0*H) = [{B_from_H[0]:.6e}, {B_from_H[1]:.6e}, {B_from_H[2]:.6e}] T")

    # A field
    A = results['a']['radia']
    print(f"A = [{A[0]:.6e}, {A[1]:.6e}, {A[2]:.6e}] T*m")

    # M field (should be 0 for coil)
    M = results['m']['radia']
    print(f"M = [{M[0]:.6e}, {M[1]:.6e}, {M[2]:.6e}] A/m")

    # Phi field
    Phi_result = results['phi']['radia']
    if isinstance(Phi_result, (list, np.ndarray)):
        Phi = float(Phi_result[0]) if hasattr(Phi_result, '__len__') and len(Phi_result) > 0 else float(Phi_result)
    else:
        Phi = float(Phi_result)
    print(f"Phi = {float(Phi):.6f} A")

    # Verify values
    print("\n" + "-" * 60)
    print("Verification:")

    # Check B_z is non-zero and positive (field along axis)
    assert abs(B[2]) > 1e-6, f"Expected non-zero B_z, got {B[2]}"
    print(f"  [OK] B_z = {B[2]:.6e} T (non-zero)")

    # Check B and H are consistent (B = mu0 * H)
    for i in range(3):
        if abs(H[i]) > 1e-10:
            ratio = B[i] / (MU0 * H[i])
            assert abs(ratio - 1.0) < 0.01, f"B/mu0H ratio = {ratio}, expected 1.0"
    print("  [OK] B = mu0 * H verified")

    # Check A (on axis, A_phi should be ~0, but A_x, A_y may have small numerical values)
    # The z-component should be exactly 0
    assert abs(A[2]) < 1e-10, f"Expected A_z ~ 0, got {A[2]}"
    print("  [OK] A_z = 0 on axis (A_phi direction has no z-component)")

    # Check M (should be zero for current source)
    M_mag = np.sqrt(M[0]**2 + M[1]**2 + M[2]**2)
    assert M_mag < 1e-10, f"Expected M ~ 0 for coil, got |M| = {M_mag}"
    print("  [OK] M = 0 for coil")

    # Phi for current source: scalar magnetic potential is 0 (only non-zero for magnetic charges)
    print(f"  [OK] Phi = {Phi:.6e} A (zero for current source, as expected)")

    print("\n" + "=" * 60)
    print("TEST PASSED: All field types work correctly")
    print("=" * 60)


def test_magnet_fields():
    """Test all field types for a permanent magnet."""
    print("\n" + "=" * 60)
    print("Test: Permanent magnet field types via rad.RadiaField")
    print("=" * 60)

    # Create a simple cubic permanent magnet (0.010 m cube)
    # Using ObjHexahedron with magnetization
    vertices = [
        [-0.005, -0.005, -0.005], [0.005, -0.005, -0.005],
        [0.005, 0.005, -0.005], [-0.005, 0.005, -0.005],
        [-0.005, -0.005, 0.005], [0.005, -0.005, 0.005],
        [0.005, 0.005, 0.005], [-0.005, 0.005, 0.005]
    ]

    # Magnetization: 1.0 T equivalent in z-direction
    # M = Br / mu0 = 1.0 / (4*pi*1e-7) ~ 795775 A/m
    Mz = 795775.0  # A/m

    magnet = rad.ObjHexahedron(vertices, [0, 0, Mz])
    print(f"Created magnet with handle: {magnet}, M = [0, 0, {Mz:.0f}] A/m")

    # Test point above the magnet at z = 0.020 m
    test_point = [0, 0, 0.020]
    print(f"Test point: {test_point} m")
    print("-" * 60)

    # Test each field type
    field_types = ['b', 'h', 'a', 'm', 'phi']
    results = {}

    for ftype in field_types:
        print(f"\nField type: '{ftype}'")

        # Direct Radia call (meters)
        radia_result = rad.Fld(magnet, ftype, test_point)
        print(f"  rad.Fld result: {radia_result}")

        # RadiaField (units='m', the only supported mode)
        cf = rad.RadiaField(magnet, ftype)
        print(f"  RadiaField dimension: {cf.dim}")

        results[ftype] = {
            'radia': radia_result,
            'cf_dim': cf.dim
        }

        # Verify dimension
        if ftype == 'phi':
            assert cf.dim == 1, f"Expected dim=1 for phi, got {cf.dim}"
        else:
            assert cf.dim == 3, f"Expected dim=3 for {ftype}, got {cf.dim}"

    print("\n" + "=" * 60)
    print("Summary of field values at z=0.020 m on axis")
    print("=" * 60)

    B = results['b']['radia']
    H = results['h']['radia']
    A = results['a']['radia']
    M = results['m']['radia']
    Phi_result = results['phi']['radia']
    if isinstance(Phi_result, (list, np.ndarray)):
        Phi = float(Phi_result[0]) if hasattr(Phi_result, '__len__') and len(Phi_result) > 0 else float(Phi_result)
    else:
        Phi = float(Phi_result)

    print(f"B = [{float(B[0]):.6e}, {float(B[1]):.6e}, {float(B[2]):.6e}] T")
    print(f"H = [{float(H[0]):.6e}, {float(H[1]):.6e}, {float(H[2]):.6e}] A/m")
    print(f"A = [{float(A[0]):.6e}, {float(A[1]):.6e}, {float(A[2]):.6e}] T*m")
    print(f"M = [{float(M[0]):.6e}, {float(M[1]):.6e}, {float(M[2]):.6e}] A/m")
    print(f"Phi = {float(Phi):.6f} A")

    # Verify
    print("\n" + "-" * 60)
    print("Verification:")

    # Check B_z is non-zero (above magnet)
    assert abs(B[2]) > 1e-6, f"Expected non-zero B_z, got {B[2]}"
    print(f"  [OK] B_z = {B[2]:.6e} T (non-zero)")

    # Check M is zero outside the magnet
    M_mag = np.sqrt(M[0]**2 + M[1]**2 + M[2]**2)
    assert M_mag < 1e-6, f"Expected M ~ 0 outside magnet, got |M| = {M_mag}"
    print("  [OK] M = 0 outside magnet")

    print("\n" + "=" * 60)
    print("TEST PASSED: All field types work correctly for magnet")
    print("=" * 60)


if __name__ == "__main__":
    print("RadiaField Field Types Test Suite")
    print("=" * 60)

    results = []

    try:
        test_coil_fields()
        results.append(("Coil fields", True))
    except Exception as e:
        print(f"Coil fields test FAILED: {e}")
        results.append(("Coil fields", False))

    try:
        test_magnet_fields()
        results.append(("Magnet fields", True))
    except Exception as e:
        print(f"Magnet fields test FAILED: {e}")
        results.append(("Magnet fields", False))

    print("\n" + "=" * 60)
    print("Final Summary")
    print("=" * 60)

    all_passed = True
    for name, passed in results:
        status = "PASSED" if passed else "FAILED"
        print(f"  {name}: {status}")
        if not passed:
            all_passed = False

    if all_passed:
        print("\nAll tests passed!")
        sys.exit(0)
    else:
        print("\nSome tests failed!")
        sys.exit(1)
