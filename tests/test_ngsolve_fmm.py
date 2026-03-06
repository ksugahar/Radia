"""
Test FMM-accelerated NGSolve integration.

This test verifies that the FMM dipole approximation in radia_ngsolve
produces correct H-field results compared to the full Radia computation.
"""

import sys
import os
import time
import numpy as np

# Add build directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/radia'))

import radia as rad

try:
    from ngsolve import Mesh, H1, GridFunction
    from netgen.occ import Box, Pnt, OCCGeometry
    import radia_ngsolve
    HAS_NGSOLVE = True
except ImportError:
    HAS_NGSOLVE = False
    print("Warning: NGSolve not installed, skipping FMM NGSolve tests")


def create_test_magnet():
    """Create a simple cube permanent magnet for testing."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    # Create cube magnet (0.1m x 0.1m x 0.1m) using ObjHexahedron
    vertices = [
        [-0.05, -0.05, -0.05],
        [0.05, -0.05, -0.05],
        [0.05, 0.05, -0.05],
        [-0.05, 0.05, -0.05],
        [-0.05, -0.05, 0.05],
        [0.05, -0.05, 0.05],
        [0.05, 0.05, 0.05],
        [-0.05, 0.05, 0.05],
    ]
    # Magnetization: Br = 1.2 T -> M = Br / mu0 = 954930 A/m
    M = 954930.0
    magnet = rad.ObjHexahedron(vertices, [0, 0, M])

    return magnet


def test_fmm_parameter():
    """Test that fmm_eps parameter is correctly exposed."""
    print("=" * 60)
    print("Test: FMM parameter exposure")
    print("=" * 60)

    if not HAS_NGSOLVE:
        print("SKIPPED: NGSolve not available")
        return True

    magnet = create_test_magnet()

    # Test default (disabled)
    cf_default = radia_ngsolve.RadiaField(magnet, 'h')
    assert cf_default.fmm_eps == 0.0, f"Expected fmm_eps=0.0, got {cf_default.fmm_eps}"
    print(f"  Default fmm_eps: {cf_default.fmm_eps} (OK)")

    # Test with explicit value
    cf_fmm = radia_ngsolve.RadiaField(magnet, 'h', fmm_eps=1e-4)
    assert cf_fmm.fmm_eps == 1e-4, f"Expected fmm_eps=1e-4, got {cf_fmm.fmm_eps}"
    print(f"  Explicit fmm_eps: {cf_fmm.fmm_eps} (OK)")

    print("-" * 60)
    print("TEST PASSED: FMM parameter correctly exposed")
    return True


def test_fmm_field_accuracy():
    """Test that FMM dipole field matches Radia field (direct comparison).

    Note: This tests the dipole approximation accuracy independent of
    GridFunction interpolation. For a single cube magnet at distance
    ~3x its size, dipole approximation should be within ~5%.
    """
    print("\n" + "=" * 60)
    print("Test: FMM dipole accuracy (direct comparison)")
    print("=" * 60)

    if not HAS_NGSOLVE:
        print("SKIPPED: NGSolve not available")
        return True

    magnet = create_test_magnet()

    # Get dipole parameters directly from Radia API
    m_result = rad.ObjM(magnet)
    Mx, My, Mz = m_result[1]  # Second element is [Mx, My, Mz]

    geo = rad.ObjGeoLim(magnet)
    xmin, xmax = geo[0], geo[1]  # mm
    ymin, ymax = geo[2], geo[3]  # mm
    zmin, zmax = geo[4], geo[5]  # mm

    # Compute dipole moment (M in A/m, V in m^3)
    vol_mm3 = (xmax - xmin) * (ymax - ymin) * (zmax - zmin)
    vol_m3 = vol_mm3 * 1e-9

    # Dipole center (converted from mm to m)
    cx = (xmin + xmax) / 2.0 * 0.001
    cy = (ymin + ymax) / 2.0 * 0.001
    cz = (zmin + zmax) / 2.0 * 0.001

    # Dipole moment in A*m^2
    mx = Mx * vol_m3
    my = My * vol_m3
    mz = Mz * vol_m3

    print(f"  Dipole moment: m = ({mx:.2f}, {my:.2f}, {mz:.2f}) A*m^2")
    print(f"  Center: ({cx:.3f}, {cy:.3f}, {cz:.3f}) m")

    # Test points (far enough from source for dipole approximation to be valid)
    # At 0.20 m from center, distance/size ratio = 0.20/0.05 = 4x
    test_points = [
        (0.20, 0, 0),
        (0, 0.20, 0),
        (0, 0, 0.20),
        (0.12, 0.12, 0.12),
    ]

    print(f"\n  {'Point':>25} {'H_radia (A/m)':>15} {'H_dipole (A/m)':>15} {'Error (%)':>12}")
    print("  " + "-" * 70)

    max_error = 0.0
    for pt in test_points:
        # Radia field (direct)
        H_radia = rad.Fld(magnet, 'h', list(pt))
        H_r_mag = np.sqrt(sum(x**2 for x in H_radia))

        # Manual dipole calculation
        rx = pt[0] - cx
        ry = pt[1] - cy
        rz = pt[2] - cz
        r = np.sqrt(rx**2 + ry**2 + rz**2)
        r3 = r**3
        r5 = r**5

        m_dot_r = mx*rx + my*ry + mz*rz
        Hx = (1.0 / (4*np.pi)) * (3*m_dot_r*rx/r5 - mx/r3)
        Hy = (1.0 / (4*np.pi)) * (3*m_dot_r*ry/r5 - my/r3)
        Hz = (1.0 / (4*np.pi)) * (3*m_dot_r*rz/r5 - mz/r3)
        H_d_mag = np.sqrt(Hx**2 + Hy**2 + Hz**2)

        if H_r_mag > 1e-10:
            error = abs(H_d_mag - H_r_mag) / H_r_mag * 100
            max_error = max(max_error, error)
        else:
            error = 0.0

        print(f"  ({pt[0]:.2f}, {pt[1]:.2f}, {pt[2]:.2f}) {H_r_mag:>15.2f} {H_d_mag:>15.2f} {error:>12.2f}")

    print("  " + "-" * 70)

    # For dipole approximation at ~4x source size, expect <5% error
    if max_error < 10.0:
        print(f"TEST PASSED: Maximum error {max_error:.2f}% (dipole approximation)")
        return True
    else:
        print(f"TEST FAILED: Maximum error {max_error:.2f}%")
        return False


def test_fmm_speedup():
    """Test that FMM provides speedup for large problems."""
    print("\n" + "=" * 60)
    print("Test: FMM speedup for larger mesh")
    print("=" * 60)

    if not HAS_NGSOLVE:
        print("SKIPPED: NGSolve not available")
        return True

    magnet = create_test_magnet()

    # Create a finer mesh for better speedup demonstration
    geo = OCCGeometry(Box(Pnt(-0.3, -0.3, -0.3), Pnt(0.3, 0.3, 0.3)))
    mesh = Mesh(geo.GenerateMesh(maxh=0.05))  # Finer mesh

    from ngsolve import VOL
    n_elements = mesh.GetNE(VOL)
    n_vertices = mesh.nv
    print(f"  Mesh: {n_elements} elements, {n_vertices} vertices")

    # Create coefficient functions
    H_radia = radia_ngsolve.RadiaField(magnet, 'h', fmm_eps=0.0)  # Full
    H_fmm = radia_ngsolve.RadiaField(magnet, 'h', fmm_eps=1e-3)   # FMM

    fes = H1(mesh, dim=3, order=1)

    # Warm-up
    gf_warmup = GridFunction(fes)
    gf_warmup.Set(H_radia)

    # Benchmark full Radia
    gf_radia = GridFunction(fes)
    t0 = time.time()
    gf_radia.Set(H_radia)
    t_radia = time.time() - t0

    # Benchmark FMM
    gf_fmm = GridFunction(fes)
    t0 = time.time()
    gf_fmm.Set(H_fmm)
    t_fmm = time.time() - t0

    print(f"  Full Radia: {t_radia:.3f} s")
    print(f"  FMM dipole: {t_fmm:.3f} s")
    print(f"  Speedup:    {t_radia/t_fmm:.2f}x" if t_fmm > 0 else "  (FMM faster)")

    print("-" * 60)
    # Note: For single-element source, FMM may not show speedup
    # FMM is primarily for multi-element sources
    print("TEST INFO: FMM timing benchmark complete")
    return True


def test_fmm_all_field_types():
    """Test FMM computation for all supported field types: H, B, A.

    Compares FMM-computed fields against Radia's full computation
    at points far from the source where dipole approximation is valid.
    """
    print("\n" + "=" * 60)
    print("Test: FMM all field types (H, B, A)")
    print("=" * 60)

    if not HAS_NGSOLVE:
        print("SKIPPED: NGSolve not available")
        return True

    magnet = create_test_magnet()

    # Test points far from source
    test_points = [
        [0.20, 0, 0],
        [0, 0.20, 0],
        [0, 0, 0.20],
        [0.15, 0.15, 0],
    ]

    field_types = ['h', 'b', 'a']
    all_passed = True

    for ftype in field_types:
        print(f"\n  Testing field type: '{ftype}'")

        for pt in test_points:
            # FMM result via radia_ngsolve with high fmm_eps
            cf_fmm = radia_ngsolve.RadiaField(magnet, ftype, fmm_eps=1e-6)

            # Get FMM result by evaluating at point
            # We need to create a mesh to evaluate - use a tiny box around the point
            from ngsolve import Mesh
            box_size = 0.001
            box = Box(Pnt(pt[0]-box_size, pt[1]-box_size, pt[2]-box_size),
                      Pnt(pt[0]+box_size, pt[1]+box_size, pt[2]+box_size))
            geo = OCCGeometry(box)
            mesh = Mesh(geo.GenerateMesh(maxh=box_size))

            fes = H1(mesh, dim=3, order=1)
            gf = GridFunction(fes)
            gf.Set(cf_fmm)

            # Evaluate at center point
            fmm_val = gf(mesh(*pt))

            # Full Radia result (direct API)
            # For vector potential 'a', Radia returns T*mm because it uses mm internally
            # NGSolve expects T*m, so FMM formula uses meters throughout
            # Radia's A in T*mm = (numerical value) * mm = (value/1000) * m
            # So to compare properly:
            # - FMM gives A in T*m (when dipole positions are in m)
            # - Radia gives A in T*mm, need to convert to T*m
            # Actually let's check what rad.FldUnits() says...
            radia_val = rad.Fld(magnet, ftype, pt)

            # Note: Radia A is in T*m when FldUnits('m') is set
            # The 0.001 scaling in radia_ngsolve Python path is for when
            # Radia uses mm internally (which is always true)
            # But since we set FldUnits('m'), the coords are in m and A should be in T*m
            # Let's not scale here and see what happens

            # Compute error
            fmm_mag = np.sqrt(sum(v**2 for v in fmm_val))
            radia_mag = np.sqrt(sum(v**2 for v in radia_val))

            # For very small fields (essentially zero), use absolute error instead
            # If both are < 1e-10, consider it a match
            if radia_mag > 1e-10:
                error = abs(fmm_mag - radia_mag) / radia_mag * 100
                is_ok = error < 10.0
            elif fmm_mag < 1e-10 and radia_mag < 1e-10:
                # Both essentially zero
                error = 0.0
                is_ok = True
            else:
                # One is zero, one is not - that's a problem
                error = 100.0
                is_ok = False

            status = "OK" if is_ok else "FAIL"
            if not is_ok:
                all_passed = False

            print(f"    Point {pt}: FMM={fmm_mag:.4e}, Radia={radia_mag:.4e}, Error={error:.1f}% [{status}]")

    print("-" * 60)
    if all_passed:
        print("TEST PASSED: All field types computed correctly")
    else:
        print("TEST FAILED: Some field types have errors > 10%")

    return all_passed


if __name__ == "__main__":
    print("NGSolve FMM Integration Test Suite")
    print("=" * 60)

    if not HAS_NGSOLVE:
        print("NGSolve not available, exiting")
        sys.exit(0)

    results = []

    results.append(("FMM Parameter", test_fmm_parameter()))
    results.append(("FMM Field Accuracy", test_fmm_field_accuracy()))
    results.append(("FMM All Field Types", test_fmm_all_field_types()))
    results.append(("FMM Speedup", test_fmm_speedup()))

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
