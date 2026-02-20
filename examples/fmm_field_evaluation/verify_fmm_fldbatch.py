"""
Verification script for FMM-accelerated FldBatch().

Compares rad.FldBatch(obj, points, method=0)  [direct, exact]
     vs  rad.FldBatch(obj, points, method=1)  [FMM dipole approx]

Tests:
  1. Far-field accuracy: dipole approx should match direct within 5%
  2. Near-field: expect larger difference (dipole approx breaks down)
  3. IMA fallback: method=1 with IMA -> silently uses direct (identical results)
  4. UtiDelAll safety: no crash after clearing all elements
  5. Timing comparison: method=1 should be faster for large N_points
"""

import sys
import os
import time
import numpy as np


sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
import radia as rad


def create_test_model(n_div=5, mu_r=1000, B_ext_z=0.1):
    """Create a multi-element soft iron cube in uniform Bz field."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    half = 0.015  # 30mm cube
    dx = 2 * half / n_div

    elements = []
    for ix in range(n_div):
        for iy in range(n_div):
            for iz in range(n_div):
                x0 = -half + ix * dx
                y0 = -half + iy * dx
                z0 = -half + iz * dx
                verts = [
                    [x0,    y0,    z0],
                    [x0+dx, y0,    z0],
                    [x0+dx, y0+dx, z0],
                    [x0,    y0+dx, z0],
                    [x0,    y0,    z0+dx],
                    [x0+dx, y0,    z0+dx],
                    [x0+dx, y0+dx, z0+dx],
                    [x0,    y0+dx, z0+dx],
                ]
                elem = rad.ObjHexahedron(verts, [0, 0, 0])
                elements.append(elem)

    container = rad.ObjCnt(elements)
    mat = rad.MatLin(mu_r - 1)  # chi = mu_r - 1
    rad.MatApl(container, mat)

    bkg = rad.ObjBckg(lambda p: [0, 0, B_ext_z])
    assembly = rad.ObjCnt([container, bkg])

    # Solve
    t0 = time.time()
    result = rad.Solve(assembly, 0.001, 100, 0)  # LU solver
    t_solve = time.time() - t0
    print(f"Solve: {n_div}^3 = {n_div**3} elements, {t_solve:.2f}s")

    return container


def generate_observation_grid(center, distances, n_per_shell=50):
    """Generate observation points at various distances from center."""
    rng = np.random.default_rng(42)
    all_points = []

    for d in distances:
        # Random points on sphere of radius d
        phi = rng.uniform(0, 2*np.pi, n_per_shell)
        cos_theta = rng.uniform(-1, 1, n_per_shell)
        sin_theta = np.sqrt(1 - cos_theta**2)

        pts = np.column_stack([
            center[0] + d * sin_theta * np.cos(phi),
            center[1] + d * sin_theta * np.sin(phi),
            center[2] + d * cos_theta
        ])
        all_points.append(pts)

    return np.vstack(all_points)


def vector_relative_error(B_ref, B_test):
    """Compute vector relative error: ||B_ref - B_test|| / ||B_ref||."""
    diff = np.linalg.norm(B_ref - B_test, axis=1)
    ref_mag = np.linalg.norm(B_ref, axis=1)
    # Avoid division by zero
    mask = ref_mag > 1e-15
    errors = np.zeros_like(ref_mag)
    errors[mask] = diff[mask] / ref_mag[mask]
    return errors


def test_far_field_accuracy():
    """Test 1: Far-field accuracy of FMM vs direct."""
    print("\n=== Test 1: Far-field accuracy ===")
    container = create_test_model(n_div=5, mu_r=1000)

    half = 0.015
    distances = [3*half, 5*half, 10*half, 20*half]
    pts = generate_observation_grid([0, 0, 0], distances, n_per_shell=100)

    # Direct computation
    t0 = time.time()
    result_direct = rad.FldBatch(container, pts, 0)
    t_direct = time.time() - t0

    # FMM computation
    t0 = time.time()
    result_fmm = rad.FldBatch(container, pts, 1)
    t_fmm = time.time() - t0

    B_direct = np.array(result_direct['B'])
    B_fmm = np.array(result_fmm['B'])

    errors = vector_relative_error(B_direct, B_fmm)

    n_per = 100
    for i, d in enumerate(distances):
        errs_shell = errors[i*n_per:(i+1)*n_per]
        mean_err = np.mean(errs_shell) * 100
        max_err = np.max(errs_shell) * 100
        status = "PASS" if max_err < 10.0 else "FAIL"
        print(f"  Distance {d/half:.0f}R: mean={mean_err:.2f}%, max={max_err:.2f}% [{status}]")

    print(f"  Timing: direct={t_direct:.4f}s, FMM={t_fmm:.4f}s, speedup={t_direct/max(t_fmm,1e-9):.2f}x")

    # Overall pass: far-field (>= 5R) should be < 5%
    far_errors = errors[n_per:]  # distances >= 5R
    overall_pass = np.max(far_errors) < 0.05
    print(f"  Overall far-field (>= 5R): max={np.max(far_errors)*100:.2f}% [{'PASS' if overall_pass else 'FAIL'}]")
    return overall_pass


def test_near_field_info():
    """Test 2: Near-field - expect larger errors (informational)."""
    print("\n=== Test 2: Near-field info ===")
    container = create_test_model(n_div=3, mu_r=100)

    half = 0.015
    distances = [0.5*half, 1.0*half, 1.5*half, 2.0*half]
    pts = generate_observation_grid([0, 0, 0], distances, n_per_shell=50)

    result_direct = rad.FldBatch(container, pts, 0)
    result_fmm = rad.FldBatch(container, pts, 1)

    B_direct = np.array(result_direct['B'])
    B_fmm = np.array(result_fmm['B'])
    errors = vector_relative_error(B_direct, B_fmm)

    n_per = 50
    for i, d in enumerate(distances):
        errs_shell = errors[i*n_per:(i+1)*n_per]
        mean_err = np.mean(errs_shell) * 100
        max_err = np.max(errs_shell) * 100
        print(f"  Distance {d/half:.1f}R: mean={mean_err:.1f}%, max={max_err:.1f}% (info only)")

    print("  (Near-field errors expected to be large - dipole approximation not valid)")
    return True  # Informational only


def test_ima_fallback():
    """Test 3: IMA fallback - method=1 with IMA should give identical results to method=0."""
    print("\n=== Test 3: IMA fallback ===")
    rad.UtiDelAll()
    rad.FldUnits('m')

    # Create quarter model
    half = 0.015
    n_div = 3
    dx = half / n_div  # quarter: 0 to half

    elements = []
    for ix in range(n_div):
        for iy in range(n_div):
            for iz in range(n_div):
                x0 = ix * dx
                y0 = iy * dx
                z0 = iz * dx
                verts = [
                    [x0,    y0,    z0],
                    [x0+dx, y0,    z0],
                    [x0+dx, y0+dx, z0],
                    [x0,    y0+dx, z0],
                    [x0,    y0,    z0+dx],
                    [x0+dx, y0,    z0+dx],
                    [x0+dx, y0+dx, z0+dx],
                    [x0,    y0+dx, z0+dx],
                ]
                elem = rad.ObjHexahedron(verts, [0, 0, 0])
                elements.append(elem)

    container = rad.ObjCnt(elements)
    mat = rad.MatLin(999)
    rad.MatApl(container, mat)

    bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])
    assembly = rad.ObjCnt([container, bkg])

    # Solve with IMA
    rad.Solve(assembly, 0.001, 100, 0, image='+x-z')

    # Observation points outside
    pts = generate_observation_grid([0.03, 0.03, 0.03], [0.05, 0.1], n_per_shell=30)

    # Use assembly handle (same as Solve) so IMA context is detected
    result_direct = rad.FldBatch(assembly, pts, 0)
    result_fmm = rad.FldBatch(assembly, pts, 1)

    B_direct = np.array(result_direct['B'])
    B_fmm = np.array(result_fmm['B'])

    # With IMA active, method=1 should silently fall back to direct
    # -> results should be IDENTICAL (not just close)
    diff = np.max(np.abs(B_direct - B_fmm))
    passed = diff < 1e-15
    print(f"  Max absolute difference: {diff:.2e} [{'PASS' if passed else 'FAIL'}]")
    if passed:
        print("  method=1 correctly fell back to direct with IMA active")
    else:
        print("  WARNING: FMM did NOT fall back with IMA - results differ!")
    return passed


def test_utidelall_safety():
    """Test 4: UtiDelAll then FldBatch should not crash."""
    print("\n=== Test 4: UtiDelAll safety ===")
    _ = create_test_model(n_div=2)
    rad.UtiDelAll()

    # FMM cache should be cleared - this should not crash
    try:
        # Create a fresh simple object
        rad.FldUnits('m')
        mag = rad.ObjRecMag([0, 0, 0], [0.02, 0.02, 0.02], [0, 0, 954930])
        pts = np.array([[0.05, 0, 0], [0, 0.05, 0], [0, 0, 0.05]])
        result = rad.FldBatch(mag, pts, 1)
        B = np.array(result['B'])
        print(f"  Field at 3 points computed successfully: |B| = {np.linalg.norm(B, axis=1)}")
        print("  PASS")
        return True
    except Exception as e:
        print(f"  FAIL: {e}")
        return False


def test_timing_comparison():
    """Test 5: Timing comparison for varying N_points."""
    print("\n=== Test 5: Timing comparison ===")
    container = create_test_model(n_div=5, mu_r=1000)

    half = 0.015
    for n_points in [100, 1000, 5000]:
        pts = generate_observation_grid([0, 0, 0], [5*half, 10*half], n_per_shell=n_points//2)

        # Warm up
        _ = rad.FldBatch(container, pts[:10], 0)

        t0 = time.time()
        _ = rad.FldBatch(container, pts, 0)
        t_direct = time.time() - t0

        t0 = time.time()
        _ = rad.FldBatch(container, pts, 1)
        t_fmm = time.time() - t0

        speedup = t_direct / max(t_fmm, 1e-9)
        print(f"  N={n_points}: direct={t_direct:.4f}s, FMM={t_fmm:.4f}s, speedup={speedup:.2f}x")

    return True  # Informational


def main():
    print("=" * 60)
    print("FMM FldBatch Verification")
    print("=" * 60)

    results = {}
    results['far_field'] = test_far_field_accuracy()
    results['near_field'] = test_near_field_info()
    results['ima_fallback'] = test_ima_fallback()
    results['utidelall'] = test_utidelall_safety()
    results['timing'] = test_timing_comparison()

    print("\n" + "=" * 60)
    print("Summary:")
    for name, passed in results.items():
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: {status}")

    all_pass = all(results.values())
    print(f"\nOverall: {'ALL PASS' if all_pass else 'SOME FAILED'}")

    rad.UtiDelAll()
    return 0 if all_pass else 1


if __name__ == '__main__':
    sys.exit(main())
