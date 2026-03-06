#!/usr/bin/env python3
"""
Test script for FldBatch batch field computation.
Verifies that batch field computation produces correct results matching single-point Fld().

Note: OpenMP parallelization is currently disabled due to Intel OpenMP runtime deadlock issues.
The batch API still provides significant speedup by reducing Python-C++ call overhead.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/radia'))
import radia as rad
import time
import numpy as np

def run_test():
    """Test FldBatch batch field computation."""
    print("=" * 70)
    print("FldBatch Batch Field Computation Test")
    print("=" * 70)

    # Clear any existing objects
    rad.UtiDelAll()
    rad.FldUnits('m')

    # Create a simple hexahedral magnet
    s = 0.05  # 5 cm cube
    vertices = [
        [-s, -s, -s], [s, -s, -s], [s, s, -s], [-s, s, -s],
        [-s, -s, s], [s, -s, s], [s, s, s], [-s, s, s]
    ]
    magnet = rad.ObjHexahedron(vertices, [0, 0, 1e6])  # 1 MA/m magnetization
    print(f"\n[Setup] Created hexahedral magnet: {magnet}")

    # Test different point counts
    test_configs = [
        10,     # Below OpenMP threshold
        100,    # At OpenMP threshold
        500,    # Above threshold
        1000,   # Large batch
        5000,   # Very large batch
    ]

    print(f"\n{'Points':<10} {'Time (ms)':<12} {'Correct':<10} {'B[0] (T)'}")
    print("-" * 60)

    all_passed = True

    for n_points in test_configs:
        # Generate random observation points outside the magnet
        np.random.seed(42)  # Reproducible results
        points = []
        for i in range(n_points):
            # Points on a sphere around the magnet at radius 0.15m
            r = 0.15
            theta = np.random.uniform(0, np.pi)
            phi = np.random.uniform(0, 2*np.pi)
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            points.append([x, y, z])

        # Run FldBatch
        t0 = time.time()
        result = rad.FldBatch(magnet, points, 0)
        t1 = time.time()
        elapsed_ms = (t1 - t0) * 1000

        # Verify result structure
        has_B = "B" in result
        has_H = "H" in result
        correct_shape = len(result["B"]) == n_points if has_B else False

        # Check first point result is reasonable
        B0 = result["B"][0] if has_B else [0, 0, 0]
        B0_str = f"[{B0[0]:.2e}, {B0[1]:.2e}, {B0[2]:.2e}]"

        # Basic sanity check: field should be non-zero
        B_mag = np.sqrt(B0[0]**2 + B0[1]**2 + B0[2]**2)
        is_correct = has_B and has_H and correct_shape and B_mag > 0

        status = "PASS" if is_correct else "FAIL"
        all_passed = all_passed and is_correct

        print(f"{n_points:<10} {elapsed_ms:<12.2f} {status:<10} {B0_str}")

    print("-" * 60)

    # Verification test: compare with single-point Fld()
    print("\n[Verification] Comparing FldBatch with single-point Fld()")
    test_points = [
        [0, 0, 0.15],
        [0.1, 0, 0.1],
        [0, 0.1, 0.1],
    ]

    result_batch = rad.FldBatch(magnet, test_points, 0)

    max_rel_error = 0.0
    for i, pt in enumerate(test_points):
        B_single = rad.Fld(magnet, 'b', pt)
        B_batch = result_batch["B"][i]

        # Calculate relative error
        B_single_mag = np.sqrt(sum(b**2 for b in B_single))
        diff = [B_batch[j] - B_single[j] for j in range(3)]
        diff_mag = np.sqrt(sum(d**2 for d in diff))

        rel_error = diff_mag / B_single_mag if B_single_mag > 0 else 0
        max_rel_error = max(max_rel_error, rel_error)

        print(f"  Point {i+1}: Fld={[f'{b:.6e}' for b in B_single]}")
        print(f"         Batch={[f'{b:.6e}' for b in B_batch]}")
        print(f"         Rel error: {rel_error:.2e}")

    accuracy_pass = max_rel_error < 1e-10
    print(f"\n  Maximum relative error: {max_rel_error:.2e}")
    print(f"  Accuracy test: {'PASS' if accuracy_pass else 'FAIL'}")
    all_passed = all_passed and accuracy_pass

    print("\n" + "=" * 70)
    if all_passed:
        print("[SUCCESS] All tests passed!")
        return 0
    else:
        print("[FAILED] Some tests failed!")
        return 1

if __name__ == "__main__":
    sys.exit(run_test())
