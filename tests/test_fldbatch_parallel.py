#!/usr/bin/env python3
"""
Test script for batch field computation via unified Fld().
Verifies that batch Fld(obj, 'b', points_Nx3) produces correct results
matching single-point Fld(obj, 'b', point).
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src'))
import radia as rad
import time
import numpy as np

def run_test():
    """Test batch field computation."""
    print("=" * 70)
    print("Batch Field Computation Test")
    print("=" * 70)

    # Clear any existing objects
    rad.UtiDelAll()

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
        r = 0.15
        theta = np.random.uniform(0, np.pi, n_points)
        phi = np.random.uniform(0, 2*np.pi, n_points)
        points = np.column_stack([
            r * np.sin(theta) * np.cos(phi),
            r * np.sin(theta) * np.sin(phi),
            r * np.cos(theta),
        ])

        # Run batch Fld
        t0 = time.time()
        B_batch = np.asarray(rad.Fld(magnet, 'b', points))
        t1 = time.time()
        elapsed_ms = (t1 - t0) * 1000

        correct_shape = B_batch.shape == (n_points, 3)
        B0 = B_batch[0]
        B0_str = f"[{B0[0]:.2e}, {B0[1]:.2e}, {B0[2]:.2e}]"
        B_mag = np.linalg.norm(B0)
        is_correct = correct_shape and B_mag > 0

        status = "PASS" if is_correct else "FAIL"
        all_passed = all_passed and is_correct

        print(f"{n_points:<10} {elapsed_ms:<12.2f} {status:<10} {B0_str}")

    print("-" * 60)

    # Verification test: compare batch with single-point
    print("\n[Verification] Comparing batch Fld with single-point Fld()")
    test_points = np.array([
        [0, 0, 0.15],
        [0.1, 0, 0.1],
        [0, 0.1, 0.1],
    ])

    B_batch = np.asarray(rad.Fld(magnet, 'b', test_points))

    max_rel_error = 0.0
    for i in range(len(test_points)):
        B_single = np.asarray(rad.Fld(magnet, 'b', test_points[i]))
        diff_mag = np.linalg.norm(B_batch[i] - B_single)
        B_single_mag = np.linalg.norm(B_single)
        rel_error = diff_mag / B_single_mag if B_single_mag > 0 else 0
        max_rel_error = max(max_rel_error, rel_error)

        print(f"  Point {i+1}: Fld={B_single}")
        print(f"         Batch={B_batch[i]}")
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
