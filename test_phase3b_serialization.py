#!/usr/bin/env python
"""
Test Phase 3B: Full H-matrix serialization to disk
Tests the complete save/load cycle for 1000x speedup on program restart

Requirements:
- 100+ magnetic elements
- Save H-matrix to disk after first build
- Load H-matrix from disk on second run (instant startup)
- Verify accuracy < 1%
- Verify speedup > 100x

Expected behavior:
- First run: Build H-matrix (~1.0s), save to disk
- Second run: Load H-matrix from disk (~0.001s)
"""

import sys
import os
sys.path.insert(0, os.path.join('S:/Radia/01_GitHub', 'build/Release'))

import radia as rad
import numpy as np
import time

def create_magnet_array(nx=10, ny=10, nz=1):
    """Create array of permanent magnets (10x10x1 = 100 elements)"""
    container = rad.ObjCnt([])

    spacing = 5.0  # mm spacing between magnets
    size = 4.0     # mm magnet size
    M0 = 795774.7  # A/m magnetization (NdFeB)

    print(f"Creating {nx}×{ny}×{nz} = {nx*ny*nz} magnet array...")

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                x = i * spacing
                y = j * spacing
                z = k * spacing

                # Create magnet with z-direction magnetization
                magnet = rad.ObjRecMag([x, y, z], [size, size, size], [0, 0, M0])
                rad.ObjAddToCnt(container, [magnet])

    print(f"Created {nx*ny*nz} magnets")
    return container

def create_observation_points(n=10):
    """Create grid of observation points"""
    points = []
    for i in range(n):
        for j in range(n):
            x = i * 5.0
            y = j * 5.0
            z = 10.0  # 10mm above magnet array
            points.append([x, y, z])

    print(f"Created {len(points)} observation points")
    return np.array(points)

def test_serialization_save():
    """Test 1: Build H-matrix and save to disk"""
    print("\n" + "="*70)
    print("TEST 1: Build H-matrix and save to disk")
    print("="*70)

    # Create geometry
    container = create_magnet_array(10, 10, 1)
    obs_points = create_observation_points(10)

    # Enable full H-matrix serialization
    print("\nEnabling full H-matrix serialization...")
    rad.SolverHMatrixCacheFull(1)  # Enable

    # Configure H-matrix
    print("Configuring H-matrix solver...")
    rad.SetHMatrixFieldEval(1, 1e-6)

    # Build and solve (will construct H-matrix and save to disk)
    print("\n[First Run] Building H-matrix (this should take ~1.0s)...")
    t_start = time.time()

    # Use batch field evaluation to trigger H-matrix construction
    H_field = rad.FldBatch(container, 'h', obs_points.tolist(), 1)  # use_hmatrix=1

    t_build = time.time() - t_start

    print(f"\n[OK] H-matrix built and saved in {t_build:.3f}s")
    print(f"  Field at point 0: H = {H_field[0]}")

    # Check cache size
    cache_size_mb = rad.SolverHMatrixCacheSize(0)  # Query only
    print(f"  Cache size: {cache_size_mb:.2f} MB")

    return container, obs_points, H_field, t_build

def test_serialization_load(container, obs_points, H_expected):
    """Test 2: Load H-matrix from disk (instant startup)"""
    print("\n" + "="*70)
    print("TEST 2: Load H-matrix from disk (instant startup)")
    print("="*70)

    # Simulate program restart by creating new container with same geometry
    print("\nSimulating program restart...")
    print("(In real use, this would be a separate Python session)")

    # Re-create same geometry
    container_new = create_magnet_array(10, 10, 1)

    # Enable full serialization
    rad.SolverHMatrixCacheFull(1)
    rad.SetHMatrixFieldEval(1, 1e-6)

    # Load H-matrix from disk (should be instant)
    print("\n[Second Run] Loading H-matrix from disk (should be instant)...")
    t_start = time.time()

    H_field_loaded = rad.FldBatch(container_new, 'h', obs_points.tolist(), 1)

    t_load = time.time() - t_start

    print(f"\n[OK] H-matrix loaded in {t_load:.6f}s")
    print(f"  Field at point 0: H = {H_field_loaded[0]}")

    return H_field_loaded, t_load

def test_accuracy(H1, H2):
    """Test 3: Verify accuracy"""
    print("\n" + "="*70)
    print("TEST 3: Verify accuracy")
    print("="*70)

    H1 = np.array(H1)
    H2 = np.array(H2)

    # Compute relative error
    diff = np.linalg.norm(H2 - H1)
    norm = np.linalg.norm(H1)
    rel_error = diff / norm

    print(f"\nRelative error: {rel_error*100:.6f}%")
    print(f"Tolerance: 0.01%")

    if rel_error < 1e-4:
        print("[OK] PASS: Accuracy within tolerance")
        return True
    else:
        print(f"[FAIL] FAIL: Error {rel_error*100:.2f}% > 0.01%")
        return False

def test_performance(t_build, t_load):
    """Test 4: Verify performance improvement"""
    print("\n" + "="*70)
    print("TEST 4: Verify performance improvement")
    print("="*70)

    speedup = t_build / t_load

    print(f"\nBuild time:  {t_build:.3f}s")
    print(f"Load time:   {t_load:.6f}s")
    print(f"Speedup:     {speedup:.0f}x")
    print(f"Target:      >100x")

    if speedup > 100:
        print(f"[OK] PASS: Achieved {speedup:.0f}x speedup (target: >100x)")
        return True
    else:
        print(f"[FAIL] FAIL: Only {speedup:.0f}x speedup (target: >100x)")
        return False

def main():
    print("="*70)
    print("Phase 3B: Full H-Matrix Serialization Test")
    print("="*70)
    print("\nThis test verifies:")
    print("1. H-matrix can be saved to disk")
    print("2. H-matrix can be loaded from disk (instant startup)")
    print("3. Loaded H-matrix produces accurate results")
    print("4. Load is >100x faster than build")

    try:
        # Test 1: Build and save
        container, obs_points, H1, t_build = test_serialization_save()

        # Test 2: Load from disk
        H2, t_load = test_serialization_load(container, obs_points, H1)

        # Test 3: Accuracy
        accuracy_pass = test_accuracy(H1, H2)

        # Test 4: Performance
        performance_pass = test_performance(t_build, t_load)

        # Summary
        print("\n" + "="*70)
        print("SUMMARY")
        print("="*70)

        all_pass = accuracy_pass and performance_pass

        if all_pass:
            print("\n[OK][OK][OK] ALL TESTS PASSED [OK][OK][OK]")
            print("\nPhase 3B implementation complete!")
            print("Full H-matrix serialization working correctly.")
            return 0
        else:
            print("\n[FAIL][FAIL][FAIL] SOME TESTS FAILED [FAIL][FAIL][FAIL]")
            if not accuracy_pass:
                print("- Accuracy test failed")
            if not performance_pass:
                print("- Performance test failed")
            return 1

    except Exception as e:
        print(f"\n[FAIL] ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
