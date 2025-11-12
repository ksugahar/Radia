"""
Test Phase 2-B: Geometry Change Detection

This test validates:
1. Geometry hash computation
2. H-matrix reuse when geometry unchanged
3. H-matrix rebuild when geometry changes
4. Adaptive parameter selection based on problem size
"""

import sys
import os
import time

# Add paths for Radia module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'dist'))

import radia as rad

print("="*80)
print("PHASE 2-B: GEOMETRY CHANGE DETECTION TEST")
print("="*80)

def test_geometry_reuse():
    """Test H-matrix reuse with unchanged geometry"""
    print("\n" + "="*80)
    print("TEST 1: H-Matrix Reuse (Geometry Unchanged)")
    print("="*80)

    rad.UtiDelAll()

    # Create geometry (N=343)
    print("\n[1] Creating geometry: 7x7x7 = 343 elements...")
    n = 7
    container = rad.ObjCnt([])

    for i in range(n):
        for j in range(n):
            for k in range(n):
                x, y, z = i*10, j*10, k*10
                block = rad.ObjRecMag([x, y, z], [8, 8, 8], [0, 0, 1])
                rad.ObjAddToCnt(container, [block])

    print(f"    Created N={n**3} elements")

    # Enable H-matrix
    rad.SolverHMatrixEnable(1, 1e-4, 30)

    # First solve - should build H-matrix
    print("\n[2] First solve (H-matrix will be built)...")
    print("-"*80)
    t0 = time.time()
    result1 = rad.Solve(container, 0.0001, 1000)
    t1 = time.time() - t0
    print("-"*80)
    print(f"    Result: {result1}")
    print(f"    Time: {t1*1000:.1f} ms")

    # Second solve - should reuse H-matrix
    print("\n[3] Second solve (H-matrix should be reused)...")
    print("-"*80)
    t0 = time.time()
    result2 = rad.Solve(container, 0.0001, 1000)
    t2 = time.time() - t0
    print("-"*80)
    print(f"    Result: {result2}")
    print(f"    Time: {t2*1000:.1f} ms")

    # Compare times
    print(f"\n[4] Performance comparison:")
    print(f"    First solve:  {t1*1000:.1f} ms (with H-matrix construction)")
    print(f"    Second solve: {t2*1000:.1f} ms (H-matrix reused)")
    if t1 > t2:
        speedup = t1 / t2
        print(f"    Speedup: {speedup:.1f}x (reuse working!)")

    return t1, t2

def test_geometry_change():
    """Test H-matrix rebuild when geometry changes"""
    print("\n" + "="*80)
    print("TEST 2: H-Matrix Rebuild (Geometry Changed)")
    print("="*80)

    rad.UtiDelAll()

    # Create initial geometry
    print("\n[1] Creating initial geometry: 7x7x7 = 343 elements...")
    n = 7
    container = rad.ObjCnt([])

    for i in range(n):
        for j in range(n):
            for k in range(n):
                x, y, z = i*10, j*10, k*10
                block = rad.ObjRecMag([x, y, z], [8, 8, 8], [0, 0, 1])
                rad.ObjAddToCnt(container, [block])

    # Enable H-matrix
    rad.SolverHMatrixEnable(1, 1e-4, 30)

    # First solve
    print("\n[2] First solve with initial geometry...")
    print("-"*80)
    t0 = time.time()
    result1 = rad.Solve(container, 0.0001, 1000)
    t1 = time.time() - t0
    print("-"*80)
    print(f"    Time: {t1*1000:.1f} ms")

    # Change geometry
    print("\n[3] Changing geometry (adding 1 element)...")
    new_block = rad.ObjRecMag([100, 100, 100], [8, 8, 8], [0, 0, 1])
    rad.ObjAddToCnt(container, [new_block])
    print(f"    New element count: {n**3 + 1}")

    # Second solve - should rebuild H-matrix
    print("\n[4] Second solve with changed geometry...")
    print("-"*80)
    t0 = time.time()
    result2 = rad.Solve(container, 0.0001, 1000)
    t2 = time.time() - t0
    print("-"*80)
    print(f"    Time: {t2*1000:.1f} ms")

    print(f"\n[5] Expected behavior:")
    print(f"    - First solve: H-matrix built")
    print(f"    - Second solve: H-matrix rebuilt (geometry changed)")
    print(f"    - Look for '[Phase 2-B] Geometry changed' message above")

def test_adaptive_parameters():
    """Test adaptive parameter selection for different problem sizes"""
    print("\n" + "="*80)
    print("TEST 3: Adaptive Parameter Selection")
    print("="*80)

    sizes = [
        (5, 125),   # Small-medium
        (8, 512),   # Large
        (10, 1000), # Very large
    ]

    for n, total in sizes:
        print(f"\n[Test N={total}] Creating {n}x{n}x{n} = {total} elements...")
        rad.UtiDelAll()

        container = rad.ObjCnt([])
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    x, y, z = i*10, j*10, k*10
                    block = rad.ObjRecMag([x, y, z], [8, 8, 8], [0, 0, 1])
                    rad.ObjAddToCnt(container, [block])

        rad.SolverHMatrixEnable(1, 0, 0)  # Use adaptive parameters (0 means auto)

        print(f"    Solving...")
        print("-"*80)
        result = rad.Solve(container, 0.0001, 1000)
        print("-"*80)
        print(f"    Look for '[Phase 2-B] H-matrix parameters' message above")

if __name__ == "__main__":
    # Test 1: Reuse with unchanged geometry
    t1, t2 = test_geometry_reuse()

    # Test 2: Rebuild when geometry changes
    test_geometry_change()

    # Test 3: Adaptive parameter selection
    test_adaptive_parameters()

    print("\n" + "="*80)
    print("PHASE 2-B VALIDATION COMPLETE")
    print("="*80)
    print("\nExpected console messages:")
    print("  [Phase 2-B] Reusing H-matrix (geometry unchanged, hash=...)")
    print("  [Phase 2-B] Geometry changed (hash: ... -> ...), rebuilding...")
    print("  [Phase 2-B] H-matrix parameters: eps=..., max_rank=... (N=...)")
    print("\nIf all these messages appeared, Phase 2-B is working correctly!")
    print("="*80)
