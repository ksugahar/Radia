"""
Simple test to verify Phase 3 disk cache functionality
"""

import sys
import os
import time

# Add paths for Radia module
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'dist'))

import radia as rad

print("="*80)
print("PHASE 3: DISK CACHE TEST")
print("="*80)

def create_soft_iron_cube(n=7):
    """Create a soft iron cube with n^3 subdivisions"""
    # Soft iron material (standard soft iron parameters)
    mat = rad.MatSatIsoFrm([1596.3, 1.1488], [133.11, 0.4268], [18.713, 0.4759])

    # Create cube
    cube = rad.ObjRecMag([0, 0, 0], [100, 100, 100], [0, 0, 0])
    rad.MatApl(cube, mat)

    # Subdivide
    rad.ObjDivMag(cube, [n, n, n])

    print(f"Created soft iron cube: {n}x{n}x{n} = {n**3} elements")
    return cube

# Test 1: First run - should build H-matrix and save to cache
print("\n" + "="*80)
print("TEST 1: First solve - should save to cache")
print("="*80)

rad.UtiDelAll()
cube1 = create_soft_iron_cube(7)
rad.SolverHMatrixEnable(1, 1e-4, 30)

print("\nRunning solver...")
print("-"*80)
t0 = time.time()
result1 = rad.Solve(cube1, 0.0001, 1000)
t1 = time.time() - t0
print("-"*80)
print(f"Solve time: {t1*1000:.1f} ms")
print(f"Result: {result1}")

# Test 2: Same geometry - should find cache hit
print("\n" + "="*80)
print("TEST 2: Same geometry (new session) - should show cache hit")
print("="*80)

rad.UtiDelAll()
cube2 = create_soft_iron_cube(7)
rad.SolverHMatrixEnable(1, 1e-4, 30)

print("\nRunning solver...")
print("-"*80)
t0 = time.time()
result2 = rad.Solve(cube2, 0.0001, 1000)
t2 = time.time() - t0
print("-"*80)
print(f"Solve time: {t2*1000:.1f} ms")
print(f"Result: {result2}")

# Test 3: Different geometry - should be cache miss
print("\n" + "="*80)
print("TEST 3: Different geometry (6x6x6) - should be cache miss")
print("="*80)

rad.UtiDelAll()
cube3 = create_soft_iron_cube(6)
rad.SolverHMatrixEnable(1, 1e-4, 30)

print("\nRunning solver...")
print("-"*80)
t0 = time.time()
result3 = rad.Solve(cube3, 0.0001, 1000)
t3 = time.time() - t0
print("-"*80)
print(f"Solve time: {t3*1000:.1f} ms")
print(f"Result: {result3}")

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print("Look for these console messages:")
print("  [Phase 3] Cache miss (new geometry, hash=...) - First run")
print("  [Phase 3] Cache hit! (hash=...) - Second run (same geometry)")
print("  [Phase 3] Saved to cache (.../.radia_cache/hmatrix_cache.bin)")
print("")
print("Cache file should exist at: ./.radia_cache/hmatrix_cache.bin")
print("="*80)
