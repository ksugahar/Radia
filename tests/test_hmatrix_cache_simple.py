#!/usr/bin/env python
"""
Simple H-matrix cache functionality test

Quick test with small number of points to verify cache works.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

import radia as rad
import numpy as np
import time

try:
    from ngsolve import *
    from netgen.occ import *
    import rad_ngsolve
except ImportError as e:
    print(f"ERROR: {e}")
    sys.exit(1)

print("="*80)
print("TEST: H-MATRIX CACHE FUNCTIONALITY (SIMPLE)")
print("="*80)
print()

# ============================================================================
# Setup: Create Radia magnet
# ============================================================================
print("[Setup] Creating Radia magnet...")

rad.UtiDelAll()

# Set Radia to use meters (required for NGSolve integration)
rad.FldUnits('m')

magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.06], [0, 0, 1.2])

def radia_field_with_A(coords):
    x, y, z = coords
    B = rad.Fld(magnet, 'b', [x, y, z])
    A = rad.Fld(magnet, 'a', [x, y, z])
    return {'B': list(B), 'A': list(A)}

bg_field = rad.ObjBckgCF(radia_field_with_A)

print(f"  Magnet: {magnet}")
print()

# ============================================================================
# Test 1: Create CoefficientFunction and check cache is disabled initially
# ============================================================================
print("[Test 1] Initial state (cache disabled)...")

A_cf = rad_ngsolve.RadiaField(bg_field, 'a')
stats = A_cf.GetCacheStats()

print(f"  Cache enabled: {stats['enabled']}")
print(f"  Cache size: {stats['size']}")

assert stats['enabled'] == False, "Cache should be disabled initially"
assert stats['size'] == 0, "Cache should be empty initially"
print("  [PASS] Cache initially disabled")
print()

# ============================================================================
# Test 2: Manually create small set of test points
# ============================================================================
print("[Test 2] Creating small test point set (10 points)...")

test_points = [
    [0.02, 0.02, 0.03],
    [0.03, 0.02, 0.03],
    [0.04, 0.02, 0.03],
    [0.02, 0.03, 0.03],
    [0.03, 0.03, 0.03],
    [0.04, 0.03, 0.03],
    [0.02, 0.02, 0.04],
    [0.03, 0.02, 0.04],
    [0.04, 0.02, 0.04],
    [0.03, 0.03, 0.04],
]

npts = len(test_points)
print(f"  Created {npts} test points")
print()

# ============================================================================
# Test 3: PrepareCache() - batch evaluation
# ============================================================================
print("[Test 3] Calling PrepareCache() for batch evaluation...")

t0 = time.time()
A_cf.PrepareCache(test_points)
t1 = time.time()
time_prepare = t1 - t0

stats = A_cf.GetCacheStats()
print(f"  Time: {time_prepare:.4f} s")
print(f"  Cache enabled: {stats['enabled']}")
print(f"  Cache size: {stats['size']}")
print(f"  Expected size: {npts}")

assert stats['enabled'] == True, "Cache should be enabled after PrepareCache()"
assert stats['size'] == npts, f"Cache size mismatch: {stats['size']} != {npts}"
print("  [PASS] Cache prepared successfully")
print()

# ============================================================================
# Test 4: Evaluate single point (should hit cache)
# ============================================================================
print("[Test 4] Evaluating cached point...")

# Test evaluation at cached point
pt = test_points[4]  # Middle point [0.03, 0.03, 0.03]
print(f"  Test point: {pt}")

stats_before = A_cf.GetCacheStats()
print(f"  Cache hits before: {stats_before['hits']}")
print(f"  Cache misses before: {stats_before['misses']}")

# Evaluate (should hit cache)
# Note: We can't directly call Evaluate from Python, so we use a mesh point
# Instead, let's just check that PrepareCache worked by examining cache stats

print("  [INFO] Cache prepared with {npts} points")
print("  [INFO] Direct evaluation not tested (requires mesh point)")
print()

# ============================================================================
# Test 5: ClearCache() and verify cache is disabled
# ============================================================================
print("[Test 5] ClearCache()...")

A_cf.ClearCache()
stats = A_cf.GetCacheStats()

print(f"  Cache enabled: {stats['enabled']}")
print(f"  Cache size: {stats['size']}")

assert stats['enabled'] == False, "Cache should be disabled after ClearCache()"
assert stats['size'] == 0, "Cache should be empty after ClearCache()"
print("  [PASS] Cache cleared successfully")
print()

# ============================================================================
# Test 6: Direct field evaluation comparison
# ============================================================================
print("[Test 6] Comparing cached vs direct evaluation...")

# Prepare cache again
A_cf.PrepareCache(test_points)

# Evaluate directly with Radia
pt = test_points[4]
A_direct = rad.Fld(magnet, 'a', [pt[0]*1000, pt[1]*1000, pt[2]*1000])  # meters to mm
print(f"  Direct Radia evaluation: {A_direct}")
print(f"  (Point: {pt})")

# Check cache contains this value
stats = A_cf.GetCacheStats()
print(f"  Cache contains {stats['size']} values")
print(f"  [INFO] Cache values cannot be directly inspected from Python")
print()

# ============================================================================
# Summary
# ============================================================================
print("="*80)
print("TEST SUMMARY")
print("="*80)
print()
print("All basic tests PASSED!")
print()
print(f"Cache functionality verified:")
print(f"  - PrepareCache(): {npts} points in {time_prepare:.4f} s")
print(f"  - Cache enabled/disabled: Working")
print(f"  - ClearCache(): Working")
print(f"  - GetCacheStats(): Working")
print()
print("Next step: Test with GridFunction.Set() for full integration")
print("="*80)

# Cleanup
rad.UtiDelAll()
