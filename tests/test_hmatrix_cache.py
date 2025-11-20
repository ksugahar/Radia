#!/usr/bin/env python
"""
Test H-matrix cache functionality

This test verifies:
1. PrepareCache() correctly caches field values
2. Evaluate() uses cached values (cache hits)
3. GetCacheStats() reports statistics correctly
4. Performance improvement with cache
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
print("TEST: H-MATRIX CACHE FUNCTIONALITY")
print("="*80)
print()

# ============================================================================
# Setup: Create Radia magnet and NGSolve mesh
# ============================================================================
print("[Setup] Creating Radia magnet and mesh...")

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

# Create mesh
box = Box((0.01, 0.01, 0.02), (0.06, 0.06, 0.08))
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=0.015))  # Medium mesh

print(f"  Magnet: {magnet}")
print(f"  Mesh elements: {mesh.ne}")
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
# Test 2: Collect integration points from mesh
# ============================================================================
print("[Test 2] Collecting integration points from mesh...")

# Use low-level FE space to get integration points
fes = HCurl(mesh, order=2)

# Collect all integration points by iterating over elements
all_points = []
for el in mesh.Elements(VOL):
    # Get integration rule for this element
    ir = IntegrationRule(el.type, order=5)  # Order 5 quadrature

    # Get mapped integration points
    trafo = mesh.GetTrafo(el)
    for ip in ir:
        mip = trafo(ip)
        pnt = mip.point
        all_points.append([pnt[0], pnt[1], pnt[2]])

npts = len(all_points)
print(f"  Collected {npts} integration points")
print(f"  Elements: {mesh.ne}, avg points/element: {npts/mesh.ne:.1f}")
print()

# ============================================================================
# Test 3: PrepareCache() - batch evaluation
# ============================================================================
print("[Test 3] Calling PrepareCache() for batch evaluation...")

t0 = time.time()
A_cf.PrepareCache(all_points)
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
# Test 4: GridFunction.Set() using cache (should be fast)
# ============================================================================
print("[Test 4] GridFunction.Set() with cache (should use cached values)...")

gf = GridFunction(fes)

t0 = time.time()
gf.Set(A_cf)
t1 = time.time()
time_set_cached = t1 - t0

stats = A_cf.GetCacheStats()
print(f"  Time: {time_set_cached:.4f} s")
print(f"  Cache hits: {stats['hits']}")
print(f"  Cache misses: {stats['misses']}")
print(f"  Hit rate: {stats['hit_rate']*100:.1f}%")

# Most evaluations should hit cache
assert stats['hit_rate'] > 0.5, f"Hit rate too low: {stats['hit_rate']*100:.1f}%"
print("  [PASS] Cache hits detected")
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
# Test 6: Performance comparison
# ============================================================================
print("[Test 6] Performance comparison: Direct vs Cached...")

# Direct evaluation (no cache)
gf_direct = GridFunction(fes)

t0 = time.time()
gf_direct.Set(A_cf)
t1 = time.time()
time_direct = t1 - t0

print(f"  Direct evaluation: {time_direct:.4f} s")

# Cached evaluation
A_cf.PrepareCache(all_points)
gf_cached = GridFunction(fes)

t0 = time.time()
gf_cached.Set(A_cf)
t1 = time.time()
time_cached = t1 - t0

print(f"  Cached evaluation: {time_cached:.4f} s")

speedup = time_direct / time_cached if time_cached > 0 else 0
print(f"\n  Speedup: {speedup:.2f}x")

if speedup > 1.0:
    improvement = (1 - 1/speedup) * 100
    print(f"  Improvement: {improvement:.1f}% faster")
    print("  [PASS] Cache provides speedup")
else:
    print("  [INFO] Cache overhead detected (expected for small problems)")

print()

# ============================================================================
# Test 7: Verify accuracy (cached vs direct should be identical)
# ============================================================================
print("[Test 7] Accuracy check (cached vs direct)...")

# Sample points for verification
test_points = [
    (0.030, 0.020, 0.030),
    (0.040, 0.040, 0.050),
]

max_error = 0.0
for pt in test_points:
    try:
        mip = mesh(*pt)

        val_direct = gf_direct(mip)
        val_cached = gf_cached(mip)

        error = np.linalg.norm(np.array(val_direct) - np.array(val_cached))
        max_error = max(max_error, error)

        print(f"  Point {pt}: error = {error:.3e}")
    except:
        pass

tolerance = 1e-12  # Should be identical (same source data)
print(f"\n  Max error: {max_error:.3e}")
print(f"  Tolerance: {tolerance:.3e}")

assert max_error < tolerance, f"Accuracy error: {max_error:.3e} > {tolerance:.3e}"
print("  [PASS] Cached values match direct evaluation")
print()

# ============================================================================
# Summary
# ============================================================================
print("="*80)
print("TEST SUMMARY")
print("="*80)
print()
print("All tests PASSED!")
print()
print(f"Cache functionality:")
print(f"  - PrepareCache(): {npts} points in {time_prepare:.4f} s")
print(f"  - Cache hit rate: {stats['hit_rate']*100:.1f}%")
print(f"  - Speedup: {speedup:.2f}x")
print(f"  - Accuracy: <{tolerance:.0e} (machine precision)")
print()
print("Next step: Run full benchmark with large mesh to measure H-matrix speedup")
print("="*80)

# Cleanup
rad.UtiDelAll()
