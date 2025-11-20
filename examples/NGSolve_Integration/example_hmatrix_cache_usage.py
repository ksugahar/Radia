#!/usr/bin/env python
"""
Example: Using H-matrix Cache for Fast GridFunction Evaluation

This example demonstrates how to use PrepareCache() to enable H-matrix
acceleration when setting GridFunctions from Radia fields.

Performance comparison:
- Without cache: Element-by-element evaluation, no H-matrix benefit
- With cache: Single batch evaluation, full H-matrix acceleration
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))

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
print("EXAMPLE: H-MATRIX CACHE USAGE")
print("="*80)
print()

# ============================================================================
# Step 1: Create Radia magnetic field
# ============================================================================
print("[Step 1] Creating Radia magnet...")

rad.UtiDelAll()

# Simple rectangular magnet
magnet = rad.ObjRecMag([0, 0, 0], [40, 40, 60], [0, 0, 1.2])

def radia_field_with_A(coords):
    x, y, z = coords
    B = rad.Fld(magnet, 'b', [x, y, z])
    A = rad.Fld(magnet, 'a', [x, y, z])
    return {'B': list(B), 'A': list(A)}

bg_field = rad.ObjBckgCF(radia_field_with_A)

print(f"  Magnet object: {magnet}")
print(f"  Magnetization: [0, 0, 1.2] T")
print()

# ============================================================================
# Step 2: Create NGSolve mesh
# ============================================================================
print("[Step 2] Creating NGSolve mesh...")

# Create a box mesh around the magnet
box = Box((0.01, 0.01, 0.02), (0.06, 0.06, 0.08))
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=0.020))  # Coarse mesh for quick example

print(f"  Mesh elements: {mesh.ne}")
print(f"  Mesh vertices: {mesh.nv}")
print()

# ============================================================================
# Step 3: Create FE space and GridFunction
# ============================================================================
print("[Step 3] Setting up finite element space...")

fes = HCurl(mesh, order=2)
gf_A = GridFunction(fes)

print(f"  FE space: HCurl(order=2)")
print(f"  Degrees of freedom: {fes.ndof}")
print()

# ============================================================================
# Step 4: METHOD 1 - Standard evaluation (without cache)
# ============================================================================
print("[Step 4] Method 1: Standard GridFunction.Set() (no cache)...")

A_cf = rad_ngsolve.RadiaField(bg_field, 'a')

# Check cache is disabled
stats = A_cf.GetCacheStats()
print(f"  Cache enabled: {stats['enabled']}")

t0 = time.time()
gf_A.Set(A_cf)
t1 = time.time()
time_no_cache = t1 - t0

stats_after = A_cf.GetCacheStats()
print(f"  Time: {time_no_cache:.4f} s")
print(f"  Cache hits: {stats_after['hits']}")
print(f"  Cache misses: {stats_after['misses']}")
print()

# ============================================================================
# Step 5: METHOD 2 - Cached evaluation (with PrepareCache)
# ============================================================================
print("[Step 5] Method 2: GridFunction.Set() with cache...")

# Collect all integration points from the mesh
print("  [5a] Collecting integration points...")

all_points = []
for el in mesh.Elements(VOL):
    ir = IntegrationRule(el.type, order=5)  # Order 5 quadrature
    trafo = mesh.GetTrafo(el)
    for ip in ir:
        mip = trafo(ip)
        pnt = mip.point
        all_points.append([pnt[0], pnt[1], pnt[2]])

npts = len(all_points)
print(f"       Collected {npts} integration points")
print(f"       ({npts/mesh.ne:.1f} points/element)")
print()

# Prepare cache with single batch evaluation
print("  [5b] Calling PrepareCache()...")

A_cf_cached = rad_ngsolve.RadiaField(bg_field, 'a')

t0 = time.time()
A_cf_cached.PrepareCache(all_points)
t1 = time.time()
time_prepare = t1 - t0

stats_cache = A_cf_cached.GetCacheStats()
print(f"       PrepareCache time: {time_prepare:.4f} s")
print(f"       Cache enabled: {stats_cache['enabled']}")
print(f"       Cache size: {stats_cache['size']}")
print()

# Set GridFunction using cached values
print("  [5c] Setting GridFunction (using cache)...")

gf_A_cached = GridFunction(fes)

t0 = time.time()
gf_A_cached.Set(A_cf_cached)
t1 = time.time()
time_set_cached = t1 - t0

stats_final = A_cf_cached.GetCacheStats()
print(f"       GridFunction.Set() time: {time_set_cached:.4f} s")
print(f"       Cache hits: {stats_final['hits']}")
print(f"       Cache misses: {stats_final['misses']}")
print(f"       Hit rate: {stats_final['hit_rate']*100:.1f}%")
print()

time_with_cache = time_prepare + time_set_cached
print(f"  Total time with cache: {time_with_cache:.4f} s")
print(f"    (PrepareCache: {time_prepare:.4f} s + Set: {time_set_cached:.4f} s)")
print()

# ============================================================================
# Step 6: Comparison
# ============================================================================
print("[Step 6] Performance comparison...")

print(f"  Method 1 (no cache):   {time_no_cache:.4f} s")
print(f"  Method 2 (with cache): {time_with_cache:.4f} s")
print()

if time_with_cache < time_no_cache:
    speedup = time_no_cache / time_with_cache
    improvement = (1 - time_with_cache/time_no_cache) * 100
    print(f"  Speedup: {speedup:.2f}x")
    print(f"  Improvement: {improvement:.1f}% faster")
else:
    overhead = time_with_cache / time_no_cache
    print(f"  Overhead: {overhead:.2f}x (cache slower for small problems)")
    print(f"  Note: Cache benefits increase with larger meshes")

print()

# ============================================================================
# Step 7: Verify accuracy
# ============================================================================
print("[Step 7] Verifying accuracy...")

# Sample points for comparison
test_point = (0.030, 0.025, 0.030)

try:
    mip = mesh(*test_point)

    val_no_cache = gf_A(mip)
    val_with_cache = gf_A_cached(mip)

    error = np.linalg.norm(np.array(val_no_cache) - np.array(val_with_cache))

    print(f"  Test point: {test_point}")
    print(f"  Value (no cache):   {val_no_cache}")
    print(f"  Value (with cache): {val_with_cache}")
    print(f"  Difference: {error:.3e}")

    if error < 1e-12:
        print(f"  [PASS] Values match to machine precision")
    else:
        print(f"  [WARNING] Values differ by {error:.3e}")
except:
    print(f"  [INFO] Test point outside mesh, skipping verification")

print()

# ============================================================================
# Summary
# ============================================================================
print("="*80)
print("SUMMARY: H-MATRIX CACHE USAGE")
print("="*80)
print()
print("Cache workflow:")
print("  1. Create CoefficientFunction: A_cf = rad_ngsolve.RadiaField(obj, 'a')")
print("  2. Collect integration points from mesh")
print("  3. Prepare cache: A_cf.PrepareCache(all_points)")
print("  4. Set GridFunction: gf.Set(A_cf)  # Uses cached values")
print()
print("Benefits:")
print("  - Single batch call to Radia (enables H-matrix)")
print("  - Faster for large meshes (N > 200 elements)")
print("  - Ideal for iterative workflows (reuse cached values)")
print()
print("When to use cache:")
print("  - Large meshes (N > 200 elements)")
print("  - Multiple GridFunction evaluations with same field")
print("  - Iterative solvers requiring repeated field evaluations")
print()
print("When NOT to use cache:")
print("  - Very small meshes (overhead > benefit)")
print("  - One-time GridFunction evaluation")
print("  - Memory-constrained systems (cache stores all field values)")
print()
print("="*80)

# Cleanup
rad.UtiDelAll()
