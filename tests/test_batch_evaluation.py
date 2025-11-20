#!/usr/bin/env python3
"""
Test script for PrepareCache() batch evaluation optimization
Verifies that H-matrix acceleration works with GridFunction.Set()
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../src/python'))

import radia as rad
import rad_ngsolve
from ngsolve import *
from netgen.occ import *
import time
import numpy as np

print("="*70)
print("Batch Evaluation Performance Test")
print("Testing PrepareCache() for H-matrix acceleration")
print("="*70)

# Create magnet array (N=125 elements)
rad.UtiDelAll()

# Set Radia to use meters (required for NGSolve integration)
rad.FldUnits('m')

n = 5
elements = []
for i in range(n):
    for j in range(n):
        for k in range(n):
            x = (i - n/2 + 0.5) * 0.02
            y = (j - n/2 + 0.5) * 0.02
            z = (k - n/2 + 0.5) * 0.02
            elem = rad.ObjRecMag([x, y, z], [0.02, 0.02, 0.02], [0, 0, 1.2])
            elements.append(elem)

magnet = rad.ObjCnt(elements)
print(f"\n[Setup] Created magnet array: {n**3} elements")

# Create mesh
box = Box((0.015, 0.015, 0.015), (0.063, 0.063, 0.063))
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.010))
print(f"[Setup] Mesh: {mesh.ne} elements, {mesh.nv} vertices")

# Enable H-matrix
rad.SetHMatrixFieldEval(1, 1e-6)
print(f"[Setup] H-matrix enabled (eps=1e-6)")

# Test 1: Standard evaluation (element-by-element)
print("\n" + "="*70)
print("TEST 1: Standard GridFunction.Set() (element-by-element)")
print("="*70)

fes = HDiv(mesh, order=2)
B_gf_standard = GridFunction(fes)
B_cf_standard = rad_ngsolve.RadiaField(magnet, 'b')

print(f"[Test 1] FE Space: HDiv order=2, {fes.ndof} DOFs")

t0 = time.time()
B_gf_standard.Set(B_cf_standard)
t_standard = time.time() - t0

print(f"[Test 1] Time: {t_standard*1000:.1f} ms")

# Test 2: Batch evaluation with PrepareCache()
print("\n" + "="*70)
print("TEST 2: Optimized GridFunction.Set() with PrepareCache()")
print("="*70)

B_gf_batch = GridFunction(fes)
B_cf_batch = rad_ngsolve.RadiaField(magnet, 'b')

print("[Test 2] Calling PrepareCache()...")
t0 = time.time()
B_cf_batch.PrepareCache(mesh)
t_cache_prep = time.time() - t0

print(f"[Test 2] PrepareCache time: {t_cache_prep*1000:.1f} ms")

print("[Test 2] Calling GridFunction.Set()...")
t0 = time.time()
B_gf_batch.Set(B_cf_batch)
t_batch_set = time.time() - t0

t_batch_total = t_cache_prep + t_batch_set

print(f"[Test 2] Set() time: {t_batch_set*1000:.1f} ms")
print(f"[Test 2] Total time: {t_batch_total*1000:.1f} ms")

print("\n[Test 2] Cache statistics:")
B_cf_batch.PrintCacheStats()

# Verify accuracy
test_points = [
    (0.030, 0.020, 0.040),
    (0.040, 0.040, 0.050),
    (0.050, 0.030, 0.060),
]

print("\n" + "="*70)
print("ACCURACY VERIFICATION")
print("="*70)

errors = []
for i, pt in enumerate(test_points):
    try:
        mip = mesh(*pt)
        B_std = np.array(B_gf_standard(mip))
        B_bat = np.array(B_gf_batch(mip))

        error = np.linalg.norm(B_bat - B_std)
        B_norm = np.linalg.norm(B_std)
        rel_error = error / B_norm * 100 if B_norm > 0 else 0
        errors.append(rel_error)

        print(f"  Point {i+1} {pt}:")
        print(f"    Standard: [{B_std[0]:.6f}, {B_std[1]:.6f}, {B_std[2]:.6f}] T")
        print(f"    Batch:    [{B_bat[0]:.6f}, {B_bat[1]:.6f}, {B_bat[2]:.6f}] T")
        print(f"    Error: {rel_error:.6f}%")
    except Exception as e:
        print(f"  Point {i+1} {pt}: [FAILED] {e}")

# Performance summary
print("\n" + "="*70)
print("PERFORMANCE SUMMARY")
print("="*70)

if t_batch_total > 0:
    speedup = t_standard / t_batch_total
else:
    speedup = 0

print(f"\n  Method                Time (ms)    Speedup")
print(f"  {'Standard':<18s}  {t_standard*1000:8.1f}       1.0x")
print(f"  {'Batch (total)':<18s}  {t_batch_total*1000:8.1f}       {speedup:.1f}x")
print(f"    PrepareCache:       {t_cache_prep*1000:8.1f}")
print(f"    Set():              {t_batch_set*1000:8.1f}")

if errors:
    mean_error = np.mean(errors)
    max_error = np.max(errors)
    print(f"\n  Accuracy:")
    print(f"    Mean error: {mean_error:.6f}%")
    print(f"    Max error:  {max_error:.6f}%")

print("\n" + "="*70)
print("TEST RESULTS")
print("="*70)

success_criteria = []

if speedup > 2.0:
    print(f"  [OK] Speedup {speedup:.1f}x > 2.0x (target achieved)")
    success_criteria.append(True)
else:
    print(f"  [WARNING] Speedup {speedup:.1f}x < 2.0x (expected >2x)")
    success_criteria.append(False)

if errors and np.mean(errors) < 1.0:
    print(f"  [OK] Mean accuracy error {np.mean(errors):.6f}% < 1.0%")
    success_criteria.append(True)
else:
    print(f"  [WARNING] Mean accuracy error {np.mean(errors):.6f}% >= 1.0%")
    success_criteria.append(False)

if all(success_criteria):
    print("\n[SUCCESS] PrepareCache() provides significant speedup with good accuracy!")
    sys.exit(0)
else:
    print("\n[INFO] Some criteria not met - review results above")
    sys.exit(1)
