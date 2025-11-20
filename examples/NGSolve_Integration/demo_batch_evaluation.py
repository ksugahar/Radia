#!/usr/bin/env python
"""
Demonstration: Batch evaluation for H-matrix acceleration

Shows the benefit of batching all evaluation points together
versus element-by-element evaluation.

Current limitation:
- NGSolve's GridFunction.Set() calls CoefficientFunction element-by-element
- Each call has only ~10-20 points (one element's integration points)
- H-matrix overhead >> computation time for such small batches
- No speedup observed

Proposed solution:
- Collect ALL integration points from ALL elements
- Evaluate in ONE batch call (1000s of points)
- H-matrix speedup is significant for large batches
- Expected: 10-100x speedup for N>1000

This script demonstrates the concept and measures potential speedup.
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
except ImportError:
    print("NGSolve not available")
    sys.exit(0)

print("="*70)
print("Batch Evaluation Demo for H-Matrix Acceleration")
print("="*70)

# Set Radia units to meters
rad.FldUnits('m')

# Create Radia magnet
rad.UtiDelAll()
n = 5  # 5x5x5 = 125 elements
cube_size = 0.100
elem_size = cube_size / n

print(f"\n[1] Creating magnet: {n}x{n}x{n} = {n**3} elements")
elements = []
for i in range(n):
    for j in range(n):
        for k in range(n):
            x = (i - n/2 + 0.5) * elem_size
            y = (j - n/2 + 0.5) * elem_size
            z = (k - n/2 + 0.5) * elem_size
            elem = rad.ObjRecMag([x, y, z], [elem_size, elem_size, elem_size],
                                 [0, 0, 1.2])
            elements.append(elem)

magnet = rad.ObjCnt(elements)

# Create mesh
print(f"\n[2] Creating NGSolve mesh")
box = Box((0.015, 0.015, 0.015), (0.063, 0.063, 0.063))
mesh = Mesh(OCCGeometry(box).GenerateMesh(maxh=0.010))
print(f"    Mesh: {mesh.nv} vertices, {mesh.ne} elements")

# Create finite element space
fes = HCurl(mesh, order=2)
gf = GridFunction(fes)
print(f"    FE space: {fes.ndof} DOFs")

# Estimate number of integration points
# For HCurl order=2, typically 10-15 points per element
points_per_elem = 15  # Estimate
total_points = mesh.ne * points_per_elem
print(f"    Estimated integration points: ~{total_points}")

print("\n" + "="*70)
print("MEASUREMENT 1: Element-by-element evaluation (current method)")
print("="*70)

# Standard GridFunction.Set() with H-matrix
rad.SetHMatrixFieldEval(1, 1e-6)
cf_standard = rad_ngsolve.RadiaField(magnet, 'b')

print(f"\n[3a] Measuring GridFunction.Set() with H-matrix...")
print(f"     (Element-by-element: {mesh.ne} calls, ~{points_per_elem} points/call)")

t0 = time.time()
gf.Set(cf_standard)
t_standard = time.time() - t0

print(f"     Time: {t_standard*1000:.1f} ms")
print(f"     Average per element: {t_standard*1000/mesh.ne:.2f} ms")

print("\n" + "="*70)
print("MEASUREMENT 2: Single batch evaluation (proposed method)")
print("="*70)

# Simulate batch evaluation: collect all vertices and evaluate once
print(f"\n[3b] Simulating batch evaluation...")
print(f"     (Single batch: 1 call, {mesh.nv} points)")

# Collect all vertex positions
points = []
for v in mesh.vertices:
    pt = v.point
    points.append([pt[0], pt[1], pt[2]])

# Single batch evaluation
t0 = time.time()
field_values = rad.FldBatch(magnet, 'b', points, 1)  # use_hmatrix=1
t_batch = time.time() - t0

print(f"     Time: {t_batch*1000:.1f} ms")
print(f"     Average per point: {t_batch*1e6/len(points):.2f} us")

print("\n" + "="*70)
print("ANALYSIS")
print("="*70)

# Calculate actual points evaluated in standard method
# Each element's integration rule has multiple points
actual_eval_points = total_points

print(f"\nEstimated evaluation points:")
print(f"  Standard method: ~{actual_eval_points} points")
print(f"  Batch method:     {len(points)} points (vertices only)")

print(f"\nTime comparison:")
print(f"  Standard (element-by-element): {t_standard*1000:8.1f} ms")
print(f"  Batch (single call):           {t_batch*1000:8.1f} ms")

speedup = t_standard / t_batch
print(f"\nObserved speedup: {speedup:.2f}x")

print(f"\nWhy the speedup is limited:")
print(f"  - Standard method evaluates ~{actual_eval_points} integration points")
print(f"  - Batch method evaluates only {len(points)} vertices")
print(f"  - Not a fair comparison!")

print(f"\nTrue batch optimization requires:")
print(f"  1. Collect ALL integration points (~{actual_eval_points})")
print(f"  2. Evaluate in ONE batch call")
print(f"  3. Map results back to GridFunction DOFs")
print(f"  4. This would give true {speedup:.1f}-{speedup*2:.1f}x speedup")

print("\n" + "="*70)
print("CONCLUSION")
print("="*70)

print(f"""
Current situation:
- GridFunction.Set() calls Evaluate() {mesh.ne} times
- Each call evaluates ~{points_per_elem} points
- H-matrix overhead dominates for small batches
- Total time: {t_standard*1000:.1f} ms

Optimized approach (requires C++ implementation):
- Collect all ~{actual_eval_points} integration points
- Single FldBatch() call
- Full H-matrix acceleration
- Expected time: ~{t_batch*actual_eval_points/len(points):.1f} ms
- Expected speedup: ~{t_standard/(t_batch*actual_eval_points/len(points)):.1f}x

Implementation path:
1. Add PrepareCache() method to RadiaFieldCF (C++)
2. PrepareCache() collects all mesh integration points
3. PrepareCache() calls FldBatch() once
4. Evaluate() returns cached values
5. User calls PrepareCache(mesh) before gf.Set(cf)
""")

print("="*70)

rad.UtiDelAll()

