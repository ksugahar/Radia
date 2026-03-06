#!/usr/bin/env python
"""
Test memory leak in radia_ngsolve with PrepareCache enabled

Tests whether PrepareCache reduces the memory leak by using batch evaluation
instead of individual Evaluate() calls.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

import gc
import psutil
import numpy as np
import radia as rad
from ngsolve import *
import radia_ngsolve

print("=" * 80)
print("radia_ngsolve Memory Test with PrepareCache")
print("=" * 80)
print()

# Create simple Radia geometry
rad.FldUnits('m')
magnet = rad.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 1.0])

# Create small NGSolve mesh (100 vertices)
from netgen.occ import *
box = Box((0.02, 0, 0), (0.05, 0.03, 0.03))
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=0.015))

print(f"Mesh: {mesh.nv} vertices, {mesh.ne} elements")
print()

# Create field object
B_cf = radia_ngsolve.RadiaField(magnet, 'b', units='m')

# **Enable PrepareCache** - this should use batch evaluation
print("[INFO] Preparing cache (batch evaluation)...")
# Extract mesh vertices as list of coordinates
points = []
for v in mesh.vertices:
    pt = mesh[v].point
    points.append([pt[0], pt[1], pt[2]])

B_cf.PrepareCache(points)
print(f"[INFO] Cache prepared with {len(points)} points")
print()

# Create GridFunction
fes = HCurl(mesh, order=1)
gf = GridFunction(fes)

print(f"FE Space: {fes.ndof} DOFs")
print()

# Measure memory usage
process = psutil.Process()

def get_memory_mb():
    """Get current memory usage in MB"""
    return process.memory_info().rss / 1024 / 1024

print("Testing memory leak with PrepareCache enabled...")
print()
print("Step    Memory (MB)   Delta (MB)   Notes")
print("-" * 60)

mem_start = get_memory_mb()
print(f"0       {mem_start:8.2f}       --         Initial")

# Perform 100 repeated Set() calls WITH PrepareCache
num_iterations = 100
mem_values = [mem_start]

for i in range(1, num_iterations + 1):
    # Force garbage collection before each iteration
    gc.collect()

    # Call GridFunction.Set() - should use cached values
    gf.Set(B_cf)

    # Measure memory after Set()
    mem_current = get_memory_mb()
    mem_delta = mem_current - mem_values[-1]
    mem_values.append(mem_current)

    # Print every 10 iterations
    if i % 10 == 0:
        total_increase = mem_current - mem_start
        print(f"{i:3d}     {mem_current:8.2f}    {mem_delta:+8.2f}    Total: +{total_increase:.2f} MB")

print()

# Final statistics
mem_end = get_memory_mb()
total_leak = mem_end - mem_start
avg_leak_per_iteration = total_leak / num_iterations

print("=" * 60)
print("Summary:")
print(f"  Mesh vertices:        {mesh.nv}")
print(f"  Iterations:           {num_iterations}")
print(f"  PrepareCache:         Enabled (batch evaluation)")
print(f"  Initial memory:       {mem_start:.2f} MB")
print(f"  Final memory:         {mem_end:.2f} MB")
print(f"  Total increase:       {total_leak:.2f} MB")
print(f"  Avg per iteration:    {avg_leak_per_iteration*1024:.1f} KB")
print()

# Compare with non-cached version
print("Comparison:")
print(f"  WITHOUT PrepareCache: ~2.6 MB leak (26.7 KB/iteration)")
print(f"  WITH PrepareCache:    {total_leak:.2f} MB leak ({avg_leak_per_iteration*1024:.1f} KB/iteration)")
if total_leak < 1.5:
    improvement = ((2.6 - total_leak) / 2.6) * 100
    print(f"  Improvement:          {improvement:.1f}% reduction in memory leak!")
print()

# Verdict
if total_leak > 10.0:  # More than 10 MB leak
    print(f"[FAIL] Significant memory leak: {total_leak:.2f} MB")
    print(f"       PrepareCache did not fix the leak")
elif total_leak > 1.0:  # 1-10 MB
    print(f"[WARNING] Possible memory leak: {total_leak:.2f} MB")
    print(f"          PrepareCache helped but leak remains")
else:
    print(f"[PASS] No significant memory leak: {total_leak:.2f} MB")
    print(f"       PrepareCache successfully prevents memory accumulation!")

print("=" * 60)
