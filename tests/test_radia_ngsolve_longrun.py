#!/usr/bin/env python
"""
Test memory leak with extended run (500 iterations)

Tests whether memory increase is:
- Linear (true leak)
- Saturating (Python memory pool initialization)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

import gc
import psutil
import radia as rad
from ngsolve import *
import radia_ngsolve

print("=" * 80)
print("radia_ngsolve Extended Memory Test (500 iterations)")
print("=" * 80)
print()

# Create simple Radia geometry
rad.FldUnits('m')
magnet = rad.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 1.0])

# Create small NGSolve mesh
from netgen.occ import *
box = Box((0.02, 0, 0), (0.05, 0.03, 0.03))
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=0.015))

print(f"Mesh: {mesh.nv} vertices, {mesh.ne} elements")
print()

# Create field object
B_cf = radia_ngsolve.RadiaField(magnet, 'b', units='m')

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

print("Testing memory with 500 iterations...")
print()
print("Step    Memory (MB)   Delta (MB)   MB/100iter   Notes")
print("-" * 70)

mem_start = get_memory_mb()
print(f"0       {mem_start:8.2f}       --         --         Initial")

# Perform 500 iterations
num_iterations = 500
mem_values = [mem_start]
mem_per_100 = []

for i in range(1, num_iterations + 1):
    gc.collect()
    gf.Set(B_cf)

    mem_current = get_memory_mb()
    mem_delta = mem_current - mem_values[-1]
    mem_values.append(mem_current)

    # Print every 50 iterations
    if i % 50 == 0:
        total_increase = mem_current - mem_start
        mb_per_100 = (total_increase / i) * 100
        mem_per_100.append(mb_per_100)
        print(f"{i:3d}     {mem_current:8.2f}    {mem_delta:+8.2f}    {mb_per_100:8.2f}   Total: +{total_increase:.2f} MB")

print()

# Final statistics
mem_end = get_memory_mb()
total_leak = mem_end - mem_start

print("=" * 70)
print("Summary:")
print(f"  Iterations:           {num_iterations}")
print(f"  Initial memory:       {mem_start:.2f} MB")
print(f"  Final memory:         {mem_end:.2f} MB")
print(f"  Total increase:       {total_leak:.2f} MB")
print()

# Analyze growth pattern
print("Growth pattern analysis:")
for i, (itr, mb) in enumerate(zip(range(50, num_iterations+1, 50), mem_per_100)):
    print(f"  At {itr:3d} iterations: {mb:.2f} MB per 100 iterations")

# Check if growth is slowing down (saturating) or constant (true leak)
if len(mem_per_100) >= 2:
    first_rate = mem_per_100[0]
    last_rate = mem_per_100[-1]
    rate_change = ((last_rate - first_rate) / first_rate) * 100

    print()
    if abs(rate_change) < 10:
        print(f"[LINEAR LEAK] Growth rate stable: {rate_change:+.1f}%")
        print(f"              This indicates a true memory leak")
    elif rate_change < -30:
        print(f"[SATURATING] Growth rate decreasing: {rate_change:+.1f}%")
        print(f"             This indicates Python memory pool initialization")
        print(f"             Not a true leak - memory will stabilize")
    else:
        print(f"[UNCLEAR] Growth rate change: {rate_change:+.1f}%")

print("=" * 70)
