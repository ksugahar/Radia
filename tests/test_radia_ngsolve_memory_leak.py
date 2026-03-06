#!/usr/bin/env python
"""
Test memory leak in radia_ngsolve Evaluate() function

Tests whether py::list coords and py::object field_result
cause memory accumulation when GridFunction.Set() is called repeatedly.
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
print("radia_ngsolve Memory Leak Test")
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

print("Testing memory leak with repeated GridFunction.Set() calls...")
print()
print("Step    Memory (MB)   Delta (MB)   Notes")
print("-" * 60)

mem_start = get_memory_mb()
print(f"0       {mem_start:8.2f}       --         Initial")

# Perform 100 repeated Set() calls
num_iterations = 100
mem_values = [mem_start]

for i in range(1, num_iterations + 1):
    # Force garbage collection before each iteration
    gc.collect()

    # Call GridFunction.Set() - this calls Evaluate() for each vertex
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
print(f"  Initial memory:       {mem_start:.2f} MB")
print(f"  Final memory:         {mem_end:.2f} MB")
print(f"  Total increase:       {total_leak:.2f} MB")
print(f"  Avg per iteration:    {avg_leak_per_iteration*1024:.1f} KB")
print()

# Estimate per-vertex leak
if mesh.nv > 0:
    avg_leak_per_vertex = (avg_leak_per_iteration * 1024 * 1024) / mesh.nv
    print(f"  Estimated per-vertex: {avg_leak_per_vertex:.1f} bytes")
    print()

# Theoretical leak estimate
theoretical_per_call = 304  # bytes (152 for coords + 152 for field_result)
theoretical_per_iteration = theoretical_per_call * mesh.nv
theoretical_total = (theoretical_per_iteration * num_iterations) / (1024 * 1024)

print(f"Theoretical leak estimate (if no cleanup):")
print(f"  Per Evaluate() call:  {theoretical_per_call} bytes")
print(f"  Per iteration:        {theoretical_per_iteration/1024:.1f} KB")
print(f"  Total expected:       {theoretical_total:.2f} MB")
print()

# Verdict
if total_leak > 10.0:  # More than 10 MB leak
    print(f"[FAIL] Significant memory leak detected: {total_leak:.2f} MB")
    print(f"       This confirms py::list and py::object accumulation")
elif total_leak > 1.0:  # 1-10 MB
    print(f"[WARNING] Possible memory leak: {total_leak:.2f} MB")
    print(f"          May be normal GC behavior or small leak")
else:
    print(f"[PASS] No significant memory leak: {total_leak:.2f} MB")
    print(f"       pybind11 reference counting appears to work correctly")

print("=" * 60)
