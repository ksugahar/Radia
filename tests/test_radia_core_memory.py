#!/usr/bin/env python
"""
Test memory behavior of direct rad.Fld() calls (without radia_ngsolve)

This test checks if memory accumulation occurs in Radia core itself,
not in the radia_ngsolve wrapper.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

import gc
import psutil
import radia as rad

print("=" * 80)
print("Radia Core Memory Test (Direct rad.Fld() calls)")
print("=" * 80)
print()

# Create simple Radia geometry
rad.FldUnits('m')
magnet = rad.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 1.0])

# Generate 20 observation points (same as NGSolve mesh test)
points = []
for i in range(20):
    x = 0.02 + i * 0.001
    y = 0.01
    z = 0.015
    points.append([x, y, z])

print(f"Number of observation points: {len(points)}")
print()

# Measure memory usage
process = psutil.Process()

def get_memory_mb():
    """Get current memory usage in MB"""
    return process.memory_info().rss / 1024 / 1024

print("Testing memory with repeated rad.Fld() calls...")
print()
print("Step    Memory (MB)   Delta (MB)   Notes")
print("-" * 60)

mem_start = get_memory_mb()
print(f"0       {mem_start:8.2f}       --         Initial")

# Perform 100 iterations, calling rad.Fld() for each point
num_iterations = 100
mem_values = [mem_start]

for i in range(1, num_iterations + 1):
    # Force garbage collection before each iteration
    gc.collect()

    # Call rad.Fld() for all points (same pattern as GridFunction.Set())
    for pt in points:
        B = rad.Fld(magnet, 'b', pt)  # Direct Radia call

    # Measure memory after all calls
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
print(f"  Observation points:   {len(points)}")
print(f"  Iterations:           {num_iterations}")
print(f"  Total Fld() calls:    {num_iterations * len(points)}")
print(f"  Initial memory:       {mem_start:.2f} MB")
print(f"  Final memory:         {mem_end:.2f} MB")
print(f"  Total increase:       {total_leak:.2f} MB")
print(f"  Avg per iteration:    {avg_leak_per_iteration*1024:.1f} KB")
print()

# Per-call estimate
if num_iterations * len(points) > 0:
    avg_leak_per_call = (total_leak * 1024 * 1024) / (num_iterations * len(points))
    print(f"  Estimated per Fld():  {avg_leak_per_call:.1f} bytes")
    print()

# Verdict
if total_leak > 10.0:  # More than 10 MB leak
    print(f"[FAIL] Significant memory leak in Radia core: {total_leak:.2f} MB")
    print(f"       This suggests the leak is in rad.Fld(), not radia_ngsolve")
elif total_leak > 1.0:  # 1-10 MB
    print(f"[WARNING] Possible memory leak in Radia core: {total_leak:.2f} MB")
    print(f"          May be normal GC behavior or small leak")
else:
    print(f"[PASS] No significant memory leak: {total_leak:.2f} MB")
    print(f"       Radia core memory management appears correct")

print("=" * 60)
