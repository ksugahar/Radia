#!/usr/bin/env python
"""
Test memory leak in Radia-only workflow (without NGSolve)

This test checks if memory leak occurs when repeatedly creating and deleting
Radia objects without involving NGSolve at all.

If this test shows NO leak, then the leak is in NGSolve integration.
If this test DOES show a leak, then the leak is in Radia core.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

import gc
import tracemalloc
import numpy as np
import radia as rad

print("=" * 80)
print("Radia-Only Memory Leak Test (No NGSolve)")
print("=" * 80)
print()

# Configuration
NUM_STEPS = 100
MOVE_DISTANCE = 0.001
NUM_FIELD_POINTS = 100  # Number of field evaluation points

print(f"Test Configuration:")
print(f"  Number of time steps: {NUM_STEPS}")
print(f"  Move distance per step: {MOVE_DISTANCE} m")
print(f"  Field evaluation points per step: {NUM_FIELD_POINTS}")
print()

# Set Radia to use meters
rad.FldUnits('m')

# Start memory tracking
print("[Test] Starting memory leak test (Radia only)...")
print()
tracemalloc.start()
gc.collect()

# Record initial memory
mem_start = tracemalloc.get_traced_memory()[0]
print(f"Initial memory: {mem_start / 1024 / 1024:.2f} MB")
print()

# Store memory usage at checkpoints
memory_samples = []
positions = []

# Define field observation points (reused across steps)
obs_points = []
for i in range(NUM_FIELD_POINTS):
    x = -0.05 + (i / NUM_FIELD_POINTS) * 0.1
    obs_points.append([x, 0.0, 0.0])

print("Time step progress:")
print("-" * 80)

magnet_size = [0.01, 0.01, 0.01]
magnetization = [0, 0, 1.0]
initial_pos = [0.0, 0.0, 0.0]

for step in range(NUM_STEPS):
    # Move magnet along x-axis
    x_pos = initial_pos[0] + step * MOVE_DISTANCE
    new_pos = [x_pos, initial_pos[1], initial_pos[2]]

    # Delete all Radia objects
    rad.UtiDelAll()

    # Create new magnet at new position
    magnet = rad.ObjRecMag(new_pos, magnet_size, magnetization)

    # Evaluate field at multiple points using Radia's Fld() function
    # This is the equivalent of what NGSolve would do internally
    B_values = []
    for pt in obs_points:
        B = rad.Fld(magnet, 'b', pt)
        B_values.append(B)

    # Record memory every 10 steps
    if step % 10 == 0:
        gc.collect()
        mem_current = tracemalloc.get_traced_memory()[0]
        mem_mb = mem_current / 1024 / 1024
        memory_samples.append(mem_mb)
        positions.append(x_pos)

        # Sample one B field value
        B_sample = B_values[0] if B_values else [0, 0, 0]

        print(f"  Step {step:3d}: x={x_pos:.4f} m, Memory={mem_mb:.2f} MB, "
              f"B_z={B_sample[2]:.6e} T")

print("-" * 80)
print()

# Final garbage collection
gc.collect()

# Get final memory
mem_end = tracemalloc.get_traced_memory()[0]
mem_peak = tracemalloc.get_traced_memory()[1]

print("Memory Statistics:")
print("-" * 80)
print(f"  Initial memory:  {mem_start / 1024 / 1024:.2f} MB")
print(f"  Final memory:    {mem_end / 1024 / 1024:.2f} MB")
print(f"  Peak memory:     {mem_peak / 1024 / 1024:.2f} MB")
print(f"  Memory increase: {(mem_end - mem_start) / 1024 / 1024:.2f} MB")
print()

# Analyze memory growth
if len(memory_samples) >= 2:
    x = np.arange(len(memory_samples))
    y = np.array(memory_samples)

    # Fit line: y = a*x + b
    coeffs = np.polyfit(x, y, 1)
    slope = coeffs[0]

    # Memory growth per step
    growth_per_step = slope / 10  # MB per step

    print("Memory Growth Analysis:")
    print("-" * 80)
    print(f"  Samples collected:     {len(memory_samples)}")
    print(f"  Memory growth trend:   {slope:.4f} MB per 10 steps")
    print(f"  Growth per step:       {growth_per_step * 1024:.2f} KB/step")
    print(f"  Projected for 1000 steps: {growth_per_step * 1000:.2f} MB")
    print()

# Stop memory tracking
tracemalloc.stop()

# Memory leak detection criteria
print("=" * 80)
print("MEMORY LEAK ASSESSMENT (Radia Only)")
print("=" * 80)
print()

# Criteria 1: Total memory increase
total_increase_mb = (mem_end - mem_start) / 1024 / 1024
THRESHOLD_TOTAL_MB = 10  # Radia-only should not increase more than 10 MB

print(f"Criterion 1: Total memory increase")
print(f"  Increase: {total_increase_mb:.2f} MB")
print(f"  Threshold: {THRESHOLD_TOTAL_MB} MB")
if total_increase_mb < THRESHOLD_TOTAL_MB:
    print(f"  Result: [PASS] Memory increase within acceptable range")
    criterion1_pass = True
else:
    print(f"  Result: [FAIL] Excessive memory increase in Radia core!")
    criterion1_pass = False
print()

# Criteria 2: Memory growth rate
if len(memory_samples) >= 2:
    THRESHOLD_GROWTH_KB = 5  # Should not grow more than 5 KB per step

    print(f"Criterion 2: Memory growth rate")
    print(f"  Growth rate: {growth_per_step * 1024:.2f} KB/step")
    print(f"  Threshold: {THRESHOLD_GROWTH_KB} KB/step")

    if growth_per_step * 1024 < THRESHOLD_GROWTH_KB:
        print(f"  Result: [PASS] No significant linear growth in Radia")
        criterion2_pass = True
    else:
        print(f"  Result: [FAIL] Linear memory growth in Radia core!")
        criterion2_pass = False
    print()
else:
    print(f"Criterion 2: Insufficient samples for growth analysis")
    criterion2_pass = True
    print()

# Overall result
print("=" * 80)
if criterion1_pass and criterion2_pass:
    print("OVERALL: [PASS] No Radia core memory leak detected")
    print()
    print("Conclusion: Memory leak is in NGSolve integration, NOT in Radia core.")
    print("The leak occurs in rad_ngsolve.cpp or NGSolve's GridFunction.Set().")
else:
    print("OVERALL: [FAIL] Memory leak detected in Radia core!")
    print()
    print("Conclusion: Memory leak is in Radia core (UtiDelAll() or field computation).")
    print("Check rad_transform_impl.cpp DeleteAllElements() implementation.")

print("=" * 80)
print()

# Cleanup
rad.UtiDelAll()

# Exit with appropriate code
if criterion1_pass and criterion2_pass:
    sys.exit(0)
else:
    sys.exit(1)
