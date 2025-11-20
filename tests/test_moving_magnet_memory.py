#!/usr/bin/env python
"""
Test for memory leaks when repeatedly updating moving magnet fields in NGSolve

This test simulates a time-stepping scenario where a magnet moves and its field
is repeatedly evaluated in NGSolve. This was reported to cause memory leaks in
previous versions.

Author: Radia development team
Date: 2025-11-21
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

import gc
import tracemalloc
import numpy as np

# Check if NGSolve is available
try:
    from ngsolve import *
    from netgen.occ import *
    import radia_ngsolve
    NGSOLVE_AVAILABLE = True
except ImportError:
    print("ERROR: NGSolve not available. This test requires NGSolve.")
    NGSOLVE_AVAILABLE = False
    sys.exit(1)

import radia as rad

print("=" * 80)
print("Memory Leak Test: Moving Magnet with NGSolve Integration")
print("=" * 80)
print()

# Configuration
NUM_STEPS = 100  # Number of time steps to simulate
MOVE_DISTANCE = 0.001  # Distance to move per step (meters)

print(f"Test Configuration:")
print(f"  Number of time steps: {NUM_STEPS}")
print(f"  Move distance per step: {MOVE_DISTANCE} m")
print()

# Set Radia to use meters
rad.FldUnits('m')

# Create NGSolve mesh (reused across all steps)
print("[Setup] Creating NGSolve mesh...")
box = Box((-0.05, -0.05, -0.05), (0.05, 0.05, 0.05))
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=0.02))
print(f"  Mesh: {mesh.nv} vertices, {mesh.ne} elements")
print()

# Create finite element space (reused)
fes = HCurl(mesh, order=2)
print(f"  FE Space: {fes.ndof} DOFs")
print()

# Initial magnet position
print("[Setup] Creating magnet...")
rad.UtiDelAll()
magnet_size = [0.01, 0.01, 0.01]  # 10mm cube
magnetization = [0, 0, 1.0]  # 1 T in z-direction
initial_pos = [0.0, 0.0, 0.0]

magnet = rad.ObjRecMag(initial_pos, magnet_size, magnetization)
print(f"  Magnet: {magnet_size} m at {initial_pos}")
print(f"  Magnetization: {magnetization} T")
print()

# Start memory tracking
print("[Test] Starting memory leak test...")
print()
tracemalloc.start()
gc.collect()  # Clean up before starting

# Record initial memory
mem_start = tracemalloc.get_traced_memory()[0]
print(f"Initial memory: {mem_start / 1024 / 1024:.2f} MB")
print()

# Store memory usage at checkpoints
memory_samples = []
positions = []

print("Time step progress:")
print("-" * 80)

for step in range(NUM_STEPS):
    # Move magnet along x-axis
    x_pos = initial_pos[0] + step * MOVE_DISTANCE
    new_pos = [x_pos, initial_pos[1], initial_pos[2]]

    # Update magnet position using transformation
    # Method 1: Delete and recreate (simpler but potentially leaky)
    rad.UtiDelAll()
    magnet = rad.ObjRecMag(new_pos, magnet_size, magnetization)

    # Create RadiaField CoefficientFunction for B field
    # This is where memory leak might occur - creating new CF each time
    B_cf = radia_ngsolve.RadiaField(magnet, 'b')

    # Create GridFunction and set field
    # This allocates memory and should be released
    gf_B = GridFunction(fes)
    gf_B.Set(B_cf)

    # Evaluate field at a test point (forces computation)
    test_point = mesh(0.02, 0.0, 0.0)
    B_val = gf_B(test_point)

    # Record memory every 10 steps
    if step % 10 == 0:
        gc.collect()  # Force garbage collection
        mem_current = tracemalloc.get_traced_memory()[0]
        mem_mb = mem_current / 1024 / 1024
        memory_samples.append(mem_mb)
        positions.append(x_pos)

        print(f"  Step {step:3d}: x={x_pos:.4f} m, Memory={mem_mb:.2f} MB, "
              f"B_z={B_val[2]:.6e} T")

    # Clean up explicitly (test if this helps)
    del gf_B
    del B_cf

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
    # Linear regression to detect memory growth trend
    x = np.arange(len(memory_samples))
    y = np.array(memory_samples)

    # Fit line: y = a*x + b
    coeffs = np.polyfit(x, y, 1)
    slope = coeffs[0]  # MB per checkpoint (10 steps)

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
print("MEMORY LEAK ASSESSMENT")
print("=" * 80)
print()

# Criteria 1: Total memory increase
total_increase_mb = (mem_end - mem_start) / 1024 / 1024
THRESHOLD_TOTAL_MB = 50  # Should not increase more than 50 MB for 100 steps

print(f"Criterion 1: Total memory increase")
print(f"  Increase: {total_increase_mb:.2f} MB")
print(f"  Threshold: {THRESHOLD_TOTAL_MB} MB")
if total_increase_mb < THRESHOLD_TOTAL_MB:
    print(f"  Result: [PASS] Memory increase within acceptable range")
    criterion1_pass = True
else:
    print(f"  Result: [FAIL] Excessive memory increase detected!")
    criterion1_pass = False
print()

# Criteria 2: Memory growth rate
if len(memory_samples) >= 2:
    THRESHOLD_GROWTH_KB = 10  # Should not grow more than 10 KB per step

    print(f"Criterion 2: Memory growth rate")
    print(f"  Growth rate: {growth_per_step * 1024:.2f} KB/step")
    print(f"  Threshold: {THRESHOLD_GROWTH_KB} KB/step")

    if growth_per_step * 1024 < THRESHOLD_GROWTH_KB:
        print(f"  Result: [PASS] No significant linear growth")
        criterion2_pass = True
    else:
        print(f"  Result: [FAIL] Linear memory growth detected!")
        criterion2_pass = False
    print()
else:
    print(f"Criterion 2: Insufficient samples for growth analysis")
    criterion2_pass = True
    print()

# Overall result
print("=" * 80)
if criterion1_pass and criterion2_pass:
    print("OVERALL: [PASS] No memory leak detected")
    print()
    print("The moving magnet simulation does not show signs of memory leaks.")
    print("Memory usage remains stable across time steps.")
else:
    print("OVERALL: [FAIL] Potential memory leak detected!")
    print()
    print("WARNING: Memory usage increases significantly with each time step.")
    print("This indicates a memory leak in the RadiaField/GridFunction workflow.")
    print()
    print("Possible causes:")
    print("  1. RadiaField CoefficientFunction not properly released")
    print("  2. GridFunction Set() operation leaks internal structures")
    print("  3. Radia objects not properly freed in rad.UtiDelAll()")
    print("  4. pybind11 Python/C++ boundary issues")
    print()
    print("Recommended actions:")
    print("  - Review radia_ngsolve.cpp for proper pybind11 lifetime management")
    print("  - Check if PrepareCache() allocations are properly freed")
    print("  - Verify Radia's UtiDelAll() clears all internal caches")

print("=" * 80)
print()

# Cleanup
rad.UtiDelAll()

# Exit with appropriate code
if criterion1_pass and criterion2_pass:
    sys.exit(0)
else:
    sys.exit(1)
