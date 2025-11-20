#!/usr/bin/env python
"""
Test Radia memory usage using process RSS (Resident Set Size)
instead of Python's tracemalloc, to catch C++ memory leaks.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

import gc
import psutil
import radia as rad

process = psutil.Process()
rad.FldUnits('m')

print("=" * 80)
print("Process Memory Test (RSS)")
print("=" * 80)
print()

def get_memory_mb():
    """Get current process memory in MB"""
    return process.memory_info().rss / 1024 / 1024

# Warm-up
rad.UtiDelAll()
magnet = rad.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 1.0])
B = rad.Fld(magnet, 'b', [0.02, 0, 0])
rad.UtiDelAll()
gc.collect()

# Baseline measurement
mem_start = get_memory_mb()
print(f"Starting memory: {mem_start:.2f} MB")
print()

# Run test with 100 field evaluations
NUM_STEPS = 50
NUM_POINTS = 100

for step in range(NUM_STEPS):
    x_pos = step * 0.001
    new_pos = [x_pos, 0.0, 0.0]

    rad.UtiDelAll()
    magnet = rad.ObjRecMag(new_pos, [0.01, 0.01, 0.01], [0, 0, 1.0])

    # Compute fields for 100 points
    for i in range(NUM_POINTS):
        pt = [0.02 + i * 0.001, 0.0, 0.0]
        B = rad.Fld(magnet, 'b', pt)

    if step % 10 == 9:
        gc.collect()
        mem = get_memory_mb()
        print(f"Step {step+1:3d}: {mem:.2f} MB (growth: {mem - mem_start:.2f} MB)")

rad.UtiDelAll()
gc.collect()

mem_end = get_memory_mb()
growth = mem_end - mem_start

print()
print("=" * 80)
print(f"Total memory growth: {growth:.2f} MB")
print(f"Per step: {growth / NUM_STEPS * 1024:.2f} KB/step")
print(f"Per field evaluation: {growth / NUM_STEPS / NUM_POINTS * 1024:.2f} KB/point")
print("=" * 80)

if growth < 1.0:
    print("\n[PASS] Memory usage is stable (< 1 MB growth)")
else:
    print(f"\n[FAIL] Memory leak detected ({growth:.2f} MB growth)")
