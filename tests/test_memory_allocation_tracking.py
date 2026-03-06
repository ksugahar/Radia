#!/usr/bin/env python
"""
Test to track memory allocation patterns in detail

This test uses tracemalloc's snapshot feature to identify exactly
which lines of code are allocating memory during field computation.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

import gc
import tracemalloc
import radia as rad

print("=" * 80)
print("Memory Allocation Tracking Test")
print("=" * 80)
print()

rad.FldUnits('m')

# Start tracemalloc
tracemalloc.start()
gc.collect()

# Initial snapshot
snapshot1 = tracemalloc.take_snapshot()

# Create magnet
rad.UtiDelAll()
magnet = rad.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 1.0])

# Compute fields for multiple points
print("Computing fields for 10 points...")
for i in range(10):
    pt = [0.02 + i * 0.001, 0.0, 0.0]
    B = rad.Fld(magnet, 'b', pt)

gc.collect()

# Second snapshot
snapshot2 = tracemalloc.take_snapshot()

# Compare snapshots to find what allocated memory
print()
print("Top 20 memory allocations:")
print("-" * 80)

top_stats = snapshot2.compare_to(snapshot1, 'lineno')

for stat in top_stats[:20]:
    print(f"{stat.size / 1024:.2f} KB - {stat.count} allocations")
    print(f"  {stat.traceback.format()[0]}")

print()
print("-" * 80)

# Check total increase
total_increase = sum(stat.size_diff for stat in top_stats if stat.size_diff > 0)
print(f"Total memory increase: {total_increase / 1024:.2f} KB")
print(f"Per field evaluation: {total_increase / 10 / 1024:.2f} KB")

print()
print("=" * 80)

tracemalloc.stop()
rad.UtiDelAll()
