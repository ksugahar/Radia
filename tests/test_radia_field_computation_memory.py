#!/usr/bin/env python
"""
Test to isolate exactly which Radia operation causes memory leak

Tests three scenarios:
1. Only create/delete magnets (no field computation)
2. Create magnets + compute fields
3. Create magnets + compute fields with many points

This will identify if leak is in object creation or field computation.
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
print("Radia Memory Leak Isolation Test")
print("=" * 80)
print()

rad.FldUnits('m')
NUM_STEPS = 50

def run_test(test_name, compute_fields=False, num_points=1):
    """Run memory test with different operations"""
    print(f"[{test_name}]")
    print("-" * 80)

    tracemalloc.start()
    gc.collect()
    mem_start = tracemalloc.get_traced_memory()[0]

    magnet_size = [0.01, 0.01, 0.01]
    magnetization = [0, 0, 1.0]

    memory_samples = []

    for step in range(NUM_STEPS):
        x_pos = step * 0.001
        new_pos = [x_pos, 0.0, 0.0]

        # Always create/delete magnet
        rad.UtiDelAll()
        magnet = rad.ObjRecMag(new_pos, magnet_size, magnetization)

        # Optionally compute fields
        if compute_fields:
            for i in range(num_points):
                pt = [0.02 + i * 0.001, 0.0, 0.0]
                B = rad.Fld(magnet, 'b', pt)

        if step % 10 == 0:
            gc.collect()
            mem_current = tracemalloc.get_traced_memory()[0]
            memory_samples.append(mem_current / 1024 / 1024)

    gc.collect()
    mem_end = tracemalloc.get_traced_memory()[0]
    tracemalloc.stop()

    mem_increase = (mem_end - mem_start) / 1024 / 1024

    # Calculate growth rate
    if len(memory_samples) >= 2:
        x = np.arange(len(memory_samples))
        y = np.array(memory_samples)
        slope = np.polyfit(x, y, 1)[0]
        growth_per_step = slope / 10 * 1024  # KB per step
    else:
        growth_per_step = 0

    print(f"  Memory increase: {mem_increase:.2f} MB")
    print(f"  Growth rate: {growth_per_step:.2f} KB/step")
    print()

    rad.UtiDelAll()
    return mem_increase, growth_per_step

# Test scenarios
print("Testing to isolate memory leak source:")
print("=" * 80)
print()

results = {}

# Test 1: Only object creation/deletion (no field computation)
results['create_only'] = run_test(
    "Test 1: Create/Delete magnet only (no field computation)",
    compute_fields=False
)

# Test 2: Object creation + 1 field evaluation per step
results['one_field'] = run_test(
    "Test 2: Create/Delete + 1 field evaluation",
    compute_fields=True,
    num_points=1
)

# Test 3: Object creation + 10 field evaluations per step
results['ten_fields'] = run_test(
    "Test 3: Create/Delete + 10 field evaluations",
    compute_fields=True,
    num_points=10
)

# Test 4: Object creation + 100 field evaluations per step
results['hundred_fields'] = run_test(
    "Test 4: Create/Delete + 100 field evaluations",
    compute_fields=True,
    num_points=100
)

# Summary
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()

print(f"{'Test':<40} {'Mem (MB)':<12} {'Rate (KB/step)':<15}")
print("-" * 80)

for test_name, (mem_increase, growth_rate) in results.items():
    print(f"{test_name:<40} {mem_increase:<12.2f} {growth_rate:<15.2f}")

print()
print("=" * 80)
print("DIAGNOSIS")
print("=" * 80)
print()

create_only_rate = results['create_only'][1]
one_field_rate = results['one_field'][1]
ten_fields_rate = results['ten_fields'][1]
hundred_fields_rate = results['hundred_fields'][1]

print(f"Object creation/deletion leak: {create_only_rate:.2f} KB/step")
print(f"Additional leak per field eval (1 pt): {one_field_rate - create_only_rate:.2f} KB/step")
print(f"Additional leak per field eval (10 pt): {(ten_fields_rate - create_only_rate) / 10:.2f} KB/point")
print(f"Additional leak per field eval (100 pt): {(hundred_fields_rate - create_only_rate) / 100:.2f} KB/point")
print()

if create_only_rate > 1.0:
    print("Primary leak: Object creation/deletion (rad.ObjRecMag / rad.UtiDelAll)")
    print("Location: rad_transform_impl.cpp DeleteAllElements() or object constructors")
elif one_field_rate > create_only_rate + 1.0:
    print("Primary leak: Field computation (rad.Fld)")
    print("Location: Field computation functions in radapl4.cpp")
else:
    print("Memory usage appears stable.")

print("=" * 80)
