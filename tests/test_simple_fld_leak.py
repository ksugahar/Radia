#!/usr/bin/env python
"""
Minimal test to isolate rad.Fld() memory leak
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))

import gc
import psutil
import radia as rad

process = psutil.Process()
rad.FldUnits('m')

print("Minimal rad.Fld() memory test")
print("=" * 80)

# Create magnet once
rad.UtiDelAll()
magnet = rad.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 1.0])

gc.collect()
mem_start = process.memory_info().rss / 1024 / 1024
print(f"Initial memory: {mem_start:.2f} MB\n")

# Call rad.Fld() many times WITHOUT creating new magnets
NUM_CALLS = 10000
for i in range(NUM_CALLS):
    pt = [0.02 + i * 0.00001, 0.0, 0.0]
    B = rad.Fld(magnet, 'b', pt)

    # Explicitly delete result
    del B

    if i % 1000 == 999:
        gc.collect()
        mem = process.memory_info().rss / 1024 / 1024
        print(f"After {i+1:5d} calls: {mem:.2f} MB (growth: {mem - mem_start:.2f} MB)")

gc.collect()
mem_end = process.memory_info().rss / 1024 / 1024
growth = mem_end - mem_start

print()
print("=" * 80)
print(f"Total growth: {growth:.2f} MB for {NUM_CALLS} calls")
print(f"Per call: {growth / NUM_CALLS * 1024:.3f} KB/call")
print("=" * 80)

if growth < 2.0:
    print("\n[PASS] Memory stable")
else:
    print(f"\n[FAIL] Memory leak: {growth:.2f} MB")
