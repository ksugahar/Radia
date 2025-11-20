#!/usr/bin/env python
"""
Diagnostic test to identify the exact source of memory leak in moving magnet scenario

This test tries different approaches to isolate the leak:
1. Reusing RadiaField CF vs creating new one each time
2. Reusing GridFunction vs creating new one each time
3. Explicitly calling ClearCache()
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

import gc
import tracemalloc
import numpy as np

try:
    from ngsolve import *
    from netgen.occ import *
    import radia_ngsolve
except ImportError:
    print("ERROR: NGSolve not available")
    sys.exit(1)

import radia as rad

print("=" * 80)
print("Memory Leak Diagnostic Test")
print("=" * 80)
print()

# Configuration
NUM_STEPS = 50
rad.FldUnits('m')

# Create mesh (reused)
box = Box((-0.05, -0.05, -0.05), (0.05, 0.05, 0.05))
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=0.02))
fes = HCurl(mesh, order=2)

print(f"Mesh: {mesh.nv} vertices, FE Space: {fes.ndof} DOFs")
print()

def run_test(test_name, reuse_cf=False, reuse_gf=False, clear_cache=False):
    """Run memory test with different strategies"""
    print(f"[{test_name}]")
    print("-" * 80)

    tracemalloc.start()
    gc.collect()
    mem_start = tracemalloc.get_traced_memory()[0]

    # Initial objects
    rad.UtiDelAll()
    magnet = rad.ObjRecMag([0, 0, 0], [0.01, 0.01, 0.01], [0, 0, 1.0])

    B_cf = None
    gf_B = None

    if reuse_cf:
        B_cf = radia_ngsolve.RadiaField(magnet, 'b')
    if reuse_gf:
        gf_B = GridFunction(fes)

    memory_samples = []

    for step in range(NUM_STEPS):
        # Move magnet
        x_pos = step * 0.001
        rad.UtiDelAll()
        magnet = rad.ObjRecMag([x_pos, 0, 0], [0.01, 0.01, 0.01], [0, 0, 1.0])

        # Create or reuse CF
        if not reuse_cf:
            B_cf = radia_ngsolve.RadiaField(magnet, 'b')

        # Create or reuse GridFunction
        if not reuse_gf:
            gf_B = GridFunction(fes)

        # Set field
        gf_B.Set(B_cf)

        # Evaluate
        test_point = mesh(0.02, 0.0, 0.0)
        B_val = gf_B(test_point)

        # Clear cache if requested
        if clear_cache and hasattr(B_cf, 'ClearCache'):
            B_cf.ClearCache()

        # Cleanup if not reusing
        if not reuse_cf:
            del B_cf
        if not reuse_gf:
            del gf_B

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
    print(f"  Result: {'PASS' if mem_increase < 25 else 'FAIL'}")
    print()

    rad.UtiDelAll()
    return mem_increase, growth_per_step

# Test scenarios
print("Testing different memory management strategies:")
print("=" * 80)
print()

results = {}

# Test 1: Baseline (create new CF and GF each time)
results['baseline'] = run_test(
    "Test 1: Baseline (new CF + new GF)",
    reuse_cf=False,
    reuse_gf=False,
    clear_cache=False
)

# Test 2: Reuse CF
results['reuse_cf'] = run_test(
    "Test 2: Reuse CoefficientFunction",
    reuse_cf=True,
    reuse_gf=False,
    clear_cache=False
)

# Test 3: Reuse GF
results['reuse_gf'] = run_test(
    "Test 3: Reuse GridFunction",
    reuse_cf=False,
    reuse_gf=True,
    clear_cache=False
)

# Test 4: Reuse both
results['reuse_both'] = run_test(
    "Test 4: Reuse CF + GF",
    reuse_cf=True,
    reuse_gf=True,
    clear_cache=False
)

# Test 5: Reuse both + clear cache
results['reuse_clear'] = run_test(
    "Test 5: Reuse CF + GF + ClearCache",
    reuse_cf=True,
    reuse_gf=True,
    clear_cache=True
)

# Summary
print("=" * 80)
print("SUMMARY")
print("=" * 80)
print()

print(f"{'Test':<40} {'Mem (MB)':<12} {'Rate (KB/step)':<15} {'Status':<10}")
print("-" * 80)

for test_name, (mem_increase, growth_rate) in results.items():
    status = "PASS" if mem_increase < 25 else "FAIL"
    print(f"{test_name:<40} {mem_increase:<12.2f} {growth_rate:<15.2f} {status:<10}")

print()
print("=" * 80)
print("DIAGNOSIS")
print("=" * 80)
print()

# Determine root cause
baseline_mem = results['baseline'][0]
reuse_cf_mem = results['reuse_cf'][0]
reuse_gf_mem = results['reuse_gf'][0]
reuse_both_mem = results['reuse_both'][0]

if reuse_cf_mem < baseline_mem * 0.5:
    print("RadiaField CoefficientFunction is a major contributor to memory leak.")
    print("Creating new CF objects accumulates memory.")
elif reuse_gf_mem < baseline_mem * 0.5:
    print("GridFunction is a major contributor to memory leak.")
    print("Creating new GridFunction objects accumulates memory.")
elif reuse_both_mem < baseline_mem * 0.5:
    print("Both CF and GF contribute to memory leak.")
    print("Reusing objects significantly reduces memory growth.")
else:
    print("Memory leak persists even when reusing CF and GF.")
    print("Likely cause: Internal cache or Radia object accumulation.")
    print("Check: PrepareCache() allocations, Radia's UtiDelAll() implementation.")

print()
print("=" * 80)
