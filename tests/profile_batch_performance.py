#!/usr/bin/env python
"""
Profile batch evaluation performance to identify bottleneck

This script measures performance at different stages:
1. Python list preparation
2. rad.Fld() batch call (Radia-side)
3. Result storage in cache (C++ side)
4. Overall PrepareCache() time
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'python'))

import radia as rad
import numpy as np
import time

print("="*80)
print("PROFILE: BATCH EVALUATION PERFORMANCE")
print("="*80)
print()

# Setup
# Set Radia to use meters for consistency
rad.FldUnits('m')
rad.UtiDelAll()
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.06], [0, 0, 1.2])

print("[Setup] Created magnet (0.04 x 0.04 x 0.06 m):", magnet)
print()

# Test different point counts
test_configs = [
    {"npts": 10, "label": "10 points"},
    {"npts": 50, "label": "50 points"},
    {"npts": 100, "label": "100 points"},
    {"npts": 500, "label": "500 points"},
    {"npts": 1000, "label": "1000 points"},
    {"npts": 2000, "label": "2000 points"},
]

print("Testing Radia.Fld() batch performance:")
print("-" * 80)
print(f"{'Points':<10} {'Prep (ms)':<12} {'Fld (ms)':<12} {'Total (ms)':<12} {'us/point':<12}")
print("-" * 80)

results = []

for config in test_configs:
    npts = config["npts"]
    label = config["label"]

    # Generate test points
    t0 = time.time()
    points_python = []
    for i in range(npts):
        x = 0.02 + (i % 10) * 0.001
        y = 0.02 + ((i // 10) % 10) * 0.001
        z = 0.03 + (i // 100) * 0.001
        points_python.append([x, y, z])  # Already in meters (rad.FldUnits('m'))
    t1 = time.time()
    time_prep = (t1 - t0) * 1000  # ms

    # Call Radia.Fld() in batch
    t0 = time.time()
    try:
        field_results = rad.Fld(magnet, 'a', points_python)
        t1 = time.time()
        time_fld = (t1 - t0) * 1000  # ms
        success = True
    except Exception as e:
        t1 = time.time()
        time_fld = (t1 - t0) * 1000  # ms
        success = False
        error_msg = str(e)

    time_total = time_prep + time_fld
    time_per_point = time_total * 1000 / npts  # us/point

    if success:
        print(f"{npts:<10} {time_prep:<12.2f} {time_fld:<12.2f} {time_total:<12.2f} {time_per_point:<12.2f}")
        results.append({
            'npts': npts,
            'time_prep': time_prep,
            'time_fld': time_fld,
            'time_total': time_total,
            'time_per_point': time_per_point,
            'success': True
        })
    else:
        print(f"{npts:<10} {time_prep:<12.2f} {time_fld:<12.2f} {time_total:<12.2f} ERROR")
        results.append({
            'npts': npts,
            'time_prep': time_prep,
            'time_fld': time_fld,
            'time_total': time_total,
            'success': False,
            'error': error_msg
        })
        print(f"         Error: {error_msg}")

print("-" * 80)
print()

# Analysis
print("[Analysis] Performance characteristics:")
print()

successful_results = [r for r in results if r['success']]

if len(successful_results) >= 2:
    # Check if time scales linearly with npts
    first = successful_results[0]
    last = successful_results[-1]

    # Skip scaling analysis if first measurement is too fast to measure accurately
    if first['time_total'] > 0.01:  # At least 0.01 ms
        scaling_factor = last['time_total'] / first['time_total']
        npts_factor = last['npts'] / first['npts']

        print(f"Scaling from {first['npts']} to {last['npts']} points:")
        print(f"  Time increased by: {scaling_factor:.2f}x")
        print(f"  Points increased by: {npts_factor:.2f}x")
        print()

        if abs(scaling_factor - npts_factor) < 0.2 * npts_factor:
            print("  Result: Linear O(N) scaling (GOOD - expected)")
        elif scaling_factor > npts_factor * 1.5:
            print("  Result: Super-linear scaling (BAD - bottleneck exists)")
            print(f"  Actual: O(N^{np.log(scaling_factor) / np.log(npts_factor):.2f})")
        else:
            print("  Result: Sub-linear scaling (EXCELLENT - benefits from batching)")
        print()
    else:
        print(f"Scaling analysis skipped: first measurement ({first['time_total']:.3f} ms) too fast to measure accurately")
        print()

# Time per point analysis
print("Time per point trend:")
for r in successful_results:
    print(f"  {r['npts']:>4} points: {r['time_per_point']:>8.2f} us/point")
print()

if len(successful_results) >= 2:
    first_tpp = successful_results[0]['time_per_point']
    last_tpp = successful_results[-1]['time_per_point']

    if last_tpp < first_tpp * 0.8:
        print("  Trend: Time/point DECREASING (GOOD - amortized overhead)")
    elif last_tpp > first_tpp * 1.2:
        print("  Trend: Time/point INCREASING (BAD - degrading performance)")
    else:
        print("  Trend: Time/point CONSTANT (expected for O(N) scaling)")
    print()

# Breakdown analysis
print("Time breakdown (for largest successful test):")
if successful_results:
    largest = successful_results[-1]
    prep_pct = largest['time_prep'] / largest['time_total'] * 100
    fld_pct = largest['time_fld'] / largest['time_total'] * 100

    print(f"  Preparation: {largest['time_prep']:.2f} ms ({prep_pct:.1f}%)")
    print(f"  Radia.Fld(): {largest['time_fld']:.2f} ms ({fld_pct:.1f}%)")
    print()

    if fld_pct > 90:
        print("  Bottleneck: Radia.Fld() dominates (Radia-side issue)")
    elif prep_pct > 50:
        print("  Bottleneck: Python list preparation (Python-side issue)")
    else:
        print("  Bottleneck: Balanced (no single dominant factor)")
    print()

# Estimate time for larger point counts
print("Projected time for larger point counts:")
if len(successful_results) >= 2:
    # Use last two points to estimate scaling
    r1 = successful_results[-2]
    r2 = successful_results[-1]

    # Calculate time per point (should be constant for linear scaling)
    avg_time_per_point = (r1['time_per_point'] + r2['time_per_point']) / 2

    for target_npts in [3000, 5000, 10000]:
        projected_time = avg_time_per_point * target_npts / 1000  # ms
        print(f"  {target_npts} points: ~{projected_time:.0f} ms ({projected_time/1000:.1f} s)")
    print()

print("="*80)
print("RECOMMENDATION")
print("="*80)
print()

if successful_results:
    largest = successful_results[-1]
    if largest['time_per_point'] < 10:
        print("Performance: EXCELLENT (< 10 us/point)")
        print("Recommendation: Safe to use batch evaluation with large point counts")
    elif largest['time_per_point'] < 100:
        print("Performance: GOOD (< 100 us/point)")
        print("Recommendation: Batch evaluation beneficial for N > 100 points")
    elif largest['time_per_point'] < 1000:
        print("Performance: MODERATE (< 1 ms/point)")
        print("Recommendation: Batch evaluation beneficial only for N > 500 points")
    else:
        print("Performance: POOR (> 1 ms/point)")
        print("Recommendation: Investigate Radia.Fld() performance issue")
        print("Possible causes:")
        print("  - H-matrix construction overhead")
        print("  - Python/C interface overhead")
        print("  - Radia internal implementation issue")
else:
    print("All tests failed - cannot evaluate performance")

print()
print("="*80)

rad.UtiDelAll()
