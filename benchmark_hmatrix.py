#!/usr/bin/env python
"""
H-Matrix Performance Benchmark

Comprehensive benchmark testing H-matrix acceleration for various problem sizes.
Tests accuracy, performance, and memory usage.

Author: Radia development team
Date: 2025-11-08
"""

import sys
sys.path.insert(0, r"S:\Radia\01_GitHub\build\Release")

import radia as rad
import numpy as np
import time

print("=" * 80)
print("H-Matrix Performance Benchmark")
print("=" * 80)

# Test configurations
test_cases = [
	{"N": 50, "M": 50, "name": "Small"},
	{"N": 100, "M": 100, "name": "Medium"},
	{"N": 200, "M": 200, "name": "Large"},
	{"N": 500, "M": 500, "name": "Very Large"},
]

results = []

for case in test_cases:
	N = case["N"]
	M = case["M"]
	name = case["name"]

	print(f"\n{'=' * 80}")
	print(f"Test Case: {name} (N={N}, M={M})")
	print(f"{'=' * 80}")

	# Create geometry
	print(f"\n[1] Creating geometry (N={N})...")
	magnets = []
	L = 8.0
	grid_size = int(np.ceil(np.sqrt(N)))

	for i in range(N):
		x = (i % grid_size) * 25.0
		y = (i // grid_size) * 25.0
		z = 0.0

		angle = i * np.pi / 20.0
		Mx = 0.8 * np.cos(angle)
		My = 0.8 * np.sin(angle)
		Mz = 0.6

		magnet = rad.ObjRecMag([x, y, z], [L, L, L], [Mx, My, Mz])
		magnets.append(magnet)

	container = rad.ObjCnt(magnets)
	print(f"  Created {N} magnets")

	# Create observation points
	print(f"[2] Creating observation points (M={M})...")
	obs_points = []
	obs_grid = int(np.ceil(np.sqrt(M)))

	for i in range(M):
		x = 125.0 + (i % obs_grid) * 10.0
		y = 125.0 + (i // obs_grid) * 10.0
		z = 30.0
		obs_points.append([x, y, z])

	print(f"  Created {M} observation points")

	# Test 1: Direct calculation
	print(f"[3] Direct calculation...")
	t0 = time.time()
	H_direct = rad.FldBatch(container, 'h', obs_points, 0)
	t_direct = time.time() - t0

	H_direct = np.array(H_direct)
	print(f"  Time: {t_direct*1000:.2f} ms")

	# Test 2: H-matrix
	print(f"[4] H-matrix acceleration...")

	# Test different eps values
	eps_values = [1e-4, 1e-6, 1e-8]
	best_eps = None
	best_accuracy = float('inf')

	for eps in eps_values:
		rad.SetHMatrixFieldEval(1, eps)

		try:
			t0 = time.time()
			H_hmat = rad.FldBatch(container, 'h', obs_points, 1)
			t_hmat = time.time() - t0

			H_hmat = np.array(H_hmat)

			# Accuracy
			diff = H_hmat - H_direct
			max_H = np.max(np.abs(H_direct))
			rms_error = np.sqrt(np.mean(diff**2))
			rel_rms = rms_error / max_H if max_H > 0 else 0

			print(f"  eps={eps:.0e}: time={t_hmat*1000:.2f}ms, rel_rms={rel_rms:.2e}")

			if rel_rms < best_accuracy:
				best_accuracy = rel_rms
				best_eps = eps
				best_time = t_hmat
				best_H = H_hmat

		except Exception as e:
			print(f"  eps={eps:.0e}: FAILED ({e})")

	# Results
	if best_eps is not None:
		speedup = t_direct / best_time if best_time > 0 else 0

		result = {
			"name": name,
			"N": N,
			"M": M,
			"t_direct": t_direct,
			"t_hmat": best_time,
			"speedup": speedup,
			"accuracy": best_accuracy,
			"eps": best_eps,
		}
		results.append(result)

		print(f"\n[Summary]")
		print(f"  Best eps: {best_eps:.0e}")
		print(f"  Direct: {t_direct*1000:.2f} ms")
		print(f"  H-matrix: {best_time*1000:.2f} ms")
		print(f"  Speedup: {speedup:.2f}x")
		print(f"  Accuracy: {best_accuracy:.2e} (relative RMS)")

	# Cleanup
	rad.UtiDelAll()

# Overall summary
print("\n" + "=" * 80)
print("Overall Summary")
print("=" * 80)

print(f"\n{'Name':<15} {'N':<6} {'M':<6} {'Direct(ms)':<12} {'H-mat(ms)':<12} {'Speedup':<10} {'Accuracy':<12} {'eps':<8}")
print("-" * 80)

for r in results:
	print(f"{r['name']:<15} {r['N']:<6} {r['M']:<6} {r['t_direct']*1000:<12.2f} {r['t_hmat']*1000:<12.2f} {r['speedup']:<10.2f} {r['accuracy']:<12.2e} {r['eps']:<8.0e}")

# Analysis
print(f"\n{'=' * 80}")
print("Analysis")
print("=" * 80)

if len(results) > 0:
	max_speedup = max(r['speedup'] for r in results)
	best_case = max(results, key=lambda r: r['speedup'])

	print(f"\nBest performance:")
	print(f"  Case: {best_case['name']} (N={best_case['N']}, M={best_case['M']})")
	print(f"  Speedup: {best_case['speedup']:.2f}x")
	print(f"  Accuracy: {best_case['accuracy']:.2e}")

	# Complexity analysis
	print(f"\nComplexity verification:")
	for r in results:
		ops_direct = r['N'] * r['M']
		ops_hmat_theory = (r['N'] + r['M']) * np.log(r['N'] + r['M'])
		theoretical_speedup = ops_direct / ops_hmat_theory
		print(f"  {r['name']}: actual={r['speedup']:.2f}x, theoretical={theoretical_speedup:.2f}x")

print("\n" + "=" * 80)
