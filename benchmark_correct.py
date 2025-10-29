#!/usr/bin/env python
"""
Correct OpenMP Benchmark for Radia
Uses batch point calculation to leverage C++ parallelization
"""

import subprocess
import sys
import os

def create_benchmark_script():
	"""Create benchmark script with proper multi-point calculations"""
	script = '''
import sys
import os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'build', 'lib', 'Release'))

import radia as rad

def benchmark_batch_calculation(grid_size, test_name):
	"""Benchmark using batch point calculation (leverages OpenMP)"""
	# Create magnet
	magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 20])
	rad.ObjSetM(magnet, [0, 0, 1200])

	# Create grid of points
	points = []
	for i in range(grid_size):
		for j in range(grid_size):
			x = -50 + (100.0 * i / (grid_size - 1))
			y = -50 + (100.0 * j / (grid_size - 1))
			z = 30.0
			points.append([x, y, z])

	# Benchmark - single batch call (this uses OpenMP internally)
	start = time.perf_counter()
	fields = rad.Fld(magnet, 'b', points)
	end = time.perf_counter()

	elapsed = end - start
	rad.UtiDelAll()

	print(f"{test_name}: {elapsed:.4f}s ({len(points)} points)")
	return elapsed, len(points)

def benchmark_3d_volume(nx, ny, nz, test_name):
	"""Benchmark 3D volume calculation"""
	# Create magnet
	magnet = rad.ObjRecMag([0, 0, 0], [15, 15, 15])
	rad.ObjSetM(magnet, [0, 0, 1500])

	# Create 3D grid
	points = []
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				x = -40 + (80.0 * i / (nx - 1))
				y = -40 + (80.0 * j / (ny - 1))
				z = 20 + (60.0 * k / (nz - 1))
				points.append([x, y, z])

	# Benchmark
	start = time.perf_counter()
	fields = rad.Fld(magnet, 'b', points)
	end = time.perf_counter()

	elapsed = end - start
	rad.UtiDelAll()

	print(f"{test_name}: {elapsed:.4f}s ({len(points)} points)")
	return elapsed, len(points)

if __name__ == '__main__':
	# Test different problem sizes
	benchmark_batch_calculation(100, "Grid 100x100")
	benchmark_batch_calculation(200, "Grid 200x200")
	benchmark_batch_calculation(300, "Grid 300x300")
	benchmark_3d_volume(40, 40, 40, "3D 40x40x40")
'''
	return script

def run_with_threads(num_threads):
	"""Run benchmark with specified thread count"""
	env = os.environ.copy()
	env['OMP_NUM_THREADS'] = str(num_threads)

	script = create_benchmark_script()
	with open('_temp_bench.py', 'w') as f:
		f.write(script)

	try:
		result = subprocess.run(
			[sys.executable, '_temp_bench.py'],
			env=env,
			capture_output=True,
			text=True,
			timeout=300
		)
		return result.stdout
	finally:
		if os.path.exists('_temp_bench.py'):
			os.remove('_temp_bench.py')

def parse_results(output):
	"""Parse timing results"""
	results = {}
	for line in output.strip().split('\n'):
		if ':' in line and 's (' in line:
			parts = line.split(':')
			name = parts[0].strip()
			time_str = parts[1].split('s')[0].strip()
			try:
				results[name] = float(time_str)
			except:
				pass
	return results

def main():
	print("=" * 75)
	print("RADIA OPENMP PARALLEL PERFORMANCE BENCHMARK")
	print("=" * 75)
	print("\nThis benchmark uses batch field calculations to properly test")
	print("the OpenMP parallelization at the C++ level.\n")

	threads = [1, 2, 4, 8]
	all_results = {}

	# Run benchmarks
	for nt in threads:
		print(f"{'='*75}")
		print(f"Testing with {nt} thread(s)...")
		print(f"{'='*75}")
		output = run_with_threads(nt)
		print(output)
		all_results[nt] = parse_results(output)

	# Results table
	print("\n" + "=" * 75)
	print("PERFORMANCE SUMMARY")
	print("=" * 75)

	tests = [
		'Grid 100x100',
		'Grid 200x200',
		'Grid 300x300',
		'3D 40x40x40'
	]

	for test_name in tests:
		print(f"\n{test_name}:")
		print("-" * 75)
		print(f"{'Threads':<10} {'Time (s)':<12} {'Speedup':<12} {'Efficiency':12}")
		print("-" * 75)

		baseline = all_results.get(1, {}).get(test_name)

		for nt in threads:
			time_val = all_results.get(nt, {}).get(test_name)

			if time_val and baseline and baseline > 0:
				speedup = baseline / time_val
				efficiency = (speedup / nt) * 100
				print(f"{nt:<10} {time_val:<12.4f} {speedup:<12.2f}x {efficiency:>11.1f}%")
			elif time_val:
				print(f"{nt:<10} {time_val:<12.4f} {'N/A':<12} {'N/A':>12}")
			else:
				print(f"{nt:<10} {'FAILED':<12} {'N/A':<12} {'N/A':>12}")

	# Final analysis
	print("\n" + "=" * 75)
	print("SCALABILITY ANALYSIS")
	print("=" * 75)

	# Calculate average speedup for each thread count
	for nt in [2, 4, 8]:
		speedups = []
		for test_name in tests:
			baseline = all_results.get(1, {}).get(test_name)
			time_nt = all_results.get(nt, {}).get(test_name)
			if baseline and time_nt and baseline > 0:
				speedups.append(baseline / time_nt)

		if speedups:
			avg_speedup = sum(speedups) / len(speedups)
			avg_efficiency = (avg_speedup / nt) * 100
			print(f"\n{nt} threads:")
			print(f"  Average speedup:    {avg_speedup:.2f}x")
			print(f"  Average efficiency: {avg_efficiency:.1f}%")

	# Highlight 8-core performance
	speedups_8 = []
	for test_name in tests:
		baseline = all_results.get(1, {}).get(test_name)
		time_8 = all_results.get(8, {}).get(test_name)
		if baseline and time_8 and baseline > 0:
			speedups_8.append(baseline / time_8)

	if speedups_8:
		avg_speedup_8 = sum(speedups_8) / len(speedups_8)
		print(f"\n" + "=" * 75)
		print(f"RESULT: 8 CPU cores provide {avg_speedup_8:.2f}x speedup on average")
		print(f"        (Parallel efficiency: {(avg_speedup_8/8)*100:.1f}%)")
		print("=" * 75)

if __name__ == '__main__':
	main()
