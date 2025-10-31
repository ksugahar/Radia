#!/usr/bin/env python
"""
OpenMP Parallelization Benchmark for Radia
Tests performance with different thread counts on large problem sizes
"""

import subprocess
import sys
import os
import time

# Add project root's build directory to path
import sys
import os
from pathlib import Path

# Find project root (works from any test subdirectory)
current_file = Path(__file__).resolve()
if 'tests' in current_file.parts:
	# Find the 'tests' directory and go up one level
	tests_index = current_file.parts.index('tests')
	project_root = Path(*current_file.parts[:tests_index])
else:
	# Fallback
	project_root = current_file.parent

# Add build directory to path
build_dir = project_root / 'build' / 'lib' / 'Release'
if build_dir.exists():
	sys.path.insert(0, str(build_dir))


def create_benchmark_script():
	"""Create a benchmark script that will be run with different thread counts"""
	script = '''
import sys
import os
import time

from pathlib import Path

# Find project root (works from any test subdirectory)
current_file = Path(__file__).resolve()
if 'tests' in current_file.parts:
	tests_index = current_file.parts.index('tests')
	project_root = Path(*current_file.parts[:tests_index])
else:
	project_root = current_file.parent

# Add build directory to path
build_dir = project_root / 'build' / 'lib' / 'Release'
if build_dir.exists():
	sys.path.insert(0, str(build_dir))


import radia as rad

def benchmark_large_grid(grid_points=100):
	"""Benchmark field computation on a large 2D grid"""
	# Create magnet
	magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 20])
	rad.ObjSetM(magnet, [0, 0, 1000])

	# Create grid
	points = []
	for i in range(grid_points):
		for j in range(grid_points):
			x = -50 + (100.0 * i / (grid_points - 1))
			y = -50 + (100.0 * j / (grid_points - 1))
			z = 30.0
			points.append([x, y, z])

	# Benchmark
	start = time.perf_counter()
	fields = []
	for point in points:
		field = rad.Fld(magnet, 'b', point)
		fields.append(field)
	end = time.perf_counter()

	rad.UtiDelAll()
	return end - start, len(points)

def benchmark_3d_volume(nx=30, ny=30, nz=30):
	"""Benchmark field computation in a 3D volume"""
	# Create magnet
	magnet = rad.ObjRecMag([0, 0, 0], [15, 15, 15])
	rad.ObjSetM(magnet, [0, 0, 1200])

	# Create 3D grid
	points = []
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				x = -40 + (80.0 * i / (nx - 1))
				y = -40 + (80.0 * j / (ny - 1))
				z = 20 + (40.0 * k / (nz - 1))
				points.append([x, y, z])

	# Benchmark
	start = time.perf_counter()
	fields = []
	for point in points:
		field = rad.Fld(magnet, 'b', point)
		fields.append(field)
	end = time.perf_counter()

	rad.UtiDelAll()
	return end - start, len(points)

if __name__ == '__main__':
	print(f"Grid 100x100 (10000 points):")
	time1, pts1 = benchmark_large_grid(100)
	print(f"  Time: {time1:.4f}s, Points: {pts1}")

	print(f"\\nGrid 150x150 (22500 points):")
	time2, pts2 = benchmark_large_grid(150)
	print(f"  Time: {time2:.4f}s, Points: {pts2}")

	print(f"\\n3D Volume 30x30x30 (27000 points):")
	time3, pts3 = benchmark_3d_volume(30, 30, 30)
	print(f"  Time: {time3:.4f}s, Points: {pts3}")
'''
	return script

def run_benchmark_with_threads(num_threads):
	"""Run benchmark with specified number of threads"""
	env = os.environ.copy()
	env['OMP_NUM_THREADS'] = str(num_threads)

	# Create temporary script
	script = create_benchmark_script()
	with open('_temp_benchmark.py', 'w') as f:
		f.write(script)

	try:
		result = subprocess.run(
			[sys.executable, '_temp_benchmark.py'],
			env=env,
			capture_output=True,
			text=True,
			timeout=300
		)
		return result.stdout
	finally:
		if os.path.exists('_temp_benchmark.py'):
			os.remove('_temp_benchmark.py')

def parse_results(output):
	"""Parse benchmark results from output"""
	results = {}
	lines = output.strip().split('\n')

	for i, line in enumerate(lines):
		if 'Grid 100x100' in line:
			if i+1 < len(lines):
				time_line = lines[i+1]
				if 'Time:' in time_line:
					time_str = time_line.split('Time:')[1].split('s')[0].strip()
					results['grid_100'] = float(time_str)
		elif 'Grid 150x150' in line:
			if i+1 < len(lines):
				time_line = lines[i+1]
				if 'Time:' in time_line:
					time_str = time_line.split('Time:')[1].split('s')[0].strip()
					results['grid_150'] = float(time_str)
		elif '3D Volume' in line:
			if i+1 < len(lines):
				time_line = lines[i+1]
				if 'Time:' in time_line:
					time_str = time_line.split('Time:')[1].split('s')[0].strip()
					results['volume_3d'] = float(time_str)

	return results

def main():
	print("=" * 70)
	print("RADIA OPENMP SCALABILITY BENCHMARK")
	print("=" * 70)
	print("\nThis benchmark tests large problem sizes with 1, 2, 4, and 8 threads\n")

	thread_counts = [1, 2, 4, 8]
	all_results = {}

	# Run benchmarks
	for num_threads in thread_counts:
		print(f"Testing with {num_threads} thread(s)...")
		output = run_benchmark_with_threads(num_threads)
		print(output)
		all_results[num_threads] = parse_results(output)
		print()

	# Display results table
	print("\n" + "=" * 70)
	print("PERFORMANCE RESULTS")
	print("=" * 70)

	test_names = [
		('grid_100', '100x100 Grid (10,000 points)'),
		('grid_150', '150x150 Grid (22,500 points)'),
		('volume_3d', '3D Volume 30x30x30 (27,000 points)')
	]

	for test_key, test_name in test_names:
		print(f"\n{test_name}:")
		print("-" * 70)
		print(f"{'Threads':<10} {'Time (s)':<15} {'Speedup':<15} {'Efficiency':<15}")
		print("-" * 70)

		baseline = all_results.get(1, {}).get(test_key)

		for num_threads in thread_counts:
			time_val = all_results.get(num_threads, {}).get(test_key)

			if time_val:
				if baseline and baseline > 0:
					speedup = baseline / time_val
					efficiency = (speedup / num_threads) * 100
					print(f"{num_threads:<10} {time_val:<15.4f} {speedup:<15.2f}x {efficiency:<15.1f}%")
				else:
					print(f"{num_threads:<10} {time_val:<15.4f} {'N/A':<15} {'N/A':<15}")
			else:
				print(f"{num_threads:<10} {'FAILED':<15} {'N/A':<15} {'N/A':<15}")

	print("\n" + "=" * 70)
	print("ANALYSIS")
	print("=" * 70)

	# Calculate average speedup for 8 threads
	speedups_8 = []
	for test_key, _ in test_names:
		baseline = all_results.get(1, {}).get(test_key)
		time_8 = all_results.get(8, {}).get(test_key)
		if baseline and time_8 and baseline > 0:
			speedups_8.append(baseline / time_8)

	if speedups_8:
		avg_speedup = sum(speedups_8) / len(speedups_8)
		avg_efficiency = (avg_speedup / 8) * 100
		print(f"\nAverage speedup with 8 threads: {avg_speedup:.2f}x")
		print(f"Average parallel efficiency: {avg_efficiency:.1f}%")
		print(f"\nThis means 8 CPU cores provide {avg_speedup:.1f}x faster computation")
		print(f"compared to single-threaded execution.")

	print("=" * 70)

if __name__ == '__main__':
	main()
