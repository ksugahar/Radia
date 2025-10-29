#!/usr/bin/env python
"""
Heavy Computational Benchmark
Tests with more complex geometry and larger point counts
"""

import subprocess
import sys
import os

def create_heavy_benchmark():
	script = '''
import sys
import os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'build', 'lib', 'Release'))

import radia as rad

def create_complex_magnet():
	"""Create a more complex magnet assembly"""
	magnets = []

	# Central magnet
	mag1 = rad.ObjRecMag([0, 0, 0], [30, 30, 30])
	rad.ObjSetM(mag1, [0, 0, 1500])
	magnets.append(mag1)

	# Surrounding magnets
	for i in range(4):
		angle = i * 90
		x = 40 * (1 if i % 2 == 0 else 0) * (1 if i < 2 else -1)
		y = 40 * (1 if i % 2 == 1 else 0) * (1 if (i==1 or i==2) else -1)

		mag = rad.ObjRecMag([x, y, 0], [15, 15, 30])
		rad.ObjSetM(mag, [0, 0, 1200])
		magnets.append(mag)

	# Create container
	return rad.ObjCnt(magnets)

def benchmark_large_grid(grid_size, z_plane, test_name):
	"""Benchmark with complex magnet and large grid"""
	magnet = create_complex_magnet()

	# Create large grid
	points = []
	for i in range(grid_size):
		for j in range(grid_size):
			x = -100 + (200.0 * i / (grid_size - 1))
			y = -100 + (200.0 * j / (grid_size - 1))
			z = z_plane
			points.append([x, y, z])

	# Benchmark
	start = time.perf_counter()
	fields = rad.Fld(magnet, 'b', points)
	end = time.perf_counter()

	elapsed = end - start
	rad.UtiDelAll()

	print(f"{test_name}: {elapsed:.4f}s ({len(points)} points)")
	return elapsed

if __name__ == '__main__':
	# Test with increasing problem sizes
	benchmark_large_grid(200, 50, "Complex 200x200")
	benchmark_large_grid(400, 50, "Complex 400x400")
	benchmark_large_grid(600, 50, "Complex 600x600")
'''
	return script

def run_with_threads(num_threads):
	env = os.environ.copy()
	env['OMP_NUM_THREADS'] = str(num_threads)

	script = create_heavy_benchmark()
	with open('_temp_heavy.py', 'w') as f:
		f.write(script)

	try:
		result = subprocess.run(
			[sys.executable, '_temp_heavy.py'],
			env=env,
			capture_output=True,
			text=True,
			timeout=600
		)
		return result.stdout
	finally:
		if os.path.exists('_temp_heavy.py'):
			os.remove('_temp_heavy.py')

def parse_results(output):
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
	print("RADIA OPENMP BENCHMARK - COMPLEX GEOMETRY")
	print("=" * 75)
	print("\nTesting with complex magnet assembly (5 magnets)\n")

	threads = [1, 2, 4, 8]
	all_results = {}

	for nt in threads:
		print(f"Testing with {nt} thread(s)...")
		output = run_with_threads(nt)
		print(output)
		all_results[nt] = parse_results(output)

	# Results
	print("\n" + "=" * 75)
	print("PERFORMANCE RESULTS")
	print("=" * 75)

	tests = ['Complex 200x200', 'Complex 400x400', 'Complex 600x600']

	for test in tests:
		print(f"\n{test}:")
		print("-" * 75)
		print(f"{'Threads':<10} {'Time (s)':<12} {'Speedup':<12} {'Efficiency':12}")
		print("-" * 75)

		baseline = all_results.get(1, {}).get(test)

		for nt in threads:
			time_val = all_results.get(nt, {}).get(test)

			if time_val and baseline and baseline > 0:
				speedup = baseline / time_val
				efficiency = (speedup / nt) * 100
				print(f"{nt:<10} {time_val:<12.4f} {speedup:<12.2f}x {efficiency:>11.1f}%")

	# Summary
	print("\n" + "=" * 75)
	speedups_8 = []
	for test in tests:
		baseline = all_results.get(1, {}).get(test)
		time_8 = all_results.get(8, {}).get(test)
		if baseline and time_8 and baseline > 0:
			speedups_8.append(baseline / time_8)

	if speedups_8:
		avg = sum(speedups_8) / len(speedups_8)
		eff = (avg / 8) * 100
		print(f"Average 8-thread speedup: {avg:.2f}x (efficiency: {eff:.1f}%)")
	print("=" * 75)

if __name__ == '__main__':
	main()
