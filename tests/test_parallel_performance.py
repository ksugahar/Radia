#!/usr/bin/env python
"""
Radia OpenMP Parallel Performance Test
Tests field computation performance with different thread counts
"""

import sys
import time
import os

# Add dist directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'build', 'lib', 'Release'))

import radia as rad

def create_test_magnet():
	"""Create a test magnet configuration"""
	# Create a rectangular block magnet
	magnet = rad.ObjRecMag([0, 0, 0], [10, 10, 10])
	rad.ObjSetM(magnet, [0, 0, 1000])  # 1000 A/m magnetization in Z
	return magnet

def test_field_computation_1d(magnet, num_points=1000):
	"""Test field computation along a line"""
	start_time = time.perf_counter()

	# Compute field along Z axis from z=0 to z=100mm
	# Use individual point calculations instead of line calculation
	fields = []
	for i in range(num_points):
		z = 100.0 * i / (num_points - 1)
		field = rad.Fld(magnet, 'b', [0, 0, z])
		fields.append(field)

	end_time = time.perf_counter()
	return end_time - start_time, len(fields)

def test_field_computation_2d(magnet, grid_size=50):
	"""Test field computation on a 2D grid"""
	start_time = time.perf_counter()

	# Create a grid of points
	points = []
	for i in range(grid_size):
		for j in range(grid_size):
			x = -25 + (50.0 * i / (grid_size - 1))
			y = -25 + (50.0 * j / (grid_size - 1))
			z = 20.0
			points.append([x, y, z])

	# Compute field at all points
	fields = []
	for point in points:
		field = rad.Fld(magnet, 'b', point)
		fields.append(field)

	end_time = time.perf_counter()
	return end_time - start_time, len(points)

def test_relaxation_performance():
	"""Test relaxation solver performance with multiple elements"""
	start_time = time.perf_counter()

	# Create multiple magnetic elements
	elements = []
	for i in range(10):
		for j in range(10):
			x = -50 + i * 10
			y = -50 + j * 10
			mag = rad.ObjRecMag([x, y, 0], [8, 8, 10])
			rad.ObjSetM(mag, [0, 0, 500])
			elements.append(mag)

	# Create container
	container = rad.ObjCnt(elements)

	# Create interaction
	interaction = rad.Solve(container, 0.0001, 1000)

	end_time = time.perf_counter()
	return end_time - start_time

def set_thread_count(num_threads):
	"""Set the number of OpenMP threads"""
	os.environ['OMP_NUM_THREADS'] = str(num_threads)
	# Note: This needs to be set before importing radia
	# So we'll need to test this by restarting Python

def run_performance_test(test_name, test_func, *args):
	"""Run a performance test and return timing"""
	print(f"  Running {test_name}...", end='', flush=True)
	try:
		result = test_func(*args)
		if isinstance(result, tuple):
			elapsed_time = result[0]
			extra_info = result[1] if len(result) > 1 else None
		else:
			elapsed_time = result
			extra_info = None

		print(f" {elapsed_time:.4f}s", end='')
		if extra_info is not None:
			print(f" ({extra_info} points)", end='')
		print()
		return elapsed_time
	except Exception as e:
		print(f" FAILED: {e}")
		return None

def main():
	print("=" * 70)
	print("RADIA OPENMP PARALLEL PERFORMANCE TEST")
	print("=" * 70)
	print()

	# Check current thread setting
	num_threads = os.environ.get('OMP_NUM_THREADS', 'not set')
	print(f"Current OMP_NUM_THREADS: {num_threads}")
	print()

	# Get Radia version
	try:
		version = rad.UtiVer()
		print(f"Radia version: {version}")
	except:
		print("Radia version: Unknown")
	print()

	print("=" * 70)
	print("TEST 1: Field Computation Along Line (1000 points)")
	print("=" * 70)
	magnet = create_test_magnet()
	time_1d = run_performance_test("1D field computation", test_field_computation_1d, magnet, 1000)
	print()

	print("=" * 70)
	print("TEST 2: Field Computation on 2D Grid (50x50 = 2500 points)")
	print("=" * 70)
	time_2d = run_performance_test("2D field computation", test_field_computation_2d, magnet, 50)
	print()

	print("=" * 70)
	print("TEST 3: Large 1D Field Computation (5000 points)")
	print("=" * 70)
	time_1d_large = run_performance_test("Large 1D computation", test_field_computation_1d, magnet, 5000)
	print()

	print("=" * 70)
	print("TEST 4: Very Large 1D Field Computation (10000 points)")
	print("=" * 70)
	time_1d_xlarge = run_performance_test("Very large 1D computation", test_field_computation_1d, magnet, 10000)
	print()

	# Summary
	print("=" * 70)
	print("PERFORMANCE SUMMARY")
	print("=" * 70)
	if time_1d:
		print(f"  1D (1000 pts):   {time_1d:.4f}s")
	if time_2d:
		print(f"  2D (2500 pts):   {time_2d:.4f}s")
	if time_1d_large:
		print(f"  1D (5000 pts):   {time_1d_large:.4f}s")
	if time_1d_xlarge:
		print(f"  1D (10000 pts):  {time_1d_xlarge:.4f}s")
	print("=" * 70)

	# Cleanup
	rad.UtiDelAll()

if __name__ == '__main__':
	main()
