#!/usr/bin/env python
"""
Benchmark Radia with different thread counts
Tests 1, 2, 4, and 8 threads to measure speedup
"""

import subprocess
import sys
import os

def run_test_with_threads(num_threads):
	"""Run the performance test with specified number of threads"""
	env = os.environ.copy()
	env['OMP_NUM_THREADS'] = str(num_threads)

	print(f"\n{'=' * 70}")
	print(f"TESTING WITH {num_threads} THREAD(S)")
	print(f"{'=' * 70}\n")

	result = subprocess.run(
		[sys.executable, 'test_parallel_performance.py'],
		env=env,
		capture_output=True,
		text=True
	)

	print(result.stdout)
	if result.stderr:
		print("STDERR:", result.stderr)

	return result.stdout

def parse_time_from_output(output, test_name):
	"""Parse execution time from test output"""
	lines = output.split('\n')
	for line in lines:
		if test_name in line and 'Running' in line:
			# Extract time from line like "  Running 1D field computation... 0.1234s"
			parts = line.split('...')
			if len(parts) > 1:
				time_part = parts[1].strip().split()[0]
				try:
					return float(time_part.replace('s', ''))
				except:
					pass
	return None

def main():
	print("=" * 70)
	print("RADIA OPENMP SCALABILITY BENCHMARK")
	print("=" * 70)
	print()
	print("This benchmark will test Radia performance with 1, 2, 4, and 8 threads")
	print("to measure the speedup from CPU parallelization.")
	print()

	# Thread counts to test
	thread_counts = [1, 2, 4, 8]

	# Store results
	results = {}

	# Run tests
	for num_threads in thread_counts:
		output = run_test_with_threads(num_threads)
		results[num_threads] = output

	# Parse and display results
	print("\n" + "=" * 70)
	print("SCALABILITY ANALYSIS")
	print("=" * 70)

	test_names = [
		("1D field computation", "1D (1000 pts)"),
		("Large 1D computation", "1D (5000 pts)"),
		("Very large 1D computation", "1D (10000 pts)"),
	]

	for test_key, display_name in test_names:
		print(f"\n{display_name}:")
		print("-" * 70)
		print(f"{'Threads':<10} {'Time (s)':<15} {'Speedup':<15} {'Efficiency':<15}")
		print("-" * 70)

		times = {}
		for num_threads in thread_counts:
			time = parse_time_from_output(results[num_threads], test_key)
			times[num_threads] = time

		baseline_time = times.get(1)

		for num_threads in thread_counts:
			time = times.get(num_threads)
			if time is not None:
				if baseline_time and baseline_time > 0:
					speedup = baseline_time / time
					efficiency = (speedup / num_threads) * 100
					print(f"{num_threads:<10} {time:<15.4f} {speedup:<15.2f}x {efficiency:<15.1f}%")
				else:
					print(f"{num_threads:<10} {time:<15.4f} {'N/A':<15} {'N/A':<15}")
			else:
				print(f"{num_threads:<10} {'FAILED':<15} {'N/A':<15} {'N/A':<15}")

	print("\n" + "=" * 70)
	print("SUMMARY")
	print("=" * 70)
	print("\nSpeedup = Time(1 thread) / Time(N threads)")
	print("Efficiency = (Speedup / N threads) * 100%")
	print("\nIdeal scaling:")
	print("  - 2 threads: 2.0x speedup (100% efficiency)")
	print("  - 4 threads: 4.0x speedup (100% efficiency)")
	print("  - 8 threads: 8.0x speedup (100% efficiency)")
	print("\nGood scaling (typical for real applications):")
	print("  - 2 threads: 1.8-2.0x speedup (90-100% efficiency)")
	print("  - 4 threads: 3.2-3.8x speedup (80-95% efficiency)")
	print("  - 8 threads: 5.0-7.0x speedup (60-90% efficiency)")
	print("=" * 70)

if __name__ == '__main__':
	main()
