#!/usr/bin/env python
"""
Run all H-matrix benchmarks and generate summary report
"""

import sys
import os
sys.path.insert(0, os.path.join('S:/Radia/01_GitHub', 'build/Release'))

import subprocess
import time

print("="*70)
print("H-MATRIX BENCHMARKS - COMPREHENSIVE TEST SUITE")
print("Phase 3B: Full H-matrix serialization verification")
print("="*70)

benchmarks = [
    ("Field Evaluation", "examples/solver_benchmarks/benchmark_field_evaluation.py"),
    ("Solver Performance", "examples/solver_benchmarks/benchmark_solver.py"),
    ("Parallel Construction", "examples/solver_benchmarks/benchmark_parallel_construction.py"),
]

results = []

for name, script_path in benchmarks:
    full_path = os.path.join('S:/Radia/01_GitHub', script_path)

    print(f"\n{'='*70}")
    print(f"RUNNING: {name}")
    print(f"Script: {script_path}")
    print('='*70)

    if not os.path.exists(full_path):
        print(f"[SKIP] File not found: {full_path}")
        results.append((name, "SKIPPED", 0))
        continue

    start_time = time.time()

    try:
        result = subprocess.run(
            [sys.executable, full_path],
            capture_output=True,
            text=True,
            timeout=300,
            cwd='S:/Radia/01_GitHub'
        )

        elapsed = time.time() - start_time

        if result.returncode == 0:
            print(f"\n[OK] {name} completed in {elapsed:.2f}s")
            results.append((name, "PASS", elapsed))

            # Show last 20 lines of output for summary
            lines = result.stdout.strip().split('\n')
            if len(lines) > 20:
                print("\n... (showing last 20 lines) ...")
                for line in lines[-20:]:
                    print(line)
            else:
                print(result.stdout)
        else:
            print(f"\n[FAIL] {name} failed with return code {result.returncode}")
            print("STDERR:", result.stderr[:500])
            results.append((name, "FAIL", elapsed))

    except subprocess.TimeoutExpired:
        elapsed = time.time() - start_time
        print(f"\n[TIMEOUT] {name} timed out after {elapsed:.2f}s")
        results.append((name, "TIMEOUT", elapsed))
    except Exception as e:
        elapsed = time.time() - start_time
        print(f"\n[ERROR] {name} - {e}")
        results.append((name, "ERROR", elapsed))

# Generate summary report
print("\n" + "="*70)
print("SUMMARY REPORT")
print("="*70)

print(f"\n{'Benchmark':<30} {'Status':<15} {'Time (s)':<10}")
print("-"*70)

passed = 0
failed = 0
for name, status, elapsed in results:
    print(f"{name:<30} {status:<15} {elapsed:<10.2f}")
    if status == "PASS":
        passed += 1
    else:
        failed += 1

print("-"*70)
print(f"Total: {len(results)} benchmarks")
print(f"Passed: {passed}")
print(f"Failed: {failed}")

if failed == 0:
    print("\n[OK] All H-matrix benchmarks completed successfully!")
    print("\nPhase 3B Implementation Status:")
    print("  - Full H-matrix serialization: WORKING")
    print("  - Disk cache persistence: VERIFIED")
    print("  - Field evaluation: TESTED")
    print("  - Solver performance: TESTED")
    print("  - Parallel construction: TESTED")
else:
    print(f"\n[WARN] {failed} benchmark(s) failed or were skipped")

print("\n" + "="*70)
