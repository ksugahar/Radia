#!/usr/bin/env python
"""
Run H-matrix examples to verify Phase 3B implementation
"""

import sys
import os
import subprocess
import time

# Add build directory to path
sys.path.insert(0, os.path.join('S:/Radia/01_GitHub', 'build/Release'))

examples_dir = 'S:/Radia/01_GitHub/examples/H-matrix'

examples = [
    'verify_field_accuracy.py',
    'benchmark_field_evaluation.py',
    'benchmark_solver.py',
]

print("="*70)
print("Running H-Matrix Examples (Phase 3B Verification)")
print("="*70)

for example in examples:
    example_path = os.path.join(examples_dir, example)

    if not os.path.exists(example_path):
        print(f"\n[SKIP] {example} - File not found")
        continue

    print(f"\n{'='*70}")
    print(f"Running: {example}")
    print('='*70)

    try:
        result = subprocess.run(
            [sys.executable, example_path],
            capture_output=True,
            text=True,
            timeout=180
        )

        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)

        if result.returncode == 0:
            print(f"\n[OK] {example} completed successfully")
        else:
            print(f"\n[FAIL] {example} failed with return code {result.returncode}")

    except subprocess.TimeoutExpired:
        print(f"\n[TIMEOUT] {example} timed out after 180 seconds")
    except Exception as e:
        print(f"\n[ERROR] {example} - {e}")

print("\n" + "="*70)
print("All examples completed")
print("="*70)
