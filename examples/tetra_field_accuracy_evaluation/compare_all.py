#!/usr/bin/env python
"""
Compare Radia Tetrahedral MSC, Hexahedral, and NGSolve Reference Solutions

This script loads results from:
1. evaluate_tetra_field.py (Radia tetra/hexa results)
2. ngsolve_reference.py (NGSolve A-formulation reference)

And creates a comprehensive comparison table.

Author: Radia Development Team
Date: 2025-12-13
"""
import os
import sys
import json
import numpy as np

print('=' * 70)
print('Comprehensive Comparison: Radia vs NGSolve')
print('=' * 70)

# =============================================================================
# Load Results
# =============================================================================
script_dir = os.path.dirname(__file__)

# Load Radia results
radia_file = os.path.join(script_dir, 'evaluation_results.json')
ngsolve_file = os.path.join(script_dir, 'ngsolve_reference_results.json')

try:
    with open(radia_file, 'r') as f:
        radia_data = json.load(f)
    print('Loaded Radia results: %s' % radia_file)
except FileNotFoundError:
    print('ERROR: Radia results not found. Run evaluate_tetra_field.py first.')
    radia_data = None

try:
    with open(ngsolve_file, 'r') as f:
        ngsolve_data = json.load(f)
    print('Loaded NGSolve results: %s' % ngsolve_file)
except FileNotFoundError:
    print('WARNING: NGSolve results not found. Run ngsolve_reference.py first.')
    ngsolve_data = None

if radia_data is None:
    print('Cannot continue without Radia results.')
    sys.exit(1)

# =============================================================================
# Extract Data
# =============================================================================
test_points = radia_data['test_points']
labels = radia_data['labels']
n_points = len(test_points)

# Radia hexahedral reference (finest mesh)
hexa_results = radia_data['hexa_results']
hexa_finest = hexa_results[-1]
B_hexa = [np.array(b) for b in hexa_finest['B_values']]

# Radia tetrahedral results
tetra_results = radia_data['tetra_results']

# NGSolve reference
if ngsolve_data:
    B_ngsolve = [np.array(b) for b in ngsolve_data['ngsolve']['B_values']]
else:
    B_ngsolve = None

# =============================================================================
# Comparison Tables
# =============================================================================
print()
print('=' * 70)
print('Comparison at Each Test Point')
print('=' * 70)

# Use finest tetrahedral mesh
if tetra_results:
    tetra_finest = tetra_results[-1]
    B_tetra = [np.array(b) for b in tetra_finest['B_values']]
else:
    tetra_finest = None
    B_tetra = None

print()
if ngsolve_data:
    print('Reference: NGSolve A-formulation (%d DOFs)' % ngsolve_data['ngsolve']['ndof'])
else:
    print('Reference: Radia Hexahedral (n_div=%d, %d elements)' %
          (hexa_finest['n_div'], hexa_finest['n_elements']))

print()

# Header
header = '%s  %s  %s' % (
    'Point'.center(25),
    '|B| Reference'.center(14),
    '|B| Tetra MSC'.center(14)
)
if ngsolve_data:
    header += '  %s' % '|B| Hexa'.center(14)
header += '  %s  %s' % ('Err Tetra%'.center(10), 'Err Hexa%'.center(10))
print(header)
print('-' * len(header))

errors_tetra = []
errors_hexa = []

for i in range(n_points):
    label = labels[i][:25].ljust(25)

    # Reference
    if B_ngsolve:
        B_ref = B_ngsolve[i]
    else:
        B_ref = B_hexa[i]

    B_ref_mag = np.linalg.norm(B_ref)

    # Tetrahedral
    if B_tetra:
        B_t = B_tetra[i]
        B_t_mag = np.linalg.norm(B_t)
        if B_ref_mag > 1e-15 and not np.isnan(B_t_mag):
            err_t = abs(B_t_mag - B_ref_mag) / B_ref_mag * 100
        else:
            err_t = np.nan
    else:
        B_t_mag = np.nan
        err_t = np.nan

    # Hexahedral
    B_h = B_hexa[i]
    B_h_mag = np.linalg.norm(B_h)
    if B_ngsolve:
        if B_ref_mag > 1e-15 and not np.isnan(B_h_mag):
            err_h = abs(B_h_mag - B_ref_mag) / B_ref_mag * 100
        else:
            err_h = np.nan
    else:
        err_h = 0.0  # Hexa is the reference

    if not np.isnan(err_t):
        errors_tetra.append(err_t)
    if not np.isnan(err_h):
        errors_hexa.append(err_h)

    # Print row
    row = '%s  %14.6e  %14.6e' % (label, B_ref_mag, B_t_mag)
    if ngsolve_data:
        row += '  %14.6e' % B_h_mag
    row += '  %10.2f  %10.2f' % (err_t, err_h)
    print(row)

# =============================================================================
# Summary Statistics
# =============================================================================
print()
print('=' * 70)
print('Summary Statistics')
print('=' * 70)

print()
print('Tetrahedral MSC Errors (vs Reference):')
if errors_tetra:
    print('  Average: %.2f%%' % np.mean(errors_tetra))
    print('  Maximum: %.2f%%' % np.max(errors_tetra))
    print('  Minimum: %.2f%%' % np.min(errors_tetra))
    print('  Std Dev: %.2f%%' % np.std(errors_tetra))
else:
    print('  No data available')

if ngsolve_data and errors_hexa:
    print()
    print('Hexahedral Errors (vs NGSolve Reference):')
    print('  Average: %.2f%%' % np.mean(errors_hexa))
    print('  Maximum: %.2f%%' % np.max(errors_hexa))
    print('  Minimum: %.2f%%' % np.min(errors_hexa))

# =============================================================================
# Convergence Analysis
# =============================================================================
print()
print('=' * 70)
print('Tetrahedral Convergence Analysis')
print('=' * 70)

print()
print('%10s  %10s  %12s  %12s  %12s' % (
    'maxh', 'Elements', 'Avg Err%', 'Max Err%', 'Time (s)'
))
print('-' * 60)

for result in tetra_results:
    maxh = result['maxh']
    n_elem = result['n_elements']
    t_time = result['time']

    B_values = [np.array(b) for b in result['B_values']]

    errors = []
    for i in range(n_points):
        if B_ngsolve:
            B_ref = B_ngsolve[i]
        else:
            B_ref = B_hexa[i]

        B_ref_mag = np.linalg.norm(B_ref)
        B_t_mag = np.linalg.norm(B_values[i])

        if B_ref_mag > 1e-15 and not np.isnan(B_t_mag):
            err = abs(B_t_mag - B_ref_mag) / B_ref_mag * 100
            errors.append(err)

    if errors:
        avg_err = np.mean(errors)
        max_err = np.max(errors)
    else:
        avg_err = np.nan
        max_err = np.nan

    print('%10.2f  %10d  %12.2f  %12.2f  %12.3f' % (
        maxh, n_elem, avg_err, max_err, t_time
    ))

# =============================================================================
# Pass/Fail Assessment
# =============================================================================
print()
print('=' * 70)
print('Assessment')
print('=' * 70)

if errors_tetra:
    avg_err = np.mean(errors_tetra)
    max_err = np.max(errors_tetra)

    print()
    if avg_err < 5.0:
        print('[PASS] Average error < 5%% (actual: %.2f%%)' % avg_err)
        print('       Tetrahedral MSC field computation is ACCURATE')
    elif avg_err < 10.0:
        print('[WARNING] Average error between 5%% and 10%% (actual: %.2f%%)' % avg_err)
        print('          Tetrahedral MSC field has MODERATE accuracy')
    elif avg_err < 20.0:
        print('[WARNING] Average error between 10%% and 20%% (actual: %.2f%%)' % avg_err)
        print('          Tetrahedral MSC field has ACCEPTABLE accuracy')
    else:
        print('[FAIL] Average error >= 20%% (actual: %.2f%%)' % avg_err)
        print('       Tetrahedral MSC field computation needs investigation')

    if max_err > 50.0:
        print()
        print('[WARNING] Maximum error is %.2f%% - check individual points' % max_err)
else:
    print()
    print('[INFO] No tetrahedral results available for assessment')

# =============================================================================
# Save Summary
# =============================================================================
summary_file = os.path.join(script_dir, 'comparison_summary.json')

summary = {
    'reference': 'NGSolve' if ngsolve_data else 'Hexahedral',
    'tetra_finest': {
        'maxh': tetra_finest['maxh'] if tetra_finest else None,
        'n_elements': tetra_finest['n_elements'] if tetra_finest else None,
        'avg_error': float(np.mean(errors_tetra)) if errors_tetra else None,
        'max_error': float(np.max(errors_tetra)) if errors_tetra else None,
    },
    'hexa_finest': {
        'n_div': hexa_finest['n_div'],
        'n_elements': hexa_finest['n_elements'],
        'avg_error': float(np.mean(errors_hexa)) if ngsolve_data and errors_hexa else 0.0,
    },
    'assessment': 'PASS' if errors_tetra and np.mean(errors_tetra) < 5.0 else 'CHECK'
}

with open(summary_file, 'w') as f:
    json.dump(summary, f, indent=2)

print()
print('Summary saved to: %s' % summary_file)
print()
print('=' * 70)
