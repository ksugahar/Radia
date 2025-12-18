#!/usr/bin/env python
"""
Hexahedral Benchmark using Unified Conditions (Linear Material)

Generates benchmark results for hexahedron_msc/{lu,bicgstab}/ directories
using Radia's ObjDivMag for hexahedral mesh generation.

This script uses the same conditions as benchmark_tetra_unified.py
for fair comparison between tetrahedral and hexahedral elements.

Usage:
    python benchmark_hexa_unified.py --lu 3 4 5 6 8 10
    python benchmark_hexa_unified.py --bicgstab 3 4 5 6 8 10
    python benchmark_hexa_unified.py 3 4 5 6 8 10  # runs both

Author: Radia Development Team
Date: 2025-12-12
"""

import sys
import os
import time
import json
import argparse

# Path setup
_build_path = os.path.join(os.path.dirname(__file__), '../../../build/Release')
_src_path = os.path.join(os.path.dirname(__file__), '../../../src/radia')
sys.path.insert(0, _build_path)
sys.path.append(_src_path)

import numpy as np
import radia as rad

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False


def get_peak_memory_mb():
    """Get peak memory usage in MB (Windows: peak_wset)"""
    if not HAS_PSUTIL:
        return None
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    if hasattr(mem_info, 'peak_wset'):
        return mem_info.peak_wset / (1024 * 1024)  # MB
    else:
        return mem_info.rss / (1024 * 1024)  # MB (fallback)

# Import unified benchmark conditions
from benchmark_conditions import (
    MU_0, CUBE_SIZE, H_EXT, B_EXT, CHI, MU_R,
    SOLVER_TOLERANCE, MAX_ITERATIONS, n_div_to_maxh, M_ANALYTICAL_Z
)


def benchmark_hexahedra(n_div, solver_method, output_dir):
    """Benchmark hexahedral mesh (ObjDivMag) with linear material."""
    rad.FldUnits('m')
    rad.UtiDelAll()

    solver_name = 'lu' if solver_method == 0 else 'bicgstab'
    n_elements = n_div ** 3

    print('=' * 70)
    print('HEXAHEDRAL MESH: n_div=%d (%d elements), solver=%s' % (n_div, n_elements, solver_name))
    print('=' * 70)

    # Create mesh
    t_mesh_start = time.time()
    cube = rad.ObjRecMag([0, 0, 0], [CUBE_SIZE, CUBE_SIZE, CUBE_SIZE], [0, 0, 0])
    rad.ObjDivMag(cube, [n_div, n_div, n_div])
    t_mesh = time.time() - t_mesh_start

    print('Generated %d hexahedral elements' % n_elements)

    # Apply linear material
    mat = rad.MatLin(CHI)
    rad.MatApl(cube, mat)

    # External field
    ext = rad.ObjBckg([0, 0, B_EXT])
    grp = rad.ObjCnt([cube, ext])

    # Solve
    print('Solving...')
    t_solve_start = time.time()
    try:
        result = rad.Solve(grp, SOLVER_TOLERANCE, MAX_ITERATIONS, solver_method)
        t_solve = time.time() - t_solve_start

        # Measure peak memory after solve
        peak_memory_mb = get_peak_memory_mb()

        # Get magnetization
        all_M = rad.ObjM(cube)
        M_list = [m[1] for m in all_M]
        M_avg_z = np.mean([m[2] for m in M_list])

        n_iter = int(result[3]) if result[3] else 0
        converged = n_iter < MAX_ITERATIONS and not np.isnan(M_avg_z)
        residual = result[0] if result[0] else 0.0
    except Exception as e:
        print('Solve failed: %s' % e)
        return None

    # Compare with analytical solution
    error_vs_analytical = abs(M_avg_z - M_ANALYTICAL_Z) / M_ANALYTICAL_Z * 100

    print('Mesh time:    %.4f s' % t_mesh)
    print('Solve time:   %.3f s' % t_solve)
    print('Iterations:   %d' % n_iter)
    print('Converged:    %s' % ('Yes' if converged else 'No'))
    print('M_avg_z:      %.0f A/m' % M_avg_z)
    print('Analytical:   %.0f A/m' % M_ANALYTICAL_Z)
    print('Error:        %.2f%%' % error_vs_analytical)
    if peak_memory_mb is not None:
        print('Peak memory:  %.1f MB' % peak_memory_mb)
    print()

    result_data = {
        'element_type': 'hexa',
        'mesh_description': 'n_div=%d' % n_div,
        'n_div': n_div,
        'n_elements': n_elements,
        'ndof': n_elements * 3,
        'H_ext': H_EXT,
        'mu_r': MU_R,
        'chi': CHI,
        't_mesh': t_mesh,
        't_solve': t_solve,
        'solver_method': solver_method,
        'solver_name': solver_name,
        'converged': converged,
        'residual': residual,
        'iterations': n_iter,
        'M_avg_z': M_avg_z,
        'M_analytical_z': M_ANALYTICAL_Z,
        'error_percent': error_vs_analytical,
        # For comparison with tetra: approximate equivalent maxh
        'equiv_maxh': n_div_to_maxh(n_div),
    }
    if peak_memory_mb is not None:
        result_data['peak_memory_mb'] = peak_memory_mb

    # Save result
    os.makedirs(output_dir, exist_ok=True)
    filename = 'hexa_n%d_results.json' % n_div
    filepath = os.path.join(output_dir, filename)

    with open(filepath, 'w') as f:
        json.dump(result_data, f, indent=2)
    print('Saved: %s' % filepath)

    return result_data


def main():
    parser = argparse.ArgumentParser(description='Hexahedral benchmark using unified conditions (linear material)')
    parser.add_argument('--lu', action='store_true', help='Use LU solver (saves to hexahedron_msc/lu/)')
    parser.add_argument('--bicgstab', action='store_true', help='Use BiCGSTAB solver (saves to hexahedron_msc/bicgstab/)')
    parser.add_argument('n_div_values', nargs='*', type=int, default=[2, 3, 4, 5, 6, 8, 10],
                       help='n_div values (subdivisions per edge, default: 2 3 4 5 6 8 10)')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # If neither --lu nor --bicgstab is specified, run both
    run_lu = args.lu or (not args.lu and not args.bicgstab)
    run_bicgstab = args.bicgstab or (not args.lu and not args.bicgstab)

    print('=' * 70)
    print('HEXAHEDRAL BENCHMARK - LINEAR MATERIAL (Unified Conditions)')
    print('=' * 70)
    print('Cube size: %.1f m' % CUBE_SIZE)
    print('mu_r: %d, chi: %d' % (MU_R, CHI))
    print('H_ext: %.0f A/m' % H_EXT)
    print('Analytical M_z: %.0f A/m' % M_ANALYTICAL_Z)
    print('n_div values: %s' % args.n_div_values)
    print()
    print('Mesh correspondence (n_div -> approx maxh):')
    for n in args.n_div_values:
        print('  n_div=%d -> maxh ~= %.3fm, elements=%d' % (n, n_div_to_maxh(n), n**3))
    print()

    results_lu = []
    results_bicgstab = []

    for n_div in args.n_div_values:
        if run_lu:
            output_dir = os.path.join(script_dir, 'hexahedron_msc', 'lu')
            r = benchmark_hexahedra(n_div, 0, output_dir)
            if r:
                results_lu.append(r)

        if run_bicgstab:
            output_dir = os.path.join(script_dir, 'hexahedron_msc', 'bicgstab')
            r = benchmark_hexahedra(n_div, 1, output_dir)
            if r:
                results_bicgstab.append(r)

    # Summary
    print('=' * 70)
    print('SUMMARY')
    print('=' * 70)

    if results_lu:
        print('\nLU Solver (hexahedron_msc/lu/):\n')
        print('%-10s %10s %10s %8s %12s %10s %8s' % ('n_div', 'Elements', 'Time (s)', 'Iter', 'M_avg_z', 'Error', 'Conv'))
        print('-' * 75)
        for r in results_lu:
            print('%-10d %10d %10.3f %8d %12.0f %9.2f%% %8s' % (
                r['n_div'], r['n_elements'], r['t_solve'],
                r['iterations'], r['M_avg_z'], r['error_percent'],
                'Yes' if r['converged'] else 'No'))

    if results_bicgstab:
        print('\nBiCGSTAB Solver (hexahedron_msc/bicgstab/):\n')
        print('%-10s %10s %10s %8s %12s %10s %8s' % ('n_div', 'Elements', 'Time (s)', 'Iter', 'M_avg_z', 'Error', 'Conv'))
        print('-' * 75)
        for r in results_bicgstab:
            print('%-10d %10d %10.3f %8d %12.0f %9.2f%% %8s' % (
                r['n_div'], r['n_elements'], r['t_solve'],
                r['iterations'], r['M_avg_z'], r['error_percent'],
                'Yes' if r['converged'] else 'No'))

    print('=' * 70)


if __name__ == '__main__':
    main()
