#!/usr/bin/env python
"""
Hexahedral Benchmark using Unified Conditions

Generates benchmark results for hexahedron/{lu,bicgstab}/ directories
using Radia's ObjDivMag for hexahedral mesh generation.

This script uses the same conditions as benchmark_tetra_netgen.py
for fair comparison between tetrahedral and hexahedral elements.

Usage:
    python benchmark_hexa_unified.py --lu 3 4 5 6
    python benchmark_hexa_unified.py --bicgstab 3 4 5 6
    python benchmark_hexa_unified.py 3 4 5 6  # runs both

Author: Radia Development Team
Date: 2025-12-06
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

# Import unified benchmark conditions
from benchmark_conditions import (
    MU_0, CUBE_SIZE, H_EXT, B_EXT, HM_DATA,
    SOLVER_TOLERANCE, MAX_ITERATIONS, n_div_to_maxh
)


def benchmark_hexahedra(n_div, solver_method, output_dir):
    """Benchmark hexahedral mesh (ObjDivMag)."""
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

    # Apply nonlinear material
    mat = rad.MatSatIsoTab(HM_DATA)
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

    print('Mesh time:    %.4f s' % t_mesh)
    print('Solve time:   %.3f s' % t_solve)
    print('Iterations:   %d' % n_iter)
    print('Converged:    %s' % ('Yes' if converged else 'No'))
    print('M_avg_z:      %.0f A/m' % M_avg_z)
    print()

    result_data = {
        'element_type': 'hexa',
        'mesh_description': 'n_div=%d' % n_div,
        'n_div': n_div,
        'n_elements': n_elements,
        'ndof': n_elements * 3,
        'H_ext': H_EXT,
        't_mesh': t_mesh,
        't_solve': t_solve,
        'solver_method': solver_method,
        'solver_name': solver_name,
        'converged': converged,
        'residual': residual,
        'nonl_iterations': n_iter,
        'M_avg_z': M_avg_z,
        # For comparison with tetra: approximate equivalent maxh
        'equiv_maxh': n_div_to_maxh(n_div),
    }

    # Save result
    os.makedirs(output_dir, exist_ok=True)
    filename = 'hexa_n%d_results.json' % n_div
    filepath = os.path.join(output_dir, filename)

    with open(filepath, 'w') as f:
        json.dump(result_data, f, indent=2)
    print('Saved: %s' % filepath)

    return result_data


def main():
    parser = argparse.ArgumentParser(description='Hexahedral benchmark using unified conditions')
    parser.add_argument('--lu', action='store_true', help='Use LU solver (saves to hexahedron/lu/)')
    parser.add_argument('--bicgstab', action='store_true', help='Use BiCGSTAB solver (saves to hexahedron/bicgstab/)')
    parser.add_argument('n_div_values', nargs='*', type=int, default=[3, 4, 5, 6],
                       help='n_div values (subdivisions per edge, default: 3 4 5 6)')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # If neither --lu nor --bicgstab is specified, run both
    run_lu = args.lu or (not args.lu and not args.bicgstab)
    run_bicgstab = args.bicgstab or (not args.lu and not args.bicgstab)

    print('=' * 70)
    print('HEXAHEDRAL BENCHMARK (Unified Conditions)')
    print('=' * 70)
    print('Cube size: %.1f m' % CUBE_SIZE)
    print('H_ext: %.0f A/m' % H_EXT)
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
            output_dir = os.path.join(script_dir, 'hexahedron', 'lu')
            r = benchmark_hexahedra(n_div, 0, output_dir)
            if r:
                results_lu.append(r)

        if run_bicgstab:
            output_dir = os.path.join(script_dir, 'hexahedron', 'bicgstab')
            r = benchmark_hexahedra(n_div, 1, output_dir)
            if r:
                results_bicgstab.append(r)

    # Summary
    print('=' * 70)
    print('SUMMARY')
    print('=' * 70)

    if results_lu:
        print('\nLU Solver (hexahedron/lu/):\n')
        print('%-10s %10s %10s %8s %12s %10s' % ('n_div', 'Elements', 'Time (s)', 'Iter', 'M_avg_z', 'Conv'))
        print('-' * 65)
        for r in results_lu:
            print('%-10d %10d %10.3f %8d %12.0f %10s' % (
                r['n_div'], r['n_elements'], r['t_solve'],
                r['nonl_iterations'], r['M_avg_z'],
                'Yes' if r['converged'] else 'No'))

    if results_bicgstab:
        print('\nBiCGSTAB Solver (hexahedron/bicgstab/):\n')
        print('%-10s %10s %10s %8s %12s %10s' % ('n_div', 'Elements', 'Time (s)', 'Iter', 'M_avg_z', 'Conv'))
        print('-' * 65)
        for r in results_bicgstab:
            print('%-10d %10d %10.3f %8d %12.0f %10s' % (
                r['n_div'], r['n_elements'], r['t_solve'],
                r['nonl_iterations'], r['M_avg_z'],
                'Yes' if r['converged'] else 'No'))

    print('=' * 70)


if __name__ == '__main__':
    main()
