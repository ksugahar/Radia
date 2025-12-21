#!/usr/bin/env python
"""
Tetrahedral Benchmark using Unified Conditions (Linear Material)

Generates benchmark results for tetrahedron_msc/{lu,bicgstab}/ directories
using Netgen mesh with various maxh values.

This script uses the same conditions as benchmark_hexa_unified.py
for fair comparison between tetrahedral and hexahedral elements.

Usage:
    python benchmark_tetra_unified.py --lu 0.5 0.4 0.35 0.3 0.25
    python benchmark_tetra_unified.py --bicgstab 0.5 0.4 0.35 0.3 0.25
    python benchmark_tetra_unified.py 0.5 0.4 0.35 0.3 0.25  # runs both

Author: Radia Development Team
Date: 2025-12-12
"""

import sys
import os
import time
import json
import argparse

# Path setup
_src_path = os.path.join(os.path.dirname(__file__), '../../../src/radia')
sys.path.insert(0, _src_path)

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
    MU_0, CUBE_SIZE, CUBE_HALF, H_EXT, B_EXT, CHI, MU_R,
    SOLVER_TOLERANCE, MAX_ITERATIONS, maxh_to_n_div, M_ANALYTICAL_Z
)


def benchmark_tetrahedra(maxh, solver_method, output_dir):
    """Benchmark tetrahedral mesh (Netgen + ObjPolyhdr) with linear material."""
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh
        from netgen_mesh_import import netgen_mesh_to_radia
    except ImportError as e:
        print('[SKIP] Tetrahedra: %s' % e)
        return None

    rad.FldUnits('m')
    rad.UtiDelAll()

    solver_name = 'lu' if solver_method == 0 else 'bicgstab'

    print('=' * 70)
    print('TETRAHEDRAL MESH: maxh=%.3fm, solver=%s' % (maxh, solver_name))
    print('=' * 70)

    # Create mesh with Netgen
    t_mesh_start = time.time()
    cube_solid = Box(Pnt(-CUBE_HALF, -CUBE_HALF, -CUBE_HALF),
                      Pnt(CUBE_HALF, CUBE_HALF, CUBE_HALF))
    cube_solid.mat('magnetic')
    geo = OCCGeometry(cube_solid)

    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)
    n_elements = mesh.ne

    cube = netgen_mesh_to_radia(mesh,
                                 material={'magnetization': [0, 0, 0]},
                                 units='m',
                                 material_filter='magnetic',
                                 verbose=False)
    t_mesh = time.time() - t_mesh_start

    print('Generated %d tetrahedral elements' % n_elements)

    # Apply linear material (mu_r - industry standard)
    mat = rad.MatLin(MU_R)
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
        'element_type': 'tetra',
        'mesh_description': 'maxh=%.2fm' % maxh,
        'maxh': maxh,
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
        # For comparison with hexa: approximate equivalent n_div
        'equiv_n_div': maxh_to_n_div(maxh),
    }
    if peak_memory_mb is not None:
        result_data['peak_memory_mb'] = peak_memory_mb

    # Save result
    os.makedirs(output_dir, exist_ok=True)
    maxh_str = ('%.2f' % maxh).replace('.', '_')
    filename = 'tetra_maxh%s_results.json' % maxh_str
    filepath = os.path.join(output_dir, filename)

    with open(filepath, 'w') as f:
        json.dump(result_data, f, indent=2)
    print('Saved: %s' % filepath)

    return result_data


def main():
    parser = argparse.ArgumentParser(description='Tetrahedral benchmark using unified conditions (linear material)')
    parser.add_argument('--lu', action='store_true', help='Use LU solver (saves to tetrahedron_msc/lu/)')
    parser.add_argument('--bicgstab', action='store_true', help='Use BiCGSTAB solver (saves to tetrahedron_msc/bicgstab/)')
    parser.add_argument('maxh_values', nargs='*', type=float, default=[0.5, 0.4, 0.35, 0.3, 0.25],
                       help='maxh values for Netgen mesh (default: 0.5 0.4 0.35 0.3 0.25)')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # If neither --lu nor --bicgstab is specified, run both
    run_lu = args.lu or (not args.lu and not args.bicgstab)
    run_bicgstab = args.bicgstab or (not args.lu and not args.bicgstab)

    print('=' * 70)
    print('TETRAHEDRAL BENCHMARK - LINEAR MATERIAL (Unified Conditions)')
    print('=' * 70)
    print('Cube size: %.1f m' % CUBE_SIZE)
    print('mu_r: %d, chi: %d' % (MU_R, CHI))
    print('H_ext: %.0f A/m' % H_EXT)
    print('Analytical M_z: %.0f A/m' % M_ANALYTICAL_Z)
    print('maxh values: %s' % args.maxh_values)
    print()
    print('Mesh correspondence (maxh -> approx n_div):')
    for maxh in args.maxh_values:
        n_div = maxh_to_n_div(maxh)
        print('  maxh=%.2fm -> n_div ~= %d, hexa_elements=%d' % (maxh, n_div, n_div**3))
    print()

    results_lu = []
    results_bicgstab = []

    for maxh in args.maxh_values:
        if run_lu:
            output_dir = os.path.join(script_dir, 'tetrahedron_msc', 'lu')
            r = benchmark_tetrahedra(maxh, 0, output_dir)
            if r:
                results_lu.append(r)

        if run_bicgstab:
            output_dir = os.path.join(script_dir, 'tetrahedron_msc', 'bicgstab')
            r = benchmark_tetrahedra(maxh, 1, output_dir)
            if r:
                results_bicgstab.append(r)

    # Summary
    print('=' * 70)
    print('SUMMARY')
    print('=' * 70)

    if results_lu:
        print('\nLU Solver (tetrahedron_msc/lu/):\n')
        print('%-12s %10s %10s %8s %12s %10s %8s' % ('maxh', 'Elements', 'Time (s)', 'Iter', 'M_avg_z', 'Error', 'Conv'))
        print('-' * 75)
        for r in results_lu:
            print('maxh=%.2fm %10d %10.3f %8d %12.0f %9.2f%% %8s' % (
                r['maxh'], r['n_elements'], r['t_solve'],
                r['iterations'], r['M_avg_z'], r['error_percent'],
                'Yes' if r['converged'] else 'No'))

    if results_bicgstab:
        print('\nBiCGSTAB Solver (tetrahedron_msc/bicgstab/):\n')
        print('%-12s %10s %10s %8s %12s %10s %8s' % ('maxh', 'Elements', 'Time (s)', 'Iter', 'M_avg_z', 'Error', 'Conv'))
        print('-' * 75)
        for r in results_bicgstab:
            print('maxh=%.2fm %10d %10.3f %8d %12.0f %9.2f%% %8s' % (
                r['maxh'], r['n_elements'], r['t_solve'],
                r['iterations'], r['M_avg_z'], r['error_percent'],
                'Yes' if r['converged'] else 'No'))

    print('=' * 70)


if __name__ == '__main__':
    main()
