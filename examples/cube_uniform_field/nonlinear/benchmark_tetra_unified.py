#!/usr/bin/env python
"""
Tetrahedral Benchmark using Unified Conditions

Generates benchmark results for tetrahedron/bicgstab/ directory
using Netgen mesh with various maxh values.

This script uses the same conditions as benchmark_hexa_unified.py
for fair comparison between tetrahedral and hexahedral elements.

NOTE: LU solver (Method 0) is NOT recommended for tetrahedral elements.
      BiCGSTAB (Method 1) provides more stable convergence.

Usage:
    python benchmark_tetra_unified.py --bicgstab 0.333 0.25 0.2 0.167
    python benchmark_tetra_unified.py 0.333 0.25 0.2 0.167  # default: bicgstab

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
    MU_0, CUBE_SIZE, CUBE_HALF, H_EXT, B_EXT, HM_DATA,
    SOLVER_TOLERANCE, MAX_ITERATIONS, maxh_to_n_div
)


def benchmark_tetrahedra(maxh, solver_method, output_dir):
    """Benchmark tetrahedral mesh (Netgen + ObjPolyhdr)."""
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh
        from netgen_mesh_import import netgen_mesh_to_radia
    except ImportError as e:
        print('[SKIP] Tetrahedra: %s' % e)
        return None

    rad.FldUnits('m')
    rad.UtiDelAll()
    rad.SolverTetraMethod(0)  # Original method

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
        'element_type': 'tetra',
        'mesh_description': 'maxh=%.3fm' % maxh,
        'maxh': maxh,
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
        # For comparison with hexa: approximate equivalent n_div
        'equiv_n_div': maxh_to_n_div(maxh),
    }

    # Save result
    os.makedirs(output_dir, exist_ok=True)
    maxh_str = ('%.3f' % maxh).replace('.', '_')
    filename = 'tetra_maxh%s_results.json' % maxh_str
    filepath = os.path.join(output_dir, filename)

    with open(filepath, 'w') as f:
        json.dump(result_data, f, indent=2)
    print('Saved: %s' % filepath)

    return result_data


def main():
    parser = argparse.ArgumentParser(description='Tetrahedral benchmark using unified conditions')
    parser.add_argument('--bicgstab', action='store_true', help='Use BiCGSTAB solver (default, saves to tetrahedron/bicgstab/)')
    parser.add_argument('maxh_values', nargs='*', type=float, default=[0.333, 0.25, 0.2, 0.167],
                       help='maxh values for Netgen mesh (default: 0.333 0.25 0.2 0.167)')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Always use BiCGSTAB (Method 1) for tetrahedra
    # LU (Method 0) is not recommended for tetrahedral elements
    run_bicgstab = True

    print('=' * 70)
    print('TETRAHEDRAL BENCHMARK (Unified Conditions)')
    print('=' * 70)
    print('Cube size: %.1f m' % CUBE_SIZE)
    print('H_ext: %.0f A/m' % H_EXT)
    print('maxh values: %s' % args.maxh_values)
    print()
    print('Mesh correspondence (maxh -> approx n_div):')
    for maxh in args.maxh_values:
        n_div = maxh_to_n_div(maxh)
        print('  maxh=%.3fm -> n_div ~= %d, hexa_elements=%d' % (maxh, n_div, n_div**3))
    print()

    results_bicgstab = []

    for maxh in args.maxh_values:
        if run_bicgstab:
            output_dir = os.path.join(script_dir, 'tetrahedron', 'bicgstab')
            r = benchmark_tetrahedra(maxh, 1, output_dir)
            if r:
                results_bicgstab.append(r)

    # Summary
    print('=' * 70)
    print('SUMMARY')
    print('=' * 70)

    if results_bicgstab:
        print('\nBiCGSTAB Solver (tetrahedron/bicgstab/):\n')
        print('%-12s %10s %10s %8s %12s %10s' % ('maxh', 'Elements', 'Time (s)', 'Iter', 'M_avg_z', 'Conv'))
        print('-' * 65)
        for r in results_bicgstab:
            print('maxh=%.3fm %10d %10.3f %8d %12.0f %10s' % (
                r['maxh'], r['n_elements'], r['t_solve'],
                r['nonl_iterations'], r['M_avg_z'],
                'Yes' if r['converged'] else 'No'))

    print('=' * 70)


if __name__ == '__main__':
    main()
