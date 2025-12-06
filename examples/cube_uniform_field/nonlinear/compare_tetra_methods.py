#!/usr/bin/env python
"""
Compare Tetrahedral Methods: Method 0 (Original/Gauss) vs Method 1 (Analytical)

This script compares the two tetrahedral field computation methods to verify
they produce equivalent results before removing Method 0 from the codebase.

Author: Radia Development Team
Date: 2025-12-06
"""

import sys
import os
import time

# Path setup
_build_path = os.path.join(os.path.dirname(__file__), '../../../build/Release')
_src_path = os.path.join(os.path.dirname(__file__), '../../../src/radia')
sys.path.insert(0, _build_path)
sys.path.append(_src_path)

import numpy as np
import radia as rad

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # T/(A/m)

# Problem parameters
CUBE_SIZE = 1.0      # 1.0 m cube
CUBE_HALF = 0.5      # half size
H_EXT = 50000.0      # External field (A/m)
B_EXT = MU_0 * H_EXT  # External B field (T)

# B-H curve data: [H (A/m), B (T)]
BH_DATA = [
    [0.0, 0.0],
    [100.0, 0.1],
    [200.0, 0.3],
    [500.0, 0.8],
    [1000.0, 1.2],
    [2000.0, 1.5],
    [5000.0, 1.7],
    [10000.0, 1.8],
    [50000.0, 2.0],
    [100000.0, 2.1],
]

# Convert B-H to H-M format: M = B/mu_0 - H
HM_DATA = [[h, b/MU_0 - h] for h, b in BH_DATA]


def benchmark_tetra_method(maxh, tetra_method, solver_method=0):
    """Benchmark tetrahedral mesh with specified method."""
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh
        from netgen_mesh_import import netgen_mesh_to_radia
    except ImportError as e:
        print('[SKIP] Tetrahedra: %s' % e)
        return None

    rad.FldUnits('m')
    rad.UtiDelAll()
    rad.SolverTetraMethod(tetra_method)

    method_name = 'Original (Gauss)' if tetra_method == 0 else 'Analytical'
    solver_name = 'LU' if solver_method == 0 else 'BiCGSTAB'

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

    # Apply nonlinear material
    mat = rad.MatSatIsoTab(HM_DATA)
    rad.MatApl(cube, mat)

    # External field
    ext = rad.ObjBckg([0, 0, B_EXT])
    grp = rad.ObjCnt([cube, ext])

    # Solve
    t_solve_start = time.time()
    try:
        result = rad.Solve(grp, 0.001, 1000, solver_method)
        t_solve = time.time() - t_solve_start

        # Get magnetization
        all_M = rad.ObjM(cube)
        M_list = [m[1] for m in all_M]
        M_avg_z = np.mean([m[2] for m in M_list])

        n_iter = int(result[3]) if result[3] else 0
        converged = n_iter < 1000 and not np.isnan(M_avg_z)
        residual = result[0] if result[0] else 0.0
    except Exception as e:
        print('Solve failed: %s' % e)
        return None

    return {
        'tetra_method': tetra_method,
        'method_name': method_name,
        'solver_name': solver_name,
        'maxh': maxh,
        'n_elements': n_elements,
        't_mesh': t_mesh,
        't_solve': t_solve,
        'converged': converged,
        'n_iter': n_iter,
        'residual': residual,
        'M_avg_z': M_avg_z,
    }


def main():
    print('=' * 70)
    print('TETRAHEDRAL METHOD COMPARISON')
    print('Method 0: Original (Gauss Integration via Polygon faces)')
    print('Method 1: Analytical (Closed-form surface charge formula)')
    print('=' * 70)
    print()

    # Test with multiple mesh sizes
    maxh_values = [0.5, 0.4, 0.35, 0.3]

    results = []

    for maxh in maxh_values:
        print('-' * 70)
        print('Testing maxh = %.2f m' % maxh)
        print('-' * 70)

        # Method 0 (Original)
        print('\n[Method 0] Original (Gauss)...')
        r0 = benchmark_tetra_method(maxh, tetra_method=0, solver_method=0)
        if r0:
            results.append(r0)
            print('  Elements: %d' % r0['n_elements'])
            print('  Time: %.3f s' % r0['t_solve'])
            print('  Iterations: %d' % r0['n_iter'])
            print('  M_avg_z: %.0f A/m' % r0['M_avg_z'])
            print('  Converged: %s' % ('Yes' if r0['converged'] else 'No'))

        # Method 1 (Analytical)
        print('\n[Method 1] Analytical...')
        r1 = benchmark_tetra_method(maxh, tetra_method=1, solver_method=0)
        if r1:
            results.append(r1)
            print('  Elements: %d' % r1['n_elements'])
            print('  Time: %.3f s' % r1['t_solve'])
            print('  Iterations: %d' % r1['n_iter'])
            print('  M_avg_z: %.0f A/m' % r1['M_avg_z'])
            print('  Converged: %s' % ('Yes' if r1['converged'] else 'No'))

        # Compare
        if r0 and r1:
            diff_M = abs(r0['M_avg_z'] - r1['M_avg_z'])
            rel_diff = diff_M / max(abs(r0['M_avg_z']), 1e-10) * 100
            speedup = r0['t_solve'] / max(r1['t_solve'], 1e-10)

            print('\n  [Comparison]')
            print('  M_avg_z difference: %.0f A/m (%.2f%%)' % (diff_M, rel_diff))
            print('  Speedup (Method 1 vs 0): %.2fx' % speedup)

            if rel_diff < 1.0:
                print('  Status: MATCH (< 1%% difference)')
            else:
                print('  Status: MISMATCH (>= 1%% difference)')

    # Summary table
    print()
    print('=' * 70)
    print('SUMMARY TABLE')
    print('=' * 70)
    print()
    print('%-8s %-12s %-10s %10s %8s %12s %10s' % (
        'maxh', 'Method', 'Elements', 'Time (s)', 'Iter', 'M_avg_z', 'Conv'))
    print('-' * 70)

    for r in results:
        print('%-8.2f %-12s %-10d %10.3f %8d %12.0f %10s' % (
            r['maxh'], r['method_name'][:12], r['n_elements'],
            r['t_solve'], r['n_iter'], r['M_avg_z'],
            'Yes' if r['converged'] else 'No'))

    print('=' * 70)

    # Final recommendation
    print()
    print('RECOMMENDATION:')
    method0_results = [r for r in results if r['tetra_method'] == 0]
    method1_results = [r for r in results if r['tetra_method'] == 1]

    if method0_results and method1_results:
        # Check if all results match within 1%
        all_match = True
        for r0 in method0_results:
            for r1 in method1_results:
                if r0['maxh'] == r1['maxh']:
                    rel_diff = abs(r0['M_avg_z'] - r1['M_avg_z']) / max(abs(r0['M_avg_z']), 1e-10) * 100
                    if rel_diff >= 1.0:
                        all_match = False
                        break

        if all_match:
            print('  Both methods produce equivalent results (< 1%% difference).')
            print('  -> Safe to remove Method 0 (Original/Gauss) from repository.')
        else:
            print('  Methods produce different results.')
            print('  -> Further investigation needed before removing Method 0.')


if __name__ == '__main__':
    main()
