#!/usr/bin/env python
"""
Tetrahedral Benchmark using Netgen mesh

Generates benchmark results for tetrahedron_msc/{lu,bicgstab}/ directories
using Netgen mesh with various maxh values.

Usage:
    python benchmark_tetra_netgen.py --lu 0.4 0.3 0.25
    python benchmark_tetra_netgen.py --bicgstab 0.4 0.3 0.25
    python benchmark_tetra_netgen.py 0.4 0.3 0.25  # runs both

Author: Radia Development Team
Date: 2025-12-05
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

    solver_name = 'lu' if solver_method == 0 else 'bicgstab'

    print('=' * 70)
    print('TETRAHEDRAL MESH: maxh=%.2fm, solver=%s' % (maxh, solver_name))
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

    print('Mesh time:    %.4f s' % t_mesh)
    print('Solve time:   %.3f s' % t_solve)
    print('Iterations:   %d' % n_iter)
    print('Converged:    %s' % ('Yes' if converged else 'No'))
    print('M_avg_z:      %.0f A/m' % M_avg_z)
    print()

    result_data = {
        'element_type': 'tetra',
        'mesh_description': 'maxh=%.2fm' % maxh,
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
    }

    # Save result
    os.makedirs(output_dir, exist_ok=True)
    maxh_str = str(maxh).replace('.', '_')
    filename = 'tetra_maxh%s_results.json' % maxh_str
    filepath = os.path.join(output_dir, filename)

    with open(filepath, 'w') as f:
        json.dump(result_data, f, indent=2)
    print('Saved: %s' % filepath)

    return result_data


def main():
    parser = argparse.ArgumentParser(description='Tetrahedral benchmark using Netgen mesh')
    parser.add_argument('--lu', action='store_true', help='Use LU solver (saves to tetrahedron_msc/lu/)')
    parser.add_argument('--bicgstab', action='store_true', help='Use BiCGSTAB solver (saves to tetrahedron_msc/bicgstab/)')
    parser.add_argument('maxh_values', nargs='*', type=float, default=[0.4, 0.3, 0.25],
                       help='maxh values for Netgen mesh (default: 0.4 0.3 0.25)')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # If neither --lu nor --bicgstab is specified, run both
    run_lu = args.lu or (not args.lu and not args.bicgstab)
    run_bicgstab = args.bicgstab or (not args.lu and not args.bicgstab)

    print('=' * 70)
    print('TETRAHEDRAL BENCHMARK (Netgen mesh)')
    print('=' * 70)
    print('Cube size: %.1f m' % CUBE_SIZE)
    print('H_ext: %.0f A/m' % H_EXT)
    print('maxh values: %s' % args.maxh_values)
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
        print('%-15s %10s %10s %8s %12s %10s' % ('maxh', 'Elements', 'Time (s)', 'Iter', 'M_avg_z', 'Conv'))
        print('-' * 70)
        for r in results_lu:
            print('%-15s %10d %10.3f %8d %12.0f %10s' % (
                r['mesh_description'], r['n_elements'], r['t_solve'],
                r['nonl_iterations'], r['M_avg_z'],
                'Yes' if r['converged'] else 'No'))

    if results_bicgstab:
        print('\nBiCGSTAB Solver (tetrahedron_msc/bicgstab/):\n')
        print('%-15s %10s %10s %8s %12s %10s' % ('maxh', 'Elements', 'Time (s)', 'Iter', 'M_avg_z', 'Conv'))
        print('-' * 70)
        for r in results_bicgstab:
            print('%-15s %10d %10.3f %8d %12.0f %10s' % (
                r['mesh_description'], r['n_elements'], r['t_solve'],
                r['nonl_iterations'], r['M_avg_z'],
                'Yes' if r['converged'] else 'No'))

    print('=' * 70)


if __name__ == '__main__':
    main()
