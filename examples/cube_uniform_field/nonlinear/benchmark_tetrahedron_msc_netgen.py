#!/usr/bin/env python
"""
Tetrahedral Benchmark using Netgen mesh (MSC method)

Generates benchmark results for tetrahedron_msc/{lu,bicgstab,hacapk}/ directories
using Netgen mesh with various maxh values.

This script uses the same parameters as ELF_MAGIC for fair comparison.

Solver types:
  lu       - Dense LU decomposition (Method 0)
  bicgstab - BiCGSTAB iterative solver (Method 1)
  hacapk   - BiCGSTAB with H-matrix acceleration (Method 2)

Usage:
    python benchmark_tetrahedron_msc_netgen.py --lu 0.4 0.2 0.15 0.10
    python benchmark_tetrahedron_msc_netgen.py --bicgstab 0.4 0.2 0.15 0.10
    python benchmark_tetrahedron_msc_netgen.py --hacapk 0.4 0.2 0.15 0.10
    python benchmark_tetrahedron_msc_netgen.py --hacapk --eps 1e-4 0.15  # custom ACA tolerance

Author: Radia Development Team
Date: 2025-12-05
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

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # T/(A/m)

# Problem parameters
CUBE_SIZE = 1.0      # 1.0 m cube
CUBE_HALF = 0.5      # half size
H_EXT = 200000.0     # External field (A/m) - matches hexahedron benchmark
B_EXT = MU_0 * H_EXT  # External B field (T)

# B-H curve data (industry standard format)
# Format: [[H (A/m), B (T)], ...]
# Radia v1.3.14+ accepts B-H directly - internal conversion: M = B/mu_0 - H
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


def benchmark_tetrahedra(maxh, solver_method, output_dir, hacapk_eps=1e-4):
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

    if solver_method == 0:
        solver_name = 'lu'
    elif solver_method == 1:
        solver_name = 'bicgstab'
    else:
        solver_name = 'hacapk'

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

    # Apply nonlinear material (B-H curve input - industry standard)
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(cube, mat)

    # External field
    ext = rad.ObjBckg([0, 0, B_EXT])
    grp = rad.ObjCnt([cube, ext])

    # Configure H-matrix if using HACApK (method=2)
    hmatrix_enabled = False
    hmatrix_stats = None
    if solver_method == 2:
        try:
            rad.SetHACApKParams(hacapk_eps, 10, 2.0)
            hmatrix_enabled = True
            print('H-matrix: Enabled (eps=%.0e, leaf_size=10, eta=2.0)' % hacapk_eps)
        except AttributeError:
            print('H-matrix: Not available (API not found)')

    # Solve with custom convergence check
    # Radia's internal convergence for ObjPolyhdr has issues, so we track M_avg_z stability
    print('Solving...')
    t_solve_start = time.time()
    try:
        # Run iterations manually to check M_avg_z convergence
        converged = False
        n_iter = 0
        prev_M_z = 0
        max_iter = 100
        tol_rel = 0.001  # 0.1% relative change in M_avg_z

        for i in range(max_iter):
            result = rad.Solve(grp, 1e-10, 1, solver_method)  # 1 iteration at a time
            n_iter += 1

            all_M = rad.ObjM(cube)
            M_list = [m[1] for m in all_M]
            M_avg_z = np.mean([m[2] for m in M_list])

            if i > 0:
                rel_change = abs(M_avg_z - prev_M_z) / max(abs(M_avg_z), 1)
                if rel_change < tol_rel:
                    converged = True
                    break
            prev_M_z = M_avg_z

        t_solve = time.time() - t_solve_start

        # Measure peak memory after solve
        peak_memory_mb = get_peak_memory_mb()

        # Get H-matrix statistics if HACApK was used
        if hmatrix_enabled:
            try:
                hmatrix_stats = rad.GetHACApKStats()
                if hmatrix_stats:
                    print('H-matrix stats:')
                    print('  Leaves: %d (low-rank: %d, dense: %d)' % (
                        hmatrix_stats['n_leaves'], hmatrix_stats['n_lowrank'], hmatrix_stats['n_dense']))
                    print('  Max rank: %d' % hmatrix_stats['max_rank'])
                    print('  Compression: %.4f' % hmatrix_stats['compression'])
                    print('  Build time: %.4f s' % hmatrix_stats['build_time'])
            except AttributeError:
                pass

        residual = result[0] if result[0] else 0.0
    except Exception as e:
        print('Solve failed: %s' % e)
        return None

    print('Mesh time:    %.4f s' % t_mesh)
    print('Solve time:   %.3f s' % t_solve)
    print('Iterations:   %d' % n_iter)
    print('Converged:    %s' % ('Yes' if converged else 'No'))
    print('M_avg_z:      %.0f A/m' % M_avg_z)
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
        't_mesh': t_mesh,
        't_solve': t_solve,
        'solver_method': solver_method,
        'solver_name': solver_name,
        'hmatrix_enabled': hmatrix_enabled,
        'hacapk_eps': hacapk_eps if solver_method == 2 else None,
        'converged': converged,
        'residual': residual,
        'nonl_iterations': n_iter,
        'M_avg_z': M_avg_z,
    }
    if peak_memory_mb is not None:
        result_data['peak_memory_mb'] = peak_memory_mb
    if hmatrix_stats is not None:
        hmatrix_data = {
            'n_lowrank': hmatrix_stats['n_lowrank'],
            'n_dense': hmatrix_stats['n_dense'],
            'max_rank': hmatrix_stats['max_rank'],
            'compression_ratio': hmatrix_stats['compression'],
            'build_time': hmatrix_stats['build_time'],
            'nlf': hmatrix_stats['n_leaves'],
        }
        # Add timing statistics if available (v1.3.16+)
        if 't_hmatrix_build' in hmatrix_stats:
            hmatrix_data['t_hmatrix_build'] = hmatrix_stats['t_hmatrix_build']
        if 't_linear_solve' in hmatrix_stats:
            hmatrix_data['t_linear_solve'] = hmatrix_stats['t_linear_solve']
        if 'linear_iterations' in hmatrix_stats:
            hmatrix_data['linear_iterations'] = hmatrix_stats['linear_iterations']
        result_data['hmatrix'] = hmatrix_data

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
    parser.add_argument('--hacapk', action='store_true', help='Use HACApK solver (saves to tetrahedron_msc/hacapk/)')
    parser.add_argument('--eps', type=float, default=1e-4,
                       help='ACA tolerance for HACApK (default: 1e-4)')
    parser.add_argument('maxh_values', nargs='*', type=float, default=[0.4, 0.2, 0.15, 0.10],
                       help='maxh values for Netgen mesh (default: 0.4 0.2 0.15 0.10)')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # If no solver specified, run only LU
    any_solver = args.lu or args.bicgstab or args.hacapk
    run_lu = args.lu or not any_solver
    run_bicgstab = args.bicgstab
    run_hacapk = args.hacapk

    print('=' * 70)
    print('TETRAHEDRAL BENCHMARK (Netgen mesh)')
    print('=' * 70)
    print('Cube size: %.1f m' % CUBE_SIZE)
    print('H_ext: %.0f A/m' % H_EXT)
    print('maxh values: %s' % args.maxh_values)
    if run_hacapk:
        print('HACApK ACA eps: %.0e' % args.eps)
    print()

    results_lu = []
    results_bicgstab = []
    results_hacapk = []

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

        if run_hacapk:
            output_dir = os.path.join(script_dir, 'tetrahedron_msc', 'hacapk')
            r = benchmark_tetrahedra(maxh, 2, output_dir, hacapk_eps=args.eps)
            if r:
                results_hacapk.append(r)

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

    if results_hacapk:
        print('\nHACApK Solver (tetrahedron_msc/hacapk/):\n')
        print('%-12s %8s %10s %8s %12s %10s %8s %10s' % (
            'maxh', 'Elements', 'Time (s)', 'Iter', 'M_avg_z', 'Compress', 'Leaves', 'Conv'))
        print('-' * 90)
        for r in results_hacapk:
            hm = r.get('hmatrix', {})
            compression = hm.get('compression_ratio', 0.0)
            n_leaves = hm.get('nlf', 0)
            print('%-12s %8d %10.3f %8d %12.0f %10.4f %8d %10s' % (
                r['mesh_description'], r['n_elements'], r['t_solve'],
                r['nonl_iterations'], r['M_avg_z'],
                compression, n_leaves,
                'Yes' if r['converged'] else 'No'))

    print('=' * 70)


if __name__ == '__main__':
    main()
