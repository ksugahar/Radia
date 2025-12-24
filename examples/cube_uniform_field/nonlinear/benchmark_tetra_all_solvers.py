#!/usr/bin/env python
"""
Tetrahedral Benchmark - All Solvers Comparison (LU, BiCGSTAB, HACApK)

Benchmarks 3DOF tetrahedral elements with all three solver methods.
Uses Netgen mesh for consistent comparison with ELF_MAGIC.

Usage:
    python benchmark_tetra_all_solvers.py 0.4 0.3 0.25
    python benchmark_tetra_all_solvers.py --lu 0.3
    python benchmark_tetra_all_solvers.py --bicgstab 0.3
    python benchmark_tetra_all_solvers.py --hacapk 0.3
"""

import sys
import os
import time
import json
import argparse
import math

# Path setup
_src_path = os.path.join(os.path.dirname(__file__), '../../../src/radia')
sys.path.insert(0, _src_path)

import radia as rad

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

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
MU_0 = 4 * math.pi * 1e-7  # T/(A/m)

# Problem parameters - SAME as ELF_MAGIC for fair comparison
CUBE_SIZE = 1.0      # 1.0 m cube
CUBE_HALF = 0.5      # half size
H_EXT = 200000.0     # External field (A/m) - same as ELF_MAGIC
B_EXT = MU_0 * H_EXT  # External B field (T)

# B-H curve data (industry standard format) - SAME as ELF_MAGIC
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
    """Benchmark tetrahedral mesh with specified solver.

    Args:
        maxh: Maximum element size for Netgen mesh
        solver_method: 0=LU, 1=BiCGSTAB, 2=HACApK
        output_dir: Directory to save results
        hacapk_eps: ACA tolerance for HACApK solver
    """
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh
        from netgen_mesh_import import netgen_mesh_to_radia
    except ImportError as e:
        print('[SKIP] Tetrahedra: %s' % e)
        return None

    rad.FldUnits('m')
    rad.UtiDelAll()

    solver_names = {0: 'lu', 1: 'bicgstab', 2: 'hacapk'}
    solver_name = solver_names.get(solver_method, 'unknown')

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

    print('Generated %d tetrahedral elements (%d DOF)' % (n_elements, n_elements * 3))

    # Apply nonlinear material (B-H curve input - industry standard)
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(cube, mat)

    # External field
    ext = rad.ObjBckg([0, 0, B_EXT])
    grp = rad.ObjCnt([cube, ext])

    # Configure HACApK if using Method 2
    hmatrix_stats = None
    if solver_method == 2:
        try:
            rad.SetHACApKParams(hacapk_eps, 10, 2.0)
            print('H-matrix: Enabled (eps=%.0e, leaf_size=10, eta=2.0)' % hacapk_eps)
        except AttributeError:
            print('H-matrix: Not available (API not found)')

    # Solve
    print('Solving...')
    t_solve_start = time.time()
    try:
        result = rad.Solve(grp, 0.001, 1000, solver_method)
        t_solve = time.time() - t_solve_start

        # Measure peak memory after solve
        peak_memory_mb = get_peak_memory_mb()

        # Get H-matrix statistics if HACApK was used
        if solver_method == 2:
            try:
                hmatrix_stats = rad.GetHACApKStats()
                if hmatrix_stats:
                    print('H-matrix stats:')
                    print('  Leaf nodes: %d (low-rank: %d, dense: %d)' % (
                        hmatrix_stats['n_leaves'], hmatrix_stats['n_lowrank'], hmatrix_stats['n_dense']))
                    print('  Max rank: %d' % hmatrix_stats['max_rank'])
                    print('  Compression: %.4f' % hmatrix_stats['compression'])
                    print('  Build time: %.4f s' % hmatrix_stats['build_time'])
            except AttributeError:
                pass

        # Get magnetization
        all_M = rad.ObjM(cube)
        M_list = [m[1] for m in all_M]
        if HAS_NUMPY:
            M_avg_z = float(np.mean([m[2] for m in M_list])) if M_list else 0.0
        else:
            M_avg_z = sum(m[2] for m in M_list) / len(M_list) if M_list else 0.0

        n_iter = int(result[3]) if result[3] else 0
        converged = n_iter < 1000
        residual = result[0] if result[0] else 0.0

    except Exception as e:
        print('Solve failed: %s' % e)
        import traceback
        traceback.print_exc()
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
        'converged': converged,
        'residual': residual,
        'nonl_iterations': n_iter,
        'M_avg_z': M_avg_z,
    }
    if peak_memory_mb is not None:
        result_data['peak_memory_mb'] = peak_memory_mb
    if hmatrix_stats is not None:
        result_data['hmatrix'] = {
            'n_lowrank': hmatrix_stats['n_lowrank'],
            'n_dense': hmatrix_stats['n_dense'],
            'max_rank': hmatrix_stats['max_rank'],
            'compression_ratio': hmatrix_stats['compression'],
            'build_time': hmatrix_stats['build_time'],
            'nlf': hmatrix_stats['n_leaves'],
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
    parser = argparse.ArgumentParser(description='Tetrahedral benchmark - all solvers')
    parser.add_argument('--lu', action='store_true', help='Use LU solver only')
    parser.add_argument('--bicgstab', action='store_true', help='Use BiCGSTAB solver only')
    parser.add_argument('--hacapk', action='store_true', help='Use HACApK solver only')
    parser.add_argument('--eps', type=float, default=1e-4,
                       help='ACA tolerance for HACApK (default: 1e-4)')
    parser.add_argument('maxh_values', nargs='*', type=float, default=[0.3],
                       help='maxh values for Netgen mesh (default: 0.3)')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # If no solver is specified, run all
    any_solver = args.lu or args.bicgstab or args.hacapk
    run_lu = args.lu or not any_solver
    run_bicgstab = args.bicgstab or not any_solver
    run_hacapk = args.hacapk or not any_solver

    print('=' * 70)
    print('TETRAHEDRAL BENCHMARK (3DOF) - All Solvers')
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
        print('\nLU Solver:\n')
        print('%-12s %10s %10s %8s %12s %10s' % ('maxh', 'Elements', 'Time (s)', 'Iter', 'M_avg_z', 'Conv'))
        print('-' * 70)
        for r in results_lu:
            print('%-12s %10d %10.3f %8d %12.0f %10s' % (
                r['mesh_description'], r['n_elements'], r['t_solve'],
                r['nonl_iterations'], r['M_avg_z'],
                'Yes' if r['converged'] else 'No'))

    if results_bicgstab:
        print('\nBiCGSTAB Solver:\n')
        print('%-12s %10s %10s %8s %12s %10s' % ('maxh', 'Elements', 'Time (s)', 'Iter', 'M_avg_z', 'Conv'))
        print('-' * 70)
        for r in results_bicgstab:
            print('%-12s %10d %10.3f %8d %12.0f %10s' % (
                r['mesh_description'], r['n_elements'], r['t_solve'],
                r['nonl_iterations'], r['M_avg_z'],
                'Yes' if r['converged'] else 'No'))

    if results_hacapk:
        print('\nHACApK Solver:\n')
        print('%-12s %10s %10s %8s %12s %10s %8s' % (
            'maxh', 'Elements', 'Time (s)', 'Iter', 'M_avg_z', 'Compress', 'Conv'))
        print('-' * 80)
        for r in results_hacapk:
            hm = r.get('hmatrix', {})
            compression = hm.get('compression_ratio', 0.0)
            print('%-12s %10d %10.3f %8d %12.0f %10.4f %8s' % (
                r['mesh_description'], r['n_elements'], r['t_solve'],
                r['nonl_iterations'], r['M_avg_z'], compression,
                'Yes' if r['converged'] else 'No'))

    print('=' * 70)


if __name__ == '__main__':
    main()
