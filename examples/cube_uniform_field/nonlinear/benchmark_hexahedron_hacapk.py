#!/usr/bin/env python
"""
Hexahedron HACApK Benchmark Script

Benchmarks hexahedral elements using 6DOF MSC with HACApK (H-matrix) acceleration.
This matches ELF_MAGIC's bicgstab_hacapk configuration for fair comparison.

Solver types:
  lu              - Dense LU decomposition (Method 0)
  bicgstab        - BiCGSTAB iterative solver (Method 1)
  hacapk          - BiCGSTAB + HACApK H-matrix (Method 2)

Usage:
    python benchmark_hexahedron_hacapk.py --lu 5
    python benchmark_hexahedron_hacapk.py --bicgstab 5
    python benchmark_hexahedron_hacapk.py --hacapk 5
    python benchmark_hexahedron_hacapk.py 3 5  # runs all solvers

Note: HACApK is slow for compact geometry (single cube).
      H-matrix acceleration only helps for spatially separated objects.
"""

import sys
import os
import time
import json
import argparse
import numpy as np

# Add Radia to path
_src_path = os.path.join(os.path.dirname(__file__), '../../../src/radia')
sys.path.insert(0, _src_path)
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


MU_0 = 4 * np.pi * 1e-7
H_EXT = 50000.0  # 50,000 A/m (same as ELF_MAGIC)

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

# Hexahedral face topology (1-indexed for Radia) - 6 DOF MSC
HEX_FACES = [
    [1, 4, 3, 2],  # Bottom face (z=0)
    [5, 6, 7, 8],  # Top face (z=1)
    [1, 2, 6, 5],  # Front face (y=0)
    [3, 4, 8, 7],  # Back face (y=1)
    [1, 5, 8, 4],  # Left face (x=0)
    [2, 3, 7, 6]   # Right face (x=1)
]


def create_hexahedron_mesh(n_div):
    """
    Create NxNxN hexahedral mesh using ObjPolyhdr (6 DOF MSC method).

    Args:
        n_div: Number of divisions per edge

    Returns:
        List of Radia object indices
    """
    cube_size = 1.0 / n_div
    elements = []

    for ix in range(n_div):
        for iy in range(n_div):
            for iz in range(n_div):
                # Corner of this sub-cube
                x0 = ix * cube_size - 0.5
                y0 = iy * cube_size - 0.5
                z0 = iz * cube_size - 0.5

                # 8 vertices of hexahedron
                vertices = [
                    [x0, y0, z0],
                    [x0 + cube_size, y0, z0],
                    [x0 + cube_size, y0 + cube_size, z0],
                    [x0, y0 + cube_size, z0],
                    [x0, y0, z0 + cube_size],
                    [x0 + cube_size, y0, z0 + cube_size],
                    [x0 + cube_size, y0 + cube_size, z0 + cube_size],
                    [x0, y0 + cube_size, z0 + cube_size]
                ]

                obj = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 0])
                elements.append(obj)

    return elements


def benchmark_hexahedron(n_div, solver_method, solver_name):
    """
    Benchmark hexahedral mesh with specified solver.

    Args:
        n_div: Number of divisions per edge
        solver_method: 0=LU, 1=BiCGSTAB, 2=HACApK
        solver_name: 'lu', 'bicgstab', 'hacapk'
    """
    rad.FldUnits('m')
    rad.UtiDelAll()

    n_elements = n_div ** 3
    ndof = n_elements * 6  # 6 DOF per hexahedron

    print('=' * 70)
    print('HEXAHEDRAL MESH: N=%d (%d elements, %d DOF), solver=%s' % (
        n_div, n_elements, ndof, solver_name))
    print('=' * 70)

    # Create mesh
    t_mesh_start = time.time()
    elements = create_hexahedron_mesh(n_div)
    t_mesh = time.time() - t_mesh_start
    print('Mesh creation: %.4f s (%d elements)' % (t_mesh, len(elements)))

    # Create container
    container = rad.ObjCnt(elements)

    # Apply nonlinear material (B-H curve input - industry standard)
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(container, mat)

    # External field
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([container, ext])

    # Solve
    print('Solving...')
    t_solve_start = time.time()
    result = rad.Solve(grp, 0.001, 1000, solver_method)
    t_solve = time.time() - t_solve_start

    # Measure peak memory after solve
    peak_memory_mb = get_peak_memory_mb()

    # Get magnetization
    all_M = rad.ObjM(container)
    M_list = [m[1] for m in all_M]
    M_avg_z = np.mean([m[2] for m in M_list]) if M_list else 0.0

    n_iter = int(result[3]) if result[3] else 0
    converged = n_iter < 1000
    residual = result[0] if result[0] else 0.0

    print('Mesh time:       %.4f s' % t_mesh)
    print('Solve time:      %.3f s' % t_solve)
    print('Iterations:      %d' % n_iter)
    print('Converged:       %s' % ('Yes' if converged else 'No'))
    print('M_avg_z:         %.0f A/m' % M_avg_z)
    if peak_memory_mb is not None:
        print('Peak memory:     %.1f MB' % peak_memory_mb)
    print()

    result_data = {
        'element_type': 'hex',
        'mesh_description': 'N=%d' % n_div,
        'n_div': n_div,
        'n_elements': len(elements),
        'ndof': ndof,
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
    return result_data


def main():
    parser = argparse.ArgumentParser(description='Hexahedron HACApK benchmark (Radia)')
    parser.add_argument('--lu', action='store_true', help='Run LU solver benchmark')
    parser.add_argument('--bicgstab', action='store_true', help='Run BiCGSTAB solver benchmark')
    parser.add_argument('--hacapk', action='store_true', help='Run HACApK solver benchmark')
    parser.add_argument('sizes', nargs='*', type=int, default=[3, 5],
                       help='Mesh sizes (N values, default: 3 5)')

    args = parser.parse_args()

    # Default to all if none specified
    any_solver = args.lu or args.bicgstab or args.hacapk
    run_lu = args.lu or not any_solver
    run_bicgstab = args.bicgstab or not any_solver
    run_hacapk = args.hacapk or not any_solver

    script_dir = os.path.dirname(os.path.abspath(__file__))

    print('=' * 70)
    print('HEXAHEDRAL HACApK BENCHMARK - Radia')
    print('=' * 70)
    print('Cube size: 1.0 m')
    print('H_ext: %.0f A/m' % H_EXT)
    print('Material: Nonlinear (BH curve)')
    print('N values: %s' % args.sizes)
    print()

    results_lu = []
    results_bicgstab = []
    results_hacapk = []

    for n in args.sizes:
        # LU Benchmark
        if run_lu:
            output_dir = os.path.join(script_dir, 'hexahedron_msc', 'lu')
            os.makedirs(output_dir, exist_ok=True)

            result = benchmark_hexahedron(n, solver_method=0, solver_name='lu')
            results_lu.append(result)

            filename = 'msc_N%d_results.json' % n
            filepath = os.path.join(output_dir, filename)
            with open(filepath, 'w') as f:
                json.dump(result, f, indent=2)
            print('Saved: %s' % filepath)

        # BiCGSTAB Benchmark
        if run_bicgstab:
            output_dir = os.path.join(script_dir, 'hexahedron_msc', 'bicgstab')
            os.makedirs(output_dir, exist_ok=True)

            result = benchmark_hexahedron(n, solver_method=1, solver_name='bicgstab')
            results_bicgstab.append(result)

            filename = 'msc_N%d_results.json' % n
            filepath = os.path.join(output_dir, filename)
            with open(filepath, 'w') as f:
                json.dump(result, f, indent=2)
            print('Saved: %s' % filepath)

        # HACApK Benchmark
        if run_hacapk:
            output_dir = os.path.join(script_dir, 'hexahedron_msc', 'hacapk')
            os.makedirs(output_dir, exist_ok=True)

            result = benchmark_hexahedron(n, solver_method=2, solver_name='hacapk')
            results_hacapk.append(result)

            filename = 'msc_N%d_results.json' % n
            filepath = os.path.join(output_dir, filename)
            with open(filepath, 'w') as f:
                json.dump(result, f, indent=2)
            print('Saved: %s' % filepath)

    # Summary
    print()
    print('=' * 70)
    print('SUMMARY')
    print('=' * 70)

    if results_lu:
        print()
        print('LU Solver (hexahedron_msc/lu/):')
        print()
        print('%-8s %10s %8s %10s %10s %12s %8s' % (
            'N', 'Elements', 'DOF', 'Time (s)', 'Nonl Iter', 'M_avg_z', 'Conv'))
        print('-' * 70)
        for r in results_lu:
            print('%-8d %10d %8d %10.3f %10d %12.0f %8s' % (
                r['n_div'], r['n_elements'], r['ndof'], r['t_solve'],
                r['nonl_iterations'], r['M_avg_z'],
                'Yes' if r['converged'] else 'No'))

    if results_bicgstab:
        print()
        print('BiCGSTAB Solver (hexahedron_msc/bicgstab/):')
        print()
        print('%-8s %10s %8s %10s %10s %12s %8s' % (
            'N', 'Elements', 'DOF', 'Time (s)', 'Nonl Iter', 'M_avg_z', 'Conv'))
        print('-' * 70)
        for r in results_bicgstab:
            print('%-8d %10d %8d %10.3f %10d %12.0f %8s' % (
                r['n_div'], r['n_elements'], r['ndof'], r['t_solve'],
                r['nonl_iterations'], r['M_avg_z'],
                'Yes' if r['converged'] else 'No'))

    if results_hacapk:
        print()
        print('HACApK Solver (hexahedron_msc/hacapk/):')
        print()
        print('%-8s %10s %8s %10s %10s %12s %8s' % (
            'N', 'Elements', 'DOF', 'Time (s)', 'Nonl Iter', 'M_avg_z', 'Conv'))
        print('-' * 70)
        for r in results_hacapk:
            print('%-8d %10d %8d %10.3f %10d %12.0f %8s' % (
                r['n_div'], r['n_elements'], r['ndof'], r['t_solve'],
                r['nonl_iterations'], r['M_avg_z'],
                'Yes' if r['converged'] else 'No'))

    print('=' * 70)


if __name__ == '__main__':
    main()
