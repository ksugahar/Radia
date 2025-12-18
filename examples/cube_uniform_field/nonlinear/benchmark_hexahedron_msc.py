#!/usr/bin/env python
"""
Hexahedron MSC (Magnetic Surface Charge) Benchmark Script

Benchmarks hexahedral elements using ObjThckPgn (thick polygon) which implements
the MSC (Magnetic Surface Charge) method - equivalent to ELF_MAGIC's MSC implementation.

This script tests three solver methods:
- LU decomposition (Method 0)
- BiCGSTAB (Method 1)
- BiCGSTAB with H-matrix acceleration (Method 1 + HMatrix enabled)

Usage:
    python benchmark_hexahedron_msc.py [--lu] [--bicgstab] [--hmatrix] [--all] [N1 N2 ...]

Examples:
    python benchmark_hexahedron_msc.py --all 10 15 20
    python benchmark_hexahedron_msc.py --lu 10 15
    python benchmark_hexahedron_msc.py --bicgstab 10 15 20
"""

import sys
import os
import time
import json
import argparse
import numpy as np

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))
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

# B-H curve (same as ELF_MAGIC)
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

# Convert to [H, M] format for Radia
HM_DATA = [[h, b/MU_0 - h] for h, b in BH_DATA]


def create_hexahedron_msc_mesh(n_div):
    """
    Create NxNxN hexahedral mesh using ObjThckPgn (MSC method).

    ObjThckPgn creates a thick polygon (extruded polygon) which uses
    the Magnetic Surface Charge approach for field computation.

    Args:
        n_div: Number of divisions per edge

    Returns:
        List of Radia object indices
    """
    cube_size = 1.0 / n_div
    half_size = cube_size / 2.0
    elements = []

    for ix in range(n_div):
        for iy in range(n_div):
            for iz in range(n_div):
                # Center of this sub-cube
                cx = (ix + 0.5) * cube_size - 0.5
                cy = (iy + 0.5) * cube_size - 0.5
                cz = (iz + 0.5) * cube_size - 0.5

                # ObjThckPgn: thick polygon (extruded rectangle)
                # Arguments: center_xy, lx, ly, lz, magnetization
                # For a cube: square base (lx x ly) extruded along z (lz)
                # We need to create it at the correct z position
                z_base = cz - half_size
                polygon_vertices = [
                    [cx - half_size, cy - half_size],
                    [cx + half_size, cy - half_size],
                    [cx + half_size, cy + half_size],
                    [cx - half_size, cy + half_size],
                ]

                # ObjThckPgn(vertices_2d, thickness_z, z_base, magnetization)
                obj = rad.ObjThckPgn(z_base, cube_size, polygon_vertices, 'z', [0, 0, 0])
                elements.append(obj)

    return elements


def benchmark_hexahedron_msc(n_div, solver_method=1, use_hmatrix=False):
    """
    Benchmark hexahedral MSC mesh with specified solver.

    Args:
        n_div: Number of divisions per edge
        solver_method: 0=LU, 1=BiCGSTAB
        use_hmatrix: Enable H-matrix acceleration (BiCGSTAB only)
    """
    rad.FldUnits('m')
    rad.UtiDelAll()

    # Determine solver name
    if solver_method == 0:
        solver_name = 'lu'
    elif use_hmatrix:
        solver_name = 'bicgstab_hmatrix'
    else:
        solver_name = 'bicgstab'

    n_elements = n_div ** 3
    print(f'Hexahedron MSC Benchmark: N={n_div} ({n_elements} elements)')
    print(f'Solver: {solver_name}')
    print('=' * 60)

    # Create mesh
    t_mesh_start = time.time()
    elements = create_hexahedron_msc_mesh(n_div)
    t_mesh = time.time() - t_mesh_start
    print(f'Mesh creation: {t_mesh:.4f} s ({len(elements)} elements)')

    # Create container
    container = rad.ObjCnt(elements)

    # Apply nonlinear material
    mat = rad.MatSatIsoTab(HM_DATA)
    rad.MatApl(container, mat)

    # External field
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([container, ext])

    # Configure H-matrix if requested
    # Note: H-matrix API may vary - check if available
    hmatrix_enabled = False
    if use_hmatrix and solver_method == 1:
        try:
            rad.SolverHMatrixEnable()
            hmatrix_enabled = True
            print('H-matrix: Enabled')
        except AttributeError:
            print('H-matrix: Not available (API not found)')

    # Solve
    print(f'Solving...')
    t_solve_start = time.time()
    result = rad.Solve(grp, 0.001, 1000, solver_method)
    t_solve = time.time() - t_solve_start

    # Measure peak memory after solve
    peak_memory_mb = get_peak_memory_mb()

    # Disable H-matrix after solve
    if hmatrix_enabled:
        try:
            rad.SolverHMatrixDisable()
        except AttributeError:
            pass

    # Get magnetization
    all_M = rad.ObjM(container)
    M_list = [m[1] for m in all_M]
    M_avg_z = np.mean([m[2] for m in M_list]) if M_list else 0.0

    n_iter = int(result[3]) if result[3] else 0
    converged = n_iter < 1000
    residual = result[0] if result[0] else 0.0

    print(f'Time: {t_solve:.3f} s')
    print(f'Iterations: {n_iter}')
    print(f'Converged: {converged}')
    print(f'M_avg_z: {M_avg_z:.0f} A/m')
    if peak_memory_mb is not None:
        print(f'Peak memory: {peak_memory_mb:.1f} MB')
    print()

    result_data = {
        'element_type': 'hexahedron_msc',
        'mesh_description': f'{n_div}x{n_div}x{n_div}',
        'n_div': n_div,
        'n_elements': len(elements),
        'ndof': len(elements) * 3,
        'H_ext': H_EXT,
        't_mesh': t_mesh,
        't_solve': t_solve,
        'solver_method': solver_method,
        'solver_name': solver_name,
        'hmatrix_enabled': hmatrix_enabled,
        'converged': converged,
        'residual': residual,
        'nonl_iterations': n_iter,
        'M_avg_z': M_avg_z,
    }
    if peak_memory_mb is not None:
        result_data['peak_memory_mb'] = peak_memory_mb
    return result_data


def main():
    parser = argparse.ArgumentParser(description='Hexahedron MSC benchmark')
    parser.add_argument('--lu', action='store_true', help='Run LU solver benchmark')
    parser.add_argument('--bicgstab', action='store_true', help='Run BiCGSTAB solver benchmark')
    parser.add_argument('--hmatrix', action='store_true', help='Run BiCGSTAB with H-matrix benchmark')
    parser.add_argument('--all', action='store_true', help='Run all solver benchmarks')
    parser.add_argument('sizes', nargs='*', type=int, default=[10, 15],
                       help='Mesh sizes (N values)')

    args = parser.parse_args()

    # Default to all if none specified
    if not (args.lu or args.bicgstab or args.hmatrix or args.all):
        args.all = True

    script_dir = os.path.dirname(os.path.abspath(__file__))

    all_results = []

    # LU Benchmark
    if args.lu or args.all:
        print('\n' + '=' * 70)
        print('LU SOLVER BENCHMARK')
        print('=' * 70 + '\n')

        output_dir = os.path.join(script_dir, 'hexahedron_msc', 'lu')
        os.makedirs(output_dir, exist_ok=True)

        for n in args.sizes:
            result = benchmark_hexahedron_msc(n, solver_method=0, use_hmatrix=False)
            all_results.append(result)

            filename = f'msc_N{n}_results.json'
            filepath = os.path.join(output_dir, filename)
            with open(filepath, 'w') as f:
                json.dump(result, f, indent=2)
            print(f'Saved: {filepath}\n')

    # BiCGSTAB Benchmark
    if args.bicgstab or args.all:
        print('\n' + '=' * 70)
        print('BiCGSTAB SOLVER BENCHMARK')
        print('=' * 70 + '\n')

        output_dir = os.path.join(script_dir, 'hexahedron_msc', 'bicgstab')
        os.makedirs(output_dir, exist_ok=True)

        for n in args.sizes:
            result = benchmark_hexahedron_msc(n, solver_method=1, use_hmatrix=False)
            all_results.append(result)

            filename = f'msc_N{n}_results.json'
            filepath = os.path.join(output_dir, filename)
            with open(filepath, 'w') as f:
                json.dump(result, f, indent=2)
            print(f'Saved: {filepath}\n')

    # BiCGSTAB + H-matrix Benchmark
    if args.hmatrix or args.all:
        print('\n' + '=' * 70)
        print('BiCGSTAB + H-MATRIX SOLVER BENCHMARK')
        print('=' * 70 + '\n')

        output_dir = os.path.join(script_dir, 'hexahedron_msc', 'hmatrix')
        os.makedirs(output_dir, exist_ok=True)

        for n in args.sizes:
            result = benchmark_hexahedron_msc(n, solver_method=1, use_hmatrix=True)
            all_results.append(result)

            filename = f'msc_N{n}_results.json'
            filepath = os.path.join(output_dir, filename)
            with open(filepath, 'w') as f:
                json.dump(result, f, indent=2)
            print(f'Saved: {filepath}\n')

    # Summary
    if all_results:
        print('\n' + '=' * 70)
        print('SUMMARY')
        print('=' * 70)
        print(f'{"Solver":<20} {"N":<5} {"Elements":<10} {"Time (s)":<12} {"Iter":<8} {"M_avg_z (A/m)":<15} {"Conv"}')
        print('-' * 80)
        for r in all_results:
            print(f'{r["solver_name"]:<20} {r["n_div"]:<5} {r["n_elements"]:<10} '
                  f'{r["t_solve"]:<12.3f} {r["nonl_iterations"]:<8} {r["M_avg_z"]:<15.0f} '
                  f'{"Yes" if r["converged"] else "No"}')


if __name__ == '__main__':
    main()
