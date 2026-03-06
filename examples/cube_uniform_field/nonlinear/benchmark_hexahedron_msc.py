#!/usr/bin/env python
"""
Hexahedron MSC (Magnetic Surface Charge) Benchmark Script

Benchmarks hexahedral elements using ObjPolyhdr with 6 DOF (surface charges).
Radia internally uses 6 DOF per hexahedron (one sigma per face), matching
ELF_MAGIC's MSC implementation for hexahedral elements.

Note: Tetrahedra use 3 DOF (Mx, My, Mz), hexahedra use 6 DOF (sigma_1..sigma_6).

This script uses the same parameters as ELF_MAGIC for fair comparison.

Solver types:
  lu       - Dense LU decomposition (Method 0)
  bicgstab - BiCGSTAB iterative solver (Method 1)

Usage:
    python benchmark_hexahedron_msc.py --lu 5 10 15 20
    python benchmark_hexahedron_msc.py --bicgstab 5 10 15 20
    python benchmark_hexahedron_msc.py --hacapk 5 10 15 20
    python benchmark_hexahedron_msc.py --hacapk --eps 1e-3 5 10 15  # custom ACA tolerance

Examples:
    python benchmark_hexahedron_msc.py --lu 5 10 15 20
    python benchmark_hexahedron_msc.py --bicgstab 5 10 15 20
    python benchmark_hexahedron_msc.py --hacapk --eps 1e-4 5 10 15  # ELF-compatible default
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
H_EXT = 200000.0  # 200,000 A/m (same as ELF_MAGIC nonlinear benchmark)

# B-H curve data (industry standard format)
# Format: [[H (A/m), B (T)], ...]
# Radia MatSatIsoTab accepts B-H directly - internal conversion: M = B/mu_0 - H
# This is the EXACT same B-H curve used in ELF benchmark
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
# Standard brick/hexahedron with 8 vertices
# Each element has 6 DOF (surface charges sigma on each face)
HEX_FACES = [
    [1, 4, 3, 2],  # Bottom face (z=0)
    [5, 6, 7, 8],  # Top face (z=1)
    [1, 2, 6, 5],  # Front face (y=0)
    [3, 4, 8, 7],  # Back face (y=1)
    [1, 5, 8, 4],  # Left face (x=0)
    [2, 3, 7, 6]   # Right face (x=1)
]


def create_hexahedron_msc_mesh(n_div):
    """
    Create NxNxN hexahedral mesh using ObjPolyhdr (6 DOF MSC method).

    ObjPolyhdr with HEX_FACES creates hexahedral elements with 6 DOF
    (surface charges on each face), matching ELF_MAGIC's implementation.

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
                # Corner of this sub-cube (minimum x, y, z)
                x0 = ix * cube_size - 0.5
                y0 = iy * cube_size - 0.5
                z0 = iz * cube_size - 0.5

                # 8 vertices of hexahedron
                # Vertex ordering: bottom face (0-3), top face (4-7)
                vertices = [
                    [x0, y0, z0],                          # v0: bottom-front-left
                    [x0 + cube_size, y0, z0],              # v1: bottom-front-right
                    [x0 + cube_size, y0 + cube_size, z0],  # v2: bottom-back-right
                    [x0, y0 + cube_size, z0],              # v3: bottom-back-left
                    [x0, y0, z0 + cube_size],              # v4: top-front-left
                    [x0 + cube_size, y0, z0 + cube_size],  # v5: top-front-right
                    [x0 + cube_size, y0 + cube_size, z0 + cube_size],  # v6: top-back-right
                    [x0, y0 + cube_size, z0 + cube_size]   # v7: top-back-left
                ]

                # ObjPolyhdr(vertices, faces, magnetization)
                obj = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 0])
                elements.append(obj)

    return elements


def benchmark_hexahedron_msc(n_div, solver_method=1, use_hmatrix=False,
                             hmat_eps=1e-4, leaf_size=10, eta=2.0, relax=0.0):
    """
    Benchmark hexahedral MSC mesh with specified solver.

    Args:
        n_div: Number of divisions per edge
        solver_method: 0=LU, 1=BiCGSTAB, 2=HACApK
        use_hmatrix: Enable H-matrix acceleration (BiCGSTAB only)
        hmat_eps: ACA tolerance for HACApK solver (default: 1e-4)
        leaf_size: Minimum cluster size in elements (default: 10)
        eta: Admissibility parameter (default: 2.0)
    """
    rad.FldUnits('m')
    rad.UtiDelAll()

    # Determine solver name
    if solver_method == 0:
        solver_name = 'lu'
    elif solver_method == 2:
        solver_name = 'hacapk'
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

    # Apply nonlinear material (B-H curve input - industry standard)
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(container, mat)

    # External field
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([container, ext])

    # Configure H-matrix if using HACApK (method=2)
    hmatrix_enabled = False
    hmatrix_stats = None
    if solver_method == 2:
        try:
            # Set HACApK parameters
            # eps: ACA tolerance (lower = more accurate, higher = faster)
            # leaf_size: minimum cluster size in elements
            # eta: admissibility parameter
            rad.SetHACApKParams(hmat_eps, leaf_size, eta)
            hmatrix_enabled = True
            print(f'H-matrix: Enabled (eps={hmat_eps:.0e}, leaf_size={leaf_size}, eta={eta})')
        except AttributeError:
            print('H-matrix: Not available (API not found)')

    # Tolerance parameters (all explicitly set for reproducibility)
    NONL_TOL = 0.001     # Nonlinear iteration tolerance (outer loop)
    BICG_TOL = 1e-4      # BiCGSTAB inner loop tolerance (Method 1 and 2)
    MAX_ITER = 100       # Maximum nonlinear iterations

    # Set BiCGSTAB inner loop tolerance explicitly (ELF-compatible default)
    rad.SetBiCGSTABTol(BICG_TOL)

    # Set under-relaxation coefficient (0.0 = full step, 0.0-1.0 = under-relaxation)
    rad.SetRelaxParam(relax)

    # Solve
    relax_str = f', relax={relax:.1f}' if relax > 0 else ''
    print(f'Solving... (nonl_tol={NONL_TOL:.0e}, bicg_tol={BICG_TOL:.0e}{relax_str})')
    t_solve_start = time.time()
    result = rad.Solve(grp, NONL_TOL, MAX_ITER, solver_method)
    t_solve = time.time() - t_solve_start

    # Measure peak memory after solve
    peak_memory_mb = get_peak_memory_mb()

    # Get H-matrix statistics if HACApK was used
    if hmatrix_enabled:
        try:
            hmatrix_stats = rad.GetHACApKStats()
            if hmatrix_stats:
                print(f'H-matrix stats:')
                print(f'  Leaf nodes: {hmatrix_stats["n_leaves"]} (low-rank: {hmatrix_stats["n_lowrank"]}, dense: {hmatrix_stats["n_dense"]})')
                print(f'  Max rank: {hmatrix_stats["max_rank"]}')
                print(f'  Compression: {hmatrix_stats["compression"]:.4f}')
                print(f'  Build time: {hmatrix_stats["build_time"]:.4f} s')
        except AttributeError:
            pass

    # Get magnetization
    all_M = rad.ObjM(container)
    M_list = [m[1] for m in all_M]
    M_avg_z = np.mean([m[2] for m in M_list]) if M_list else 0.0

    n_iter = int(result[3]) if result[3] else 0
    converged = n_iter < MAX_ITER
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
        'ndof': len(elements) * 6,  # 6 DOF per hexahedron
        'H_ext': H_EXT,
        't_mesh': t_mesh,
        't_solve': t_solve,
        'solver_method': solver_method,
        'solver_name': solver_name,
        # Tolerance parameters (all explicitly recorded)
        'nonl_tol': NONL_TOL,        # Nonlinear iteration tolerance (outer loop)
        'bicg_tol': BICG_TOL,        # BiCGSTAB inner loop tolerance
        # H-matrix parameters (Method 2 only)
        'hmat_eps': hmat_eps if solver_method == 2 else None,
        'leaf_size': leaf_size if solver_method == 2 else None,
        'eta': eta if solver_method == 2 else None,
        'hmatrix_enabled': hmatrix_enabled,
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
    return result_data


def main():
    parser = argparse.ArgumentParser(description='Hexahedron MSC benchmark (Radia)')
    parser.add_argument('--lu', action='store_true', help='Run LU solver benchmark')
    parser.add_argument('--bicgstab', action='store_true', help='Run BiCGSTAB solver benchmark')
    parser.add_argument('--hacapk', action='store_true', help='Run HACApK (H-matrix) solver benchmark')
    parser.add_argument('--eps', type=float, default=1e-4,
                       help='ACA tolerance for HACApK (default: 1e-4, lower=accurate, higher=fast)')
    parser.add_argument('--leaf_size', type=int, default=10,
                       help='Minimum cluster size in elements (default: 10)')
    parser.add_argument('--eta', type=float, default=2.0,
                       help='Admissibility parameter (default: 2.0)')
    parser.add_argument('--relax', type=float, default=0.0,
                       help='Under-relaxation coefficient (default: 0.0 = full step, 0.0-1.0)')
    parser.add_argument('sizes', nargs='*', type=int, default=[5, 10, 15, 20],
                       help='Mesh sizes (N values, default: 5 10 15 20)')

    args = parser.parse_args()

    # Default to LU only if none specified
    any_solver = args.lu or args.bicgstab or args.hacapk
    run_lu = args.lu or not any_solver
    run_bicgstab = args.bicgstab
    run_hacapk = args.hacapk

    script_dir = os.path.dirname(os.path.abspath(__file__))

    print('=' * 70)
    print('HEXAHEDRAL BENCHMARK (MSC) - Radia')
    print('=' * 70)
    print('Cube size: 1.0 m')
    print('H_ext: %.0f A/m' % H_EXT)
    print('N values: %s' % args.sizes)
    if args.relax > 0:
        print('Relaxation: %.2f (under-relaxation)' % args.relax)
    if run_hacapk:
        print('HACApK: eps=%.0e, leaf_size=%d, eta=%.1f' % (args.eps, args.leaf_size, args.eta))
    print()

    results_lu = []
    results_bicgstab = []
    results_hacapk = []

    for n in args.sizes:
        # LU Benchmark
        if run_lu:
            print('\n' + '=' * 70)
            print('LU SOLVER: N=%d' % n)
            print('=' * 70 + '\n')

            output_dir = os.path.join(script_dir, 'hexahedron_msc', 'lu')
            os.makedirs(output_dir, exist_ok=True)

            result = benchmark_hexahedron_msc(n, solver_method=0, use_hmatrix=False, relax=args.relax)
            results_lu.append(result)

            filename = 'msc_N%d_results.json' % n
            filepath = os.path.join(output_dir, filename)
            with open(filepath, 'w') as f:
                json.dump(result, f, indent=2)
            print('Saved: %s\n' % filepath)

        # BiCGSTAB Benchmark
        if run_bicgstab:
            print('\n' + '=' * 70)
            print('BiCGSTAB SOLVER: N=%d' % n)
            print('=' * 70 + '\n')

            output_dir = os.path.join(script_dir, 'hexahedron_msc', 'bicgstab')
            os.makedirs(output_dir, exist_ok=True)

            result = benchmark_hexahedron_msc(n, solver_method=1, use_hmatrix=False, relax=args.relax)
            results_bicgstab.append(result)

            filename = 'msc_N%d_results.json' % n
            filepath = os.path.join(output_dir, filename)
            with open(filepath, 'w') as f:
                json.dump(result, f, indent=2)
            print('Saved: %s\n' % filepath)

        # HACApK Benchmark (Method 2)
        if run_hacapk:
            print('\n' + '=' * 70)
            print('HACApK SOLVER: N=%d (eps=%.0e, leaf=%d, eta=%.1f)' % (n, args.eps, args.leaf_size, args.eta))
            print('=' * 70 + '\n')

            output_dir = os.path.join(script_dir, 'hexahedron_msc', 'hacapk')
            os.makedirs(output_dir, exist_ok=True)

            result = benchmark_hexahedron_msc(n, solver_method=2, use_hmatrix=False,
                                              hmat_eps=args.eps, leaf_size=args.leaf_size, eta=args.eta, relax=args.relax)
            results_hacapk.append(result)

            filename = 'msc_N%d_results.json' % n
            filepath = os.path.join(output_dir, filename)
            with open(filepath, 'w') as f:
                json.dump(result, f, indent=2)
            print('Saved: %s\n' % filepath)

    # Summary
    print('\n' + '=' * 70)
    print('SUMMARY')
    print('=' * 70)

    if results_lu:
        print('\nLU Solver (hexahedron_msc/lu/):\n')
        print('%-10s %10s %10s %10s %12s %10s' % ('N', 'Elements', 'Time (s)', 'Nonl Iter', 'M_avg_z', 'Conv'))
        print('-' * 70)
        for r in results_lu:
            print('%-10d %10d %10.3f %10d %12.0f %10s' % (
                r['n_div'], r['n_elements'], r['t_solve'],
                r['nonl_iterations'], r['M_avg_z'],
                'Yes' if r['converged'] else 'No'))

    if results_bicgstab:
        print('\nBiCGSTAB Solver (hexahedron_msc/bicgstab/):\n')
        print('%-10s %10s %10s %10s %12s %10s' % ('N', 'Elements', 'Time (s)', 'Nonl Iter', 'M_avg_z', 'Conv'))
        print('-' * 70)
        for r in results_bicgstab:
            print('%-10d %10d %10.3f %10d %12.0f %10s' % (
                r['n_div'], r['n_elements'], r['t_solve'],
                r['nonl_iterations'], r['M_avg_z'],
                'Yes' if r['converged'] else 'No'))

    if results_hacapk:
        print('\nHACApK Solver (hexahedron_msc/hacapk/):\n')
        print('%-6s %8s %10s %8s %10s %10s %8s %10s' % (
            'N', 'Elements', 'Time (s)', 'Nonl It', 'M_avg_z', 'Compress', 'Leaves', 'Conv'))
        print('-' * 85)
        for r in results_hacapk:
            hm = r.get('hmatrix', {})
            compression = hm.get('compression_ratio', 0.0)
            n_leaves = hm.get('nlf', 0)
            print('%-6d %8d %10.3f %8d %10.0f %10.4f %8d %10s' % (
                r['n_div'], r['n_elements'], r['t_solve'],
                r['nonl_iterations'], r['M_avg_z'],
                compression, n_leaves,
                'Yes' if r['converged'] else 'No'))

    print('=' * 70)


if __name__ == '__main__':
    main()
