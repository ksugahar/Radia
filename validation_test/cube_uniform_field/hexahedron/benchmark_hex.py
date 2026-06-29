#!/usr/bin/env python
"""
Hexahedral Benchmark for Radia

Generates benchmark results for linear/ and nonlinear/ subdirectories
using hexahedral cube mesh with various N divisions.

Solver types:
  lu       - Dense LU decomposition (Method 0)
  bicgstab - BiCGSTAB iterative solver (Method 1)
  hacapk   - BiCGSTAB with H-matrix acceleration (Method 2)

Usage:
    python benchmark_hex.py --lu 5 10 15
    python benchmark_hex.py --bicgstab 5 10 15
    python benchmark_hex.py --hacapk 5 10 15 20
    python benchmark_hex.py --linear --lu 5 10 15
    python benchmark_hex.py --nonlinear --bicgstab 5 10 15
    python benchmark_hex.py 5 10 15  # runs lu, nonlinear by default
"""

import sys
import os
import time
import argparse
import subprocess
import json

# Add parent directory for benchmark_common
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src'))

import radia as rad
from benchmark_common import (
    run_benchmark, print_summary, generate_hex_mesh,
    CUBE_SIZE, H_EXT, MU_R, M_ANALYTICAL_Z
)


def benchmark_hexahedra(n_div, solver_type, output_dir, is_linear=False,
                        hmat_eps=1e-4, bicg_tol=1e-4, nonl_tol=0.001,
                        hmat_leaf_size=10, hmat_eta=2.0):
    """Benchmark hexahedral mesh.

    Args:
        n_div: Number of divisions per cube edge
        solver_type: 'lu', 'bicgstab', or 'hacapk'
        output_dir: Directory to save results
        is_linear: True for linear material, False for nonlinear
        hmat_eps: ACA tolerance for H-matrix
        bicg_tol: BiCGSTAB convergence tolerance
        nonl_tol: Nonlinear iteration tolerance
        hmat_leaf_size: H-matrix leaf size
        hmat_eta: H-matrix admissibility parameter
    """
    rad.UtiDelAll()

    n_elements = n_div ** 3
    ndof = n_elements * 6
    material_type = 'linear' if is_linear else 'nonlinear'

    print('=' * 70)
    print('HEXAHEDRAL MESH: N=%d (%d elements), solver=%s, material=%s' % (
        n_div, n_elements, solver_type, material_type))
    print('=' * 70)

    # Generate mesh
    t_mesh_start = time.time()
    mesh = generate_hex_mesh(n_div, CUBE_SIZE)

    # Create Radia hexahedra
    hex_objs = []
    for verts in mesh:
        obj = rad.ObjHexahedron(verts, [0, 0, 0])
        hex_objs.append(obj)

    container = rad.ObjCnt(hex_objs)
    t_mesh = time.time() - t_mesh_start

    print('Generated %d hexahedral elements' % n_elements)

    # Run benchmark
    result = run_benchmark(
        radia_obj=container,
        n_elements=n_elements,
        solver_type=solver_type,
        output_dir=output_dir,
        element_type='hex',
        mesh_description='N=%d' % n_div,
        t_mesh=t_mesh,
        is_linear=is_linear,
        nonl_tol=nonl_tol,
        bicg_tol=bicg_tol,
        hmat_eps=hmat_eps,
        hmat_leaf_size=hmat_leaf_size,
        hmat_eta=hmat_eta,
        extra_data={'n_div': n_div}
    )

    return result


def run_single_benchmark(n_div, solver_type, script_dir, args):
    """Run a single benchmark in a subprocess for accurate memory measurement."""
    cmd = [
        sys.executable, __file__,
        '--single', str(n_div), solver_type,
        '--hmat_eps', str(args.hmat_eps),
        '--bicg_tol', str(args.bicg_tol),
        '--nonl_tol', str(args.nonl_tol),
        '--hmat_leaf_size', str(args.hmat_leaf_size),
        '--hmat_eta', str(args.hmat_eta),
    ]
    if args.linear:
        cmd.append('--linear')
    else:
        cmd.append('--nonlinear')

    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        return None

    # Read result from JSON file
    material_dir = 'linear' if args.linear else 'nonlinear'
    output_dir = os.path.join(script_dir, material_dir, solver_type)
    filename = 'hex_N%d_results.json' % n_div
    filepath = os.path.join(output_dir, filename)
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            return json.load(f)
    return None


def main():
    parser = argparse.ArgumentParser(description='Hexahedral benchmark (Radia)')
    parser.add_argument('--lu', action='store_true', help='Use LU solver')
    parser.add_argument('--bicgstab', action='store_true', help='Use BiCGSTAB solver')
    parser.add_argument('--hacapk', action='store_true', help='Use HACApK solver')
    parser.add_argument('--linear', action='store_true', help='Use linear material (mu_r=1000)')
    parser.add_argument('--nonlinear', action='store_true', help='Use nonlinear material (BH curve)')
    parser.add_argument('--hmat_eps', type=float, default=1e-4,
                       help='ACA tolerance for H-matrix (default: 1e-4)')
    parser.add_argument('--bicg_tol', type=float, default=1e-4,
                       help='BiCGSTAB convergence tolerance (default: 1e-4)')
    parser.add_argument('--nonl_tol', type=float, default=0.001,
                       help='Nonlinear iteration tolerance (default: 0.001)')
    parser.add_argument('--hmat_leaf_size', type=int, default=10,
                       help='H-matrix leaf size (default: 10)')
    parser.add_argument('--hmat_eta', type=float, default=2.0,
                       help='H-matrix admissibility eta (default: 2.0)')
    parser.add_argument('--single', nargs=2, metavar=('N', 'SOLVER'),
                       help='Run single benchmark (internal use)')
    parser.add_argument('n_values', nargs='*', type=int, default=[5, 10, 15],
                       help='N values for mesh divisions (default: 5 10 15)')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Default to nonlinear if neither specified
    if not args.linear and not args.nonlinear:
        args.nonlinear = True

    # Single benchmark mode (called by subprocess)
    if args.single:
        n_div = int(args.single[0])
        solver_type = args.single[1]
        material_dir = 'linear' if args.linear else 'nonlinear'
        output_dir = os.path.join(script_dir, material_dir, solver_type)
        benchmark_hexahedra(n_div, solver_type, output_dir,
                           is_linear=args.linear,
                           hmat_eps=args.hmat_eps, bicg_tol=args.bicg_tol,
                           nonl_tol=args.nonl_tol,
                           hmat_leaf_size=args.hmat_leaf_size, hmat_eta=args.hmat_eta)
        return

    # If no solver is specified, run lu only
    any_solver = args.lu or args.bicgstab or args.hacapk
    run_lu = args.lu or not any_solver
    run_bicgstab = args.bicgstab
    run_hacapk = args.hacapk

    material_type = 'LINEAR' if args.linear else 'NONLINEAR'

    print('=' * 70)
    print('HEXAHEDRAL BENCHMARK - %s MATERIAL (Radia)' % material_type)
    print('=' * 70)
    print('Cube size: %.1f m' % CUBE_SIZE)
    print('H_ext: %.0f A/m' % H_EXT)
    if args.linear:
        print('mu_r: %d' % MU_R)
        print('M_analytical: %.2f A/m' % M_ANALYTICAL_Z)
    print('N values: %s' % args.n_values)
    print()

    results_lu = []
    results_bicgstab = []
    results_hacapk = []

    for n_div in args.n_values:
        if run_lu:
            r = run_single_benchmark(n_div, 'lu', script_dir, args)
            if r:
                results_lu.append(r)

        if run_bicgstab:
            r = run_single_benchmark(n_div, 'bicgstab', script_dir, args)
            if r:
                results_bicgstab.append(r)

        if run_hacapk:
            r = run_single_benchmark(n_div, 'hacapk', script_dir, args)
            if r:
                results_hacapk.append(r)

    # Summary
    print('=' * 70)
    print('SUMMARY')
    print('=' * 70)

    print_summary(results_lu, 'LU', 'lu')
    print_summary(results_bicgstab, 'BiCGSTAB', 'bicgstab')
    print_summary(results_hacapk, 'HACApK', 'hacapk')

    print('=' * 70)


if __name__ == '__main__':
    main()
