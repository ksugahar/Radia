#!/usr/bin/env python
"""
Hexahedral Benchmark for Radia

Generates benchmark results for {lu,bicgstab,hacapk}/ directories
using hexahedral cube mesh with various N divisions.

This structure matches the ELF benchmark for consistent organization.

Each benchmark is run in a separate subprocess to ensure accurate memory
measurement (memory is not shared between runs).

Solver types:
  lu       - Dense LU decomposition (Method 0)
  bicgstab - BiCGSTAB iterative solver (Method 1)
  hacapk   - BiCGSTAB with H-matrix acceleration (Method 2)

Usage:
    python benchmark_hex.py --lu 5 10
    python benchmark_hex.py --bicgstab 5 10
    python benchmark_hex.py --hacapk 5 10
    python benchmark_hex.py --hacapk --hmat_eps 1e-4 5
    python benchmark_hex.py 5 10  # runs lu only
"""

import sys
import os
import time
import argparse
import subprocess
import json

import numpy as np

# Add paths for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../../src/radia'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import radia as rad

from benchmark_common import (
    run_nonlinear_benchmark, print_summary,
    CUBE_SIZE, H_EXT, HEX_FACES, get_peak_memory_mb
)


def generate_cube_mesh(n_div: int, size: float = 1.0):
    """Generate cubic mesh with n_div divisions per edge"""
    vertices_list = []
    dx = size / n_div
    offset = size / 2

    for iz in range(n_div):
        for iy in range(n_div):
            for ix in range(n_div):
                x0 = ix * dx - offset
                y0 = iy * dx - offset
                z0 = iz * dx - offset
                verts = [
                    [x0, y0, z0],
                    [x0 + dx, y0, z0],
                    [x0 + dx, y0 + dx, z0],
                    [x0, y0 + dx, z0],
                    [x0, y0, z0 + dx],
                    [x0 + dx, y0, z0 + dx],
                    [x0 + dx, y0 + dx, z0 + dx],
                    [x0, y0 + dx, z0 + dx]
                ]
                vertices_list.append(verts)

    return vertices_list


def benchmark_hexahedra(n_div, solver_type, output_dir, hmat_eps=1e-4,
                        bicg_tol=1e-4, nonl_tol=0.001):
    """Benchmark hexahedral mesh with Radia.

    Args:
        n_div: Number of divisions per cube edge
        solver_type: 'lu', 'bicgstab', or 'hacapk'
        output_dir: Directory to save results
        hmat_eps: ACA tolerance for H-matrix (only used with hacapk)
        bicg_tol: BiCGSTAB convergence tolerance
        nonl_tol: Nonlinear iteration convergence tolerance
    """
    rad.FldUnits('m')
    rad.UtiDelAll()

    n_elements = n_div ** 3

    print('=' * 70)
    print('HEXAHEDRAL MESH: N=%d (%d elements), solver=%s' % (n_div, n_elements, solver_type))
    print('=' * 70)

    # Generate mesh
    t_mesh_start = time.time()
    mesh = generate_cube_mesh(n_div, size=CUBE_SIZE)

    # Create Radia polyhedra
    hex_objs = []
    for verts in mesh:
        obj = rad.ObjPolyhdr(verts, HEX_FACES, [0, 0, 0])
        hex_objs.append(obj)

    container = rad.ObjCnt(hex_objs)
    t_mesh = time.time() - t_mesh_start

    print('Generated %d hexahedral elements' % n_elements)

    if not hex_objs:
        print('ERROR: No hexahedra created!')
        return None

    # Run benchmark using common function
    result = run_nonlinear_benchmark(
        radia_obj=container,
        n_elements=n_elements,
        solver_type=solver_type,
        output_dir=output_dir,
        element_type='hex',
        mesh_description='N=%d' % n_div,
        t_mesh=t_mesh,
        nonl_tol=nonl_tol,
        bicg_tol=bicg_tol,
        hmat_eps=hmat_eps,
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
    ]

    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        return None

    # Read result from JSON file
    output_dir = os.path.join(script_dir, solver_type)
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
    parser.add_argument('--hmat_eps', type=float, default=1e-4,
                       help='ACA tolerance for H-matrix (default: 1e-4)')
    parser.add_argument('--bicg_tol', type=float, default=1e-4,
                       help='BiCGSTAB convergence tolerance (default: 1e-4)')
    parser.add_argument('--nonl_tol', type=float, default=0.001,
                       help='Nonlinear iteration tolerance (default: 0.001)')
    parser.add_argument('--single', nargs=2, metavar=('N', 'SOLVER'),
                       help='Run single benchmark (internal use)')
    parser.add_argument('n_values', nargs='*', type=int, default=[5, 10],
                       help='N values for mesh divisions (default: 5 10)')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Single benchmark mode (called by subprocess)
    if args.single:
        n_div = int(args.single[0])
        solver_type = args.single[1]
        output_dir = os.path.join(script_dir, solver_type)
        benchmark_hexahedra(n_div, solver_type, output_dir,
                           hmat_eps=args.hmat_eps, bicg_tol=args.bicg_tol,
                           nonl_tol=args.nonl_tol)
        return

    # If no solver is specified, run lu only
    any_solver = args.lu or args.bicgstab or args.hacapk
    run_lu = args.lu or not any_solver
    run_bicgstab = args.bicgstab
    run_hacapk = args.hacapk

    print('=' * 70)
    print('HEXAHEDRAL BENCHMARK - Radia')
    print('=' * 70)
    print('Cube size: %.1f m' % CUBE_SIZE)
    print('H_ext: %.0f A/m' % H_EXT)
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
