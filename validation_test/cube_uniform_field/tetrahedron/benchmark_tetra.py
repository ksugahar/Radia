#!/usr/bin/env python
"""
Tetrahedral Benchmark using Netgen mesh for Radia

Generates benchmark results for linear/ and nonlinear/ subdirectories
using Netgen mesh with various maxh values.

Solver types:
  lu       - Dense LU decomposition (Method 0)
  bicgstab - BiCGSTAB iterative solver (Method 1)
  hacapk   - BiCGSTAB with H-matrix acceleration (Method 2)

Usage:
    python benchmark_tetra.py --lu 0.4 0.3 0.25
    python benchmark_tetra.py --bicgstab 0.4 0.3 0.25
    python benchmark_tetra.py --hacapk 0.15
    python benchmark_tetra.py --linear --lu 0.3 0.25 0.20
    python benchmark_tetra.py --nonlinear --bicgstab 0.3 0.25
    python benchmark_tetra.py 0.4 0.3 0.25  # runs lu, nonlinear by default
"""

import sys
import os
import time
import argparse
import subprocess
import json

# Add parent directory for benchmark_common
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad
from benchmark_common import (
    run_benchmark, print_summary,
    CUBE_SIZE, CUBE_HALF, H_EXT, MU_R, M_ANALYTICAL_Z
)


def benchmark_tetrahedra(maxh, solver_type, output_dir, is_linear=False,
                         hmat_eps=1e-4, bicg_tol=1e-4, nonl_tol=0.001,
                         hmat_leaf_size=10, hmat_eta=2.0):
    """Benchmark tetrahedral mesh (Netgen + Radia).

    Args:
        maxh: Maximum element size for Netgen mesh
        solver_type: 'lu', 'bicgstab', or 'hacapk'
        output_dir: Directory to save results
        is_linear: True for linear material, False for nonlinear
        hmat_eps: ACA tolerance for H-matrix
        bicg_tol: BiCGSTAB convergence tolerance
        nonl_tol: Nonlinear iteration convergence tolerance
        hmat_leaf_size: H-matrix leaf size
        hmat_eta: H-matrix admissibility parameter
    """
    try:
        from netgen.occ import Box, Pnt, OCCGeometry
        from ngsolve import Mesh
        from ngsolve import TaskManager
        from netgen_mesh_import import netgen_mesh_to_radia
    except ImportError as e:
        print('[SKIP] Tetrahedra: %s' % e)
        return None

    rad.UtiDelAll()

    material_type = 'linear' if is_linear else 'nonlinear'

    print('=' * 70)
    print('TETRAHEDRAL MESH: maxh=%.2fm, solver=%s, material=%s' % (
        maxh, solver_type, material_type))
    print('=' * 70)

    # Create mesh with Netgen
    t_mesh_start = time.time()
    cube_solid = Box(Pnt(-CUBE_HALF, -CUBE_HALF, -CUBE_HALF),
                      Pnt(CUBE_HALF, CUBE_HALF, CUBE_HALF))
    cube_solid.mat('magnetic')
    geo = OCCGeometry(cube_solid)

    with TaskManager():
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

        # Run benchmark
        result = run_benchmark(
            radia_obj=cube,
            n_elements=n_elements,
            solver_type=solver_type,
            output_dir=output_dir,
            element_type='tetra',
            mesh_description='maxh=%.2fm' % maxh,
            t_mesh=t_mesh,
            is_linear=is_linear,
            nonl_tol=nonl_tol,
            bicg_tol=bicg_tol,
            hmat_eps=hmat_eps,
            hmat_leaf_size=hmat_leaf_size,
            hmat_eta=hmat_eta,
            extra_data={'maxh': maxh}
        )

        return result


def run_single_benchmark(maxh, solver_type, script_dir, args):
    """Run a single benchmark in a subprocess for accurate memory measurement."""
    cmd = [
        sys.executable, __file__,
        '--single', str(maxh), solver_type,
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
    maxh_str = ('%.2fm' % maxh).replace('.', '_')
    filename = 'tetra_maxh%s_results.json' % maxh_str
    filepath = os.path.join(output_dir, filename)
    if os.path.exists(filepath):
        with open(filepath, 'r') as f:
            return json.load(f)
    return None


def main():
    parser = argparse.ArgumentParser(description='Tetrahedral benchmark using Netgen mesh (Radia)')
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
    parser.add_argument('--single', nargs=2, metavar=('MAXH', 'SOLVER'),
                       help='Run single benchmark (internal use)')
    parser.add_argument('maxh_values', nargs='*', type=float, default=[0.4, 0.3, 0.25],
                       help='maxh values for Netgen mesh (default: 0.4 0.3 0.25)')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Default to nonlinear if neither specified
    if not args.linear and not args.nonlinear:
        args.nonlinear = True

    # Single benchmark mode (called by subprocess)
    if args.single:
        maxh = float(args.single[0])
        solver_type = args.single[1]
        material_dir = 'linear' if args.linear else 'nonlinear'
        output_dir = os.path.join(script_dir, material_dir, solver_type)
        benchmark_tetrahedra(maxh, solver_type, output_dir,
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
    print('TETRAHEDRAL BENCHMARK - %s MATERIAL (Netgen mesh, Radia)' % material_type)
    print('=' * 70)
    print('Cube size: %.1f m' % CUBE_SIZE)
    print('H_ext: %.0f A/m' % H_EXT)
    if args.linear:
        print('mu_r: %d' % MU_R)
        print('M_analytical: %.2f A/m' % M_ANALYTICAL_Z)
    print('maxh values: %s' % args.maxh_values)
    print()

    results_lu = []
    results_bicgstab = []
    results_hacapk = []

    for maxh in args.maxh_values:
        if run_lu:
            r = run_single_benchmark(maxh, 'lu', script_dir, args)
            if r:
                results_lu.append(r)

        if run_bicgstab:
            r = run_single_benchmark(maxh, 'bicgstab', script_dir, args)
            if r:
                results_bicgstab.append(r)

        if run_hacapk:
            r = run_single_benchmark(maxh, 'hacapk', script_dir, args)
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
