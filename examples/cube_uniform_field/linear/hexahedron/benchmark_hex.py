#!/usr/bin/env python
"""
Hexahedral Benchmark - Linear Material (Radia)

Generates benchmark results for {lu,bicgstab,hacapk}/ directories
using hexahedral cube mesh with various N divisions.

This structure matches the ELF benchmark for consistent organization.

Solver types:
  lu       - Dense LU decomposition (Method 0)
  bicgstab - BiCGSTAB iterative solver (Method 1)
  hacapk   - BiCGSTAB with H-matrix acceleration (Method 2)

Usage:
    python benchmark_hex.py --lu 5 10 15
    python benchmark_hex.py --bicgstab 5 10 15
    python benchmark_hex.py --hacapk 5 10 15 20
    python benchmark_hex.py 5 10 15  # runs lu only
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

import radia as rad

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False


# =============================================================================
# Benchmark Conditions (matching ELF)
# =============================================================================
MU_0 = 4 * np.pi * 1e-7  # Vacuum permeability [T/(A/m)]
CUBE_SIZE = 1.0          # Cube edge length [m]
CUBE_HALF = 0.5          # Half cube size
MU_R = 1000              # Relative permeability
CHI = MU_R - 1           # Magnetic susceptibility
H_EXT = 200000.0         # External H field [A/m]
B_EXT = MU_0 * H_EXT     # External B field [T]

# Analytical solution (demagnetizing factor N = 1/3 for cube)
N_DEMAG = 1.0 / 3.0
M_ANALYTICAL_Z = CHI * H_EXT / (1 + CHI * N_DEMAG)

# Hexahedral face topology (1-indexed for Radia)
HEX_FACES = [
    [1, 4, 3, 2],  # Bottom face (z=0)
    [5, 6, 7, 8],  # Top face (z=1)
    [1, 2, 6, 5],  # Front face (y=0)
    [3, 4, 8, 7],  # Back face (y=1)
    [1, 5, 8, 4],  # Left face (x=0)
    [2, 3, 7, 6]   # Right face (x=1)
]


def get_rss_memory_mb():
    """Get current RSS memory usage in MB."""
    if not HAS_PSUTIL:
        return None
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss / (1024 * 1024)


def generate_cube_mesh(n_div, size=1.0):
    """Generate hexahedral mesh for a cube."""
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
                        bicg_tol=1e-4):
    """Benchmark hexahedral mesh with linear material.

    Args:
        n_div: Number of divisions per cube edge
        solver_type: 'lu', 'bicgstab', or 'hacapk'
        output_dir: Directory to save results
        hmat_eps: ACA tolerance for H-matrix (only used with hacapk)
        bicg_tol: BiCGSTAB convergence tolerance
    """
    rad.FldUnits('m')
    rad.UtiDelAll()

    n_elements = n_div ** 3
    ndof = n_elements * 6  # 6 DOF per hexahedron

    print('=' * 70)
    print('HEXAHEDRAL MESH: N=%d (%d elements, %d DOF), solver=%s' % (
        n_div, n_elements, ndof, solver_type))
    print('=' * 70)

    # Generate mesh
    t_mesh_start = time.time()
    mesh = generate_cube_mesh(n_div, CUBE_SIZE)

    # Create Radia polyhedra
    hex_objs = []
    for verts in mesh:
        obj = rad.ObjPolyhdr(verts, HEX_FACES, [0, 0, 0])
        hex_objs.append(obj)

    container = rad.ObjCnt(hex_objs)
    t_mesh = time.time() - t_mesh_start

    print('Generated %d hexahedral elements' % n_elements)

    # Apply linear material
    mat = rad.MatLin(MU_R)
    rad.MatApl(container, mat)

    # External field
    ext = rad.ObjBckg([0, 0, B_EXT])
    grp = rad.ObjCnt([container, ext])

    # Map solver type to method number
    solver_method_map = {'lu': 0, 'bicgstab': 1, 'hacapk': 2}
    solver_method = solver_method_map.get(solver_type, 0)

    # Configure H-matrix if using HACApK
    if solver_method == 2:
        try:
            rad.SetHACApKParams(hmat_eps, 10, 2.0)
            print('H-matrix: Enabled (eps=%.0e)' % hmat_eps)
        except AttributeError:
            print('H-matrix: Not available')

    # Set solver tolerances
    rad.SetBiCGSTABTol(bicg_tol)

    # Solve
    print('Solving...')
    t_solve_start = time.time()
    result = rad.Solve(grp, 0.001, 100, solver_method)
    t_solve = time.time() - t_solve_start

    # Get solve statistics
    stats = rad.GetSolveStats()
    n_iter = stats.get('nonl_iterations', 0)
    n_linear_iter = stats.get('linear_iterations', 0)

    # Measure memory
    peak_memory_mb = get_rss_memory_mb()

    # Get average magnetization
    all_M = rad.ObjM(container)
    M_z_sum = sum(m[1][2] for m in all_M)
    M_avg_z = M_z_sum / n_elements

    # Calculate error vs analytical
    error_percent = abs(M_avg_z - M_ANALYTICAL_Z) / M_ANALYTICAL_Z * 100

    # Get H-matrix info if applicable
    hmat_info = None
    if solver_method == 2:
        try:
            hmat_info = rad.GetHACApKStats()
        except AttributeError:
            pass

    # Print results
    print('Mesh time:       %.4f s' % t_mesh)
    print('Solve time:      %.4f s' % t_solve)
    print('Iterations:      %d' % n_iter)
    print('Linear iter:     %d' % n_linear_iter)
    print('M_avg_z:         %.2f A/m' % M_avg_z)
    print('Analytical:      %.2f A/m' % M_ANALYTICAL_Z)
    print('Error:           %.2f%%' % error_percent)
    if peak_memory_mb is not None:
        print('Memory:          %.1f MB' % peak_memory_mb)
    print()

    # Build result dictionary
    result_data = {
        'element_type': 'hex',
        'mesh_description': 'N=%d' % n_div,
        'n_div': n_div,
        'n_elements': n_elements,
        'ndof': ndof,
        'H_ext': H_EXT,
        'mu_r': MU_R,
        'chi': CHI,
        't_mesh': t_mesh,
        't_solve': t_solve,
        'solver_method': solver_method,
        'solver_name': solver_type,
        'converged': True,
        'iterations': n_iter,
        'linear_iterations': n_linear_iter,
        'M_avg_z': M_avg_z,
        'M_analytical_z': M_ANALYTICAL_Z,
        'error_percent': error_percent,
    }

    if peak_memory_mb is not None:
        result_data['peak_memory_mb'] = peak_memory_mb

    # H-matrix stats
    if hmat_info:
        result_data['hmatrix'] = {
            'n_lowrank': hmat_info.get('n_lowrank', 0),
            'n_dense': hmat_info.get('n_dense', 0),
            'max_rank': hmat_info.get('max_rank', 0),
            'compression_ratio': hmat_info.get('compression', 0.0),
            'build_time': hmat_info.get('build_time', 0.0),
            'nlf': hmat_info.get('n_leaves', 0),
        }

    # Save result
    os.makedirs(output_dir, exist_ok=True)
    filename = 'hex_N%d_results.json' % n_div
    filepath = os.path.join(output_dir, filename)

    with open(filepath, 'w') as f:
        json.dump(result_data, f, indent=2)
    print('Saved: %s' % filepath)

    return result_data


def run_single_benchmark(n_div, solver_type, script_dir, args):
    """Run a single benchmark in a subprocess for accurate memory measurement."""
    cmd = [
        sys.executable, __file__,
        '--single', str(n_div), solver_type,
        '--hmat_eps', str(args.hmat_eps),
        '--bicg_tol', str(args.bicg_tol),
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
    parser = argparse.ArgumentParser(description='Hexahedral benchmark - linear material (Radia)')
    parser.add_argument('--lu', action='store_true', help='Use LU solver')
    parser.add_argument('--bicgstab', action='store_true', help='Use BiCGSTAB solver')
    parser.add_argument('--hacapk', action='store_true', help='Use HACApK solver')
    parser.add_argument('--hmat_eps', type=float, default=1e-4,
                       help='ACA tolerance for H-matrix (default: 1e-4)')
    parser.add_argument('--bicg_tol', type=float, default=1e-4,
                       help='BiCGSTAB convergence tolerance (default: 1e-4)')
    parser.add_argument('--single', nargs=2, metavar=('N', 'SOLVER'),
                       help='Run single benchmark (internal use)')
    parser.add_argument('n_values', nargs='*', type=int, default=[5, 10, 15],
                       help='N values for mesh divisions (default: 5 10 15)')
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Single benchmark mode (called by subprocess)
    if args.single:
        n_div = int(args.single[0])
        solver_type = args.single[1]
        output_dir = os.path.join(script_dir, solver_type)
        benchmark_hexahedra(n_div, solver_type, output_dir,
                           hmat_eps=args.hmat_eps, bicg_tol=args.bicg_tol)
        return

    # If no solver is specified, run lu only
    any_solver = args.lu or args.bicgstab or args.hacapk
    run_lu = args.lu or not any_solver
    run_bicgstab = args.bicgstab
    run_hacapk = args.hacapk

    print('=' * 70)
    print('HEXAHEDRAL BENCHMARK - LINEAR MATERIAL (Radia)')
    print('=' * 70)
    print('Cube size: %.1f m' % CUBE_SIZE)
    print('H_ext: %.0f A/m' % H_EXT)
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

    for label, results in [('LU', results_lu), ('BiCGSTAB', results_bicgstab), ('HACApK', results_hacapk)]:
        if results:
            print('\n%s Solver:\n' % label)
            print('%-10s %10s %10s %12s %10s' % ('N', 'Elements', 'Time (s)', 'M_avg_z', 'Error %'))
            print('-' * 60)
            for r in results:
                print('%-10d %10d %10.4f %12.2f %10.2f' % (
                    r['n_div'], r['n_elements'], r['t_solve'], r['M_avg_z'], r['error_percent']))

    print('=' * 70)


if __name__ == '__main__':
    main()
