#!/usr/bin/env python
"""
Hexahedral Benchmark - Linear Material (Unified Conditions)

Benchmarks 6DOF MSC hexahedral elements with linear material,
comparing LU and BiCGSTAB solvers.

Author: Radia Development Team
Date: 2025-12-20
"""

import sys
import os
import json
import time
import numpy as np

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))
import radia as rad

# Import unified benchmark conditions
from benchmark_conditions import (
    MU_0, CUBE_SIZE, CUBE_HALF, H_EXT, B_EXT, CHI, MU_R,
    SOLVER_TOLERANCE, MAX_ITERATIONS, M_ANALYTICAL_Z, STANDARD_N_DIVS
)

# Hexahedral face topology (1-indexed for Radia)
HEX_FACES = [
    [1, 4, 3, 2],  # Bottom face (z=0)
    [5, 6, 7, 8],  # Top face (z=1)
    [1, 2, 6, 5],  # Front face (y=0)
    [3, 4, 8, 7],  # Back face (y=1)
    [1, 5, 8, 4],  # Left face (x=0)
    [2, 3, 7, 6]   # Right face (x=1)
]


def create_hex_mesh(n_div):
    """Create hexahedral mesh using ObjPolyhdr."""
    rad.FldUnits('m')
    rad.UtiDelAll()

    cube_size_elem = CUBE_SIZE / n_div
    elements = []

    for ix in range(n_div):
        for iy in range(n_div):
            for iz in range(n_div):
                x0 = ix * cube_size_elem - CUBE_HALF
                y0 = iy * cube_size_elem - CUBE_HALF
                z0 = iz * cube_size_elem - CUBE_HALF

                vertices = [
                    [x0, y0, z0],
                    [x0 + cube_size_elem, y0, z0],
                    [x0 + cube_size_elem, y0 + cube_size_elem, z0],
                    [x0, y0 + cube_size_elem, z0],
                    [x0, y0, z0 + cube_size_elem],
                    [x0 + cube_size_elem, y0, z0 + cube_size_elem],
                    [x0 + cube_size_elem, y0 + cube_size_elem, z0 + cube_size_elem],
                    [x0, y0 + cube_size_elem, z0 + cube_size_elem]
                ]

                obj = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 0])
                elements.append(obj)

    cube = rad.ObjCnt(elements)
    mat = rad.MatLin(MU_R)  # mu_r (industry standard)
    rad.MatApl(cube, mat)

    ext = rad.ObjBckg([0, 0, B_EXT])
    grp = rad.ObjCnt([cube, ext])

    return grp, cube, len(elements)


def get_avg_magnetization(cube):
    """Get average z-magnetization."""
    all_M = rad.ObjM(cube)
    # Handle single vs multiple elements
    if isinstance(all_M[0][0], (int, float)):
        return all_M[1][2]
    else:
        M_z_values = [m[1][2] for m in all_M]
        return np.mean(M_z_values)


def run_benchmark(n_div, solver_method, solver_name):
    """Run single benchmark case."""
    print('=' * 70)
    print('HEXAHEDRAL MESH: N=%d, solver=%s' % (n_div, solver_name))
    print('=' * 70)

    # Create mesh
    t0 = time.time()
    grp, cube, n_elements = create_hex_mesh(n_div)
    t_mesh = time.time() - t0

    print('Generated %d hexahedral elements' % n_elements)
    print('DOF: %d' % (n_elements * 6))

    # Solve
    print('Solving...')
    t0 = time.time()
    result = rad.Solve(grp, SOLVER_TOLERANCE, MAX_ITERATIONS, solver_method)
    t_solve = time.time() - t0

    # Get results
    M_avg_z = get_avg_magnetization(cube)
    n_iter = int(result[3])
    converged = result[0] == 0.0

    # Calculate error vs analytical
    error_pct = abs(M_avg_z - M_ANALYTICAL_Z) / M_ANALYTICAL_Z * 100

    # Print results
    print('Mesh time:    %.4f s' % t_mesh)
    print('Solve time:   %.3f s' % t_solve)
    print('Iterations:   %d' % n_iter)
    print('Converged:    %s' % ('Yes' if converged else 'No'))
    print('M_avg_z:      %.0f A/m' % M_avg_z)
    print('Analytical:   %.0f A/m' % M_ANALYTICAL_Z)
    print('Error:        %.2f%%' % error_pct)

    # Save results
    results = {
        'n_div': n_div,
        'n_elements': n_elements,
        'ndof': n_elements * 6,
        'solver': solver_name,
        't_mesh': t_mesh,
        't_solve': t_solve,
        'iterations': n_iter,
        'converged': converged,
        'M_avg_z': M_avg_z,
        'M_analytical': M_ANALYTICAL_Z,
        'error_pct': error_pct
    }

    # Create output directory
    out_dir = os.path.join(os.path.dirname(__file__), 'hexahedron_msc', solver_name)
    os.makedirs(out_dir, exist_ok=True)

    # Save JSON
    filename = 'hexa_n%d_results.json' % n_div
    filepath = os.path.join(out_dir, filename)
    with open(filepath, 'w') as f:
        json.dump(results, f, indent=2)
    print()
    print('Saved: %s' % filepath)

    return results


def main():
    print('=' * 70)
    print('HEXAHEDRAL BENCHMARK - LINEAR MATERIAL (Unified Conditions)')
    print('=' * 70)
    print('Cube size: %.1f m' % CUBE_SIZE)
    print('mu_r: %d, chi: %d' % (MU_R, CHI))
    print('H_ext: %.0f A/m' % H_EXT)
    print('Analytical M_z: %.0f A/m' % M_ANALYTICAL_Z)
    print('N values: %s' % STANDARD_N_DIVS)
    print()

    all_results = []

    # Parse command line arguments for solver selection
    solvers = [(0, 'lu'), (1, 'bicgstab')]
    if '--lu' in sys.argv:
        solvers = [(0, 'lu')]
    elif '--bicgstab' in sys.argv:
        solvers = [(1, 'bicgstab')]
    elif '--hacapk' in sys.argv:
        solvers = [(2, 'hacapk')]
    elif '--all' in sys.argv:
        solvers = [(0, 'lu'), (1, 'bicgstab'), (2, 'hacapk')]

    for n_div in STANDARD_N_DIVS:
        for solver_method, solver_name in solvers:
            try:
                result = run_benchmark(n_div, solver_method, solver_name)
                all_results.append(result)
            except Exception as e:
                print('ERROR: %s' % e)

    # Print summary
    print()
    print('=' * 70)
    print('SUMMARY')
    print('=' * 70)
    print('%5s  %6s  %10s  %10s  %8s  %10s  %8s' %
          ('N', 'Elems', 'DOF', 'Solver', 't_solve', 'M_avg_z', 'Error%'))
    print('-' * 70)

    for r in all_results:
        print('%5d  %6d  %10d  %10s  %8.3f  %10.0f  %8.2f' % (
            r['n_div'], r['n_elements'], r['ndof'], r['solver'],
            r['t_solve'], r['M_avg_z'], r['error_pct']))

    return 0


if __name__ == '__main__':
    sys.exit(main())
