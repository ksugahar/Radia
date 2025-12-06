#!/usr/bin/env python
"""
Dense BiCGSTAB Benchmark for Linear Material

Benchmarks the Dense BiCGSTAB solver (Method 1) for linear magnetic materials.
No H-matrix acceleration - pure dense matrix operations.

Problem Setup:
- Geometry: 0.1m x 0.1m x 0.1m cube centered at origin
- Material: Linear isotropic, mu_r = 100 (chi = 99)
- External Field: B0 = 1.0 T (uniform along z-axis)
- Solver: Dense BiCGSTAB (Method 1)

Author: Radia Development Team
Date: 2025-12-05
"""
import sys
import os
import time
import json
import argparse

# Add build/Release for radia.pyd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))
import radia as rad
import numpy as np

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # T/(A/m)

# Problem parameters
CUBE_SIZE = 0.1   # 0.1 m cube (100 mm)
MU_R = 100        # Relative permeability
CHI = MU_R - 1    # Magnetic susceptibility = 99
B0_Z = 1.0        # External field B0 = 1 T

# Evaluation points (external, in meters)
EVAL_POINTS = [
    [0.0, 0.0, 0.06],   # z = 60mm
    [0.0, 0.0, 0.08],   # z = 80mm
    [0.0, 0.0, 0.10],   # z = 100mm
    [0.06, 0.0, 0.0],   # x = 60mm
    [0.08, 0.0, 0.0],   # x = 80mm
]


def run_benchmark(n_div, method=1, verbose=True):
    """
    Run Dense BiCGSTAB benchmark for given mesh subdivision.

    Parameters:
    - n_div: Number of subdivisions per axis (total elements = n_div^3)
    - method: Solver method (0=LU, 1=BiCGSTAB)
    - verbose: Print progress messages

    Returns: dict with benchmark results
    """
    rad.FldUnits('m')
    rad.UtiDelAll()

    n_elem = n_div ** 3
    if verbose:
        print(f"\n{'='*60}")
        print(f"Dense BiCGSTAB Benchmark: {n_div}x{n_div}x{n_div} = {n_elem} elements")
        print(f"{'='*60}")

    # Create mesh
    t_mesh_start = time.time()
    cube = rad.ObjRecMag([0, 0, 0], [CUBE_SIZE, CUBE_SIZE, CUBE_SIZE], [0, 0, 0])
    rad.ObjDivMag(cube, [n_div, n_div, n_div])
    t_mesh = time.time() - t_mesh_start

    if verbose:
        print(f"Mesh creation: {t_mesh:.4f} s")

    # Apply linear material
    mat = rad.MatLin(CHI)  # chi = mu_r - 1
    rad.MatApl(cube, mat)

    # External field
    ext_field = rad.ObjBckg([0, 0, B0_Z])
    grp = rad.ObjCnt([cube, ext_field])

    # Solve
    solver_name = "BiCGSTAB" if method == 1 else "LU"
    if verbose:
        print(f"Solving with {solver_name} (Method {method})...")

    t_solve_start = time.time()
    result = rad.Solve(grp, 0.0001, 1000, method)
    t_solve = time.time() - t_solve_start

    # Extract result
    residual = result[0] if result[0] else 0.0
    n_iter = int(result[3]) if result[3] else 0
    converged = n_iter < 1000

    if verbose:
        print(f"Solve time: {t_solve:.3f} s")
        print(f"Iterations: {n_iter}")
        print(f"Converged: {converged}")

    # Get magnetization
    all_M = rad.ObjM(cube)
    M_list = [m[1] for m in all_M]
    M_avg_z = np.mean([m[2] for m in M_list])

    if verbose:
        print(f"M_avg_z: {M_avg_z:.0f} A/m")

    # Evaluate external field
    B_ext = []
    for pt in EVAL_POINTS:
        B = rad.Fld(grp, 'b', pt)
        B_ext.append({
            'point': pt,
            'Bz': B[2],
            'Bz_pert_mT': (B[2] - B0_Z) * 1000  # Perturbation in mT
        })

    if verbose:
        print(f"\nExternal Field Evaluation:")
        for b in B_ext:
            print(f"  {b['point']}: Bz_pert = {b['Bz_pert_mT']:.3f} mT")

    return {
        'solver': f'dense_{solver_name.lower()}',
        'method': method,
        'n_div': n_div,
        'n_elements': n_elem,
        'ndof': n_elem * 3,
        'mu_r': MU_R,
        'chi': CHI,
        'B0_z': B0_Z,
        't_mesh': t_mesh,
        't_solve': t_solve,
        'converged': converged,
        'residual': residual,
        'iterations': n_iter,
        'M_avg_z': M_avg_z,
        'B_external': B_ext
    }


def main():
    parser = argparse.ArgumentParser(description='Dense BiCGSTAB Benchmark for Linear Material')
    parser.add_argument('n_div', type=int, nargs='?', default=8,
                        help='Mesh subdivision per axis (default: 8)')
    parser.add_argument('--method', type=int, choices=[0, 1], default=1,
                        help='Solver method: 0=LU, 1=BiCGSTAB (default: 1)')
    parser.add_argument('--all', action='store_true',
                        help='Run all standard mesh sizes (4, 6, 8, 10)')
    parser.add_argument('--output', type=str, default=None,
                        help='Output JSON file (default: auto-generated)')
    args = parser.parse_args()

    if args.all:
        mesh_sizes = [4, 6, 8, 10]
    else:
        mesh_sizes = [args.n_div]

    results = []
    for n in mesh_sizes:
        result = run_benchmark(n, method=args.method, verbose=True)
        results.append(result)

    # Save results
    if args.output:
        output_file = args.output
    else:
        solver_name = 'bicgstab' if args.method == 1 else 'lu'
        if args.all:
            output_file = f'radia_dense_{solver_name}_linear_results.json'
        else:
            output_file = f'radia_dense_{solver_name}_linear_N{args.n_div}_results.json'

    output_path = os.path.join(os.path.dirname(__file__), output_file)
    with open(output_path, 'w') as f:
        json.dump(results if len(results) > 1 else results[0], f, indent=2)
    print(f"\nResults saved to: {output_file}")

    # Summary table
    if len(results) > 1:
        print(f"\n{'='*70}")
        print("SUMMARY: Dense BiCGSTAB Linear Material Benchmark")
        print(f"{'='*70}")
        print(f"{'N':>4} {'Elements':>10} {'DOF':>10} {'Time':>10} {'Iter':>8} {'M_avg_z':>12}")
        print("-" * 70)
        for r in results:
            print(f"{r['n_div']:>4} {r['n_elements']:>10} {r['ndof']:>10} "
                  f"{r['t_solve']:>9.3f}s {r['iterations']:>8} {r['M_avg_z']:>11.0f}")


if __name__ == "__main__":
    main()
