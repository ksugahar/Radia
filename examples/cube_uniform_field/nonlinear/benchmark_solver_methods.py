#!/usr/bin/env python
"""
Radia Solver Benchmark: LU vs BiCGSTAB vs BiCGSTAB+HACApK

Compares three solver methods on nonlinear magnetic cube problem.
Outputs JSON logs compatible with ELF_MAGIC benchmark format.

Problem Setup:
- Cube: 1.0m x 1.0m x 1.0m, centered at origin
- Material: Nonlinear BH curve (soft iron)
- External field: H_ext = 50,000 A/m (along z-axis)

Author: Radia Development Team
Date: 2025-12-04
"""
import sys
import os
import json
import time
import argparse

# Add build/Release first for radia.pyd
_build_path = os.path.join(os.path.dirname(__file__), '../../../build/Release')
_src_path = os.path.join(os.path.dirname(__file__), '../../../src/python')
sys.path.insert(0, _build_path)
sys.path.append(_src_path)

import numpy as np
import radia as rad

# Physical constants
MU_0 = 4 * np.pi * 1e-7  # T/(A/m)

# Problem parameters (match ELF_MAGIC benchmark)
CUBE_SIZE = 1.0      # 1.0 m cube
H_EXT = 50000.0      # 50,000 A/m external field

# B-H curve (soft iron with saturation)
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

# Convert to [H, M] format for MatSatIsoTab
HM_DATA = [[h, b/MU_0 - h] for h, b in BH_DATA]


def run_benchmark(n_div, method, use_hmatrix=False, hmat_eps=1e-4):
    """
    Run benchmark with specified solver method.

    Parameters:
    - n_div: Number of divisions per axis
    - method: Solver method (0=LU, 1=BiCGSTAB)
    - use_hmatrix: Enable H-matrix acceleration
    - hmat_eps: H-matrix tolerance

    Returns:
    - Dictionary with timing and results
    """
    rad.FldUnits('m')
    rad.UtiDelAll()

    # Create mesh
    t_mesh_start = time.time()
    cube = rad.ObjRecMag([0, 0, 0], [CUBE_SIZE, CUBE_SIZE, CUBE_SIZE], [0, 0, 0])
    rad.ObjDivMag(cube, [n_div, n_div, n_div])
    t_mesh = time.time() - t_mesh_start

    n_elements = n_div ** 3
    ndof = n_elements * 3

    # Apply nonlinear material
    mat = rad.MatSatIsoTab(HM_DATA)
    rad.MatApl(cube, mat)

    # External field (H = 50,000 A/m -> B = mu_0 * H = 0.0628 T for air)
    # ObjBckg takes B in Tesla
    B_ext = MU_0 * H_EXT
    background = rad.ObjBckg([0.0, 0.0, B_ext])
    container = rad.ObjCnt([cube, background])

    # Enable/disable H-matrix
    if use_hmatrix:
        rad.SolverHMatrixEnable(1, hmat_eps)
    else:
        rad.SolverHMatrixDisable()

    # Solve
    t_solve_start = time.time()
    res = rad.Solve(container, 0.001, 1000, method)
    t_solve = time.time() - t_solve_start

    # Get magnetization at center
    M_center = rad.ObjM(cube)
    M_center_z = M_center[0][1][2] if M_center else 0.0

    # Get average magnetization
    M_avg_x = 0.0
    M_avg_y = 0.0
    M_avg_z = 0.0

    # For hexahedral mesh, get average M
    try:
        # ObjM returns list of [[center], [M]] for each sub-element
        all_M = rad.ObjM(cube)
        if all_M and len(all_M) > 0:
            M_list = [m[1] for m in all_M]  # Extract M vectors
            M_avg_x = np.mean([m[0] for m in M_list])
            M_avg_y = np.mean([m[1] for m in M_list])
            M_avg_z = np.mean([m[2] for m in M_list])
    except:
        pass

    # Check convergence
    converged = res[3] < 1000  # Less than max iterations
    residual = res[0] if res[0] else 0.0
    n_iter = int(res[3]) if res[3] else 0

    result = {
        'n_div': n_div,
        'n_elements': n_elements,
        'ndof': ndof,
        'H_ext': H_EXT,
        't_mesh': t_mesh,
        't_solve': t_solve,
        'converged': converged,
        'residual': residual,
        'nonl_iterations': n_iter,
        'M_avg_x': M_avg_x,
        'M_avg_y': M_avg_y,
        'M_avg_z': M_avg_z,
        'M_center_z': M_center_z,
    }

    if use_hmatrix:
        result['hmat_eps'] = hmat_eps
        result['solver'] = 'bicgstab_hacapk'
    elif method == 0:
        result['solver'] = 'lu'
    else:
        result['solver'] = 'bicgstab'

    return result


def save_json(result, filename):
    """Save result to JSON file."""
    with open(filename, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {filename}")


def main():
    parser = argparse.ArgumentParser(description='Radia Solver Benchmark')
    parser.add_argument('n_div', type=int, nargs='?', default=10,
                        help='Number of divisions per axis (default: 10)')
    parser.add_argument('--all', action='store_true',
                        help='Run all three methods')
    parser.add_argument('--lu', action='store_true',
                        help='Run LU solver only')
    parser.add_argument('--bicgstab', action='store_true',
                        help='Run BiCGSTAB solver only')
    parser.add_argument('--hacapk', action='store_true',
                        help='Run BiCGSTAB+HACApK solver only')
    args = parser.parse_args()

    n_div = args.n_div

    # Default: run all methods
    run_lu = args.all or args.lu or (not args.bicgstab and not args.hacapk)
    run_bicgstab = args.all or args.bicgstab
    run_hacapk = args.all or args.hacapk

    print("=" * 70)
    print(f"Radia Solver Benchmark: N = {n_div} ({n_div**3} elements)")
    print("=" * 70)
    print(f"Cube: {CUBE_SIZE} m, H_ext = {H_EXT} A/m")
    print()

    output_dir = os.path.dirname(__file__)
    results = {}

    # 1. LU solver (Method 0)
    if run_lu:
        print(f"[1/3] LU Direct Solver (Method 0)...")
        try:
            r = run_benchmark(n_div, method=0, use_hmatrix=False)
            results['lu'] = r
            filename = os.path.join(output_dir, f'radia_lu_N{n_div}_results.json')
            save_json(r, filename)
            print(f"  Time: {r['t_solve']:.3f}s, Iter: {r['nonl_iterations']}, Mz_avg: {r['M_avg_z']:.0f} A/m")
        except Exception as e:
            print(f"  ERROR: {e}")

    # 2. BiCGSTAB solver (Method 1)
    if run_bicgstab:
        print(f"\n[2/3] BiCGSTAB Iterative Solver (Method 1)...")
        try:
            r = run_benchmark(n_div, method=1, use_hmatrix=False)
            results['bicgstab'] = r
            filename = os.path.join(output_dir, f'radia_bicgstab_N{n_div}_results.json')
            save_json(r, filename)
            print(f"  Time: {r['t_solve']:.3f}s, Iter: {r['nonl_iterations']}, Mz_avg: {r['M_avg_z']:.0f} A/m")
        except Exception as e:
            print(f"  ERROR: {e}")

    # 3. BiCGSTAB + HACApK (Method 1 with H-matrix)
    if run_hacapk:
        print(f"\n[3/3] BiCGSTAB + HACApK H-matrix (Method 1)...")
        try:
            r = run_benchmark(n_div, method=1, use_hmatrix=True, hmat_eps=1e-4)
            results['hacapk'] = r
            filename = os.path.join(output_dir, f'radia_hacapk_N{n_div}_results.json')
            save_json(r, filename)
            print(f"  Time: {r['t_solve']:.3f}s, Iter: {r['nonl_iterations']}, Mz_avg: {r['M_avg_z']:.0f} A/m")
        except Exception as e:
            print(f"  ERROR: {e}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print(f"\n{'Solver':<20} {'Time [s]':>12} {'Iterations':>12} {'Mz_avg [A/m]':>15}")
    print("-" * 60)

    for name, r in results.items():
        print(f"{name:<20} {r['t_solve']:>12.3f} {r['nonl_iterations']:>12} {r['M_avg_z']:>15.0f}")

    # Speedup comparison
    if 'lu' in results and 'hacapk' in results:
        speedup = results['lu']['t_solve'] / results['hacapk']['t_solve']
        print(f"\nHACApK speedup vs LU: {speedup:.2f}x")

    if 'bicgstab' in results and 'hacapk' in results:
        speedup = results['bicgstab']['t_solve'] / results['hacapk']['t_solve']
        print(f"HACApK speedup vs BiCGSTAB: {speedup:.2f}x")

    return results


if __name__ == "__main__":
    results = main()
