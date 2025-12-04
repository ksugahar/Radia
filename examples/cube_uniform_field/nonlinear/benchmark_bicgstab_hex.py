"""
Benchmark BiCGSTAB and BiCGSTAB + H-matrix solvers for hexahedral elements

This script benchmarks Radia's BiCGSTAB solver (Method 1) with and without
H-matrix acceleration for nonlinear material problems.

Output: JSON files for each configuration
"""

import sys
import os
import time
import json
import argparse

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))
import radia as rad
import numpy as np

MU_0 = 4 * np.pi * 1e-7
H_EXT = 50000.0  # 50,000 A/m external field

# B-H curve (same as ELF_MAGIC benchmark)
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

# Convert B-H to H-M format for Radia
HM_DATA = [[h, b / MU_0 - h] for h, b in BH_DATA]


def benchmark_bicgstab(n_div, use_hmatrix=False, hmat_eps=1e-4):
    """
    Run BiCGSTAB benchmark for hexahedral mesh.

    Args:
        n_div: Number of divisions per axis
        use_hmatrix: Enable H-matrix acceleration
        hmat_eps: H-matrix tolerance (only used if use_hmatrix=True)

    Returns:
        dict: Benchmark results
    """
    rad.FldUnits('m')
    rad.UtiDelAll()

    solver_name = f"bicgstab_hmatrix" if use_hmatrix else "bicgstab"
    n_elements = n_div ** 3
    ndof = n_elements * 3

    print(f"\n{'=' * 70}")
    print(f"Benchmark: {solver_name.upper()}, N={n_div} ({n_elements} elements, {ndof} DOF)")
    print(f"{'=' * 70}")

    # Create hexahedral mesh
    t_mesh_start = time.time()
    cube = rad.ObjRecMag([0, 0, 0], [1.0, 1.0, 1.0], [0, 0, 0])
    rad.ObjDivMag(cube, [n_div, n_div, n_div])
    t_mesh = time.time() - t_mesh_start
    print(f"Mesh creation: {t_mesh:.4f} s")

    # Apply nonlinear material
    mat = rad.MatSatIsoTab(HM_DATA)
    rad.MatApl(cube, mat)

    # External field
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([cube, ext])

    # Configure H-matrix
    if use_hmatrix:
        rad.SolverHMatrixEnable(1, hmat_eps)  # (enable=1, eps=hmat_eps)
        print(f"H-matrix enabled (eps={hmat_eps})")
    else:
        rad.SolverHMatrixDisable()
        print("H-matrix disabled (dense matrix)")

    # Solve with BiCGSTAB (Method 1)
    print("Solving with BiCGSTAB (Method 1)...")
    t_solve_start = time.time()
    result = rad.Solve(grp, 0.001, 1000, 1)  # Method 1 = BiCGSTAB
    t_solve = time.time() - t_solve_start

    # Extract results
    # result = [MaxModM, MaxModH, MaxMisfitAbs, Iterations]
    n_iter = int(result[3]) if len(result) > 3 and result[3] else 0
    converged = n_iter < 1000
    residual = result[0] if len(result) > 0 else 0.0

    # Get magnetization
    all_M = rad.ObjM(cube)
    M_list = [m[1] for m in all_M]
    M_avg_z = np.mean([m[2] for m in M_list])
    M_center_z = all_M[0][1][2] if all_M else 0.0

    print(f"\nResults:")
    print(f"  Time: {t_solve:.3f} s")
    print(f"  Iterations: {n_iter}")
    print(f"  Converged: {converged}")
    print(f"  Residual: {residual:.6e}")
    print(f"  M_avg_z: {M_avg_z:.1f} A/m")
    print(f"  M_center_z: {M_center_z:.1f} A/m")

    return {
        "solver": solver_name,
        "n_div": n_div,
        "n_elements": n_elements,
        "ndof": ndof,
        "hmat_eps": hmat_eps if use_hmatrix else None,
        "H_ext": H_EXT,
        "t_mesh": t_mesh,
        "t_solve": t_solve,
        "converged": converged,
        "residual": residual,
        "nonl_iterations": n_iter,
        "M_avg_z": M_avg_z,
        "M_center_z": M_center_z,
    }


def main():
    parser = argparse.ArgumentParser(description='Benchmark BiCGSTAB solver for hexahedral mesh')
    parser.add_argument('n_div', type=int, nargs='?', default=10,
                        help='Number of divisions per axis (default: 10)')
    parser.add_argument('--all', action='store_true',
                        help='Run all benchmarks (N=10, 15, 20)')
    parser.add_argument('--hmatrix-only', action='store_true',
                        help='Only run H-matrix benchmarks')
    parser.add_argument('--dense-only', action='store_true',
                        help='Only run dense matrix benchmarks')
    args = parser.parse_args()

    print("=" * 70)
    print("Radia BiCGSTAB Benchmark for Hexahedral Elements")
    print("=" * 70)
    print(f"External field H_ext = {H_EXT} A/m")
    print(f"Material: Nonlinear B-H curve (saturation ~2.1 T)")

    # Determine which N values to run
    if args.all:
        n_divs = [10, 15, 20]
    else:
        n_divs = [args.n_div]

    results = []
    output_dir = os.path.dirname(__file__)

    for n_div in n_divs:
        # Dense BiCGSTAB (no H-matrix)
        if not args.hmatrix_only:
            result_dense = benchmark_bicgstab(n_div, use_hmatrix=False)
            results.append(result_dense)

            # Save individual JSON
            json_file = os.path.join(output_dir, f"radia_bicgstab_N{n_div}_results.json")
            with open(json_file, 'w') as f:
                json.dump(result_dense, f, indent=2)
            print(f"Saved: {json_file}")

        # BiCGSTAB + H-matrix
        if not args.dense_only:
            result_hmat = benchmark_bicgstab(n_div, use_hmatrix=True, hmat_eps=1e-4)
            results.append(result_hmat)

            # Save individual JSON
            json_file = os.path.join(output_dir, f"radia_bicgstab_hmatrix_N{n_div}_results.json")
            with open(json_file, 'w') as f:
                json.dump(result_hmat, f, indent=2)
            print(f"Saved: {json_file}")

    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"{'Solver':<20} {'N':<5} {'Elements':<10} {'Time (s)':<12} {'Iter':<8} {'M_avg_z':<12}")
    print("-" * 70)
    for r in results:
        print(f"{r['solver']:<20} {r['n_div']:<5} {r['n_elements']:<10} {r['t_solve']:<12.3f} {r['nonl_iterations']:<8} {r['M_avg_z']:<12.1f}")

    # Save combined summary
    summary_file = os.path.join(output_dir, "radia_bicgstab_benchmark_summary.json")
    with open(summary_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSummary saved: {summary_file}")


if __name__ == '__main__':
    main()
