#!/usr/bin/env python
"""
Debug script for Hex N=10 HACApK comparison with ELF_MAGIC.

Runs detailed comparison and outputs per-iteration info.
"""

import sys
import os
import time
import math

_src_path = os.path.join(os.path.dirname(__file__), '../../../src/radia')
sys.path.insert(0, _src_path)

import radia as rad

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

MU_0 = 4 * math.pi * 1e-7
H_EXT = 200000.0
N_DIV = 10
HACAPK_EPS = 1e-6

# B-H curve data (same as ELF_MAGIC)
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
    """Create NxNxN hexahedral mesh."""
    cube_size = 1.0 / n_div
    elements = []

    for ix in range(n_div):
        for iy in range(n_div):
            for iz in range(n_div):
                x0 = ix * cube_size - 0.5
                y0 = iy * cube_size - 0.5
                z0 = iz * cube_size - 0.5

                vertices = [
                    [x0, y0, z0],
                    [x0 + cube_size, y0, z0],
                    [x0 + cube_size, y0 + cube_size, z0],
                    [x0, y0 + cube_size, z0],
                    [x0, y0, z0 + cube_size],
                    [x0 + cube_size, y0, z0 + cube_size],
                    [x0 + cube_size, y0 + cube_size, z0 + cube_size],
                    [x0, y0 + cube_size, z0 + cube_size]
                ]

                obj = rad.ObjPolyhdr(vertices, HEX_FACES, [0, 0, 0])
                elements.append(obj)

    return elements


def get_M_avg_z(container):
    """Get average Mz."""
    all_M = rad.ObjM(container)
    M_list = [m[1] for m in all_M]
    if HAS_NUMPY:
        return float(np.mean([m[2] for m in M_list]))
    else:
        return sum(m[2] for m in M_list) / len(M_list)


def run_solver(solver_method, solver_name, n_iterations=1):
    """Run solver and return per-iteration info."""
    rad.FldUnits('m')
    rad.UtiDelAll()

    print(f'\n{"="*70}')
    print(f'{solver_name} Solver - N={N_DIV} ({N_DIV**3} elements, {N_DIV**3 * 6} DOF)')
    print(f'{"="*70}')

    # Create mesh
    t0 = time.time()
    elements = create_hex_mesh(N_DIV)
    container = rad.ObjCnt(elements)
    t_mesh = time.time() - t0
    print(f'Mesh creation: {t_mesh:.4f} s')

    # Apply nonlinear material
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(container, mat)

    # External field
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([container, ext])

    # Configure HACApK if Method 2
    if solver_method == 2:
        rad.SetHACApKParams(HACAPK_EPS, 10, 2.0)
        print(f'HACApK: eps={HACAPK_EPS:.0e}, leaf_size=10, eta=2.0')

    # Run iterations one at a time
    print(f'\nRunning {n_iterations} iterations (1 at a time):')
    print(f'{"Iter":>4} {"Time (s)":>10} {"Residual":>12} {"M_avg_z":>12} {"dM_z":>12}')
    print('-' * 60)

    total_time = 0.0
    prev_M_z = 0.0
    results = []

    for i in range(n_iterations):
        t_start = time.time()
        res = rad.Solve(grp, 0.001, 1, solver_method)  # 1 iteration
        t_iter = time.time() - t_start
        total_time += t_iter

        M_avg_z = get_M_avg_z(container)
        dM_z = M_avg_z - prev_M_z
        residual = res[0] if res[0] else 0.0

        print(f'{i+1:4d} {t_iter:10.4f} {residual:12.4e} {M_avg_z:12.0f} {dM_z:12.0f}')

        results.append({
            'iter': i + 1,
            't_iter': t_iter,
            'residual': residual,
            'M_avg_z': M_avg_z,
            'dM_z': dM_z
        })

        prev_M_z = M_avg_z

        # Check convergence
        if abs(dM_z) < 100 and i > 0:  # Less than 100 A/m change
            print(f'Converged at iteration {i+1}')
            break

    print(f'\nTotal time: {total_time:.3f} s')
    print(f'Final M_avg_z: {M_avg_z:.0f} A/m')

    # Get H-matrix stats if available
    if solver_method == 2:
        try:
            stats = rad.GetHACApKStats()
            if stats:
                print(f'\nH-matrix stats:')
                print(f'  Leaves: {stats["n_leaves"]} (low-rank: {stats["n_lowrank"]}, dense: {stats["n_dense"]})')
                print(f'  Max rank: {stats["max_rank"]}')
                print(f'  Compression: {stats["compression"]:.4f}')
                print(f'  Build time: {stats["build_time"]:.4f} s')
        except:
            pass

    return results


def run_solver_batch(solver_method, solver_name, max_iter=1000):
    """Run solver in batch mode (all iterations at once)."""
    rad.FldUnits('m')
    rad.UtiDelAll()

    print(f'\n{"="*70}')
    print(f'{solver_name} Solver (BATCH) - N={N_DIV} ({N_DIV**3} elements, {N_DIV**3 * 6} DOF)')
    print(f'{"="*70}')

    # Create mesh
    t0 = time.time()
    elements = create_hex_mesh(N_DIV)
    container = rad.ObjCnt(elements)
    t_mesh = time.time() - t0
    print(f'Mesh creation: {t_mesh:.4f} s')

    # Apply nonlinear material
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(container, mat)

    # External field
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([container, ext])

    # Configure HACApK if Method 2
    if solver_method == 2:
        rad.SetHACApKParams(HACAPK_EPS, 10, 2.0)
        print(f'HACApK: eps={HACAPK_EPS:.0e}, leaf_size=10, eta=2.0')

    # Run solver in batch
    t_start = time.time()
    res = rad.Solve(grp, 0.001, max_iter, solver_method)
    t_solve = time.time() - t_start

    M_avg_z = get_M_avg_z(container)
    residual = res[0] if res[0] else 0.0
    n_iter = int(res[3]) if res[3] else 0

    print(f'Solve time: {t_solve:.3f} s')
    print(f'Iterations: {n_iter}')
    print(f'Residual: {residual:.4e}')
    print(f'M_avg_z: {M_avg_z:.0f} A/m')

    # Get H-matrix stats if available
    if solver_method == 2:
        try:
            stats = rad.GetHACApKStats()
            if stats:
                print(f'\nH-matrix stats:')
                print(f'  Leaves: {stats["n_leaves"]} (low-rank: {stats["n_lowrank"]}, dense: {stats["n_dense"]})')
                print(f'  Max rank: {stats["max_rank"]}')
                print(f'  Compression: {stats["compression"]:.4f}')
                print(f'  Build time: {stats["build_time"]:.4f} s')
        except:
            pass

    return {
        't_solve': t_solve,
        'n_iter': n_iter,
        'residual': residual,
        'M_avg_z': M_avg_z
    }


def main():
    print('=' * 70)
    print('HEX N=10 DEBUG COMPARISON')
    print('=' * 70)
    print(f'H_ext = {H_EXT} A/m')
    print(f'N = {N_DIV} ({N_DIV**3} elements, {N_DIV**3 * 6} DOF)')
    print(f'HACApK eps = {HACAPK_EPS:.0e}')

    # Run in BATCH mode (all iterations at once)
    result_lu = run_solver_batch(0, 'LU')
    result_bicgstab = run_solver_batch(1, 'BiCGSTAB')
    result_hacapk = run_solver_batch(2, 'HACApK')

    # Summary
    print('\n' + '=' * 70)
    print('SUMMARY (BATCH MODE)')
    print('=' * 70)

    print(f'{"Solver":<12} {"Time (s)":>10} {"Iter":>6} {"M_avg_z":>12} {"Residual":>12}')
    print('-' * 60)
    print(f'{"LU":<12} {result_lu["t_solve"]:>10.3f} {result_lu["n_iter"]:>6} {result_lu["M_avg_z"]:>12.0f} {result_lu["residual"]:>12.4e}')
    print(f'{"BiCGSTAB":<12} {result_bicgstab["t_solve"]:>10.3f} {result_bicgstab["n_iter"]:>6} {result_bicgstab["M_avg_z"]:>12.0f} {result_bicgstab["residual"]:>12.4e}')
    print(f'{"HACApK":<12} {result_hacapk["t_solve"]:>10.3f} {result_hacapk["n_iter"]:>6} {result_hacapk["M_avg_z"]:>12.0f} {result_hacapk["residual"]:>12.4e}')

    # Differences
    print('\nDifferences from LU:')
    diff_bicgstab = abs(result_bicgstab['M_avg_z'] - result_lu['M_avg_z'])
    diff_hacapk = abs(result_hacapk['M_avg_z'] - result_lu['M_avg_z'])
    print(f'  BiCGSTAB: {diff_bicgstab:.0f} A/m ({diff_bicgstab/result_lu["M_avg_z"]*100:.3f}%)')
    print(f'  HACApK:   {diff_hacapk:.0f} A/m ({diff_hacapk/result_lu["M_avg_z"]*100:.3f}%)')


if __name__ == '__main__':
    main()
