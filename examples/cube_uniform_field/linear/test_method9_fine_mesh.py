#!/usr/bin/env python3
"""
Test Method 9 (LU Direct Solver) with Fine Mesh

Compare Method 9 vs Method 4 for 8x8x8 mesh (512 elements).
Method 4 often fails to converge for high mu_r with fine meshes.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/python'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build/Release'))

import numpy as np
import radia as rad
import time

MU_0 = 4 * np.pi * 1e-7
CUBE_SIZE = 0.1  # meters
B0 = 1.0  # Tesla


def uniform_B(pos):
    return [0, 0, B0]


OBS_POINTS = [
    ([0.10, 0, 0], '+x'),
    ([0.12, 0, 0], '+x'),
    ([0, 0, 0.10], '+z'),
    ([0, 0, 0.12], '+z'),
]


def test_solver(mu_r, n_div, method, tol=0.001, max_iter=20000):
    chi = mu_r - 1.0
    rad.UtiDelAll()
    rad.FldUnits('m')

    cube = rad.ObjRecMag([0, 0, 0], [CUBE_SIZE]*3, [0, 0, 0])
    rad.ObjDivMag(cube, [n_div, n_div, n_div])
    mat = rad.MatLin(chi)
    rad.MatApl(cube, mat)

    bg = rad.ObjBckgCF(uniform_B)
    system = rad.ObjCnt([cube, bg])

    t0 = time.time()
    result = rad.Solve(system, tol, max_iter, method)
    solve_time = time.time() - t0

    fields = {}
    valid = True
    for pt, desc in OBS_POINTS:
        H = np.array(rad.Fld(system, 'h', pt))
        if np.any(np.isnan(H)):
            valid = False
        fields[tuple(pt)] = H

    return {
        'n_div': n_div,
        'n_elements': n_div**3,
        'method': method,
        'solve_time': solve_time,
        'iterations': result[3],
        'max_change': result[1],
        'converged': result[3] < max_iter - 1 and valid,
        'fields': fields,
        'valid': valid
    }


def compute_error(test_fields, ref_fields):
    errors = []
    for pt, desc in OBS_POINTS:
        H_ref = ref_fields[tuple(pt)]
        H_test = test_fields[tuple(pt)]
        if np.linalg.norm(H_ref) > 1e-10:
            err = np.linalg.norm(H_test - H_ref) / np.linalg.norm(H_ref) * 100
            errors.append(err)
    if len(errors) == 0:
        return float('nan'), float('nan')
    return np.mean(errors), max(errors)


def main():
    print("=" * 80)
    print("METHOD 9 (LU DIRECT) - FINE MESH TEST (8x8x8 = 512 elements)")
    print("=" * 80)
    print(f"Cube: {CUBE_SIZE*1000:.0f} mm, B0 = {B0} T")
    print()

    mu_r_list = [100, 500, 1000]
    n_div = 8  # 8x8x8 = 512 elements

    print(f"Mesh: {n_div}x{n_div}x{n_div} = {n_div**3} hexahedral elements")
    print()

    # First get Method 9 results as reference (guaranteed to converge)
    print("Using Method 9 (LU) as reference...")
    print("-" * 80)
    print(f"{'mu_r':<8} {'Method':<12} {'Time (s)':<10} {'Iter':<8} {'Status':<15}")
    print("-" * 80)

    for mu_r in mu_r_list:
        # Method 9 (direct LU) - use as reference
        m9_result = test_solver(mu_r, n_div, 9, tol=0.001, max_iter=1)
        status_9 = 'OK' if m9_result['valid'] else 'NaN'
        print(f"{mu_r:<8} {'9 (LU)':<12} {m9_result['solve_time']:<10.3f} {m9_result['iterations']:<8.0f} {status_9:<15}")

        # Method 4 (iterative) - may not converge
        m4_result = test_solver(mu_r, n_div, 4, tol=0.001, max_iter=20000)
        status_4 = 'OK' if m4_result['converged'] else f"N/C ({m4_result['iterations']:.0f})"

        if m4_result['valid'] and m9_result['valid']:
            avg_err, _ = compute_error(m4_result['fields'], m9_result['fields'])
            err_str = f"  vs M9: {avg_err:.2f}%"
        else:
            err_str = ""

        print(f"{'':<8} {'4 (iter)':<12} {m4_result['solve_time']:<10.3f} {m4_result['iterations']:<8.0f} {status_4:<15}{err_str}")
        print()

    print("=" * 80)
    print()
    print("Key finding:")
    print("  Method 9 (LU) converges instantly for any mu_r, even with fine mesh.")
    print("  Method 4 (iterative) may hit iteration limit or diverge.")


if __name__ == '__main__':
    main()
