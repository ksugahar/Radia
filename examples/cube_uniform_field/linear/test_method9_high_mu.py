#!/usr/bin/env python3
"""
Test Method 9 (LU Direct Solver) for High Permeability Materials

Compare Method 9 (direct) vs Method 4 (iterative) for high mu_r values.
Method 9 should be stable for any permeability, similar to ELF/MAGIC.
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
    """Background field callback."""
    return [0, 0, B0]


OBS_POINTS = [
    ([0.10, 0, 0], '+x, r=0.10m'),
    ([0.12, 0, 0], '+x, r=0.12m'),
    ([0.15, 0, 0], '+x, r=0.15m'),
    ([0, 0, 0.10], '+z, r=0.10m'),
    ([0, 0, 0.12], '+z, r=0.12m'),
    ([0, 0, 0.15], '+z, r=0.15m'),
]


def test_solver(mu_r, n_div, method, tol=0.001, max_iter=20000):
    """Test solver with given method."""
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
    # rad.Solve takes: object, precision, max_iter, method
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
    """Compute average and max error vs reference."""
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
    print("METHOD 9 (LU DIRECT) vs METHOD 4 (ITERATIVE) COMPARISON")
    print("=" * 80)
    print(f"Cube: {CUBE_SIZE*1000:.0f} mm, B0 = {B0} T")
    print()

    mu_r_list = [100, 500, 1000, 2000, 4000]
    n_div = 4  # Use 4x4x4 = 64 elements

    print(f"Mesh: {n_div}x{n_div}x{n_div} = {n_div**3} hexahedral elements")
    print()
    print("-" * 80)
    print(f"{'mu_r':<8} {'Method':<10} {'Time (s)':<10} {'Iter':<8} {'Status':<15} {'Avg Err (%)':<12}")
    print("-" * 80)

    for mu_r in mu_r_list:
        # Method 4 (iterative) as reference
        ref_result = test_solver(mu_r, n_div, 4, tol=0.001, max_iter=50000)

        status_4 = 'OK' if ref_result['converged'] else f"NOT CONV ({ref_result['max_change']:.1e})"
        print(f"{mu_r:<8} {'4 (iter)':<10} {ref_result['solve_time']:<10.3f} {ref_result['iterations']:<8} {status_4:<15} {'(ref)':<12}")

        # Method 9 (direct LU)
        m9_result = test_solver(mu_r, n_div, 9, tol=0.001, max_iter=1)

        if ref_result['valid'] and m9_result['valid']:
            avg_err, max_err = compute_error(m9_result['fields'], ref_result['fields'])
            err_str = f"{avg_err:.2f}%"
        else:
            err_str = 'N/A' if not m9_result['valid'] else '(ref N/A)'

        status_9 = 'OK' if m9_result['valid'] else 'NaN'
        print(f"{'':<8} {'9 (LU)':<10} {m9_result['solve_time']:<10.3f} {m9_result['iterations']:<8} {status_9:<15} {err_str:<12}")
        print()

    print("=" * 80)
    print()
    print("Notes:")
    print("  - Method 4: Iterative relaxation (default Radia method)")
    print("  - Method 9: Direct LU solver (ELF/MAGIC-like)")
    print("  - Method 9 should converge in 1 'iteration' for linear materials")
    print("  - Method 9 should be stable for any permeability")


if __name__ == '__main__':
    main()
