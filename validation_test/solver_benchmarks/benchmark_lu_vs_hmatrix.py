"""
Benchmark: LU Decomposition vs BiCGSTAB vs HACApK (H-matrix)

Compares three solver methods:
- Method 0 (LU): Direct solver, O(N^3) complexity
- Method 1 (BiCGSTAB): Iterative solver, O(N^2) per iteration
- Method 2 (HACApK): H-matrix accelerated BiCGSTAB, O(N^2 log N) per iteration

Evaluates:
- Computational time scaling with problem size
- Solution accuracy consistency
- Memory usage (estimated)
"""

import sys
import os
import time
import gc
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
import radia as rad

# Set unit system to meters
mm = 1e-3  # 1 mm in meters


def hex_vertices(cx, cy, cz, dx, dy, dz):
    """Generate hexahedron vertices from center and dimensions."""
    hx, hy, hz = dx/2, dy/2, dz/2
    return [
        [cx-hx, cy-hy, cz-hz], [cx+hx, cy-hy, cz-hz],
        [cx+hx, cy+hy, cz-hz], [cx-hx, cy+hy, cz-hz],
        [cx-hx, cy-hy, cz+hz], [cx+hx, cy-hy, cz+hz],
        [cx+hx, cy+hy, cz+hz], [cx-hx, cy+hy, cz+hz]
    ]


def create_geometry(n, size, B_bg):
    """Create hexahedral mesh with material and background field.

    Args:
        B_bg: Background field in Tesla (ObjBckg callback returns B, not H).
    """
    elem_size = size / n
    mat = rad.MatLin(1000.0)  # mu_r = 1000

    elements = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                x = -size/2 + (i + 0.5) * elem_size
                y = -size/2 + (j + 0.5) * elem_size
                z = -size/2 + (k + 0.5) * elem_size
                vertices = hex_vertices(x, y, z, elem_size, elem_size, elem_size)
                elem = rad.ObjHexahedron(vertices, [0, 0, 0])
                rad.MatApl(elem, mat)
                elements.append(elem)

    bg_field = rad.ObjBckg(lambda p: B_bg)
    return rad.ObjCnt(elements + [bg_field])


def benchmark_solver(n, method, tol=1e-4, max_iter=1000):
    """
    Benchmark solver with specified method.
    Returns timing and solution.
    """
    rad.UtiDelAll()
    gc.collect()

    size = 40.0 * mm  # 40mm cube
    B_bg = [0, 0, 1.0]  # Background field in Tesla

    container = create_geometry(n, size, B_bg)

    # Configure HACApK if method 2
    if method == 2:
        rad.SolverConfig(hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)

    # Solve
    start_total = time.perf_counter()
    res = rad.Solve(container, tol, max_iter, method)
    total_time = time.perf_counter() - start_total

    # Get field at evaluation points
    B_center = rad.Fld(container, 'b', [0, 0, 0])
    B_outside = rad.Fld(container, 'b', [0, 0, 60 * mm])

    iterations = int(res[1]) if isinstance(res, (list, tuple)) and len(res) > 1 else 1

    return {
        'total_time': total_time,
        'iterations': iterations,
        'Bz_center': B_center[2],
        'Bz_outside': B_outside[2],
    }


def estimate_memory(n_elem, method):
    """Estimate memory usage in MB."""
    n_dof = n_elem * 6  # 6 face-charge DOF per hexahedral element
    if method in [0, 1]:  # LU or BiCGSTAB (dense matrix)
        return n_dof * n_dof * 8 / (1024 * 1024)
    elif method == 2:  # HACApK (H-matrix compression)
        # H-matrix: approximately 30-50% compression
        return n_dof * n_dof * 8 * 0.5 / (1024 * 1024)
    return 0


def run_benchmark():
    print('=' * 80)
    print('Benchmark: LU vs BiCGSTAB vs HACApK')
    print('=' * 80)
    print()
    print('Problem: Soft iron cube (mu_r=1000) in uniform B0=1T field')
    print('Evaluation: Bz at center and at z=60mm (outside)')
    print()

    # Test configurations (increasing problem sizes)
    test_sizes = [3, 4, 5, 6, 7, 8]  # n x n x n elements

    results = []

    for n in test_sizes:
        n_elem = n ** 3
        n_dof = n_elem * 6

        print(f'\n{"="*80}')
        print(f'Size: {n}x{n}x{n} ({n_elem} elements, {n_dof} DOF)')
        print(f'{"="*80}')

        row = {'n': n, 'n_elem': n_elem, 'n_dof': n_dof}

        # Method 0: LU
        print('  Running Method 0 (LU)...', end=' ', flush=True)
        try:
            lu_result = benchmark_solver(n, method=0)
            row['lu_time'] = lu_result['total_time']
            row['lu_Bz_outside'] = lu_result['Bz_outside']
            row['lu_memory'] = estimate_memory(n_elem, 0)
            print(f'Done ({lu_result["total_time"]:.3f}s)')
        except Exception as e:
            print(f'Failed: {e}')
            row['lu_time'] = float('nan')
            row['lu_Bz_outside'] = float('nan')
            row['lu_memory'] = 0

        # Method 1: BiCGSTAB
        print('  Running Method 1 (BiCGSTAB)...', end=' ', flush=True)
        try:
            bicg_result = benchmark_solver(n, method=1)
            row['bicg_time'] = bicg_result['total_time']
            row['bicg_iters'] = bicg_result['iterations']
            row['bicg_Bz_outside'] = bicg_result['Bz_outside']
            row['bicg_memory'] = estimate_memory(n_elem, 1)
            print(f'Done ({bicg_result["total_time"]:.3f}s, {bicg_result["iterations"]} iters)')
        except Exception as e:
            print(f'Failed: {e}')
            row['bicg_time'] = float('nan')
            row['bicg_iters'] = 0
            row['bicg_Bz_outside'] = float('nan')
            row['bicg_memory'] = 0

        # Method 2: HACApK
        print('  Running Method 2 (HACApK)...', end=' ', flush=True)
        try:
            hmat_result = benchmark_solver(n, method=2)
            row['hmat_time'] = hmat_result['total_time']
            row['hmat_iters'] = hmat_result['iterations']
            row['hmat_Bz_outside'] = hmat_result['Bz_outside']
            row['hmat_memory'] = estimate_memory(n_elem, 2)
            print(f'Done ({hmat_result["total_time"]:.3f}s, {hmat_result["iterations"]} iters)')
        except Exception as e:
            print(f'Failed: {e}')
            row['hmat_time'] = float('nan')
            row['hmat_iters'] = 0
            row['hmat_Bz_outside'] = float('nan')
            row['hmat_memory'] = 0

        results.append(row)

        # Print immediate comparison
        if not any(np.isnan([row.get('lu_time', np.nan),
                             row.get('bicg_time', np.nan),
                             row.get('hmat_time', np.nan)])):
            print(f'\n  Timing comparison:')
            print(f'    LU:      {row["lu_time"]:8.4f}s')
            print(f'    BiCGSTAB: {row["bicg_time"]:8.4f}s ({row["bicg_iters"]} iters)')
            print(f'    HACApK:  {row["hmat_time"]:8.4f}s ({row["hmat_iters"]} iters)')

            if row['lu_time'] > 0:
                print(f'\n  Speedup vs LU:')
                print(f'    BiCGSTAB: {row["lu_time"]/row["bicg_time"]:.2f}x')
                print(f'    HACApK:  {row["lu_time"]/row["hmat_time"]:.2f}x')

    # Summary table
    print('\n')
    print('=' * 80)
    print('SUMMARY TABLE')
    print('=' * 80)
    print()

    # Timing table
    print('Computation Time (seconds):')
    print('-' * 80)
    print(f'{"N_elem":<8} {"N_DOF":<8} {"LU":<12} {"BiCGSTAB":<12} {"HACApK":<12} {"LU/BiCG":<10} {"LU/HACApK":<10}')
    print('-' * 80)

    for row in results:
        n_elem = row['n_elem']
        n_dof = row['n_dof']
        lu_t = row.get('lu_time', float('nan'))
        bicg_t = row.get('bicg_time', float('nan'))
        hmat_t = row.get('hmat_time', float('nan'))

        lu_str = f'{lu_t:.4f}' if not np.isnan(lu_t) else 'N/A'
        bicg_str = f'{bicg_t:.4f}' if not np.isnan(bicg_t) else 'N/A'
        hmat_str = f'{hmat_t:.4f}' if not np.isnan(hmat_t) else 'N/A'

        ratio1 = f'{lu_t/bicg_t:.2f}x' if (not np.isnan(lu_t) and not np.isnan(bicg_t) and bicg_t > 0) else 'N/A'
        ratio2 = f'{lu_t/hmat_t:.2f}x' if (not np.isnan(lu_t) and not np.isnan(hmat_t) and hmat_t > 0) else 'N/A'

        print(f'{n_elem:<8} {n_dof:<8} {lu_str:<12} {bicg_str:<12} {hmat_str:<12} {ratio1:<10} {ratio2:<10}')

    # Accuracy table
    print()
    print('Solution Accuracy (Bz at z=60mm):')
    print('-' * 80)
    print(f'{"N_elem":<8} {"LU (T)":<16} {"BiCGSTAB (T)":<16} {"HACApK (T)":<16} {"BiCG err%":<10} {"HACApK err%":<10}')
    print('-' * 80)

    for row in results:
        n_elem = row['n_elem']
        lu_Bz = row.get('lu_Bz_outside', float('nan'))
        bicg_Bz = row.get('bicg_Bz_outside', float('nan'))
        hmat_Bz = row.get('hmat_Bz_outside', float('nan'))

        bicg_err = abs(bicg_Bz - lu_Bz) / abs(lu_Bz) * 100 if (not np.isnan(lu_Bz) and not np.isnan(bicg_Bz) and lu_Bz != 0) else float('nan')
        hmat_err = abs(hmat_Bz - lu_Bz) / abs(lu_Bz) * 100 if (not np.isnan(lu_Bz) and not np.isnan(hmat_Bz) and lu_Bz != 0) else float('nan')

        lu_str = f'{lu_Bz:.10f}' if not np.isnan(lu_Bz) else 'N/A'
        bicg_str = f'{bicg_Bz:.10f}' if not np.isnan(bicg_Bz) else 'N/A'
        hmat_str = f'{hmat_Bz:.10f}' if not np.isnan(hmat_Bz) else 'N/A'
        bicg_err_str = f'{bicg_err:.6f}' if not np.isnan(bicg_err) else 'N/A'
        hmat_err_str = f'{hmat_err:.6f}' if not np.isnan(hmat_err) else 'N/A'

        print(f'{n_elem:<8} {lu_str:<16} {bicg_str:<16} {hmat_str:<16} {bicg_err_str:<10} {hmat_err_str:<10}')

    print()
    print('=' * 80)
    print('Benchmark complete')
    print('=' * 80)


if __name__ == '__main__':
    run_benchmark()
