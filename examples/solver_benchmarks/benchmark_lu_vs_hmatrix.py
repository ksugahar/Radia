"""
Benchmark: LU Decomposition vs H-matrix BiCGSTAB

Compares:
- Method 9 (LU): Matrix assembly time + LU solve time
- Method 10 (BiCGSTAB) + H-matrix: H-matrix build time + iterative solve time

Evaluates:
- Computational time scaling with problem size
- Solution accuracy consistency
- Memory usage (estimated)

Author: Claude Code
Date: 2025-12-04
"""

import sys
import os
import time
import gc

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
import radia as rad


def hex_vertices(cx, cy, cz, dx, dy, dz):
    """Generate hexahedron vertices from center and dimensions."""
    hx, hy, hz = dx/2, dy/2, dz/2
    return [
        [cx-hx, cy-hy, cz-hz], [cx+hx, cy-hy, cz-hz],
        [cx+hx, cy+hy, cz-hz], [cx-hx, cy+hy, cz-hz],
        [cx-hx, cy-hy, cz+hz], [cx+hx, cy-hy, cz+hz],
        [cx+hx, cy+hy, cz+hz], [cx-hx, cy+hy, cz+hz]
    ]


def benchmark_lu(subdiv, tol=1e-4, max_iter=1):
    """
    Benchmark LU decomposition (Method 9).
    Returns timing breakdown and solution.
    """
    rad.UtiDelAll()
    rad.SolverHMatrixDisable()
    gc.collect()

    # Create geometry: 40x40x40 cube centered at origin
    vertices = hex_vertices(0, 0, 0, 40, 40, 40)
    cube = rad.ObjHexahedron(vertices, [0, 0, 0])
    rad.ObjDivMag(cube, subdiv)

    # Apply linear material (mu_r = 1000)
    mat = rad.MatLin(1000.0)  # relative permeability
    rad.MatApl(cube, mat)

    # Add external field
    ext_field = rad.ObjBckg([0, 0, 1.0])
    grp = rad.ObjCnt([cube, ext_field])

    # Solve with Method 9 (LU)
    # Note: Radia combines matrix assembly and solve internally
    start_total = time.perf_counter()
    res = rad.Solve(grp, tol, max_iter, 9)
    total_time = time.perf_counter() - start_total

    # Get field at evaluation points
    B_center = rad.Fld(grp, 'b', [0, 0, 0])
    B_outside = rad.Fld(grp, 'b', [0, 0, 60])

    return {
        'total_time': total_time,
        'iterations': 1,  # LU is direct solver
        'Bz_center': B_center[2],
        'Bz_outside': B_outside[2],
        'max_abs_M': res[0],
        'max_abs_H': res[1],
    }


def benchmark_bicgstab(subdiv, use_hmatrix=False, tol=1e-4, max_iter=1000):
    """
    Benchmark BiCGSTAB (Method 10) with or without H-matrix.
    Returns timing breakdown and solution.
    """
    rad.UtiDelAll()
    gc.collect()

    # Configure H-matrix
    if use_hmatrix:
        rad.SolverHMatrixEnable()
    else:
        rad.SolverHMatrixDisable()

    # Create geometry: 40x40x40 cube centered at origin
    vertices = hex_vertices(0, 0, 0, 40, 40, 40)
    cube = rad.ObjHexahedron(vertices, [0, 0, 0])
    rad.ObjDivMag(cube, subdiv)

    # Apply linear material (mu_r = 1000)
    mat = rad.MatLin(1000.0)  # relative permeability
    rad.MatApl(cube, mat)

    # Add external field
    ext_field = rad.ObjBckg([0, 0, 1.0])
    grp = rad.ObjCnt([cube, ext_field])

    # Solve with Method 10 (BiCGSTAB)
    start_total = time.perf_counter()
    res = rad.Solve(grp, tol, max_iter, 10)
    total_time = time.perf_counter() - start_total

    # Get field at evaluation points
    B_center = rad.Fld(grp, 'b', [0, 0, 0])
    B_outside = rad.Fld(grp, 'b', [0, 0, 60])

    return {
        'total_time': total_time,
        'iterations': int(res[3]) if len(res) > 3 else 0,
        'Bz_center': B_center[2],
        'Bz_outside': B_outside[2],
        'max_abs_M': res[0],
        'max_abs_H': res[1],
    }


def estimate_memory(n_elem, method):
    """Estimate memory usage in MB."""
    n_dof = n_elem * 3
    if method == 'lu':
        # Full matrix: n_dof x n_dof x 8 bytes (double)
        return n_dof * n_dof * 8 / (1024 * 1024)
    elif method == 'bicgstab':
        # Same as LU (dense matrix)
        return n_dof * n_dof * 8 / (1024 * 1024)
    elif method == 'hmatrix':
        # H-matrix: approximately 30-50% compression for MMM kernel
        # Conservative estimate: 60% of full matrix
        return n_dof * n_dof * 8 * 0.6 / (1024 * 1024)
    return 0


def run_benchmark():
    print('=' * 80)
    print('Benchmark: LU Decomposition vs H-matrix BiCGSTAB')
    print('=' * 80)
    print()
    print('Problem: Soft iron cube (mu_r=1000) in uniform B0=1T field')
    print('Evaluation: Bz at center and at z=60mm (outside)')
    print()

    # Test configurations (increasing problem sizes)
    subdivisions = [
        [3, 3, 3],    # 27 elements
        [4, 4, 4],    # 64 elements
        [5, 5, 5],    # 125 elements
        [6, 6, 6],    # 216 elements
        [7, 7, 7],    # 343 elements
        [8, 8, 8],    # 512 elements
        [9, 9, 9],    # 729 elements
        [10, 10, 10], # 1000 elements
        [12, 12, 12], # 1728 elements
        [14, 14, 14], # 2744 elements
    ]

    results = []

    for subdiv in subdivisions:
        n_elem = subdiv[0] * subdiv[1] * subdiv[2]
        n_dof = n_elem * 3

        print(f'\n{"="*80}')
        print(f'Subdivision: {subdiv[0]}x{subdiv[1]}x{subdiv[2]} ({n_elem} elements, {n_dof} DOF)')
        print(f'{"="*80}')

        row = {'n_elem': n_elem, 'n_dof': n_dof, 'subdiv': subdiv}

        # Method 9: LU decomposition
        print('  Running LU decomposition (Method 9)...', end=' ', flush=True)
        try:
            lu_result = benchmark_lu(subdiv)
            row['lu_time'] = lu_result['total_time']
            row['lu_Bz_center'] = lu_result['Bz_center']
            row['lu_Bz_outside'] = lu_result['Bz_outside']
            row['lu_memory'] = estimate_memory(n_elem, 'lu')
            print(f'Done ({lu_result["total_time"]:.3f}s)')
        except Exception as e:
            print(f'Failed: {e}')
            row['lu_time'] = float('nan')
            row['lu_Bz_center'] = float('nan')
            row['lu_Bz_outside'] = float('nan')
            row['lu_memory'] = 0

        # Method 10: BiCGSTAB (no H-matrix)
        print('  Running BiCGSTAB (Method 10)...', end=' ', flush=True)
        try:
            bicg_result = benchmark_bicgstab(subdiv, use_hmatrix=False)
            row['bicg_time'] = bicg_result['total_time']
            row['bicg_iters'] = bicg_result['iterations']
            row['bicg_Bz_center'] = bicg_result['Bz_center']
            row['bicg_Bz_outside'] = bicg_result['Bz_outside']
            row['bicg_memory'] = estimate_memory(n_elem, 'bicgstab')
            print(f'Done ({bicg_result["total_time"]:.3f}s, {bicg_result["iterations"]} iters)')
        except Exception as e:
            print(f'Failed: {e}')
            row['bicg_time'] = float('nan')
            row['bicg_iters'] = 0
            row['bicg_Bz_center'] = float('nan')
            row['bicg_Bz_outside'] = float('nan')
            row['bicg_memory'] = 0

        # Method 10 + H-matrix
        print('  Running BiCGSTAB + H-matrix (Method 10)...', end=' ', flush=True)
        try:
            hmat_result = benchmark_bicgstab(subdiv, use_hmatrix=True)
            row['hmat_time'] = hmat_result['total_time']
            row['hmat_iters'] = hmat_result['iterations']
            row['hmat_Bz_center'] = hmat_result['Bz_center']
            row['hmat_Bz_outside'] = hmat_result['Bz_outside']
            row['hmat_memory'] = estimate_memory(n_elem, 'hmatrix')
            print(f'Done ({hmat_result["total_time"]:.3f}s, {hmat_result["iterations"]} iters)')
        except Exception as e:
            print(f'Failed: {e}')
            row['hmat_time'] = float('nan')
            row['hmat_iters'] = 0
            row['hmat_Bz_center'] = float('nan')
            row['hmat_Bz_outside'] = float('nan')
            row['hmat_memory'] = 0

        results.append(row)

        # Print immediate comparison
        if not any(x != x for x in [row.get('lu_time', float('nan')),
                                     row.get('bicg_time', float('nan')),
                                     row.get('hmat_time', float('nan'))]):
            print(f'\n  Timing comparison:')
            print(f'    LU:             {row["lu_time"]:8.4f}s')
            print(f'    BiCGSTAB:       {row["bicg_time"]:8.4f}s ({row["bicg_iters"]} iters)')
            print(f'    BiCGSTAB+Hmat:  {row["hmat_time"]:8.4f}s ({row["hmat_iters"]} iters)')

            if row['lu_time'] > 0:
                print(f'\n  Speedup vs LU:')
                print(f'    BiCGSTAB:       {row["lu_time"]/row["bicg_time"]:.2f}x')
                print(f'    BiCGSTAB+Hmat:  {row["lu_time"]/row["hmat_time"]:.2f}x')

            # Accuracy check
            lu_Bz = row['lu_Bz_outside']
            bicg_diff = abs(row['bicg_Bz_outside'] - lu_Bz) / abs(lu_Bz) * 100 if lu_Bz != 0 else 0
            hmat_diff = abs(row['hmat_Bz_outside'] - lu_Bz) / abs(lu_Bz) * 100 if lu_Bz != 0 else 0
            print(f'\n  Accuracy (Bz at z=60mm, vs LU reference):')
            print(f'    LU:             {lu_Bz:.10f} T (reference)')
            print(f'    BiCGSTAB:       {row["bicg_Bz_outside"]:.10f} T ({bicg_diff:.4f}% diff)')
            print(f'    BiCGSTAB+Hmat:  {row["hmat_Bz_outside"]:.10f} T ({hmat_diff:.4f}% diff)')

    # Summary table
    print('\n')
    print('=' * 80)
    print('SUMMARY TABLE')
    print('=' * 80)
    print()

    # Timing table
    print('Computation Time (seconds):')
    print('-' * 80)
    print(f'{"N_elem":<8} {"N_DOF":<8} {"LU":<12} {"BiCGSTAB":<12} {"Hmat+BiCG":<12} {"LU/BiCG":<10} {"LU/Hmat":<10}')
    print('-' * 80)

    for row in results:
        n_elem = row['n_elem']
        n_dof = row['n_dof']
        lu_t = row.get('lu_time', float('nan'))
        bicg_t = row.get('bicg_time', float('nan'))
        hmat_t = row.get('hmat_time', float('nan'))

        lu_str = f'{lu_t:.4f}' if lu_t == lu_t else 'N/A'
        bicg_str = f'{bicg_t:.4f}' if bicg_t == bicg_t else 'N/A'
        hmat_str = f'{hmat_t:.4f}' if hmat_t == hmat_t else 'N/A'

        ratio1 = f'{lu_t/bicg_t:.2f}x' if (lu_t == lu_t and bicg_t == bicg_t and bicg_t > 0) else 'N/A'
        ratio2 = f'{lu_t/hmat_t:.2f}x' if (lu_t == lu_t and hmat_t == hmat_t and hmat_t > 0) else 'N/A'

        print(f'{n_elem:<8} {n_dof:<8} {lu_str:<12} {bicg_str:<12} {hmat_str:<12} {ratio1:<10} {ratio2:<10}')

    # Iteration count table
    print()
    print('BiCGSTAB Iteration Counts:')
    print('-' * 50)
    print(f'{"N_elem":<10} {"BiCGSTAB":<15} {"Hmat+BiCGSTAB":<15}')
    print('-' * 50)

    for row in results:
        n_elem = row['n_elem']
        bicg_iters = row.get('bicg_iters', 0)
        hmat_iters = row.get('hmat_iters', 0)
        print(f'{n_elem:<10} {bicg_iters:<15} {hmat_iters:<15}')

    # Accuracy table
    print()
    print('Solution Accuracy (Bz at z=60mm):')
    print('-' * 80)
    print(f'{"N_elem":<8} {"LU (T)":<16} {"BiCGSTAB (T)":<16} {"Hmat (T)":<16} {"BiCG err%":<10} {"Hmat err%":<10}')
    print('-' * 80)

    for row in results:
        n_elem = row['n_elem']
        lu_Bz = row.get('lu_Bz_outside', float('nan'))
        bicg_Bz = row.get('bicg_Bz_outside', float('nan'))
        hmat_Bz = row.get('hmat_Bz_outside', float('nan'))

        bicg_err = abs(bicg_Bz - lu_Bz) / abs(lu_Bz) * 100 if (lu_Bz == lu_Bz and bicg_Bz == bicg_Bz and lu_Bz != 0) else float('nan')
        hmat_err = abs(hmat_Bz - lu_Bz) / abs(lu_Bz) * 100 if (lu_Bz == lu_Bz and hmat_Bz == hmat_Bz and lu_Bz != 0) else float('nan')

        lu_str = f'{lu_Bz:.10f}' if lu_Bz == lu_Bz else 'N/A'
        bicg_str = f'{bicg_Bz:.10f}' if bicg_Bz == bicg_Bz else 'N/A'
        hmat_str = f'{hmat_Bz:.10f}' if hmat_Bz == hmat_Bz else 'N/A'
        bicg_err_str = f'{bicg_err:.6f}' if bicg_err == bicg_err else 'N/A'
        hmat_err_str = f'{hmat_err:.6f}' if hmat_err == hmat_err else 'N/A'

        print(f'{n_elem:<8} {lu_str:<16} {bicg_str:<16} {hmat_str:<16} {bicg_err_str:<10} {hmat_err_str:<10}')

    # Memory estimate table
    print()
    print('Estimated Memory Usage (MB):')
    print('-' * 50)
    print(f'{"N_elem":<10} {"LU/BiCGSTAB":<15} {"H-matrix":<15}')
    print('-' * 50)

    for row in results:
        n_elem = row['n_elem']
        dense_mem = row.get('lu_memory', 0)
        hmat_mem = row.get('hmat_memory', 0)
        print(f'{n_elem:<10} {dense_mem:<15.2f} {hmat_mem:<15.2f}')

    print()
    print('=' * 80)
    print('Benchmark complete')
    print('=' * 80)

    # Save results to CSV
    csv_file = os.path.join(os.path.dirname(__file__), 'benchmark_lu_vs_hmatrix_results.csv')
    with open(csv_file, 'w') as f:
        f.write('n_elem,n_dof,lu_time,bicg_time,hmat_time,bicg_iters,hmat_iters,')
        f.write('lu_Bz_outside,bicg_Bz_outside,hmat_Bz_outside,bicg_err_pct,hmat_err_pct\n')
        for row in results:
            lu_Bz = row.get('lu_Bz_outside', float('nan'))
            bicg_Bz = row.get('bicg_Bz_outside', float('nan'))
            hmat_Bz = row.get('hmat_Bz_outside', float('nan'))
            bicg_err = abs(bicg_Bz - lu_Bz) / abs(lu_Bz) * 100 if (lu_Bz == lu_Bz and bicg_Bz == bicg_Bz and lu_Bz != 0) else float('nan')
            hmat_err = abs(hmat_Bz - lu_Bz) / abs(lu_Bz) * 100 if (lu_Bz == lu_Bz and hmat_Bz == hmat_Bz and lu_Bz != 0) else float('nan')

            f.write(f'{row["n_elem"]},{row["n_dof"]},')
            f.write(f'{row.get("lu_time", "")},{row.get("bicg_time", "")},{row.get("hmat_time", "")},')
            f.write(f'{row.get("bicg_iters", "")},{row.get("hmat_iters", "")},')
            f.write(f'{lu_Bz},{bicg_Bz},{hmat_Bz},{bicg_err},{hmat_err}\n')

    print(f'\nResults saved to: {csv_file}')


if __name__ == '__main__':
    run_benchmark()
