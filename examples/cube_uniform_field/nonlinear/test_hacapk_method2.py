"""
Test HACApK (Method 2) with hexahedron MSC benchmark.

This script tests the H-matrix accelerated BiCGSTAB solver (Method 2)
and compares results with LU solver (Method 0).
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

import radia as rad
import time
import json

# Test parameters
N = 5  # 5x5x5 = 125 elements
H_EXT = 200000.0  # A/m

# BH curve (soft iron)
BH_CURVE = [
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

def create_hexahedron_mesh(n_div, size=1.0):
    """Create hexahedral mesh using ObjPolyhdr with 6DOF MSC."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    half = size / 2.0
    dx = size / n_div

    elements = []
    for iz in range(n_div):
        for iy in range(n_div):
            for ix in range(n_div):
                x0 = -half + ix * dx
                y0 = -half + iy * dx
                z0 = -half + iz * dx
                x1 = x0 + dx
                y1 = y0 + dx
                z1 = z0 + dx

                # Hexahedron vertices (Radia ordering)
                vertices = [
                    [x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                    [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]
                ]

                # 6 quadrilateral faces (1-indexed for Radia)
                faces = [
                    [1, 4, 3, 2],  # -z face
                    [5, 6, 7, 8],  # +z face
                    [1, 2, 6, 5],  # -y face
                    [3, 4, 8, 7],  # +y face
                    [1, 5, 8, 4],  # -x face
                    [2, 3, 7, 6],  # +x face
                ]

                elem = rad.ObjPolyhdr(vertices, faces, [0.0, 0.0, 0.0])
                elements.append(elem)

    grp = rad.ObjCnt(elements)
    return grp, len(elements)


def run_benchmark(solver_method, solver_name):
    """Run benchmark with specified solver method."""
    print(f"\n{'='*60}")
    print(f"Solver: {solver_name} (Method {solver_method})")
    print(f"{'='*60}")

    # Create mesh
    t0 = time.time()
    container, n_elem = create_hexahedron_mesh(N)
    t_mesh = time.time() - t0

    # Apply nonlinear material
    mat = rad.MatSatIsoTab(BH_CURVE)
    rad.MatApl(container, mat)

    # Apply external field using ObjBckg (like benchmark_hexahedron_msc.py)
    MU_0 = 4 * 3.14159265359 * 1e-7
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([container, ext])

    # Solve
    print(f"Solving with {n_elem} elements, {n_elem * 6} DOF...")
    t0 = time.time()
    try:
        res = rad.Solve(grp, 0.0001, 100, solver_method)
        t_solve = time.time() - t0

        # Get results
        converged = res[0] < 0.0001
        residual = res[0]
        iterations = int(res[3])

        # Get average magnetization from container (not grp which includes background)
        M_all = rad.ObjM(container)
        M_z_sum = sum(m[1][2] for m in M_all)
        M_avg_z = M_z_sum / len(M_all)

        print(f"  Time: {t_solve:.3f} s")
        print(f"  Iterations: {iterations}")
        print(f"  Residual: {residual:.6e}")
        print(f"  Converged: {converged}")
        print(f"  M_avg_z: {M_avg_z:.1f} A/m")

        return {
            'solver_method': solver_method,
            'solver_name': solver_name,
            'n_elements': n_elem,
            'ndof': n_elem * 6,
            'H_ext': H_EXT,
            't_mesh': t_mesh,
            't_solve': t_solve,
            'iterations': iterations,
            'residual': residual,
            'converged': converged,
            'M_avg_z': M_avg_z,
        }
    except Exception as e:
        print(f"  Error: {e}")
        return {
            'solver_method': solver_method,
            'solver_name': solver_name,
            'error': str(e),
        }


def main():
    print("HACApK (Method 2) Benchmark Test")
    print(f"N = {N}, Elements = {N**3}, DOF = {N**3 * 6}")
    print(f"H_ext = {H_EXT} A/m")

    results = {}

    # Test LU solver (Method 0) as reference
    results['lu'] = run_benchmark(0, 'LU')

    # Test BiCGSTAB solver (Method 1)
    results['bicgstab'] = run_benchmark(1, 'BiCGSTAB')

    # Test HACApK solver (Method 2)
    results['hacapk'] = run_benchmark(2, 'HACApK')

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)

    for name, r in results.items():
        if 'error' in r:
            print(f"{name}: ERROR - {r['error']}")
        else:
            print(f"{name}: {r['iterations']} iter, {r['t_solve']:.3f}s, M_avg_z = {r['M_avg_z']:.1f} A/m")

    # Compare with reference
    if 'error' not in results['lu'] and 'error' not in results.get('hacapk', {'error': 'not run'}):
        lu_M = results['lu']['M_avg_z']
        hacapk_M = results['hacapk']['M_avg_z']
        diff_pct = abs(hacapk_M - lu_M) / lu_M * 100
        print(f"\nHACApK vs LU difference: {diff_pct:.4f}%")

    # Save results
    output_file = os.path.join(os.path.dirname(__file__), 'hacapk_test_results.json')
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {output_file}")


if __name__ == '__main__':
    main()
