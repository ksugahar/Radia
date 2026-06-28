"""
H-Matrix Solver Performance Benchmark

Compares three solver methods:
1. Method 0 (LU): Direct solver
2. Method 1 (BiCGSTAB): Iterative solver
3. Method 2 (HACApK): H-matrix accelerated solver

Demonstrates the speedup and memory reduction from H-matrix acceleration.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build/Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import radia as rad
import time
import numpy as np

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


def create_magnet(n_per_side):
    """
    Create a cubic magnet subdivided into n x n x n elements.

    Args:
        n_per_side: Number of subdivisions per side

    Returns:
        Radia object ID
    """
    size = 20.0 * mm  # 20mm cube
    elem_size = size / n_per_side

    print(f"Creating {n_per_side}x{n_per_side}x{n_per_side} = {n_per_side**3} elements...")

    mat = rad.MatSatIsoFrm([2000, 2], [0.1, 2], [0.1, 2])

    elements = []
    start_time = time.time()

    for i in range(n_per_side):
        for j in range(n_per_side):
            for k in range(n_per_side):
                # Element center
                x = (i - n_per_side/2 + 0.5) * elem_size
                y = (j - n_per_side/2 + 0.5) * elem_size
                z = (k - n_per_side/2 + 0.5) * elem_size

                # Create element
                vertices = hex_vertices(x, y, z, elem_size, elem_size, elem_size)
                block = rad.ObjHexahedron(vertices, [0, 0, 1])
                rad.MatApl(block, mat)
                elements.append(block)

    container = rad.ObjCnt(elements)

    creation_time = time.time() - start_time
    print(f"  Magnet created in {creation_time:.3f} s")

    return container


def benchmark_solver(magnet, method, precision=0.0001, max_iter=1000):
    """
    Benchmark solver with specified method.
    """
    method_names = {0: "LU", 1: "BiCGSTAB", 2: "HACApK"}
    method_name = method_names.get(method, f"Method {method}")

    print("\n" + "="*80)
    print(f"SOLVER: {method_name} (Method {method})")
    print("="*80)

    # Configure HACApK if method 2
    if method == 2:
        rad.SolverConfig(hacapk_eps=1e-4, hacapk_leaf=10, hacapk_eta=2.0)

    # Solve
    start_time = time.time()
    result = rad.Solve(magnet, precision, max_iter, method)
    solve_time = time.time() - start_time

    iterations = int(result[1]) if isinstance(result, (list, tuple)) and len(result) > 1 else 1

    print(f"\nSolver result: {result}")
    print(f"  Solving time: {solve_time*1000:.1f} ms")
    print(f"  Iterations: {iterations}")

    # Check field at center
    B_center = rad.Fld(magnet, 'b', [0, 0, 0])
    print(f"  B at center:  [{B_center[0]:.6f}, {B_center[1]:.6f}, {B_center[2]:.6f}] T")

    return {
        'time': solve_time,
        'iterations': iterations,
        'result': result,
        'B_center': B_center
    }


def print_comparison(results):
    """
    Print comparison between solver methods.
    """
    print("\n" + "="*80)
    print("COMPARISON")
    print("="*80)

    # Reference is LU (Method 0)
    ref_time = results[0]['time'] if results[0]['time'] > 0 else 1e-10
    ref_B = results[0]['B_center']

    print(f"\n{'Method':<15} {'Time (ms)':<12} {'Iterations':<12} {'B_z (T)':<15} {'Speedup':<10}")
    print("-" * 70)

    for i, r in enumerate(results):
        method_names = ["LU", "BiCGSTAB", "HACApK"]
        speedup = ref_time / r['time'] if r['time'] > 0 else 0

        print(f"{method_names[i]:<15} {r['time']*1000:<12.1f} {r['iterations']:<12} "
              f"{r['B_center'][2]:<15.6f} {speedup:<10.2f}x")

    print("-" * 70)

    # Field accuracy check (vs LU reference)
    for i, r in enumerate(results):
        if i == 0:
            continue
        B_diff = [abs(ref_B[j] - r['B_center'][j]) for j in range(3)]
        B_rel_error = max(B_diff) / abs(ref_B[2]) * 100 if ref_B[2] != 0 else 0
        method_names = ["LU", "BiCGSTAB", "HACApK"]
        print(f"  {method_names[i]} vs LU error: {B_rel_error:.4f}%")


def main():
    """
    Main benchmark routine.
    """
    print("="*80)
    print("SOLVER BENCHMARK")
    print("="*80)
    print("\nThis benchmark compares three solver methods:")
    print("  Method 0: LU (direct solver)")
    print("  Method 1: BiCGSTAB (iterative solver)")
    print("  Method 2: HACApK (H-matrix accelerated)")
    print()

    # Test configuration
    n_per_side = 5  # 5x5x5 = 125 elements
    precision = 0.0001
    max_iter = 1000

    print(f"Configuration:")
    print(f"  Elements:     {n_per_side}x{n_per_side}x{n_per_side} = {n_per_side**3}")
    print(f"  Precision:    {precision}")
    print(f"  Max iter:     {max_iter}")
    print()

    results = []

    # Method 0: LU
    rad.UtiDelAll()
    magnet = create_magnet(n_per_side)
    results.append(benchmark_solver(magnet, method=0, precision=precision, max_iter=max_iter))

    # Method 1: BiCGSTAB
    rad.UtiDelAll()
    magnet = create_magnet(n_per_side)
    results.append(benchmark_solver(magnet, method=1, precision=precision, max_iter=max_iter))

    # Method 2: HACApK
    try:
        rad.UtiDelAll()
        magnet = create_magnet(n_per_side)
        results.append(benchmark_solver(magnet, method=2, precision=precision, max_iter=max_iter))
    except Exception as e:
        print(f"\nHACApK not available: {e}")
        results.append({'time': float('inf'), 'iterations': 0, 'result': None, 'B_center': [0, 0, 0]})

    # Print comparison
    print_comparison(results)

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print("\nSolver Selection Guidelines:")
    print("  - Method 0 (LU): Best for small problems (N < 200)")
    print("  - Method 1 (BiCGSTAB): Good for medium problems (200 < N < 2000)")
    print("  - Method 2 (HACApK): Best for large problems (N > 1000)")
    print()
    print("Usage:")
    print("  rad.Solve(obj, precision, max_iterations, method)")
    print("  rad.Solve(container, 0.0001, 1000, 2)  # HACApK for large problem")
    print()

    rad.UtiDelAll()


if __name__ == "__main__":
    main()
