#!/usr/bin/env python
"""
FMM Field Calculation Benchmark for Radia

Benchmarks field computation performance for Radia.
Currently measures direct field computation; FMM acceleration to be added.

This script is designed to match ELF_MAGIC's fmm_field_benchmark for comparison.

Usage:
    python benchmark_fmm_field.py --hex 10
    python benchmark_fmm_field.py --tetra 0.15
    python benchmark_fmm_field.py --hex 10 --n_grid 30

Author: Radia Development Team
Date: 2025-12-30
"""

import sys
import os
import time
import json
import argparse
import numpy as np

# Path setup
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))
import radia as rad

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False

try:
    import netgen.occ as occ
    from ngsolve import Mesh
    HAS_NETGEN = True
except ImportError:
    HAS_NETGEN = False


# =============================================================================
# Benchmark Conditions (same as ELF)
# =============================================================================
MU_0 = 4 * np.pi * 1e-7  # Vacuum permeability [T/(A/m)]
CUBE_SIZE = 1.0          # Cube edge length [m]
MU_R = 1000              # Relative permeability
CHI = MU_R - 1           # Magnetic susceptibility
H_EXT = 200000.0         # External H field [A/m]

# Linear BH curve: B = mu_r * mu_0 * H
BH_DATA = [
    [0.0, 0.0],
    [1000000.0, MU_R * MU_0 * 1000000.0]
]



def get_memory_mb():
    """Get current memory usage in MB."""
    if not HAS_PSUTIL:
        return None
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 * 1024)


def create_hex_mesh(n_div, size=1.0):
    """Create hexahedral mesh for a cube using Radia ObjPolyhdr."""
    rad.UtiDelAll()
    rad.FldUnits('m')

    dx = size / n_div
    offset = size / 2
    elements = []

    for iz in range(n_div):
        for iy in range(n_div):
            for ix in range(n_div):
                x0 = ix * dx - offset
                y0 = iy * dx - offset
                z0 = iz * dx - offset

                vertices = [
                    [x0, y0, z0],
                    [x0 + dx, y0, z0],
                    [x0 + dx, y0 + dx, z0],
                    [x0, y0 + dx, z0],
                    [x0, y0, z0 + dx],
                    [x0 + dx, y0, z0 + dx],
                    [x0 + dx, y0 + dx, z0 + dx],
                    [x0, y0 + dx, z0 + dx]
                ]

                elem = rad.ObjHexahedron(vertices, [0, 0, 0])
                elements.append(elem)

    return elements


def create_tetra_mesh(maxh, size=1.0):
    """Create tetrahedral mesh for a cube using Netgen."""
    if not HAS_NETGEN:
        raise RuntimeError("Netgen is required for tetrahedral mesh generation")

    rad.UtiDelAll()
    rad.FldUnits('m')

    half = size / 2
    box = occ.Box(occ.Pnt(-half, -half, -half), occ.Pnt(half, half, half))
    geo = occ.OCCGeometry(box)
    ngmesh = geo.GenerateMesh(maxh=maxh)
    mesh = Mesh(ngmesh)

    # Extract vertices
    vertices_np = np.array([[v.point[i] for i in range(3)] for v in mesh.vertices])

    # Create tetrahedral elements
    elements = []
    for el in mesh.Elements3D():
        node_ids = [v.nr for v in el.vertices]
        verts = [list(vertices_np[i]) for i in node_ids]
        elem = rad.ObjTetrahedron(verts, [0, 0, 0])
        elements.append(elem)

    return elements


def generate_field_grid(center, extent, n_grid):
    """Generate a 3D grid of field evaluation points."""
    x = np.linspace(center[0] - extent/2, center[0] + extent/2, n_grid)
    y = np.linspace(center[1] - extent/2, center[1] + extent/2, n_grid)
    z = np.linspace(center[2] - extent/2, center[2] + extent/2, n_grid)

    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
    points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])

    return points


def classify_points_simple(points, cube_half=0.5, margin=0.5):
    """
    Simple point classification (inside/near/far).

    Returns mask of points that are outside the cube + margin zone.
    """
    exclusion_half = cube_half + margin
    inside_mask = (
        (np.abs(points[:, 0]) < exclusion_half) &
        (np.abs(points[:, 1]) < exclusion_half) &
        (np.abs(points[:, 2]) < exclusion_half)
    )
    return ~inside_mask


def compute_field_batch(container, points):
    """
    Compute B and H field at multiple points.

    Args:
        container: Radia container object
        points: numpy array of shape (N, 3)

    Returns:
        B: numpy array of shape (N, 3)
        H: numpy array of shape (N, 3)
    """
    n_points = len(points)
    B = np.zeros((n_points, 3))
    H = np.zeros((n_points, 3))

    for i, pt in enumerate(points):
        B[i] = rad.Fld(container, 'b', list(pt))
        H[i] = rad.Fld(container, 'h', list(pt))

    return B, H


def run_benchmark(mesh_type, mesh_param, n_grid):
    """Run field computation benchmark.

    Args:
        mesh_type: 'hex' or 'tetra'
        mesh_param: n_div for hex, maxh for tetra
        n_grid: Number of grid points per axis
    """
    if mesh_type == 'hex':
        mesh_desc = f"Hex N={mesh_param}"
    else:
        mesh_desc = f"Tetra maxh={mesh_param}"

    print(f"\n{'='*70}")
    print(f"Benchmark: {mesh_desc}, Grid={n_grid}x{n_grid}x{n_grid} ({n_grid**3} points)")
    print(f"{'='*70}")

    # Create mesh
    t_start = time.perf_counter()
    if mesh_type == 'hex':
        elements = create_hex_mesh(mesh_param, CUBE_SIZE)
    else:
        elements = create_tetra_mesh(mesh_param, CUBE_SIZE)

    n_elements = len(elements)

    # Create container
    container = rad.ObjCnt(elements)

    # Apply material
    mat = rad.MatSatIsoTab(BH_DATA)
    rad.MatApl(container, mat)

    # Add external field (B_ext = mu0 * H_ext)
    B_ext = MU_0 * H_EXT
    ext = rad.ObjBckg([0, 0, B_ext])
    grp = rad.ObjCnt([container, ext])

    t_mesh = time.perf_counter() - t_start

    # Get DOF
    ndof = n_elements * (6 if mesh_type == 'hex' else 3)
    print(f"Mesh: {n_elements} elements, {ndof} DOF (created in {t_mesh:.3f}s)")

    # Solve
    MAX_ITER = 100
    t_start = time.perf_counter()
    result = rad.Solve(grp, 0.001, MAX_ITER, 1)  # BiCGSTAB
    t_solve = time.perf_counter() - t_start

    # Get solve statistics
    stats = rad.GetSolveStats()
    n_iter = stats.get('nonl_iterations', 0)
    converged = n_iter < MAX_ITER
    residual = result[0] if result[0] else 0.0
    print(f"Solve: {t_solve:.3f}s (converged={converged}, iterations={n_iter}, residual={residual:.2e})")

    # Generate field evaluation grid
    grid_extent = 3.0
    points = generate_field_grid((0.0, 0.0, 0.0), grid_extent, n_grid)
    n_total_points = len(points)
    print(f"Field grid: {n_total_points} points total")

    # Filter out points too close to the cube
    min_distance = 0.5
    outside_mask = classify_points_simple(points, cube_half=0.5, margin=min_distance)
    outside_points = points[outside_mask]
    n_outside = len(outside_points)
    print(f"Points outside exclusion zone (>{min_distance}m from cube): {n_outside}")

    # Benchmark direct computation
    print("\n--- Direct Computation ---")
    t_start = time.perf_counter()
    B_direct, H_direct = compute_field_batch(container, outside_points)
    t_direct = time.perf_counter() - t_start

    points_per_sec = n_outside / t_direct
    print(f"Time: {t_direct:.3f}s ({points_per_sec:.0f} points/sec)")

    # Field statistics
    B_mag = np.linalg.norm(B_direct, axis=1)
    H_mag = np.linalg.norm(H_direct, axis=1)

    print(f"\n--- Field Statistics ---")
    print(f"|B| range: {np.min(B_mag):.4f} - {np.max(B_mag):.4f} T")
    print(f"|H| range: {np.min(H_mag):.0f} - {np.max(H_mag):.0f} A/m")

    # Prepare results
    results = {
        'mesh_type': mesh_type,
        'mesh_param': mesh_param,
        'mesh_desc': mesh_desc,
        'n_elements': n_elements,
        'ndof': ndof,
        'n_grid': n_grid,
        'n_points': n_outside,
        'H_ext': H_EXT,
        'chi': CHI,
        't_mesh': t_mesh,
        't_solve': t_solve,
        't_direct': t_direct,
        'direct_points_per_sec': points_per_sec,
        'B_max': float(np.max(B_mag)),
        'B_min': float(np.min(B_mag)),
        'H_max': float(np.max(H_mag)),
        'H_min': float(np.min(H_mag)),
        'converged': converged,
        'iterations': n_iter,
        'residual': residual,
    }

    # Memory info
    if HAS_PSUTIL:
        results['peak_memory_mb'] = get_memory_mb()

    return results


def main():
    parser = argparse.ArgumentParser(description='Field Computation Benchmark')
    parser.add_argument('--hex', type=int, nargs='*',
                        help='Hexahedral mesh: N divisions (e.g., --hex 5 10 15)')
    parser.add_argument('--tetra', type=float, nargs='*',
                        help='Tetrahedral mesh: maxh values (e.g., --tetra 0.2 0.15 0.1)')
    parser.add_argument('--n_grid', type=int, default=20,
                        help='Number of grid points per axis for field evaluation')
    parser.add_argument('--output', type=str, default='fmm_benchmark_results.json',
                        help='Output JSON file')
    args = parser.parse_args()

    # Default: hex N=5, 10 if nothing specified
    if args.hex is None and args.tetra is None:
        args.hex = [5, 10]

    print("="*70)
    print("FIELD COMPUTATION BENCHMARK (Radia)")
    print("="*70)
    if args.hex:
        print(f"Hexahedral meshes: N = {args.hex}")
    if args.tetra:
        print(f"Tetrahedral meshes: maxh = {args.tetra}")
    print(f"Grid size: {args.n_grid}x{args.n_grid}x{args.n_grid}")

    all_results = []

    # Run hexahedral benchmarks
    if args.hex:
        for n_div in args.hex:
            try:
                result = run_benchmark('hex', n_div, args.n_grid)
                all_results.append(result)
            except Exception as e:
                print(f"Error for Hex N={n_div}: {e}")
                import traceback
                traceback.print_exc()

    # Run tetrahedral benchmarks
    if args.tetra:
        if not HAS_NETGEN:
            print("Warning: Netgen not available, skipping tetrahedral benchmarks")
        else:
            for maxh in args.tetra:
                try:
                    result = run_benchmark('tetra', maxh, args.n_grid)
                    all_results.append(result)
                except Exception as e:
                    print(f"Error for Tetra maxh={maxh}: {e}")
                    import traceback
                    traceback.print_exc()

    # Save results
    output_path = os.path.join(os.path.dirname(__file__), args.output)
    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nResults saved to: {output_path}")

    # Summary table
    if all_results:
        print("\n" + "="*70)
        print("SUMMARY")
        print("="*70)
        print(f"{'Mesh':>15} {'Elements':>10} {'Points':>10} {'Direct':>10} {'pts/sec':>12}")
        print("-"*60)
        for r in all_results:
            print(f"{r['mesh_desc']:>15} {r['n_elements']:>10} {r['n_points']:>10} "
                  f"{r['t_direct']:>9.3f}s {r['direct_points_per_sec']:>11.0f}")


if __name__ == '__main__':
    main()
