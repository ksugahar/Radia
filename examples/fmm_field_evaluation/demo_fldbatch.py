#!/usr/bin/env python
"""
Demo: Batch Field Evaluation with rad.Fld()

This example demonstrates how to use rad.Fld() with batch points to efficiently
compute magnetic fields at many observation points. Passing an (N,3) array is
significantly faster than calling rad.Fld() in a loop because:

1. Single Python->C++ call overhead instead of N calls
2. TaskManager parallelization across evaluation points

Use Case:
- Netgen/NGSolve tetrahedral mesh -> Radia solve -> air field distribution
- Visualization grid (1000+ points)
- Trajectory field computation

Author: Radia Development Team
Date: 2025-12-30
"""

import sys
import os
import time
import numpy as np

# Add Radia path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import radia as rad


def create_magnetized_cube(size=0.1, n_div=3, magnetization=[0, 0, 1e6]):
    """Create a subdivided cube with uniform magnetization.

    Args:
        size: Cube side length in meters
        n_div: Number of divisions per side
        magnetization: [Mx, My, Mz] in A/m

    Returns:
        container: Radia container object with n_div^3 hexahedral elements
    """
    rad.UtiDelAll()

    half = size / 2
    step = size / n_div

    hex_objs = []
    for ix in range(n_div):
        for iy in range(n_div):
            for iz in range(n_div):
                x0 = -half + ix * step
                y0 = -half + iy * step
                z0 = -half + iz * step
                x1 = x0 + step
                y1 = y0 + step
                z1 = z0 + step

                vertices = [
                    [x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                    [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]
                ]
                obj = rad.ObjHexahedron(vertices, magnetization)
                hex_objs.append(obj)

    container = rad.ObjCnt(hex_objs)
    return container


def create_observation_grid(center, extent, n_points_per_axis):
    """Create a 3D grid of observation points.

    Args:
        center: [x, y, z] center of grid
        extent: Half-size of grid in each direction
        n_points_per_axis: Number of points per axis

    Returns:
        points: numpy array of shape (N, 3)
        shape: (nx, ny, nz) tuple for reshaping
    """
    cx, cy, cz = center
    x = np.linspace(cx - extent, cx + extent, n_points_per_axis)
    y = np.linspace(cy - extent, cy + extent, n_points_per_axis)
    z = np.linspace(cz - extent, cz + extent, n_points_per_axis)

    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
    points = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=1)

    return points, (n_points_per_axis, n_points_per_axis, n_points_per_axis)


def main():
    print("=" * 70)
    print("Batch Field Evaluation Demo")
    print("=" * 70)

    # Create magnetized cube
    size = 0.1  # 10 cm cube
    n_div = 3   # 3x3x3 = 27 elements
    M = [0, 0, 1e6]  # 1 MA/m in z-direction

    print(f"\n1. Creating magnetized cube:")
    print(f"   Size: {size*1000:.0f} mm")
    print(f"   Divisions: {n_div} x {n_div} x {n_div} = {n_div**3} elements")
    print(f"   Magnetization: {M} A/m")

    t0 = time.time()
    container = create_magnetized_cube(size, n_div, M)
    t_create = time.time() - t0
    print(f"   Creation time: {t_create:.3f} s")

    # Create observation grid (in air, outside the magnet)
    center = [0, 0, 0.15]  # 15 cm above center (in air)
    extent = 0.1  # +/- 10 cm
    n_pts = 10    # 10x10x10 = 1000 points

    points, shape = create_observation_grid(center, extent, n_pts)
    n_total = len(points)

    print(f"\n2. Creating observation grid:")
    print(f"   Center: {center} m")
    print(f"   Extent: +/- {extent*1000:.0f} mm")
    print(f"   Grid: {n_pts} x {n_pts} x {n_pts} = {n_total} points")

    # Method 1: Batch Fld (recommended for many points)
    print(f"\n3. Computing B field with batch Fld()...")
    t0 = time.time()
    B_batch = np.asarray(rad.Fld(container, 'b', points))
    t_batch = time.time() - t0

    print(f"   Time: {t_batch:.3f} s")
    print(f"   Points/sec: {n_total/t_batch:.0f}")

    # Method 2: Loop with single-point Fld (for comparison)
    print(f"\n4. Computing field with Fld loop (for comparison)...")
    t0 = time.time()
    B_loop = []
    for pt in points[:100]:  # Only first 100 points for speed
        B = rad.Fld(container, 'b', pt)
        B_loop.append(B)
    t_loop = time.time() - t0
    t_loop_estimated = t_loop * n_total / 100

    print(f"   Time for 100 points: {t_loop:.3f} s")
    print(f"   Estimated time for {n_total} points: {t_loop_estimated:.3f} s")
    print(f"   Speedup: {t_loop_estimated/t_batch:.1f}x")

    # Verify results match
    B_loop = np.array(B_loop)
    B_batch_100 = B_batch[:100]
    max_diff = np.max(np.abs(B_batch_100 - B_loop))
    print(f"\n5. Verification:")
    print(f"   Max difference (first 100 pts): {max_diff:.2e} T")

    # Field statistics
    B_mag = np.linalg.norm(B_batch, axis=1)
    print(f"\n6. Field statistics (all {n_total} points):")
    print(f"   |B| min: {B_mag.min():.6e} T")
    print(f"   |B| max: {B_mag.max():.6e} T")
    print(f"   |B| mean: {B_mag.mean():.6e} T")

    # Bz at center of grid
    center_idx = n_total // 2
    B_center = B_batch[center_idx]
    print(f"\n7. Field at grid center ({center} m):")
    print(f"   B = [{B_center[0]:.6e}, {B_center[1]:.6e}, {B_center[2]:.6e}] T")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Batch Fld() is {t_loop_estimated/t_batch:.1f}x faster than single-point loop")
    print(f"for {n_total} evaluation points.")
    print()
    print("Use batch Fld(obj, 'b', points_Nx3) when:")
    print("  - Computing field at many observation points (>100)")
    print("  - Visualization grids")
    print("  - Trajectory field evaluation")
    print("=" * 70)


if __name__ == '__main__':
    main()
