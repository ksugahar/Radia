#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo: Export Radia magnetic field to VTS (Structured Grid) format

This example demonstrates how to:
1. Create a permanent magnet using ObjRecMag
2. Export magnetic field on a 3D grid to VTS format using C++ API
3. Export with Python API for flexibility
4. Visualize in ParaView

The VTS format is ideal for structured 3D grids and is efficient for
visualization of B and H fields in ParaView.

Author: Radia Development Team
Date: 2026-01-09
"""

import sys
import os

# Add package directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'radia'))

import radia as rad
from radia import RadiaVTKOutput, export_field_grid_vts


def main():
    print("=" * 60)
    print("Radia VTS Export Demo")
    print("=" * 60)

    # Set units to meters (required for consistent results)
    rad.FldUnits('m')

    # Create a permanent magnet
    # NdFeB magnet: Br = 1.2 T -> Mr = Br/mu_0 = 954930 A/m
    print("\n1. Creating permanent magnet...")

    # Rectangular magnet: 40mm x 40mm x 20mm
    magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])

    print(f"   Magnet handle: {magnet}")
    print("   Size: 40mm x 40mm x 20mm")
    print("   Magnetization: 954930 A/m (1.2 T equivalent)")

    # Define observation grid
    # Grid covers the region above the magnet (z = 20mm to 100mm)
    grid_params = {
        'x_range': [-0.08, 0.08],   # -80mm to +80mm
        'y_range': [-0.08, 0.08],   # -80mm to +80mm
        'z_range': [0.03, 0.12],    # 30mm to 120mm (above magnet)
        'nx': 33,                    # Grid points in x
        'ny': 33,                    # Grid points in y
        'nz': 19                     # Grid points in z
    }

    total_points = grid_params['nx'] * grid_params['ny'] * grid_params['nz']
    print(f"\n2. Exporting field grid to VTS (C++ API)...")
    print(f"   Grid: {grid_params['nx']}x{grid_params['ny']}x{grid_params['nz']} = {total_points} points")
    print(f"   X range: [{grid_params['x_range'][0]*1000:.0f}, {grid_params['x_range'][1]*1000:.0f}] mm")
    print(f"   Y range: [{grid_params['y_range'][0]*1000:.0f}, {grid_params['y_range'][1]*1000:.0f}] mm")
    print(f"   Z range: [{grid_params['z_range'][0]*1000:.0f}, {grid_params['z_range'][1]*1000:.0f}] mm")

    # Method 1: Using C++ API (recommended for performance)
    rad.FldVTS(
        magnet,
        'magnet_field_cpp',
        grid_params['x_range'],
        grid_params['y_range'],
        grid_params['z_range'],
        grid_params['nx'],
        grid_params['ny'],
        grid_params['nz'],
        1,  # include_B
        1   # include_H
    )
    print("   -> magnet_field_cpp.vts (C++ export)")

    # Method 2: Using Python convenience function
    print("\n3. Exporting field grid to VTS (Python API)...")
    export_field_grid_vts(
        magnet,
        grid_params,
        'magnet_field_py',
        coefs=['B', 'B_magnitude', 'Bz']
    )
    print("   -> magnet_field_py.vts (Python export)")

    # Method 3: Using RadiaVTKOutput class (more flexible)
    print("\n4. Exporting with custom coefficients...")
    vtk = RadiaVTKOutput(
        obj=magnet,
        coefs=['B', 'H', 'B_magnitude', 'H_magnitude'],
        filename='magnet_field_full',
        grid_params=grid_params,
        floatsize='double'
    )
    vtk.Do()
    print("   -> magnet_field_full.vts (class-based export)")

    # Summary
    print("\n" + "=" * 60)
    print("Export complete!")
    print("\nGenerated files:")
    print("  - magnet_field_cpp.vts    (VTS - C++ API, B and H fields)")
    print("  - magnet_field_py.vts     (VTS - Python API, B field)")
    print("  - magnet_field_full.vts   (VTS - Class API, B and H fields)")
    print("\nTo view in ParaView:")
    print("  1. Open ParaView")
    print("  2. File -> Open -> magnet_field_cpp.vts")
    print("  3. Apply, then use Glyph filter for vector arrows")
    print("  4. Use Contour filter for B_magnitude isosurfaces")
    print("=" * 60)


if __name__ == '__main__':
    main()
