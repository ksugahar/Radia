#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Demo: Export Radia magnetic field to VTS (Structured Grid) format

This example demonstrates how to:
1. Create a permanent magnet using ObjRecMag
2. Export magnetic field on a 3D grid to VTS format using rad.FldVTS()
3. Visualize in ParaView

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


def main():
    print("=" * 60)
    print("Radia VTS Export Demo")
    print("=" * 60)

    # Set units to meters (required for consistent results)

    # Create a permanent magnet
    # NdFeB magnet: Br = 1.2 T -> Mr = Br/mu_0 = 954930 A/m
    print("\n1. Creating permanent magnet...")

    # Rectangular magnet: 40mm x 40mm x 20mm
    magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])

    print(f"   Magnet handle: {magnet}")
    print("   Size: 40mm x 40mm x 20mm")
    print("   Magnetization: 954930 A/m (1.2 T equivalent)")

    # Define observation grid
    x_range = [-0.08, 0.08]   # -80mm to +80mm
    y_range = [-0.08, 0.08]   # -80mm to +80mm
    z_range = [0.03, 0.12]    # 30mm to 120mm (above magnet)
    nx, ny, nz = 33, 33, 19

    total_points = nx * ny * nz
    print(f"\n2. Exporting field grid to VTS...")
    print(f"   Grid: {nx}x{ny}x{nz} = {total_points} points")
    print(f"   X range: [{x_range[0]*1000:.0f}, {x_range[1]*1000:.0f}] mm")
    print(f"   Y range: [{y_range[0]*1000:.0f}, {y_range[1]*1000:.0f}] mm")
    print(f"   Z range: [{z_range[0]*1000:.0f}, {z_range[1]*1000:.0f}] mm")

    # Export using rad.FldVTS() - C++ implementation
    # Arguments: obj, filename, x_range, y_range, z_range, nx, ny, nz, include_B, include_H
    rad.FldVTS(
        magnet,
        'magnet_field.vts',
        x_range,
        y_range,
        z_range,
        nx, ny, nz,
        1,  # include_B
        1   # include_H
    )
    print("   -> magnet_field.vts")

    # Summary
    print("\n" + "=" * 60)
    print("Export complete!")
    print("\nGenerated file:")
    print("  - magnet_field.vts (VTS format, B and H fields)")
    print("\nTo view in ParaView:")
    print("  1. Open ParaView")
    print("  2. File -> Open -> magnet_field.vts")
    print("  3. Apply, then use Glyph filter for vector arrows")
    print("  4. Color by 'B_magnitude' for visualization")
    print("=" * 60)


if __name__ == '__main__':
    main()
