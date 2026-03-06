#!/usr/bin/env python
"""
Netgen GUI Visualization Demo for Radia-NGSolve

Demonstrates:
- OCC shape visualization (exact geometry)
- Mesh visualization with quality check
- Lightweight native GUI (Tcl/Tk)
- Integrated mesh generation workflow

Advantages over webgui:
- ✅ Lightweight (no browser overhead)
- ✅ Native GUI (Tcl/Tk, fast and stable)
- ✅ Integrated workflow (shape → mesh → check)
- ✅ Better for shape/mesh verification
- ✅ Works from .py script (Jupyter not required)

Use case: Geometry verification, mesh quality check, CAD import verification

Requirements:
- NGSolve/Netgen with GUI support
- Python script execution (Jupyter not required)

Note: For field visualization, use webgui instead.
      netgen.gui is specialized for geometry/mesh, not field data.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import radia as rad
from netgen.occ import Box, Cylinder, Sphere, Pnt, OCCGeometry
from netgen.meshing import MeshingParameters
from netgen.gui import StartGUI


def radia_recmag_to_occ(center, dimensions):
    """
    Convert Radia ObjRecMag to OCC Box (exact geometry).

    Args:
        center: [x, y, z] in meters
        dimensions: [dx, dy, dz] in meters

    Returns:
        OCC Box shape
    """
    x0, y0, z0 = center
    dx, dy, dz = dimensions

    p1 = Pnt(x0 - dx/2, y0 - dy/2, z0 - dz/2)
    p2 = Pnt(x0 + dx/2, y0 + dy/2, z0 + dz/2)

    return Box(p1, p2)


def main():
    print("="*60)
    print("Netgen GUI Visualization Demo")
    print("="*60)

    rad.FldUnits('m')

    # Create Radia magnet
    print("\n[1] Creating Radia magnet...")
    center = [0, 0, 0]
    dimensions = [0.04, 0.04, 0.02]  # 40x40x20 mm
    magnetization = [0, 0, 954930]   # 1.2T NdFeB

    magnet = rad.ObjRecMag(center, dimensions, magnetization)
    print(f"    Magnet: {dimensions[0]*1000}x{dimensions[1]*1000}x{dimensions[2]*1000} mm, Br=1.2T")

    # Convert to OCC shape (exact geometry)
    print("\n[2] Converting to OCC shape...")
    occ_magnet = radia_recmag_to_occ(center, dimensions)
    print("    OCC Box created (exact geometry)")

    # Create OCC geometry
    print("\n[3] Creating OCCGeometry...")
    geo = OCCGeometry(occ_magnet)

    # Start Netgen GUI
    print("\n[4] Starting Netgen GUI...")
    print("    Netgen GUI window will open")
    print("\n    GUI操作:")
    print("    - File > Geometry で形状確認")
    print("    - Mesh > Generate Mesh でメッシュ生成")
    print("    - View > Mesh Quality でメッシュ品質確認")
    print("    - マウスドラッグで回転、ホイールでズーム")

    # Draw geometry in GUI
    StartGUI()
    geo.Draw()

    print("\n[5] Generating mesh...")
    mp = MeshingParameters(maxh=0.005)  # 5mm max element size
    mesh = geo.GenerateMesh(mp)
    print(f"    Mesh: {mesh.nv} vertices, {mesh.ne} elements")

    # Draw mesh in GUI
    print("\n[6] Displaying mesh in GUI...")
    mesh.Draw()

    print("\n[OK] Netgen GUI visualization complete")
    print("\nGUI features:")
    print("  - 3D rotation with mouse drag")
    print("  - Zoom with mouse wheel")
    print("  - Mesh quality visualization")
    print("  - Face/Edge highlighting")
    print("\nTo close: Close Netgen GUI window")


if __name__ == '__main__':
    main()
