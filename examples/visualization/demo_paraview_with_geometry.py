#!/usr/bin/env python
"""
ParaView Visualization with Accurate Geometry

Problem: rad.FldVTS() approximates magnet shape with grid
Solution: Export geometry and field separately, overlay in ParaView

Workflow:
1. Export Radia geometry to STL/STEP (exact shape)
2. Export field to VTS (grid approximation)
3. Overlay in ParaView (shape + field)

Requirements:
- radia
- netgen/ngsolve (for OCC export)
- ParaView
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import radia as rad
from netgen.occ import Box, Pnt, OCCGeometry
from netgen.meshing import MeshingParameters


def export_radia_geometry_to_stl(radia_obj, stl_file):
    """
    Export Radia analytical object to STL with exact geometry.

    Note: This function converts ObjRecMag to OCC Box, then to STL.
          The OCC representation is exact (no approximation).
    """
    print(f"Exporting geometry to STL: {stl_file}")

    # Get magnet parameters (assuming ObjRecMag)
    # For now, manually create OCC Box matching the magnet
    # TODO: Implement rad.ExportOCC() for automatic conversion

    # Example: ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])
    # → OCC Box from [-0.02, -0.02, -0.01] to [0.02, 0.02, 0.01]

    box = Box(Pnt(-0.02, -0.02, -0.01), Pnt(0.02, 0.02, 0.01))
    geo = OCCGeometry(box)

    # Generate surface mesh (fine for accurate representation)
    mp = MeshingParameters(maxh=0.002)  # 2mm max element size
    mesh = geo.GenerateMesh(mp)

    # Export to STL
    mesh.Export(stl_file, "STL Format")
    print(f"[OK] STL exported: {stl_file}")


def main():
    print("="*60)
    print("ParaView Visualization with Accurate Geometry")
    print("="*60)

    rad.FldUnits('m')

    # Create magnet
    print("\n[1] Creating magnet...")
    magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])
    print("    Magnet: 40x40x20 mm, Br=1.2T (NdFeB)")

    # Export geometry to STL (exact shape)
    print("\n[2] Exporting geometry to STL...")
    stl_file = 'magnet_geometry.stl'
    export_radia_geometry_to_stl(magnet, stl_file)

    # Export field to VTS (grid approximation)
    print("\n[3] Exporting field to VTS...")
    vts_file = 'field_data.vts'
    rad.FldVTS(magnet, vts_file,
               [-0.1, 0.1], [-0.1, 0.1], [0.02, 0.15],
               41, 41, 27, 1, 0, 1.0)
    print(f"[OK] VTS exported: {vts_file}")

    # Instructions for ParaView
    print("\n" + "="*60)
    print("ParaView Overlay Instructions")
    print("="*60)
    print("\n1. Open ParaView")
    print("\n2. Load magnet geometry (exact shape):")
    print(f"   File > Open > {stl_file}")
    print("   Apply")
    print("   Display:")
    print("     - Representation: Surface")
    print("     - Opacity: 0.5 (semi-transparent)")
    print("     - Color: Gray or Red")

    print("\n3. Load field data (grid approximation):")
    print(f"   File > Open > {vts_file}")
    print("   Apply")
    print("   Display:")
    print("     - Coloring: B_magnitude")
    print("     - Colormap: coolwarm or rainbow")

    print("\n4. Add field vectors (optional):")
    print("   Select VTS data")
    print("   Filters > Glyph")
    print("   Glyph Type: Arrow")
    print("   Orientation Array: B_field")
    print("   Scale Array: B_magnitude")
    print("   Scale Factor: 0.02")

    print("\n5. Result:")
    print("   - Exact magnet geometry (STL, no approximation)")
    print("   - Field distribution overlay (VTS)")
    print("   - Publication-quality figure")

    print("\n[OK] Files ready for ParaView")


if __name__ == '__main__':
    main()
