#!/usr/bin/env python
"""
Demonstration: Netgen automatically generates surface elements

Key point: When you generate a mesh with Netgen/NGSolve, surface elements
           are AUTOMATICALLY created. You don't need to do anything special.

Result: The generated .vol file can be displayed in Netgen GUI without issues.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

from netgen.occ import Box, Pnt, OCCGeometry
from netgen.meshing import MeshingParameters
from netgen.gui import StartGUI


def main():
    print("="*60)
    print("Demonstration: Netgen Automatic Surface Element Generation")
    print("="*60)

    # Step 1: Create OCC geometry
    print("\n[1] Creating OCC geometry (Box)...")
    box = Box(Pnt(-0.05, -0.05, -0.05), Pnt(0.05, 0.05, 0.05))
    geo = OCCGeometry(box)
    print("    Geometry created")

    # Step 2: Generate mesh
    print("\n[2] Generating mesh...")
    mp = MeshingParameters(maxh=0.01)  # 10mm max element size
    mesh = geo.GenerateMesh(mp)
    print(f"    Mesh generated:")
    print(f"      Vertices: {mesh.nv}")
    print(f"      Volume elements: {mesh.ne}")
    print(f"      Surface elements: {mesh.nse}")  # <- Automatically created!

    # Step 3: Verify surface elements
    print("\n[3] Verifying surface elements...")
    if mesh.nse > 0:
        print("    ✅ Surface elements present!")
        print(f"       Count: {mesh.nse}")
        print("       Netgen GUI will display this mesh correctly")

        # Analyze surface element types
        surf_types = {}
        for el in mesh.Elements2D():
            el_type = len(el.vertices)
            surf_types[el_type] = surf_types.get(el_type, 0) + 1

        print("\n    Surface element types:")
        type_names = {3: 'Triangle', 4: 'Quad'}
        for el_type, count in sorted(surf_types.items()):
            name = type_names.get(el_type, f'Unknown({el_type})')
            print(f"      {name}: {count}")
    else:
        print("    ❌ No surface elements (unexpected!)")

    # Step 4: Save to .vol
    print("\n[4] Saving to .vol file...")
    vol_file = 'test_mesh_with_surface.vol'
    mesh.Save(vol_file)
    print(f"    Saved: {vol_file}")
    print("    This .vol file includes surface elements")
    print("    Can be opened in Netgen GUI by double-clicking")

    # Step 5: Display in Netgen GUI
    print("\n[5] Opening in Netgen GUI...")
    print("    Netgen GUI window will open")
    StartGUI()
    mesh.Draw()

    print("\n[OK] Demonstration complete")
    print("\nKey takeaway:")
    print("  When you generate a mesh with Netgen/NGSolve,")
    print("  surface elements are AUTOMATICALLY created.")
    print("  You don't need to worry about surface elements.")


if __name__ == '__main__':
    main()
