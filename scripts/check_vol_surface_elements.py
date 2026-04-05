#!/usr/bin/env python
"""
Check if .vol file has surface elements for Netgen GUI display

Purpose:
    Verify if a .vol file can be properly displayed in Netgen GUI.
    Netgen GUI primarily displays surface elements, not volume elements.

Usage:
    python check_vol_surface_elements.py mesh.vol

Output:
    - Number of vertices
    - Number of volume elements (tet, hex, wedge, pyramid)
    - Number of surface elements (tri, quad)
    - Netgen GUI display compatibility

Author: Radia-NGSolve Project
"""

import sys
import os


def check_vol_file(vol_file):
    """Check .vol file for surface elements."""
    from netgen.meshing import Mesh

    if not os.path.exists(vol_file):
        print(f"Error: File not found: {vol_file}")
        return False

    print(f"Analyzing: {vol_file}")
    print("=" * 60)

    # Load mesh
    mesh = Mesh()
    try:
        mesh.Load(vol_file)
    except Exception as e:
        print(f"Error: Failed to load mesh: {e}")
        return False

    # Get statistics
    nv = mesh.nv
    ne = mesh.ne       # Volume elements
    nse = mesh.nse     # Surface elements

    print(f"\nMesh Statistics:")
    print(f"  Vertices:        {nv}")
    print(f"  Volume elements: {ne}")
    print(f"  Surface elements: {nse}")

    # Analyze volume element types
    if ne > 0:
        print(f"\n  Volume element types:")
        vol_types = {}
        type_names = {4: 'Tet', 5: 'Pyramid', 6: 'Wedge', 8: 'Hex'}
        for el in mesh.Elements3D():
            el_type = len(el.vertices)
            vol_types[el_type] = vol_types.get(el_type, 0) + 1
        for el_type, count in sorted(vol_types.items()):
            name = type_names.get(el_type, f'Unknown({el_type})')
            print(f"    {name}: {count}")

    # Analyze surface element types
    if nse > 0:
        print(f"\n  Surface element types:")
        surf_types = {}
        type_names = {3: 'Triangle', 4: 'Quad'}
        for el in mesh.Elements2D():
            el_type = len(el.vertices)
            surf_types[el_type] = surf_types.get(el_type, 0) + 1
        for el_type, count in sorted(surf_types.items()):
            name = type_names.get(el_type, f'Unknown({el_type})')
            print(f"    {name}: {count}")

    # Display compatibility
    print("\n" + "=" * 60)
    print("Display Compatibility:")
    print("=" * 60)

    if nse > 0:
        print("\n✅ Netgen GUI: COMPATIBLE")
        print(f"   - Surface elements present: {nse}")
        print("   - Mesh will be displayed as surface")
        print("   - Recommended viewer: Netgen GUI")
    else:
        print("\n⚠️  Netgen GUI: LIMITED")
        print("   - No surface elements")
        print("   - Only volume elements present")
        print("   - Netgen GUI may show nothing or require clipping plane")
        print("\n   Recommended alternatives:")
        print("   - ParaView: Can display volume mesh with slice/clip")
        print("   - PyVista: Can display volume mesh with slice/clip")

    if ne > 0 and nse == 0:
        print("\n   To generate surface elements:")
        print("   1. Use NGSolve to regenerate mesh with surface")
        print("   2. Or use ParaView/PyVista for volume visualization")

    return True


def main():
    if len(sys.argv) < 2:
        print("Usage: python check_vol_surface_elements.py <file.vol>")
        print("\nCheck if .vol file has surface elements for Netgen GUI display.")
        sys.exit(1)

    vol_file = sys.argv[1]
    check_vol_file(vol_file)


if __name__ == '__main__':
    main()
