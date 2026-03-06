#!/usr/bin/env python
"""
View Netgen .vol mesh (recommended methods).

Usage:
    python view_vol.py mesh.vol                # NGSolve webgui (default)
    python view_vol.py mesh.vol --netgen       # Netgen GUI
    python view_vol.py mesh.vol --vtk          # Export to VTK
    python view_vol.py mesh.vol --info         # Print info only

Recommended: Use NGSolve webgui (browser-based, interactive)
"""

import sys
import os
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))


def view_webgui(vol_file):
    """NGSolve webgui (RECOMMENDED)."""
    from ngsolve import Mesh
    from ngsolve.webgui import Draw

    mesh = Mesh(vol_file)
    print(f"Mesh: {vol_file}")
    print(f"  Elements: {mesh.ne}, Vertices: {mesh.nv}")
    print(f"  Materials: {mesh.GetMaterials()}")
    print("\nOpening webgui (browser)...")

    Draw(mesh, name='mesh')
    input("Press Enter to close...")


def view_netgen_gui(vol_file):
    """Netgen native GUI."""
    try:
        import netgen.gui
        from netgen.meshing import Mesh as NetgenMesh

        ngmesh = NetgenMesh()
        ngmesh.Load(vol_file)
        print(f"Loaded: {vol_file}")
        print("Opening Netgen GUI...")

        netgen.gui.StartGUI()
        input("Press Enter to close...")

    except ImportError:
        print("Netgen GUI not available. Use webgui instead:")
        print(f"  python view_vol.py {vol_file}")


def export_vtk(vol_file):
    """Export to VTK for ParaView."""
    from ngsolve import Mesh, VTKOutput

    mesh = Mesh(vol_file)
    vtk_file = str(Path(vol_file).with_suffix(''))

    print(f"Exporting: {vtk_file}.vtk")
    vtk = VTKOutput(mesh, coefs=[], names=[], filename=vtk_file)
    vtk.Do()

    print(f"[OK] Exported: {vtk_file}.vtk")
    print(f"View: paraview {vtk_file}.vtk")


def print_info(vol_file):
    """Print mesh information."""
    from ngsolve import Mesh, VOL
    import numpy as np

    mesh = Mesh(vol_file)

    print(f"Mesh: {vol_file}")
    print(f"  Elements: {mesh.ne}")
    print(f"  Vertices: {mesh.nv}")
    print(f"  Materials: {mesh.GetMaterials()}")

    # Element types
    types = {}
    for el in mesh.Elements(VOL):
        t = str(el.type)
        types[t] = types.get(t, 0) + 1

    print("  Element types:")
    for t, c in types.items():
        print(f"    {t}: {c}")

    # Bounding box
    pts = np.array([[v.point[i] for i in range(3)] for v in mesh.vertices])
    for i, axis in enumerate(['X', 'Y', 'Z']):
        vmin, vmax = pts[:, i].min(), pts[:, i].max()
        print(f"  {axis}: [{vmin:.6f}, {vmax:.6f}] (size: {vmax-vmin:.6f})")


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='View .vol mesh')
    parser.add_argument('vol_file')
    parser.add_argument('--netgen', action='store_true', help='Netgen GUI')
    parser.add_argument('--vtk', action='store_true', help='Export VTK')
    parser.add_argument('--info', action='store_true', help='Info only')

    args = parser.parse_args()

    if not os.path.exists(args.vol_file):
        print(f"Error: {args.vol_file} not found")
        sys.exit(1)

    if args.info:
        print_info(args.vol_file)
    elif args.netgen:
        view_netgen_gui(args.vol_file)
    elif args.vtk:
        export_vtk(args.vol_file)
    else:
        view_webgui(args.vol_file)
