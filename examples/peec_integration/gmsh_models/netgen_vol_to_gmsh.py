#!/usr/bin/env python
"""
Convert Netgen .vol mesh to Gmsh .msh format for visualization.

Usage:
    python vol_to_gmsh.py input.vol [output.msh]

Features:
- Converts Netgen volume mesh (.vol) to Gmsh format (.msh)
- Preserves element types (hex, tet, wedge)
- Supports NGSolve mesh with materials
- Can open in Gmsh for visualization

Requirements:
- NGSolve / Netgen
- Gmsh (for viewing)
"""

import sys
import os
from pathlib import Path

# Add Radia to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))

def convert_vol_to_gmsh(vol_file, gmsh_file=None, format='gmsh2'):
    """
    Convert Netgen .vol to Gmsh .msh format.

    Parameters:
    -----------
    vol_file : str
        Input .vol file path
    gmsh_file : str, optional
        Output .msh file path (default: same as input with .msh extension)
    format : str
        'gmsh' or 'gmsh2' (default: gmsh2 for modern Gmsh)
    """
    from netgen.meshing import Mesh as NetgenMesh

    # Load Netgen mesh
    print(f"Loading: {vol_file}")
    ngmesh = NetgenMesh()
    ngmesh.Load(vol_file)

    # Get mesh info
    nv = ngmesh.GetNV()
    ne = ngmesh.GetNE()
    nse = ngmesh.GetNSE()

    print(f"  Vertices: {nv}")
    print(f"  Volume elements: {ne}")
    print(f"  Surface elements: {nse}")

    # Determine output filename
    if gmsh_file is None:
        vol_path = Path(vol_file)
        gmsh_file = str(vol_path.with_suffix('.msh'))

    # Export to Gmsh format
    print(f"Exporting to Gmsh format: {gmsh_file}")

    # Note: Netgen Export() uses format string to determine output format
    # For Gmsh: use 'Gmsh Format' or 'Gmsh2 Format'
    if format == 'gmsh2':
        format_str = 'Gmsh2 Format'
    else:
        format_str = 'Gmsh Format'

    ngmesh.Export(gmsh_file, format_str)

    print(f"✓ Conversion complete: {gmsh_file}")
    print(f"\nTo view in Gmsh:")
    print(f"  gmsh {gmsh_file}")

    return gmsh_file


def convert_vol_to_vtk(vol_file, vtk_file=None):
    """
    Convert Netgen .vol to VTK format (alternative to Gmsh).

    Uses NGSolve VTKOutput for conversion.
    """
    from ngsolve import Mesh, VTKOutput

    print(f"Loading NGSolve mesh: {vol_file}")
    mesh = Mesh(vol_file)

    # Determine output filename
    if vtk_file is None:
        vol_path = Path(vol_file)
        vtk_file = str(vol_path.with_suffix(''))  # VTKOutput adds .vtk

    print(f"Exporting to VTK: {vtk_file}.vtk")

    # Export using VTKOutput
    vtk = VTKOutput(mesh, coefs=[], names=[], filename=vtk_file)
    vtk.Do()

    print(f"✓ VTK export complete: {vtk_file}.vtk")
    print(f"\nTo view in ParaView:")
    print(f"  paraview {vtk_file}.vtk")

    return vtk_file + '.vtk'


def view_in_netgen(vol_file):
    """
    Open .vol file in Netgen GUI (if available).
    """
    try:
        from netgen.meshing import Mesh as NetgenMesh
        import netgen.gui

        print(f"Opening in Netgen GUI: {vol_file}")
        ngmesh = NetgenMesh()
        ngmesh.Load(vol_file)

        # This requires Netgen GUI support
        netgen.gui.StartGUI()
        netgen.gui.SetMesh(ngmesh)

    except ImportError:
        print("Netgen GUI not available.")
        print("Use Gmsh or ParaView instead:")
        print(f"  python vol_to_gmsh.py {vol_file}")


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Convert Netgen .vol to Gmsh .msh or VTK format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert to Gmsh (default)
  python vol_to_gmsh.py mesh.vol

  # Convert to Gmsh with custom output name
  python vol_to_gmsh.py mesh.vol output.msh

  # Convert to VTK instead
  python vol_to_gmsh.py mesh.vol --vtk

  # View in Gmsh after conversion
  python vol_to_gmsh.py mesh.vol --view
        """
    )

    parser.add_argument('vol_file', help='Input .vol file')
    parser.add_argument('output_file', nargs='?', help='Output .msh or .vtk file (optional)')
    parser.add_argument('--vtk', action='store_true', help='Export to VTK instead of Gmsh')
    parser.add_argument('--format', choices=['gmsh', 'gmsh2'], default='gmsh2',
                        help='Gmsh format version (default: gmsh2)')
    parser.add_argument('--view', action='store_true', help='Open in Gmsh after conversion')

    args = parser.parse_args()

    # Check input file exists
    if not os.path.exists(args.vol_file):
        print(f"Error: File not found: {args.vol_file}")
        sys.exit(1)

    # Convert
    if args.vtk:
        output_file = convert_vol_to_vtk(args.vol_file, args.output_file)
    else:
        output_file = convert_vol_to_gmsh(args.vol_file, args.output_file, args.format)

    # View in Gmsh if requested
    if args.view and not args.vtk:
        print(f"\nLaunching Gmsh...")
        import subprocess
        try:
            subprocess.run(['gmsh', output_file])
        except FileNotFoundError:
            print("Error: Gmsh not found in PATH")
            print("Install Gmsh: https://gmsh.info/")


if __name__ == '__main__':
    main()
