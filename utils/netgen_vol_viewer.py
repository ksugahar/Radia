#!/usr/bin/env python
"""
Netgen .vol File Viewer (Windows File Association Wrapper)

Purpose:
    Open .vol files in Netgen GUI by double-clicking in Windows Explorer

Usage:
    1. Direct execution:
       python netgen_vol_viewer.py mesh.vol

    2. Windows file association:
       - Right-click .vol file
       - "Open with" → "Choose another app"
       - "More apps" → "Look for another app on this PC"
       - Select: python.exe (or pythonw.exe for no console)
       - Check "Always use this app to open .vol files"
       - Command line: python.exe "path\to\netgen_vol_viewer.py" "%1"

    3. Registry (advanced):
       HKEY_CLASSES_ROOT\.vol
         (Default) = "NetgenMeshFile"
       HKEY_CLASSES_ROOT\NetgenMeshFile\shell\open\command
         (Default) = "pythonw.exe" "C:\path\to\netgen_vol_viewer.py" "%1"

Requirements:
    - NGSolve/Netgen with GUI support
    - Python in PATH

Author: Radia-NGSolve Project
"""

import sys
import os


def main():
    if len(sys.argv) < 2:
        print("Usage: python netgen_vol_viewer.py <file.vol>")
        print("\nThis script is designed to be used as Windows file association handler.")
        print("To open a .vol file:")
        print("  python netgen_vol_viewer.py mesh.vol")
        sys.exit(1)

    vol_file = sys.argv[1]

    # Check file exists
    if not os.path.exists(vol_file):
        print(f"Error: File not found: {vol_file}")
        sys.exit(1)

    # Import Netgen
    try:
        from netgen.meshing import Mesh
        from netgen.gui import StartGUI
    except ImportError as e:
        print(f"Error: Cannot import Netgen: {e}")
        print("\nPlease install NGSolve:")
        print("  pip install ngsolve")
        sys.exit(1)

    # Load mesh
    print(f"Loading: {vol_file}")
    mesh = Mesh()

    try:
        mesh.Load(vol_file)
    except Exception as e:
        print(f"Error: Failed to load mesh: {e}")
        sys.exit(1)

    # Start GUI and display
    print(f"Opening Netgen GUI...")
    print(f"  Vertices: {mesh.nv}")
    print(f"  Volume elements: {mesh.ne}")
    print(f"  Surface elements: {mesh.nse}")

    StartGUI()
    mesh.Draw()

    print("\nNetgen GUI opened.")
    print("To close: Close Netgen GUI window")


if __name__ == '__main__':
    main()
