#!/usr/bin/env python
"""
Generate Netgen mesh from Cubit journal file

This script reads York.jou (Cubit journal file) and generates Netgen mesh directly:
- Direct export to Netgen (in-memory, no intermediate files)
- York.vtk (VTK format) - for ParaView visualization

Workflow: Cubit geometry -> export_netgen() -> Netgen mesh -> Radia

Requirements:
- Coreform Cubit 2025.3 (or compatible version)
- cubit_mesh_export module from S:\CoreformCubit\01_GitHub
- NGSolve (for Mesh wrapper)

Usage:
    python york_cubit_mesh.py

Input:
    York.jou - Cubit journal file with yoke geometry definition

Output:
    Netgen mesh object (returned for use in main_simulation_workflow.py)
    York.vtk - VTK legacy format for ParaView
"""

import os
import sys

# Add Cubit Python API to path
CUBIT_PATH = "C:/Program Files/Coreform Cubit 2025.3/bin"
CUBIT_EXPORT_PATH = "S:/CoreformCubit/01_GitHub"

if os.path.exists(CUBIT_PATH):
    sys.path.insert(0, CUBIT_PATH)
else:
    print(f"[WARNING] Cubit not found at: {CUBIT_PATH}")
    print("          Please adjust CUBIT_PATH in this script")
    sys.exit(1)

if os.path.exists(CUBIT_EXPORT_PATH):
    sys.path.insert(0, CUBIT_EXPORT_PATH)
else:
    print(f"[WARNING] cubit_mesh_export not found at: {CUBIT_EXPORT_PATH}")
    sys.exit(1)

try:
    import cubit
except ImportError:
    print("[ERROR] Failed to import cubit module")
    print("        Make sure Coreform Cubit is installed")
    sys.exit(1)

try:
    import cubit_mesh_export
except ImportError:
    print("[ERROR] Failed to import cubit_mesh_export module")
    print("        This module should be at S:\\CoreformCubit\\01_GitHub")
    sys.exit(1)


def generate_yoke_mesh():
    """
    Generate yoke mesh from Cubit journal file.

    Returns:
        Netgen mesh object (for use with NGSolve.Mesh())
    """
    # Initialize Cubit in batch mode (no GUI)
    print("Initializing Cubit...")
    cubit.init(['cubit', '-nojournal', '-batch'])

    # Read and execute Cubit journal file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    journal_file = os.path.join(script_dir, 'York.jou')

    if not os.path.exists(journal_file):
        print(f"[ERROR] Journal file not found: {journal_file}")
        return None

    print(f"Reading journal file: {journal_file}")
    with open(journal_file, 'r', encoding='utf8') as fid:
        lines = fid.readlines()
        for line in lines:
            cubit.cmd(line)

    print(f"  Executed {len(lines)} commands")

    # Export VTK for visualization
    vtk_file = os.path.join(script_dir, 'York.vtk')
    print(f"\nExporting VTK: {vtk_file}")
    cubit_mesh_export.export_vtk(cubit, vtk_file, ORDER="1st")

    # Export directly to Netgen (no intermediate file)
    print("Exporting to Netgen mesh (direct, in-memory)...")
    ngmesh = cubit_mesh_export.export_netgen(cubit)

    print("\n[SUCCESS] Mesh generation complete!")
    print(f"          Visualize with: paraview York.vtk")

    return ngmesh


if __name__ == '__main__':
    ngmesh = generate_yoke_mesh()
    if ngmesh:
        print(f"\nNetgen mesh created: {type(ngmesh)}")
        print("Use this mesh in main_simulation_workflow.py")
