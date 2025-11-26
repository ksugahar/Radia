#!/usr/bin/env python
"""
Generate Nastran mesh from Cubit journal file

This script reads York.jou (Cubit journal file) and generates mesh files:
- York.bdf (Nastran format) - for Radia simulation
- York.msh (Gmsh format) - optional
- York.vtk (VTK format) - for ParaView visualization

Requirements:
- Coreform Cubit 2025.3 (or compatible version)
- cubit_mesh_export module (from Coreform Cubit installation)

Usage:
    python york_cubit_mesh.py

Input:
    York.jou - Cubit journal file with yoke geometry definition

Output:
    York.bdf - Nastran bulk data file (288 elements, 569 vertices)
    York.msh - Gmsh version 2 format
    York.vtk - VTK legacy format for ParaView
"""

import os
import sys

# Add Cubit Python API to path
# Adjust this path if Cubit is installed in a different location
CUBIT_PATH = "C:/Program Files/Coreform Cubit 2025.3/bin"
if os.path.exists(CUBIT_PATH):
    sys.path.append(CUBIT_PATH)
else:
    print(f"[WARNING] Cubit not found at: {CUBIT_PATH}")
    print("          Please adjust CUBIT_PATH in this script")
    print("          or install Coreform Cubit")
    sys.exit(1)

try:
    import cubit
except ImportError:
    print("[ERROR] Failed to import cubit module")
    print("        Make sure Coreform Cubit is installed")
    print(f"        and CUBIT_PATH is correct: {CUBIT_PATH}")
    sys.exit(1)

try:
    import cubit_mesh_export
except ImportError:
    print("[ERROR] Failed to import cubit_mesh_export module")
    print("        This module should be included with Coreform Cubit")
    sys.exit(1)

# Initialize Cubit in batch mode (no GUI)
print("Initializing Cubit...")
cubit.init(['cubit', '-nojournal', '-batch'])

# Read and execute Cubit journal file
JOURNAL_FILE = 'York.jou'
if not os.path.exists(JOURNAL_FILE):
    print(f"[ERROR] Journal file not found: {JOURNAL_FILE}")
    sys.exit(1)

print(f"Reading journal file: {JOURNAL_FILE}")
with open(JOURNAL_FILE, 'r', encoding='utf8') as fid:
    strLines = fid.readlines()
    for n, line in enumerate(strLines):
        cubit.cmd(line)

print(f"  Executed {len(strLines)} commands")

# Export mesh in multiple formats
FileName = 'York'

print(f"\nExporting mesh files:")
print(f"  {FileName}.msh (Gmsh v2)")
cubit_mesh_export.export_Gmsh_ver2(cubit, FileName + '.msh')

print(f"  {FileName}.bdf (Nastran)")
cubit_mesh_export.export_Nastran(cubit, FileName + '.bdf', DIM='3D', PYRAM=False)

print(f"  {FileName}.vtk (VTK)")
cubit_mesh_export.export_vtk(cubit, FileName + '.vtk', ORDER="1st")

print("\n[SUCCESS] Mesh generation complete!")
print(f"          Use {FileName}.bdf in main_simulation_workflow.py")
print(f"          Visualize with: paraview {FileName}.vtk")

