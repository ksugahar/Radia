#!/usr/bin/env python
"""
Verify Quarter Model Mesh

Loads yoke_quarter.vol and displays in Netgen GUI.
"""

import sys
import os

work_dir = os.path.dirname(os.path.abspath(__file__))

# Use custom NGSolve with GUI support
NGSOLVE_PATH = r'S:\NGSolve\01_GitHub\install_ksugahar\lib\site-packages'
sys.path.insert(0, NGSOLVE_PATH)

from netgen.meshing import Mesh as NetgenMesh
from ngsolve import Mesh, Draw
import netgen.gui

print("=" * 60)
print("VERIFY QUARTER MODEL MESH")
print("=" * 60)

# Load mesh
vol_file = os.path.join(work_dir, 'yoke_quarter.vol')
if not os.path.exists(vol_file):
    print(f"[ERROR] {vol_file} not found!")
    print("Run: python generate_quarter_mesh.py first")
    sys.exit(1)

print(f"Loading: {vol_file}")

ngmesh = NetgenMesh()
ngmesh.Load(vol_file)
mesh = Mesh(ngmesh)

print(f"Elements: {mesh.ne}")
print(f"Vertices: {mesh.nv}")
print(f"Materials: {mesh.GetMaterials()}")

# Draw mesh
Draw(mesh)
print()
print("GUI opened - showing 1/4 mesh")
print("Close window to exit")
