#!/usr/bin/env python
"""
Generate Quarter Model Mesh: C-Type Yoke (1/4 Model)

Based on ELF_MAGIC 6x6x6 Trelis.jou for uniform meshing.
Creates 1/4 model with controlled curve intervals.

Reference: S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\Cubit\6x6x6\Trelis.jou

Output:
- yoke_quarter.vol : Netgen hex mesh (1/4 model)
- yoke_quarter.vtk : VTK mesh for visualization

Run: python generate_quarter_mesh.py
"""

import sys
import os

# ============================================================
# Paths (auto-detect Cubit via install_panels.find_cubit_bin)
# ============================================================
work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(work_dir))
sys.path.insert(0, os.path.join(repo_root, 'src', 'radia'))
from install_panels import find_cubit_bin
_cubit_path = find_cubit_bin()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.insert(0, _cubit_path)

from ngsolve import Mesh, VTKOutput, Draw
import netgen.gui
import cubit
import cubit_mesh_export

print("=" * 60)
print("QUARTER MODEL MESH: C-TYPE YOKE (6x6x6 Style)")
print("=" * 60)
print("Reference: ELF_MAGIC 6x6x6 Trelis.jou")
print("Symmetry planes: X=0 and Z=0")
print()

# ============================================================
# Step 1: Create quarter geometry (same as full model)
# ============================================================
print("Step 1: Create quarter geometry")

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")

# --- Create base geometry (Trelis.jou lines 2-6) ---
cubit.cmd("brick x 62.5 y 105 z 25")
cubit.cmd("brick x 80 y 50 z 25")
cubit.cmd("move Volume 2 location 71.25 -27.5 0")
cubit.cmd("brick x 40 y 100 z 25")
cubit.cmd("move Volume 3 location 131.25 -2.5 0")

# --- Chamfer on pole piece (line 8) ---
cubit.cmd("modify curve 35 26 36 chamfer radius 8")

# --- Webcuts for divisions (lines 9-11) ---
cubit.cmd("webcut volume 3 with plane yplane offset 39.5")
cubit.cmd("webcut volume 3 with plane yplane offset 27.5")
cubit.cmd("webcut volume 1 3 with plane yplane offset -2.5")

# --- Imprint and merge quarter geometry (lines 13-17) ---
cubit.cmd("imprint volume 7 3 2 1 6")
cubit.cmd("merge volume 7 3 2 1 6")
cubit.cmd("imprint volume 5 4")
cubit.cmd("merge volume 5 4")

print(f"  Total volumes: {cubit.get_volume_count()}")

# ============================================================
# Step 2: Mesh with controlled curve intervals (6x6x6 style)
# ============================================================
print("\nStep 2: Mesh with controlled curve intervals")

# Primary curves: interval 12 (lines 19-21)
cubit.cmd("curve 62 41 42 60 47 32 interval 12")
cubit.cmd("curve 62 41 42 60 47 32 scheme equal")
cubit.cmd("mesh curve 62 41 42 60 47 32")

# Secondary curves: interval 8 (lines 22-24)
cubit.cmd("curve 61 40 38 46 45 63 interval 8")
cubit.cmd("curve 61 40 38 46 45 63 scheme equal")
cubit.cmd("mesh curve 61 40 38 46 45 63")

# Tertiary curves: interval 4 (lines 25-27)
cubit.cmd("curve 32 37 64 46 45 63 interval 4")
cubit.cmd("curve 32 37 64 46 45 63 scheme equal")
cubit.cmd("mesh curve 32 37 64 46 45 63")

# Mesh pole piece volumes (line 28-29)
cubit.cmd("volume 5 4 size auto factor 10")
cubit.cmd("mesh volume 5 4")

# Remaining curves: interval 3 (lines 31-33)
cubit.cmd("curve 55 85 34 11 69 9 12 23 96 interval 3")
cubit.cmd("curve 55 85 34 11 69 9 12 23 96 scheme equal")
cubit.cmd("mesh curve 55 85 34 11 69 9 12 23 96")

# More curves: interval 6 (lines 34-36)
cubit.cmd("curve 55 85 34 88 90 13 15 72 11 69 9 74 72 82 83 12 23 80 81 8 2 70 4 16 14 20 interval 6")
cubit.cmd("curve 55 85 34 88 90 13 15 72 11 69 9 74 72 82 83 12 23 80 81 8 2 70 4 16 14 20 scheme equal")
cubit.cmd("mesh curve 55 85 34 88 90 13 15 72 11 69 9 74 72 82 83 12 23 80 81 8 2 70 4 16 14 20")

# Mesh remaining volumes (lines 37-38)
cubit.cmd("volume 1 2 3 6 7 size auto factor 7")
cubit.cmd("mesh volume 1 2 3 6 7")

# Final imprint/merge (lines 40-41)
cubit.cmd("imprint volume all")
cubit.cmd("merge volume all")

print(f"  Hex elements: {cubit.get_hex_count()}")
print(f"  Nodes: {cubit.get_node_count()}")

# --- Define blocks ---
cubit.cmd("block 1 add hex all")
cubit.cmd('block 1 name "yoke"')
cubit.cmd("block 2 add face all")
cubit.cmd('block 2 name "boundary"')

# ============================================================
# Step 3: Transform to final orientation (NO reflections for 1/4)
# ============================================================
# Note: Full model uses two reflections (lines 43-47 in Trelis.jou)
# For 1/4 model, we skip reflections and only apply transformations
print("\nStep 3: Transform to final orientation (1/4 model, no reflections)")

cubit.cmd("move Volume all x -131.25 y -52.5 z 12.5 include_merged")
cubit.cmd("rotate Volume all angle 90 about X include_merged")
cubit.cmd("rotate Volume all angle -90 about Z include_merged")
cubit.cmd("renumber node all in Curve all start_id 1 uniqueids")

print(f"  Transformation applied (moved and rotated)")

# ============================================================
# Step 4: Export hex mesh to Netgen
# ============================================================
print("\nStep 4: Export hex mesh to Netgen")

mesh = cubit_mesh_export.export_NGSolveCurvedMesh(cubit, order=1)

# Scale to meters via ngmesh
ngmesh = mesh.ngmesh
ngmesh.Scale(0.001)
mesh = Mesh(ngmesh)

print(f"  Hex elements: {mesh.ne}")
print(f"  Vertices: {mesh.nv}")
print(f"  Materials: {mesh.GetMaterials()}")

# Save mesh
vol_file = os.path.join(work_dir, 'yoke_quarter.vol')
ngmesh.Save(vol_file)
print(f"  Saved: {vol_file}")

# Save VTK
vtk_file = os.path.join(work_dir, 'yoke_quarter')
vtk = VTKOutput(mesh, coefs=[], names=[], filename=vtk_file)
vtk.Do()
print(f"  Saved: {vtk_file}.vtk")

# ============================================================
# Step 5: Bounding box
# ============================================================
print("\nStep 5: Bounding box (mm)")

import numpy as np
vertices = []
for v in mesh.vertices:
    pt = v.point
    vertices.append([pt[0], pt[1], pt[2]])

vertices = np.array(vertices)

x_min, x_max = vertices[:, 0].min(), vertices[:, 0].max()
y_min, y_max = vertices[:, 1].min(), vertices[:, 1].max()
z_min, z_max = vertices[:, 2].min(), vertices[:, 2].max()

print(f"  X: [{x_min*1000:.1f}, {x_max*1000:.1f}] (size: {(x_max-x_min)*1000:.1f})")
print(f"  Y: [{y_min*1000:.1f}, {y_max*1000:.1f}] (size: {(y_max-y_min)*1000:.1f})")
print(f"  Z: [{z_min*1000:.1f}, {z_max*1000:.1f}] (size: {(z_max-z_min)*1000:.1f})")

# ============================================================
# Step 6: Visualize in Netgen GUI
# ============================================================
print("\nStep 6: Visualize in Netgen GUI")

Draw(mesh)
print("  GUI opened - showing 1/4 mesh (R3856 style)")

# ============================================================
# Summary
# ============================================================
print()
print("=" * 60)
print("QUARTER MODEL MESH COMPLETE (6x6x6 Style)")
print("=" * 60)
print()
print(f"Quarter model: {mesh.ne} hex elements")
print(f"Full model would have: ~{mesh.ne * 4} hex elements")
print()
print("Mesh intervals (6x6x6 style):")
print("  - Primary curves: 12 intervals")
print("  - Secondary curves: 8 intervals")
print("  - Tertiary curves: 4 intervals")
print("  - Auto factor: 10 (pole), 7 (body)")
print()
print("Symmetry planes:")
print("  X=0: TrfPlSym([0,0,0], [1,0,0])")
print("  Z=0: TrfPlSym([0,0,0], [0,0,1])")
print()
print("Output files:")
print(f"  - yoke_quarter.vol : Netgen hex mesh")
print(f"  - yoke_quarter.vtk : VTK mesh")
print("=" * 60)
