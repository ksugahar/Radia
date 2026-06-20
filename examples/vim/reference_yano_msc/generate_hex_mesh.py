#!/usr/bin/env python
"""
Mesh Generation: C-Type Yoke

Creates C-type yoke geometry and generates meshes for both Radia and NGSolve comparison.
Direct translation of ELF_MAGIC Cubit journal file (Trelis.jou).

Workflow:
1. Cubit: Create C-type yoke geometry (from Trelis.jou)
   - Quarter model: main leg + yoke back + pole piece with chamfer
   - Hex mesh with controlled intervals
   - Reflect twice to create full symmetric model
   - Transform to final orientation
2. Cubit: Save database (.cub)
3. Cubit: Export STEP (geometry for NGSolve tet mesh)
4. cubit_mesh_export: Export hex mesh to Netgen (for Radia)

Output:
- yoke.cub  : Cubit database (for editing)
- yoke.vol  : Netgen hex mesh (for Radia via run_simulation.py)
- yoke.step : STEP geometry (for NGSolve tet mesh comparison)

Run: python generate_mesh.py

Requirements:
- Coreform Cubit 2025.8+
- cubit_mesh_export module (via radia package)
- NGSolve / Netgen

Reference: S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\Cubit\6x6x6\Trelis.jou
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
# cubit_mesh_export is imported from src/radia

from ngsolve import Mesh, VTKOutput
import cubit
from cubit_mesh_export import extract_curved_mesh

# ============================================================
# Parameters (mm for Cubit)
# ============================================================
# Mesh resolution: 'coarse' (~80 hex), 'medium' (~600 hex), 'fine' (~4800 hex)
# Reference: ELF_MAGIC Cubit models 1x1x1, 3x3x3, 6x6x6
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--resolution', '-r', choices=['coarse', 'medium', 'fine'],
                    default='coarse', help='Mesh resolution (default: coarse)')
args = parser.parse_args()
MESH_RESOLUTION = args.resolution

# Mesh interval settings based on resolution
if MESH_RESOLUTION == 'coarse':
    # ~80 hex elements (1x1x1 equivalent)
    INTERVAL_MAIN = 2
    INTERVAL_SEC = 1
    INTERVAL_THIRD = 1
    SIZE_FACTOR = 10
elif MESH_RESOLUTION == 'medium':
    # ~600 hex elements (3x3x3 equivalent)
    INTERVAL_MAIN = 6
    INTERVAL_SEC = 4
    INTERVAL_THIRD = 2
    SIZE_FACTOR = 8
else:  # fine
    # ~4800 hex elements (6x6x6 equivalent)
    INTERVAL_MAIN = 12
    INTERVAL_SEC = 8
    INTERVAL_THIRD = 4
    SIZE_FACTOR = 7

# C-type yoke geometry (from ELF_MAGIC Cubit model)
#
#           ┌─────┐
#           │pole │
#     ┌─────┴─────┴─────┐
#     │   yoke back     │
#     └──┬──────────┬───┘
#        │          │
#        │   leg    │  (magnetic gap in center)
#        │          │
#     ┌──┴──────────┴───┐
#     │   yoke back     │
#     └─────┬─────┬─────┘
#           │pole │
#           └─────┘
#
# Quarter model dimensions (mm), then reflected twice:
# - Main leg: 62.5 x 105 x 25 at origin
# - Yoke back: 80 x 50 x 25 at (71.25, -27.5, 0)
# - Pole piece: 40 x 100 x 25 at (131.25, -2.5, 0)

print("=" * 60)
print("MESH GENERATION: C-TYPE YOKE")
print("=" * 60)
print(f"Resolution: {MESH_RESOLUTION}")
print("Reference: ELF_MAGIC model_C-Type Cubit")
print()

# ============================================================
# Step 1: Create C-type yoke geometry and hex mesh in Cubit
# ============================================================
# Direct translation from:
# S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\Cubit\6x6x6\Trelis.jou
print("Step 1: Create C-type yoke geometry and hex mesh in Cubit")

cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("reset")

# --- Trelis.jou lines 2-6: Create quarter geometry ---
cubit.cmd("brick x 62.5 y 105 z 25")
cubit.cmd("brick x 80 y 50 z 25")
cubit.cmd("move Volume 2 location 71.25 -27.5 0")
cubit.cmd("brick x 40 y 100 z 25")
cubit.cmd("move Volume 3 location 131.25 -2.5 0")

# --- Trelis.jou line 8: Chamfer on pole piece ---
cubit.cmd("modify curve 35 26 36 chamfer radius 8")

# --- Trelis.jou lines 9-11: Webcuts ---
cubit.cmd("webcut volume 3 with plane yplane offset 39.5")
cubit.cmd("webcut volume 3 with plane yplane offset 27.5")
cubit.cmd("webcut volume 1 3 with plane yplane offset -2.5")

# --- Trelis.jou lines 13-17: Imprint/merge quarter ---
cubit.cmd("imprint volume 7 3 2 1 6")
cubit.cmd("merge volume 7 3 2 1 6")
cubit.cmd("imprint volume 5 4")
cubit.cmd("merge volume 5 4")

print(f"  Created quarter geometry")
print(f"  Volumes: {cubit.get_volume_count()}")

# --- Hex meshing with parameterized intervals ---
print(f"  Meshing quarter geometry (hex, resolution={MESH_RESOLUTION})...")
cubit.cmd(f"curve 62 41 42 60 47 32 interval {INTERVAL_MAIN}")
cubit.cmd("curve 62 41 42 60 47 32 scheme equal")
cubit.cmd("mesh curve 62 41 42 60 47 32")
cubit.cmd(f"curve 61 40 38 46 45 63 interval {INTERVAL_SEC}")
cubit.cmd("curve 61 40 38 46 45 63 scheme equal")
cubit.cmd("mesh curve 61 40 38 46 45 63")
cubit.cmd(f"curve 32 37 64 46 45 63 interval {INTERVAL_THIRD}")
cubit.cmd("curve 32 37 64 46 45 63 scheme equal")
cubit.cmd("mesh curve 32 37 64 46 45 63")
cubit.cmd(f"volume 5 4 size auto factor {SIZE_FACTOR}")
cubit.cmd("mesh volume 5 4")

cubit.cmd(f"curve 55 85 34 11 69 9 12 23 96 interval {max(1, INTERVAL_THIRD//2)}")
cubit.cmd("curve 55 85 34 11 69 9 12 23 96 scheme equal")
cubit.cmd("mesh curve 55 85 34 11 69 9 12 23 96")
cubit.cmd(f"curve 55 85 34 88 90 13 15 72 11 69 9 74 72 82 83 12 23 80 81 8 2 70 4 16 14 20 interval {INTERVAL_SEC//2}")
cubit.cmd("curve 55 85 34 88 90 13 15 72 11 69 9 74 72 82 83 12 23 80 81 8 2 70 4 16 14 20 scheme equal")
cubit.cmd("mesh curve 55 85 34 88 90 13 15 72 11 69 9 74 72 82 83 12 23 80 81 8 2 70 4 16 14 20")
cubit.cmd(f"volume 1 2 3 6 7 size auto factor {SIZE_FACTOR}")
cubit.cmd("mesh volume 1 2 3 6 7")

# --- Trelis.jou lines 40-41: Imprint/merge ---
cubit.cmd("imprint volume all")
cubit.cmd("merge volume all")

# --- Trelis.jou lines 43-47: Reflections (copies mesh) ---
cubit.cmd("Volume all copy reflect vertex 55 53")
cubit.cmd("Volume all copy reflect vertex 3 47")
cubit.cmd("imprint volume all")
cubit.cmd("merge volume all")

print(f"  After reflections: {cubit.get_volume_count()} volume(s)")

# --- Trelis.jou lines 49-51: Transformations ---
cubit.cmd("move Volume all x -131.25 y -52.5 z 12.5 include_merged")
cubit.cmd("rotate Volume all angle 90 about X include_merged")
cubit.cmd("rotate Volume all angle -90 about Z include_merged")

# --- Trelis.jou line 53: Renumber nodes ---
cubit.cmd("renumber node all in Curve all start_id 1 uniqueids")

print(f"  Hex elements: {cubit.get_hex_count()}")
print(f"  Nodes: {cubit.get_node_count()}")

# --- Define blocks for Netgen export ---
cubit.cmd("block 1 add hex all")
cubit.cmd('block 1 name "yoke"')
cubit.cmd("block 2 add face all")
cubit.cmd('block 2 name "boundary"')

print(f"  Block 1 (yoke): {cubit.get_hex_count()} hexes")
print(f"  Block 2 (boundary): {len(cubit.get_block_faces(2))} faces")

# ============================================================
# Step 2: Save Cubit database (.cub) - in mm
# ============================================================
print("\nStep 2: Save Cubit database (mm)")

cub_file = os.path.join(work_dir, "yoke.cub")
cubit.cmd(f'save as "{cub_file}" overwrite')
print(f"  Saved: {cub_file}")

# ============================================================
# Step 3: Export STEP (BEFORE unite - preserves all volumes)
# ============================================================
print("\nStep 3: Export STEP (for OCCGeometry)")

step_file = os.path.join(work_dir, "yoke.step")
cubit.cmd(f'export step "{step_file}" overwrite')
print(f"  Exported: {step_file}")
print(f"  Volumes: {cubit.get_volume_count()}")

# ============================================================
# Step 4: Export hex mesh to Netgen (for Radia)
# ============================================================
print("\nStep 4: Export hex mesh to Netgen (for Radia)")

with TaskManager():
    mesh = Mesh(extract_curved_mesh(cubit, order=1))
    ngmesh = mesh.ngmesh

    print(f"  Hex elements: {mesh.ne}")
    print(f"  Vertices: {mesh.nv}")
    print(f"  Materials: {mesh.GetMaterials()}")

    # Scale mesh to meters
    ngmesh.Scale(0.001)

    # Save Netgen mesh (.vol)
    vol_file = os.path.join(work_dir, 'yoke.vol')
    ngmesh.Save(vol_file)
    print(f"  Saved: {vol_file}")

    # Save VTK mesh for visualization
    vtk_file = os.path.join(work_dir, 'yoke_mesh')
    mesh_scaled = Mesh(ngmesh)  # Reload mesh after scaling
    vtk = VTKOutput(mesh_scaled, coefs=[], names=[], filename=vtk_file)
    vtk.Do()
    print(f"  Saved: {vtk_file}.vtk")

    # ============================================================
    # Step 5: Verify hex mesh and bounding box
    # ============================================================
    print("\nStep 5: Verify hex mesh and bounding box")

    import numpy as np

    # Element type verification
    from ngsolve import VOL
    from ngsolve import TaskManager
    element_types = {}
    for el in mesh_scaled.Elements(VOL):
        et = str(el.type)
        element_types[et] = element_types.get(et, 0) + 1

    for et, count in element_types.items():
        print(f"  {et}: {count}")

    if 'HEX' in str(element_types) or 'HEXAHEDRON' in str(element_types):
        print("  [OK] Hexahedral mesh confirmed")
    else:
        print("  [WARNING] No hexahedral elements found!")

    # Bounding box (meters)
    vertices = []
    for v in mesh_scaled.vertices:
        pt = v.point
        vertices.append([pt[0], pt[1], pt[2]])

    vertices = np.array(vertices)

    x_min, x_max = vertices[:, 0].min(), vertices[:, 0].max()
    y_min, y_max = vertices[:, 1].min(), vertices[:, 1].max()
    z_min, z_max = vertices[:, 2].min(), vertices[:, 2].max()

    print(f"\n  Bounding box (mm):")
    print(f"    X: [{x_min*1000:.1f}, {x_max*1000:.1f}] (size: {(x_max-x_min)*1000:.1f})")
    print(f"    Y: [{y_min*1000:.1f}, {y_max*1000:.1f}] (size: {(y_max-y_min)*1000:.1f})")
    print(f"    Z: [{z_min*1000:.1f}, {z_max*1000:.1f}] (size: {(z_max-z_min)*1000:.1f})")

    # Update mesh variable for interactive mode
    mesh = mesh_scaled

    # ============================================================
    # Summary
    # ============================================================
    print()
    print("=" * 60)
    print("MESH GENERATION COMPLETE")
    print("=" * 60)
    print()
    print("Output files:")
    print(f"  - yoke.cub      : Cubit database (for editing)")
    print(f"  - yoke.vol      : Netgen hex mesh (for Radia)")
    print(f"  - yoke_mesh.vtk : VTK mesh (for ParaView)")
    print(f"  - yoke.step     : STEP geometry (for NGSolve tet mesh)")
    print()
    print("Next steps:")
    print("  Radia:   python solve_radia.py")
    print("  NGSolve: python solve_ngsolve.py")
    print("=" * 60)

    # ============================================================
    # Interactive mode: Variables available for inspection
    # ============================================================
    # When running with `py -i generate_mesh.py`:
    #   mesh       : NGSolve Mesh object (scaled to meters)
    #   ngmesh     : Netgen mesh object
    #   cubit      : Cubit module (for further inspection)
    #
    # Example commands in interactive mode:
    #   mesh.ne                    # Number of elements
    #   mesh.nv                    # Number of vertices
    #   mesh.GetMaterials()        # Material names
    #   list(mesh.Elements3D())    # List of 3D elements
    #   for el in mesh.Elements3D(): print(el.type)  # Element types

    print()
    print("Interactive mode (py -i):")
    print("  mesh.ne              # Number of hex elements")
    print("  mesh.nv              # Number of vertices")
    print("  mesh.GetMaterials()  # Material names ['yoke']")
    print("  cubit.get_hex_count()  # Cubit hex count")
