#!/usr/bin/env python
"""
Electromagnet Simulation: C-Type Yoke with Racetrack Coil

Loads mesh from generate_mesh.py and runs Radia magnetostatic simulation.

Workflow:
1. Load Netgen mesh (yoke.vol)
2. Convert to Radia geometry
3. Create racetrack coil
4. Solve magnetostatics
5. Export field distribution (VTS)

Input:
- yoke.vol : Netgen mesh file (from generate_mesh.py)

Output:
- field_distribution.vts : Magnetic field data (for ParaView)

Run: python run_simulation.py

Requirements:
- NGSolve / Netgen
- Radia
- yoke.vol (run generate_mesh.py first)

Reference: S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=10000\ELF_MMB8T_EIEM0_R288
"""

import sys
import os

# ============================================================
# Paths
# ============================================================
work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(work_dir))
# Import from src/radia package (contains latest .pyd from BuildMSVC.ps1)
sys.path.insert(0, os.path.join(repo_root, 'src/radia'))

from netgen.meshing import Mesh as NetgenMesh
from ngsolve import Mesh
import radia as rad
from netgen_mesh_import import netgen_mesh_to_radia

# ============================================================
# Parameters (SI units: meters)
# ============================================================
# Coil parameters from ELF_magic.mei:
#   AA O 0 131.25 0
#   AA SQRING 60 72.5 105 35 5 3
#
# Original position: (0, 131.25, 0) mm in Cubit/ELF coordinates (before transformation)
#
# Cubit coordinate transformations (from Trelis.jou):
#   1. move Volume all x -131.25 y -52.5 z 12.5
#   2. rotate Volume all angle 90 about X
#   3. rotate Volume all angle -90 about Z
#
# Coil center transformation:
#   (0, 131.25, 0) -> translate -> (-131.25, 78.75, 12.5)
#                  -> rot 90 X   -> (-131.25, -12.5, 78.75)
#                  -> rot -90 Z  -> (-12.5, 131.25, 78.75)
#
# Final coil center: (-12.5, 131.25, 78.75) mm = (-0.0125, 0.13125, 0.07875) m
#
# ELF SQRING parameters (from IEmesh.def):
#   lx=60:   窓X方向の大きさ (inner opening X) [mm]
#   ly=72.5: 窓Y方向の大きさ (inner opening Y) [mm]
#   lz=105:  導線Z方向の大きさ (coil height in Z) [mm]
#   lt=35:   導線の大きさ (conductor thickness) [mm]
#   r=5:     角部の半径 (corner radius) [mm]
#   Mr=3:    角部の分割数 (corner divisions)

# Coil position after coordinate transformation
COIL_CENTER = [-0.0125, 0.13125, 0.07875]  # m (transformed from ELF_magic.mei)

# Racetrack coil dimensions (from ELF_magic.mei SQRING 60 72.5 105 35 5 3)
# Corner radii: r_inner = 5mm, r_outer = r + lt = 40mm
COIL_RADII = [0.005, 0.040]            # [inner, outer] corner radii (m)
# Straight section half-lengths: lx/2 = 30mm, ly/2 = 36.25mm
COIL_STRAIGHT = [0.030, 0.03625]       # [lx/2, ly/2] half-lengths (m)
# Coil height: lz = 105mm
COIL_HEIGHT = 0.105                     # m (105mm coil height)
# Note: Radia uses 'auto' mode for arc precision (no manual segmentation needed)
COIL_CURRENT = -2000                    # A

# Field grid (m) - covers the C-type yoke volume
# After transformation: yoke is ~262.5mm wide, ~105mm high, ~50mm deep
GRID_X = [-0.08, 0.08]
GRID_Y = [-0.15, 0.15]
GRID_Z = [-0.08, 0.20]
GRID_NX, GRID_NY, GRID_NZ = 17, 31, 29

print("=" * 60)
print("ELECTROMAGNET SIMULATION: C-TYPE YOKE")
print("=" * 60)
print("Unit system: SI (meters)")
print()

# Set Radia unit system
rad.FldUnits('m')

# ============================================================
# Step 1: Load Netgen mesh
# ============================================================
print("Step 1: Load Netgen mesh")

vol_file = os.path.join(work_dir, 'yoke.vol')

if not os.path.exists(vol_file):
    print(f"  [ERROR] Mesh file not found: {vol_file}")
    print("  Please run generate_mesh.py first")
    sys.exit(1)

ngmesh = NetgenMesh()
ngmesh.Load(vol_file)
mesh = Mesh(ngmesh)

print(f"  Loaded: {vol_file}")
print(f"  Elements: {mesh.ne}")
print(f"  Vertices: {mesh.nv}")

# ============================================================
# Step 2: Convert to Radia geometry
# ============================================================
print("\nStep 2: Convert mesh to Radia (MSC hexahedra)")

yoke = netgen_mesh_to_radia(mesh,
                             material={'magnetization': [0, 0, 0]},
                             units='m',
                             material_filter='yoke',
                             allow_hex=True)

if yoke is None:
    print("  [ERROR] Failed to create Radia geometry")
    sys.exit(1)

# Apply soft iron material (linear, mu_r = 1000 for testing)
# For production: use rad.MatSatIsoFrm or rad.MatSatIsoTab for nonlinear B-H curve
yoke_mat = rad.MatLin(1000)  # mu_r = 1000 (linear iron)
rad.MatApl(yoke, yoke_mat)
print(f"  Yoke ID: {yoke}")
print(f"  Material: Linear iron (mu_r = 1000)")

# ============================================================
# Step 3: Create racetrack coil
# ============================================================
print("\nStep 3: Create racetrack coil")

coil_width = COIL_RADII[1] - COIL_RADII[0]
coil_cs = COIL_HEIGHT * coil_width
j_density = COIL_CURRENT / coil_cs

coil = rad.ObjRaceTrk(
    COIL_CENTER,
    COIL_RADII,
    COIL_STRAIGHT,
    COIL_HEIGHT,
    1,           # nseg (ignored in 'auto' mode)
    j_density,
    'auto'       # Use arc precision from FldCmpCrt
)

print(f"  Coil ID: {coil}")
print(f"  Current: {COIL_CURRENT} A")
print(f"  Current density: {j_density:.1f} A/m^2")

# ============================================================
# Step 4: Solve magnetostatics
# ============================================================
print("\nStep 4: Solve magnetostatics")

g = rad.ObjCnt([coil, yoke])
print(f"  Combined model ID: {g}")

import time
t_start = time.time()
# Method 1 = BiCGSTAB (faster for large problems)
# Method 0 = LU (accurate but slow for large problems)
# Method 2 = HACApK (H-matrix accelerated, best for very large problems)
result = rad.Solve(g, 0.001, 1000, 1)  # BiCGSTAB solver
t_solve = time.time() - t_start
print(f"  [OK] Solution converged in {t_solve:.1f}s")

# ============================================================
# Step 5: Export field distribution
# ============================================================
print("\nStep 5: Export field distribution (VTS)")

vts_file = os.path.join(work_dir, 'field_distribution.vts')
rad.FldVTS(g, vts_file, GRID_X, GRID_Y, GRID_Z,
           GRID_NX, GRID_NY, GRID_NZ, 1, 0, 1.0)
print(f"  Exported: {vts_file}")

# ============================================================
# Summary
# ============================================================
print()
print("=" * 60)
print("SIMULATION COMPLETE")
print("=" * 60)
print()
print("Output files:")
print(f"  - field_distribution.vts : Magnetic field data")
print()
print("Visualize:")
print(f"  paraview {os.path.join(work_dir, 'yoke_mesh.vtk')} {vts_file}")
print("=" * 60)
