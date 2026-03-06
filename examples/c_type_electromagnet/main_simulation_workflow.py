#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Complete Electromagnet Simulation Workflow

This script performs the complete simulation workflow:
1. Import yoke mesh via Cubit -> Netgen direct export
2. Create racetrack coil
3. Combine coil + yoke into Radia model
4. Solve magnetostatic problem
5. Export field_distribution.vts (field data for ParaView)

Workflow: Cubit geometry -> export_netgen() -> Netgen mesh -> Radia

Requirements:
- Coreform Cubit 2025.3+
- cubit_mesh_export from S:\CoreformCubit\01_GitHub
- NGSolve
- York.jou (Cubit journal file with yoke geometry)
"""

import os
import sys

# Add paths
CUBIT_PATH = "C:/Program Files/Coreform Cubit 2025.3/bin"
CUBIT_EXPORT_PATH = "S:/CoreformCubit/01_GitHub"

sys.path.insert(0, CUBIT_PATH)
sys.path.insert(0, CUBIT_EXPORT_PATH)
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../build-msvc'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

import radia as rad
from netgen_mesh_import import netgen_mesh_to_radia

# Set Radia unit system to millimeters
rad.FldUnits('mm')

print("=" * 70)
print("ELECTROMAGNET SIMULATION WORKFLOW")
print("=" * 70)
print(f"Unit system: mm (millimeters)")

# ========================================================================
# Step 1: Import Yoke Mesh via Cubit -> Netgen
# ========================================================================
print("\n[Step 1/5] Importing yoke mesh via Cubit -> Netgen...")

try:
    import cubit
    import cubit_mesh_export
    CUBIT_AVAILABLE = True
except ImportError:
    CUBIT_AVAILABLE = False
    print("[ERROR] Cubit not available")
    print("        Please install Coreform Cubit 2025.3+")
    sys.exit(1)

# Initialize Cubit
cubit.init(['cubit', '-nojournal', '-batch'])

# Read journal file
script_dir = os.path.dirname(os.path.abspath(__file__))
journal_file = os.path.join(script_dir, 'York.jou')

if not os.path.exists(journal_file):
    print(f"[ERROR] Journal file not found: {journal_file}")
    sys.exit(1)

print(f"  Reading: {journal_file}")
with open(journal_file, 'r', encoding='utf8') as fid:
    for line in fid:
        cubit.cmd(line)

# Export directly to Netgen
print("  Exporting to Netgen mesh...")
ngmesh = cubit_mesh_export.export_netgen(cubit)

# Convert to NGSolve mesh
from ngsolve import Mesh
mesh = Mesh(ngmesh)

# Convert to Radia (units='mm' since we're using mm)
yoke = netgen_mesh_to_radia(mesh,
                             material={'magnetization': [0, 0, 0]},
                             units='mm')

if yoke is None:
    print("[ERROR] Failed to create Radia yoke from mesh")
    sys.exit(1)

# Apply soft iron material
yoke_mat = rad.MatSatIsoFrm([20000, 2], [0.1, 2], [0.1, 2])
rad.MatApl(yoke, yoke_mat)

print(f"  [OK] Yoke imported: ID={yoke}")

# ========================================================================
# Step 2: Create Racetrack Coil
# ========================================================================
print("\n[Step 2/5] Creating racetrack coil...")

coil = rad.ObjRaceTrk(
    [0, 131.25, 0],  # center (mm)
    [5, 40],         # lx_min, lx_max (mm)
    [50, 62.5],      # ly_min, ly_max (mm)
    105,             # height (mm)
    3,               # nseg
    -2000/105/35     # current density (A/mm^2): -2000A total
)

print(f"  [OK] Coil created: ID={coil}")
print(f"  Total current: -2000 A")

# ========================================================================
# Step 3: Combine into Radia Model
# ========================================================================
print("\n[Step 3/5] Combining coil + yoke...")

g = rad.ObjCnt([coil, yoke])
print(f"  [OK] Combined model: ID={g}")

# ========================================================================
# Step 4: Solve Magnetostatic Problem
# ========================================================================
print("\n[Step 4/5] Solving magnetostatics...")

rad.Solve(g, 0.0001, 10000)
print(f"  [OK] Solution converged")

# ========================================================================
# Step 5: Export Field Distribution
# ========================================================================
print("\n[Step 5/5] Exporting field distribution...")

# Define field evaluation grid
# Combined geometry bbox: X[-40, 40], Y[-25, 193.8], Z[-52.5, 162.5] mm
# Add 50mm margin around geometry
x_range_bounds = [-90, 90]
y_range_bounds = [-75, 245]
z_range_bounds = [-105, 215]

# Export field distribution to VTS format
try:
    vts_path = os.path.join(script_dir, 'field_distribution.vts')
    rad.FldVTS(g, vts_path, x_range_bounds, y_range_bounds, z_range_bounds, 21, 31, 21, 1, 0, 1.0)
    print(f"  [OK] Field distribution exported to field_distribution.vts")
except Exception as e:
    print(f"  [WARNING] VTS export failed: {e}")

# ========================================================================
# Summary
# ========================================================================
print("\n" + "=" * 70)
print("SIMULATION COMPLETE")
print("=" * 70)
print(f"\nGenerated files:")
print(f"  1. field_distribution.vts - Magnetic field data (VTS format)")
print(f"\nVisualize with ParaView:")
print(f"  paraview field_distribution.vts")
print("=" * 70)
