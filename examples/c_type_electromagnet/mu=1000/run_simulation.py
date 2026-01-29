#!/usr/bin/env python
"""
Electromagnet Simulation: mu_r = 1000 case

Uses full model hex mesh from Cubit (yoke.vol).
Reference: ELF_MAGIC Radia_1x1x1 (Bz = 51.5 mT at origin)
"""

import sys
import os

# Paths
work_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(work_dir)
repo_root = os.path.dirname(os.path.dirname(parent_dir))
sys.path.insert(0, os.path.join(repo_root, 'src'))

from netgen.meshing import Mesh as NetgenMesh
from ngsolve import Mesh
import radia as rad
from radia.netgen_mesh_import import netgen_mesh_to_radia
import numpy as np
import time

# ============================================================
# Parameters
# ============================================================
MU_R = 1000  # Relative permeability

# Coil parameters (from verified Radia_1x1x1/Radia.py)
# radObjRaceTrk[{0, 131.25, 0}, {5, 40}, {50, 62.5}, 105, 3, current/105/35]
COIL_CENTER = [0, 0.13125, 0]  # m = [0, 131.25, 0] mm
COIL_RADII = [0.005, 0.040]   # m = [5, 40] mm (corner radii)
COIL_STRAIGHT = [0.050, 0.0625]  # m = [50, 62.5] mm (half-lengths)
COIL_HEIGHT = 0.105  # m = 105 mm
COIL_CURRENT = -2000  # A (verified value)

# Field grid
GRID_X = [-0.05, 0.05]
GRID_Y = [-0.05, 0.05]
GRID_Z = [0, 0.10]
GRID_NX, GRID_NY, GRID_NZ = 21, 21, 21

print("=" * 60)
print(f"ELECTROMAGNET SIMULATION: mu_r = {MU_R}")
print("=" * 60)

rad.FldUnits('m')

# ============================================================
# Load mesh
# ============================================================
print("\nStep 1: Load mesh")

vol_file = os.path.join(work_dir, 'yoke.vol')

if not os.path.exists(vol_file):
    print(f"  [ERROR] Mesh file not found: {vol_file}")
    print(f"  Run: python generate_mesh.py")
    sys.exit(1)

ngmesh = NetgenMesh()
ngmesh.Load(vol_file)
mesh = Mesh(ngmesh)

print(f"  Loaded: {vol_file}")
print(f"  Elements: {mesh.ne} (full model)")

# ============================================================
# Convert to Radia
# ============================================================
print("\nStep 2: Convert to Radia")

yoke = netgen_mesh_to_radia(mesh,
                             material={'magnetization': [0, 0, 0]},
                             units='m',
                             allow_hex=True)

if yoke is None:
    print("  [ERROR] Failed to create Radia geometry")
    sys.exit(1)

# Full model - no symmetry transformation needed
print(f"  Yoke ID: {yoke} (full model, {mesh.ne} hex elements)")

# Apply material
yoke_mat = rad.MatLin(MU_R)
rad.MatApl(yoke, yoke_mat)
print(f"  Material: mu_r = {MU_R}")

# ============================================================
# Create coil
# ============================================================
print("\nStep 3: Create racetrack coil")

# Coil cross-section: height * width (width = outer_radius - inner_radius)
coil_width = COIL_RADII[1] - COIL_RADII[0]  # 35 mm = 0.035 m
coil_cs = COIL_HEIGHT * coil_width
j_density = COIL_CURRENT / coil_cs

coil = rad.ObjRaceTrk(
    COIL_CENTER,
    COIL_RADII,
    COIL_STRAIGHT,
    COIL_HEIGHT,
    3,           # nseg (from verified script)
    'man',       # manual mode
    'z',         # axis
    j_density
)

print(f"  Coil ID: {coil}")
print(f"  Current: {COIL_CURRENT} A")
print(f"  Current density: {j_density:.1f} A/m^2")

# ============================================================
# Solve
# ============================================================
print("\nStep 4: Solve")

g = rad.ObjCnt([coil, yoke])

# Set solver tolerances (matching benchmark setup)
rad.SetBiCGSTABTol(1e-4)
rad.SetRelaxParam(0.0)

t_start = time.time()
result = rad.Solve(g, 0.0001, 100, 0)  # LU direct solver
t_solve = time.time() - t_start

# Get solve statistics
stats = rad.GetSolveStats()
nonl_iter = stats.get('nonl_iterations', 0)
linear_iter = stats.get('linear_iterations', 0)

print(f"  Residual: {result[0]:.6f}")
print(f"  Nonlinear iterations: {nonl_iter}")
print(f"  Linear iterations: {linear_iter}")
print(f"  Time: {t_solve:.2f} s")

# ============================================================
# Field at gap center
# ============================================================
print("\nStep 5: Field at observation points")

obs_points = [
    [0, 0, 0],
    [0, 0, 0.01],
    [0, 0, 0.02],
    [0.01, 0, 0.01],
]

results = []
for pt in obs_points:
    B = rad.Fld(g, 'b', pt)
    B_mag = np.linalg.norm(B)
    results.append({'point': pt, 'B': B, 'B_mag': B_mag})
    print(f"  B at {pt}: ({B[0]*1000:.4f}, {B[1]*1000:.4f}, {B[2]*1000:.4f}) mT, |B| = {B_mag*1000:.4f} mT")

# ============================================================
# Export VTS (TODO: fix GIL issue)
# ============================================================
# print("\nStep 6: Export field")
# vts_file = os.path.join(work_dir, 'field_distribution.vts')
# rad.FldVTS(g, vts_file, GRID_X, GRID_Y, GRID_Z,
#            GRID_NX, GRID_NY, GRID_NZ, 1, 0, 1.0)
# print(f"  Exported: {vts_file}")

# ============================================================
# Save results
# ============================================================
result_file = os.path.join(work_dir, 'results.txt')
with open(result_file, 'w') as f:
    f.write(f"Electromagnet Simulation Results\n")
    f.write(f"================================\n\n")
    f.write(f"Parameters:\n")
    f.write(f"  mu_r = {MU_R}\n")
    f.write(f"  Coil current = {COIL_CURRENT} A\n")
    f.write(f"  Mesh elements = {mesh.ne}\n\n")
    f.write(f"Solver:\n")
    f.write(f"  Residual = {result[0]:.6f}\n")
    f.write(f"  Iterations = {result[1]}\n")
    f.write(f"  Time = {t_solve:.2f} s\n\n")
    f.write(f"Field at observation points:\n")
    for r in results:
        f.write(f"  {r['point']}: |B| = {r['B_mag']*1000:.4f} mT\n")

print(f"  Results: {result_file}")

print("\n" + "=" * 60)
print("SIMULATION COMPLETE")
print("=" * 60)
