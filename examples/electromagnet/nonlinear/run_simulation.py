#!/usr/bin/env python
"""
Electromagnet Simulation: Nonlinear B-H curve (20000 AT)

Uses B-H curve from ELF_MAGIC reference data.
Reference: S:/ELF_MAGIC/2020_03_07_CEFC_2020/model_C-Type/nonlinear_20000AT
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
# B-H Curve Data (from ELF_MAGIC)
# Format: [H (A/m), B (T)]
# ============================================================
BH_DATA = [
    [0, 0],
    [13.898, 0.22706],
    [27.797, 0.45413],
    [42.397, 0.68119],
    [61.416, 0.90826],
    [82.382, 1.1353],
    [144.669, 1.3624],
    [897.76, 1.5894],
    [4581.74, 1.8124],
    [17736.2, 2.01],
    [41339.3, 2.1332],
    [68321.8, 2.2],
    [95685.5, 2.2548],
    [123000, 2.2999],
    [151000, 2.3425],
    [179000, 2.3788],
    [207000, 2.415],
    [235000, 2.4513],
    [263000, 2.4875],
    [290000, 2.5238],
    [318000, 2.56],
]

# ============================================================
# Parameters
# ============================================================
COIL_CURRENT = 20000  # A (Ampere-turns, 20000 AT)

# Coil parameters (simplified for quarter model)
COIL_CENTER = [0, 0, 0.05]  # m
COIL_RADII = [0.02, 0.025]  # m
COIL_HEIGHT = 0.01  # m

print("=" * 60)
print(f"ELECTROMAGNET SIMULATION: Nonlinear B-H curve")
print(f"  Coil current: {COIL_CURRENT} AT")
print("=" * 60)

rad.FldUnits('m')

# ============================================================
# Load mesh
# ============================================================
print("\nStep 1: Load mesh")

vol_file = os.path.join(parent_dir, 'yoke_1x1x1_quarter.vol')

if not os.path.exists(vol_file):
    print(f"  [ERROR] Mesh file not found: {vol_file}")
    sys.exit(1)

ngmesh = NetgenMesh()
ngmesh.Load(vol_file)
mesh = Mesh(ngmesh)

print(f"  Loaded: {vol_file}")
print(f"  Elements: {mesh.ne}")

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

# Apply nonlinear material (B-H curve)
yoke_mat = rad.MatSatIsoTab(BH_DATA)
rad.MatApl(yoke, yoke_mat)
print(f"  Yoke ID: {yoke}")
print(f"  Material: Nonlinear B-H curve ({len(BH_DATA)} points)")

# ============================================================
# Create coil
# ============================================================
print("\nStep 3: Create coil")

coil_width = COIL_RADII[1] - COIL_RADII[0]
coil_cs = COIL_HEIGHT * coil_width
j_density = COIL_CURRENT / coil_cs

coil = rad.ObjArcCur(
    COIL_CENTER,
    COIL_RADII,
    [0, 2*np.pi],
    COIL_HEIGHT,
    36,
    'man',
    'z',
    j_density
)

print(f"  Coil ID: {coil}")
print(f"  Current: {COIL_CURRENT} AT")
print(f"  Current density: {j_density:.0f} A/m^2")

# ============================================================
# Solve
# ============================================================
print("\nStep 4: Solve (nonlinear)")

g = rad.ObjCnt([coil, yoke])

t_start = time.time()
# Use LU solver (Method 0) for nonlinear problems - more stable
result = rad.Solve(g, 0.001, 1000, 0)
t_solve = time.time() - t_start

print(f"  Converged: residual = {result[0]:.6f}")
print(f"  Iterations: {result[1]}")
print(f"  Time: {t_solve:.2f} s")

# ============================================================
# Field at observation points
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
# Save results
# ============================================================
result_file = os.path.join(work_dir, 'results.txt')
with open(result_file, 'w') as f:
    f.write(f"Electromagnet Simulation Results (Nonlinear)\n")
    f.write(f"=============================================\n\n")
    f.write(f"Parameters:\n")
    f.write(f"  Material: Nonlinear B-H curve ({len(BH_DATA)} points)\n")
    f.write(f"  Coil current = {COIL_CURRENT} AT\n")
    f.write(f"  Mesh elements = {mesh.ne}\n\n")
    f.write(f"Solver:\n")
    f.write(f"  Residual = {result[0]:.6f}\n")
    f.write(f"  Iterations = {result[1]}\n")
    f.write(f"  Time = {t_solve:.2f} s\n\n")
    f.write(f"Field at observation points:\n")
    for r in results:
        f.write(f"  {r['point']}: |B| = {r['B_mag']*1000:.4f} mT\n")
    f.write(f"\nB-H curve data:\n")
    f.write(f"  H (A/m)     B (T)\n")
    for h, b in BH_DATA:
        f.write(f"  {h:10.1f}  {b:.4f}\n")

print(f"  Results: {result_file}")

print("\n" + "=" * 60)
print("SIMULATION COMPLETE")
print("=" * 60)
