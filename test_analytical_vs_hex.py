"""
Compare ANALYTICAL tetrahedral method with STANDARD hexahedral method
"""
import os
os.environ['RADIA_TETRA_METHOD'] = 'ANALYTICAL'

import sys
sys.path.insert(0, r'S:/Radia/01_GitHub/build/Release')

from numpy import *
import radia as rad

print("="*70)
print("Compare: ANALYTICAL tetrahedron vs STANDARD hexahedron")
print("="*70)

rad.FldUnits('m')

# Material
mu_r = 100
ksi = mu_r - 1

# Background field
mu0 = 4*pi*1e-7
H_ext = 1000.0  # A/m
H_s_tesla = mu0 * H_ext

# Evaluation point
eval_pt = [0.05, 0.05, 0.2]

print(f"\nMaterial: mu_r = {mu_r}")
print(f"Background: H_ext = {H_ext} A/m (z-direction)")
print(f"Evaluation point: {eval_pt}")

# Reference: Hexahedral mesh (single cube, no subdivision)
print("\n" + "-"*70)
print("[Reference] Single Hexahedron (STANDARD method)")
print("-"*70)

cube_hex = rad.ObjRecMag([0, 0, 0], [0.1, 0.1, 0.1], [0, 0, 0])
mat = rad.MatLin([ksi, ksi], [0, 0, 1])
rad.MatApl(cube_hex, mat)

def uniform_field(pos):
    return [0.0, 0.0, H_s_tesla]

bg_hex = rad.ObjBckgCF(uniform_field)
container_hex = rad.ObjCnt([cube_hex, bg_hex])

result_hex = rad.Solve(container_hex, 0.0001, 1000)
print(f"Solve result: {result_hex}")

H_hex = rad.Fld(container_hex, 'h', eval_pt)
H_hex_mag = (H_hex[0]**2 + H_hex[1]**2 + H_hex[2]**2)**0.5

print(f"H = [{H_hex[0]:.4f}, {H_hex[1]:.4f}, {H_hex[2]:.4f}] A/m")
print(f"|H| = {H_hex_mag:.4f} A/m")

# Test: Single Tetrahedron (ANALYTICAL method)
print("\n" + "-"*70)
print("[Test] Single Tetrahedron (ANALYTICAL method)")
print("-"*70)

# Simple tetrahedron (corner of cube)
vertices = [
    [0, 0, 0],
    [0.1, 0, 0],
    [0, 0.1, 0],
    [0, 0, 0.1]
]

# Faces (1-indexed)
TETRA_FACES = [
    [1, 3, 2],  # Bottom (z=0)
    [1, 2, 4],  # Front (y=0)
    [1, 4, 3],  # Left (x=0)
    [2, 3, 4]   # Slanted face
]

tetra = rad.ObjPolyhdr(vertices, TETRA_FACES, [0, 0, 0])
rad.MatApl(tetra, mat)
bg_tetra = rad.ObjBckgCF(uniform_field)
container_tetra = rad.ObjCnt([tetra, bg_tetra])

result_tetra = rad.Solve(container_tetra, 0.0001, 1000)
print(f"Solve result: {result_tetra}")

H_tetra = rad.Fld(container_tetra, 'h', eval_pt)
H_tetra_mag = (H_tetra[0]**2 + H_tetra[1]**2 + H_tetra[2]**2)**0.5

print(f"H = [{H_tetra[0]:.4f}, {H_tetra[1]:.4f}, {H_tetra[2]:.4f}] A/m")
print(f"|H| = {H_tetra_mag:.4f} A/m")

# Compare
print("\n" + "="*70)
print("Comparison")
print("="*70)

error = abs(H_tetra_mag - H_hex_mag) / H_hex_mag * 100
print(f"Error: {error:.1f}%")

if error < 20:
    print("[PASS] Error < 20%")
else:
    print("[FAIL] Error >= 20%")

print("\n" + "="*70)
