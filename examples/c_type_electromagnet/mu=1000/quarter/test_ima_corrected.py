#!/usr/bin/env python
"""
Test IMA with corrected understanding of symmetry.

Key insight:
Geometric mirroring FLIPS the magnetization direction.

IMA symmetry:
  - (antisymmetric): FLIPS M component again -> undoes geometric flip -> SAME M direction
  + (symmetric): PRESERVES M -> keeps geometric flip -> OPPOSITE M direction

So for uniform external field (both elements magnetized in same direction):
  N-S | N-S = opposite poles at boundary = 異極対称 = IMA -
  N-S | S-N = same poles at boundary = 同極対称 = IMA +

For Bz field with z-mirror, both cubes have M in +z:
  -> We want N-S | N-S (異極対称)
  -> Correct IMA: -z (antisymmetric)
"""

import sys
import os
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(work_dir)
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(parent_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))

import radia as rad

def create_cube_vertices(cx, cy, cz, size):
    s = size / 2
    return [
        [cx - s, cy - s, cz - s],
        [cx + s, cy - s, cz - s],
        [cx + s, cy + s, cz - s],
        [cx - s, cy + s, cz - s],
        [cx - s, cy - s, cz + s],
        [cx + s, cy - s, cz + s],
        [cx + s, cy + s, cz + s],
        [cx - s, cy + s, cz + s],
    ]

MU_R = 1000
cube_size = 0.02
y_center = 0.1

print("="*70)
print("IMA Test with Corrected Symmetry Understanding")
print("="*70)

print("""
Geometric mirror flips M direction.
IMA - (antisymmetric) flips M again -> same M direction -> 異極対称
IMA + (symmetric) preserves M -> opposite M directions -> 同極対称

For Bz field: both cubes have M in +z
  -> Want S|N at boundary (異極対称)
  -> Use -z (antisymmetric)
""")

# Vertices
v1 = create_cube_vertices(0, y_center, cube_size/2, cube_size)   # z > 0
v2 = create_cube_vertices(0, y_center, -cube_size/2, cube_size)  # z < 0

# ============================================================================
# Full model
# ============================================================================
print("="*70)
print("[1] Full Model (2 cubes)")
print("="*70)

rad.UtiDelAll()

c1 = rad.ObjHexahedron(v1, [0, 0, 0])
c2 = rad.ObjHexahedron(v2, [0, 0, 0])
mat = rad.MatLin(MU_R)
rad.MatApl(c1, mat)
rad.MatApl(c2, mat)
container = rad.ObjCnt([c1, c2])
bkg = rad.ObjBckg(lambda p: [0, 0, 0.1])
model = rad.ObjCnt([container, bkg])

result = rad.Solve(model, 1e-10, 100, 0)
print(f"Solve: {result}")

M1 = rad.Fld(c1, 'm', [0, y_center, cube_size/2])
M2 = rad.Fld(c2, 'm', [0, y_center, -cube_size/2])
print(f"Cube 1 M: {M1}")
print(f"Cube 2 M: {M2}")
print(f"Both Mz positive -> same direction -> 異極対称 -> use -z")

obs = [[0, y_center, 0], [0, y_center, 0.03], [0.03, y_center, 0]]
print("\nFields:")
B_full = {}
for pt in obs:
    B = rad.Fld(model, 'b', pt)
    B_full[tuple(pt)] = B[2]
    print(f"  {pt}: Bz = {B[2]*1000:.4f} mT")

# ============================================================================
# IMA with -z (antisymmetric) - CORRECT for this case
# ============================================================================
print("\n" + "="*70)
print("[2] IMA -z (antisymmetric) = 異極対称 - SHOULD BE CORRECT")
print("="*70)

rad.UtiDelAll()

c1_ima = rad.ObjHexahedron(v1, [0, 0, 0])
mat_ima = rad.MatLin(MU_R)
rad.MatApl(c1_ima, mat_ima)
container_ima = rad.ObjCnt([c1_ima])
bkg_ima = rad.ObjBckg(lambda p: [0, 0, 0.1])
model_ima = rad.ObjCnt([container_ima, bkg_ima])

result_ima = rad.Solve(model_ima, 1e-10, 100, 0, image='-z')
print(f"Solve: {result_ima}")

M_ima = rad.Fld(c1_ima, 'm', [0, y_center, cube_size/2])
print(f"Cube 1 M (IMA): {M_ima}")

print("\nFields:")
all_pass = True
for pt in obs:
    B = rad.Fld(model_ima, 'b', pt)
    ratio = B[2] / B_full[tuple(pt)] if abs(B_full[tuple(pt)]) > 1e-12 else 0
    status = "PASS" if abs(ratio - 1.0) < 0.001 else "FAIL"
    if status == "FAIL":
        all_pass = False
    print(f"  {pt}: Bz = {B[2]*1000:.4f} mT, ratio = {ratio:.6f} [{status}]")

# ============================================================================
# IMA with +z (symmetric) - WRONG for this case
# ============================================================================
print("\n" + "="*70)
print("[3] IMA +z (symmetric) = 同極対称 - SHOULD BE WRONG")
print("="*70)

rad.UtiDelAll()

c1_ima2 = rad.ObjHexahedron(v1, [0, 0, 0])
mat_ima2 = rad.MatLin(MU_R)
rad.MatApl(c1_ima2, mat_ima2)
container_ima2 = rad.ObjCnt([c1_ima2])
bkg_ima2 = rad.ObjBckg(lambda p: [0, 0, 0.1])
model_ima2 = rad.ObjCnt([container_ima2, bkg_ima2])

result_ima2 = rad.Solve(model_ima2, 1e-10, 100, 0, image='+z')
print(f"Solve: {result_ima2}")

print("\nFields:")
for pt in obs:
    B = rad.Fld(model_ima2, 'b', pt)
    ratio = B[2] / B_full[tuple(pt)] if abs(B_full[tuple(pt)]) > 1e-12 else 0
    print(f"  {pt}: Bz = {B[2]*1000:.4f} mT, ratio = {ratio:.6f}")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)

if all_pass:
    print("IMA -z (異極対称) gives correct results!")
else:
    print("IMA -z still has errors - need to debug implementation")
