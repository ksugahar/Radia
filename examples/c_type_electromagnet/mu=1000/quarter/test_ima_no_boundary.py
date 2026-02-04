#!/usr/bin/env python
"""
Test IMA with NO boundary elements.

Place cube AWAY from the symmetry plane to avoid the known limitation.
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
gap = 0.005  # 5mm gap from z=0

print("="*70)
print("IMA Test WITHOUT Boundary Elements")
print("="*70)

print(f"""
Cube is placed with 5mm gap from z=0:
  - Cube 1: z from {gap} to {gap + cube_size}
  - Cube 2 (mirror): z from {-gap - cube_size} to {-gap}

NO faces are on the symmetry plane z=0.
""")

# Vertices - cube is OFFSET from z=0
v1 = create_cube_vertices(0, y_center, gap + cube_size/2, cube_size)
v2 = create_cube_vertices(0, y_center, -(gap + cube_size/2), cube_size)

# ============================================================================
# Full model
# ============================================================================
print("="*70)
print("[1] Full Model (2 cubes with gap)")
print("="*70)

rad.UtiDelAll()
rad.FldUnits('m')

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

M1 = rad.Fld(c1, 'm', [0, y_center, gap + cube_size/2])
M2 = rad.Fld(c2, 'm', [0, y_center, -(gap + cube_size/2)])
print(f"Cube 1 M: {np.array(M1)}")
print(f"Cube 2 M: {np.array(M2)}")

obs = [
    [0, y_center, 0],           # In the gap (z=0)
    [0, y_center, 0.03],        # Off to z+
    [0.03, y_center, 0],        # Off to x+
    [0, y_center, gap/2],       # In gap but off-center
]
print("\nFields:")
B_full = {}
for pt in obs:
    B = rad.Fld(model, 'b', pt)
    B_full[tuple(pt)] = B[2]
    print(f"  {pt}: Bz = {B[2]*1000:.6f} mT")

# ============================================================================
# IMA with -z (antisymmetric)
# ============================================================================
print("\n" + "="*70)
print("[2] IMA -z (antisymmetric)")
print("="*70)

rad.UtiDelAll()
rad.FldUnits('m')

c1_ima = rad.ObjHexahedron(v1, [0, 0, 0])
mat_ima = rad.MatLin(MU_R)
rad.MatApl(c1_ima, mat_ima)
container_ima = rad.ObjCnt([c1_ima])
bkg_ima = rad.ObjBckg(lambda p: [0, 0, 0.1])
model_ima = rad.ObjCnt([container_ima, bkg_ima])

result_ima = rad.Solve(model_ima, 1e-10, 100, 0, image='-z')
print(f"Solve: {result_ima}")

print("\nFields:")
all_pass = True
for pt in obs:
    B = rad.Fld(model_ima, 'b', pt)
    ratio = B[2] / B_full[tuple(pt)] if abs(B_full[tuple(pt)]) > 1e-12 else 0
    status = "PASS" if abs(ratio - 1.0) < 0.001 else "FAIL"
    if status == "FAIL":
        all_pass = False
    print(f"  {pt}: Bz = {B[2]*1000:.6f} mT, ratio = {ratio:.6f} [{status}]")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)

if all_pass:
    print("\n*** PASS: IMA -z works correctly for NON-boundary elements ***")
else:
    print("\n*** FAIL: IMA has issues even for non-boundary elements ***")
