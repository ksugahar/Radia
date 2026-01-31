#!/usr/bin/env python
"""
Test IMA with a single hexahedral element.

This helps isolate the IMA field computation issue.
"""

import sys
import os
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(work_dir)
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(parent_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))

import radia as rad

print("=" * 70)
print("Single Element IMA Test")
print("=" * 70)

# Simple cube at (0.05, 0, 0) - in x > 0 region
# Size: 0.02 x 0.02 x 0.02 m
center = [0.05, 0.0, 0.0]
size = 0.02
half = size / 2

# Vertices in standard order
vertices = [
    [center[0] - half, center[1] - half, center[2] - half],  # 0
    [center[0] + half, center[1] - half, center[2] - half],  # 1
    [center[0] + half, center[1] + half, center[2] - half],  # 2
    [center[0] - half, center[1] + half, center[2] - half],  # 3
    [center[0] - half, center[1] - half, center[2] + half],  # 4
    [center[0] + half, center[1] - half, center[2] + half],  # 5
    [center[0] + half, center[1] + half, center[2] + half],  # 6
    [center[0] - half, center[1] + half, center[2] + half],  # 7
]

# Mirrored vertices (x -> -x)
vertices_mirror = [
    [-v[0], v[1], v[2]] for v in vertices
]

# Background field (y-directed)
B_ext_y = 0.1  # 100 mT
MU_R = 1000

# Observation point (on x=0 plane)
obs = [0.0, 0.0, 0.0]

print(f"Cube center: {center}")
print(f"Background B_y: {B_ext_y*1e3:.1f} mT")
print(f"Observation point: {obs}")

# Test 1: Full model (2 cubes - original + mirror)
print("\n--- Test 1: Full model (2 cubes) ---")

rad.UtiDelAll()
rad.FldUnits('m')

hex1 = rad.ObjHexahedron(vertices, [0, 0, 0])
hex2 = rad.ObjHexahedron(vertices_mirror, [0, 0, 0])

mat1 = rad.MatLin(MU_R)
mat2 = rad.MatLin(MU_R)
rad.MatApl(hex1, mat1)
rad.MatApl(hex2, mat2)

iron_full = rad.ObjCnt([hex1, hex2])
bkg = rad.ObjBckg(lambda p: [0, B_ext_y, 0])
container_full = rad.ObjCnt([iron_full, bkg])

result_full = rad.Solve(container_full, 0.0001, 100, 0)
print(f"Solve: {result_full}")

B_full = rad.Fld(container_full, 'b', obs)
print(f"B_full at obs: {[b*1e3 for b in B_full]} mT")

# Get field from each cube separately
B_hex1 = rad.Fld(hex1, 'b', obs)
B_hex2 = rad.Fld(hex2, 'b', obs)
print(f"B_hex1 at obs: {[b*1e3 for b in B_hex1]} mT")
print(f"B_hex2 at obs: {[b*1e3 for b in B_hex2]} mT")
print(f"Sum: {[(B_hex1[i]+B_hex2[i])*1e3 for i in range(3)]} mT")

# Test 2: Single cube only (no mirror)
print("\n--- Test 2: Single cube only ---")

rad.UtiDelAll()
rad.FldUnits('m')

hex_single = rad.ObjHexahedron(vertices, [0, 0, 0])
mat_single = rad.MatLin(MU_R)
rad.MatApl(hex_single, mat_single)

bkg_single = rad.ObjBckg(lambda p: [0, B_ext_y, 0])
container_single = rad.ObjCnt([hex_single, bkg_single])

result_single = rad.Solve(container_single, 0.0001, 100, 0)
print(f"Solve: {result_single}")

B_single = rad.Fld(container_single, 'b', obs)
print(f"B_single at obs: {[b*1e3 for b in B_single]} mT")

# Skip sigma retrieval (API not available)

# Test 3: IMA model
print("\n--- Test 3: IMA model (1 cube with +x mirror) ---")

rad.UtiDelAll()
rad.FldUnits('m')

hex_ima = rad.ObjHexahedron(vertices, [0, 0, 0])
mat_ima = rad.MatLin(MU_R)
rad.MatApl(hex_ima, mat_ima)

bkg_ima = rad.ObjBckg(lambda p: [0, B_ext_y, 0])
container_ima = rad.ObjCnt([hex_ima, bkg_ima])

result_ima = rad.Solve(container_ima, 0.0001, 100, 0, image='+x')
print(f"Solve: {result_ima}")

# Field without IMA context (direct from element)
B_ima_direct = rad.Fld(hex_ima, 'b', obs)
print(f"B_ima_direct at obs: {[b*1e3 for b in B_ima_direct]} mT")

# Field with IMA context
B_ima = rad.Fld(container_ima, 'b', obs)
print(f"B_ima at obs: {[b*1e3 for b in B_ima]} mT")

# Skip sigma retrieval (API not available)

# Summary
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)
print(f"Full model (2 cubes):      B = {[b*1e3 for b in B_full]} mT")
print(f"Single cube:               B = {[b*1e3 for b in B_single]} mT")
print(f"IMA direct (1 cube):       B = {[b*1e3 for b in B_ima_direct]} mT")
print(f"IMA with mirror (1 cube):  B = {[b*1e3 for b in B_ima]} mT")

# Check field symmetry
print("\nExpected behavior:")
print("  IMA with mirror should match Full model")
