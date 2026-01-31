#!/usr/bin/env python
"""
Test IMA by comparing field results, not raw matrix values.

The key test: Does the x-mirror IMA model give the same field as the full model?
"""

import sys
import os
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(work_dir)
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(parent_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))

import radia as rad

ELF_X_MIRROR = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\x-mirror"

MU_R = 1000
scale = 0.001

# Coil parameters
N_TURNS = 1000
CURRENT = 2.0  # A
COIL_AT = N_TURNS * CURRENT  # 2000 AT

# Gap center - in y direction based on geometry
GAP_CENTER = [0.0, 0.055, 0.0275]  # in meters


def load_elf_geometry(path):
    nodes = {}
    elements = []
    with open(os.path.join(path, "ELF_magic.meg"), 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('MGR1'):
                parts = line.split()
                node_id = int(parts[1])
                x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                nodes[node_id] = np.array([x, y, z])
            elif line.startswith('MMB8T'):
                parts = line.split()
                elem_id = int(parts[1])
                node_ids = [int(parts[i]) for i in range(4, 12)]
                elements.append((elem_id, node_ids))
    return nodes, elements


print("=" * 70)
print("Test IMA Field Result")
print("=" * 70)

# Load ELF x-mirror geometry
nodes, elements = load_elf_geometry(ELF_X_MIRROR)
n_elem = len(elements)
print(f"X-mirror model: {n_elem} elements")

# Calculate background field
H_ext = COIL_AT / (2 * 0.055)  # Rough approximation
B_ext_y = 4e-7 * np.pi * H_ext
print(f"External field approximation: B_y = {B_ext_y*1e3:.2f} mT")

# Test 1: Full model (explicit duplication, no IMA)
print("\n--- Test 1: Full model (explicit duplication) ---")

rad.UtiDelAll()
rad.FldUnits('m')

# Create iron elements for x-mirror (original)
hex_objects_orig = []
for elem_id, node_ids in elements:
    verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects_orig.append(hex_obj)

# Create mirrored iron elements (x -> -x)
hex_objects_mirror = []
for elem_id, node_ids in elements:
    verts_mirror = [[-nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
                     for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts_mirror, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects_mirror.append(hex_obj)

iron_full = rad.ObjCnt(hex_objects_orig + hex_objects_mirror)
print(f"Full model: {len(hex_objects_orig) + len(hex_objects_mirror)} elements")

# Apply uniform background field
bkg_full = rad.ObjBckg(lambda p: [0, B_ext_y, 0])
container_full = rad.ObjCnt([iron_full, bkg_full])

print("Solving full model...")
result_full = rad.Solve(container_full, 0.0001, 1000, 0)  # LU solver
print(f"Solve result: {result_full}")

B_full = rad.Fld(container_full, 'b', GAP_CENTER)
print(f"B at gap center: Bx={B_full[0]*1e3:.2f} mT, By={B_full[1]*1e3:.2f} mT, Bz={B_full[2]*1e3:.2f} mT")

# Test 2: IMA model (x-mirror symmetry)
print("\n--- Test 2: IMA model (x-mirror symmetry) ---")

rad.UtiDelAll()
rad.FldUnits('m')

hex_objects_ima = []
for elem_id, node_ids in elements:
    verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects_ima.append(hex_obj)

iron_ima = rad.ObjCnt(hex_objects_ima)

# Apply background field
bkg_ima = rad.ObjBckg(lambda p: [0, B_ext_y, 0])
container_ima = rad.ObjCnt([iron_ima, bkg_ima])

print("Solving IMA model with x-mirror symmetry...")
result_ima = rad.Solve(container_ima, 0.0001, 1000, 0, image='+x')
print(f"Solve result: {result_ima}")

B_ima = rad.Fld(container_ima, 'b', GAP_CENTER)
print(f"B at gap center: Bx={B_ima[0]*1e3:.2f} mT, By={B_ima[1]*1e3:.2f} mT, Bz={B_ima[2]*1e3:.2f} mT")

# Summary
print("\n" + "=" * 70)
print("Summary")
print("=" * 70)
print(f"Full model:  By = {B_full[1]*1e3:.2f} mT, Bz = {B_full[2]*1e3:.2f} mT")
print(f"IMA model:   By = {B_ima[1]*1e3:.2f} mT, Bz = {B_ima[2]*1e3:.2f} mT")

# Compute difference
diff_y = abs(B_full[1] - B_ima[1])
diff_z = abs(B_full[2] - B_ima[2])
rel_err_y = diff_y / abs(B_full[1]) * 100 if abs(B_full[1]) > 1e-10 else 0
rel_err_z = diff_z / abs(B_full[2]) * 100 if abs(B_full[2]) > 1e-10 else 0

print(f"\nDifference By: {diff_y*1e3:.2f} mT ({rel_err_y:.2f}%)")
print(f"Difference Bz: {diff_z*1e3:.2f} mT ({rel_err_z:.2f}%)")

if max(rel_err_y, rel_err_z) < 1.0:
    print("\n*** IMA PASSED: Error < 1% ***")
else:
    print(f"\n*** IMA FAILED: Error = {max(rel_err_y, rel_err_z):.2f}% ***")
