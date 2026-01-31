#!/usr/bin/env python
"""
Debug kernel-based IMA implementation.

Compare IMA matrix elements from:
1. Radia kernel-based IMA (new implementation)
2. ELF reference matrix

The key insight from user:
- meg2vol with kkh ordering should give identical DOF ordering
- Mirroring should be done at the kernel level (not virtual elements/DOFs)
"""

import sys
import os
import numpy as np
import struct

work_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(work_dir)
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(parent_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))

import radia as rad

ELF_X_MIRROR = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\x-mirror"

MU_R = 1000
scale = 0.001


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


def load_elf_matrix(path, n_dof):
    with open(os.path.join(path, "ELF_magic.mat"), 'rb') as f:
        data = f.read()
    mat = np.zeros((n_dof, n_dof))
    offset = 0
    for i in range(n_dof):
        rec_len = struct.unpack('<i', data[offset:offset+4])[0]
        row_data = np.frombuffer(data[offset+4:offset+4+rec_len], dtype=np.float64)
        mat[i, :] = row_data
        offset += 4 + rec_len + 4
    return mat


print("=" * 70)
print("Debug Kernel-Based IMA Implementation")
print("=" * 70)

# Load ELF x-mirror geometry
nodes, elements = load_elf_geometry(ELF_X_MIRROR)
n_elem = len(elements)
n_dof = n_elem * 6
print(f"X-mirror model: {n_elem} elements, {n_dof} DOF")

# Load ELF IMA matrix (already includes mirror contributions)
elf_mat = load_elf_matrix(ELF_X_MIRROR, n_dof)
print(f"ELF matrix shape: {elf_mat.shape}")

# Build Radia model with IMA
rad.UtiDelAll()
rad.FldUnits('m')

hex_objects = []
for elem_id, node_ids in elements:
    verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects.append(hex_obj)

container = rad.ObjCnt(hex_objects)

# Build IMA matrix
handle = rad.BuildMatrix(container, image='+x')
radia_mat, radia_rhs = rad.GetInteractMatrix(handle)
print(f"Radia IMA matrix shape: {radia_mat.shape}")

# Compare matrices
print("\n" + "=" * 70)
print("Matrix Comparison")
print("=" * 70)

# Overall statistics
diff = radia_mat - elf_mat
max_diff = np.abs(diff).max()
rel_diff = np.abs(diff) / (np.abs(elf_mat) + 1e-30)
max_rel_diff = rel_diff.max()
mean_rel_diff = rel_diff.mean()

print(f"Max absolute difference: {max_diff:.6e}")
print(f"Max relative difference: {max_rel_diff*100:.2f}%")
print(f"Mean relative difference: {mean_rel_diff*100:.2f}%")

# Compare diagonal blocks (self-interaction with IMA contribution)
print("\n--- Diagonal Block Comparison (Element 0) ---")
elf_block = elf_mat[0:6, 0:6]
radia_block = radia_mat[0:6, 0:6]

print("ELF diagonal block [0,0]:")
print(elf_block)
print("\nRadia diagonal block [0,0]:")
print(radia_block)
print(f"\nDifference:")
print(radia_block - elf_block)

# Compare off-diagonal blocks
print("\n--- Off-Diagonal Block Comparison (Element 0 -> Element 1) ---")
elf_block_01 = elf_mat[0:6, 6:12]
radia_block_01 = radia_mat[0:6, 6:12]

print("ELF off-diagonal block [0,1]:")
print(elf_block_01)
print("\nRadia off-diagonal block [0,1]:")
print(radia_block_01)
print(f"\nDifference:")
print(radia_block_01 - elf_block_01)

# Test: What if we DON'T apply IMA in Radia?
print("\n" + "=" * 70)
print("Compare with NO IMA (direct interaction only)")
print("=" * 70)

rad.UtiDelAll()
rad.FldUnits('m')

hex_objects2 = []
for elem_id, node_ids in elements:
    verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects2.append(hex_obj)

container2 = rad.ObjCnt(hex_objects2)
handle2 = rad.BuildMatrix(container2)  # No IMA
radia_mat_no_ima, _ = rad.GetInteractMatrix(handle2)

print("Radia diagonal block [0,0] (no IMA):")
print(radia_mat_no_ima[0:6, 0:6])

print("\nDifference (Radia IMA - Radia no IMA):")
ima_contribution = radia_mat[0:6, 0:6] - radia_mat_no_ima[0:6, 0:6]
print(ima_contribution)

print("\nExpected IMA contribution from ELF:")
# ELF IMA = direct + mirror contribution
# If we assume ELF matrix IS the IMA matrix, then:
# mirror_contribution = ELF - direct
# But we don't have ELF direct matrix, so compare with Radia no-IMA
expected_mirror = elf_mat[0:6, 0:6] - radia_mat_no_ima[0:6, 0:6]
print(expected_mirror)

print("\nDifference (Radia mirror - Expected mirror):")
print(ima_contribution - expected_mirror)

# Analyze row by row
print("\n" + "=" * 70)
print("Row-by-row analysis of diagonal block")
print("=" * 70)

for i in range(6):
    print(f"\nRow {i} (face {i}):")
    print(f"  ELF:    {elf_block[i, :]}")
    print(f"  Radia:  {radia_block[i, :]}")
    print(f"  NoIMA:  {radia_mat_no_ima[i, :6]}")
    print(f"  IMA contrib (Radia): {ima_contribution[i, :]}")
    print(f"  IMA contrib (expected): {expected_mirror[i, :]}")
