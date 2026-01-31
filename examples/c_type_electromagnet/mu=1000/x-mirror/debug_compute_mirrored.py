#!/usr/bin/env python
"""
Debug Compute6x6BlockMirrored: Compare Radia's raw mirrored block with ELF.
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

ELF_FULL = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\full"
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
print("Debug Compute6x6BlockMirrored")
print("=" * 70)

# Load ELF models
nodes_full, elements_full = load_elf_geometry(ELF_FULL)
nodes_xm, elements_xm = load_elf_geometry(ELF_X_MIRROR)

n_elem_full = len(elements_full)
n_elem_xm = len(elements_xm)
print(f"Full model: {n_elem_full} elements")
print(f"X-mirror model: {n_elem_xm} elements")

elf_mat_full = load_elf_matrix(ELF_FULL, n_elem_full * 6)
elf_mat_xm = load_elf_matrix(ELF_X_MIRROR, n_elem_xm * 6)

# Build Radia model WITHOUT IMA first
rad.UtiDelAll()
rad.FldUnits('m')

hex_objects = []
for elem_id, node_ids in elements_xm:
    verts = [[nodes_xm[nid][0] * scale, nodes_xm[nid][1] * scale, nodes_xm[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects.append(hex_obj)

container = rad.ObjCnt(hex_objects)

# Build normal matrix (no IMA)
handle_normal = rad.BuildMatrix(container)
radia_mat_normal, _ = rad.GetInteractMatrix(handle_normal)

# Build IMA matrix
rad.UtiDelAll()
rad.FldUnits('m')

hex_objects2 = []
for elem_id, node_ids in elements_xm:
    verts = [[nodes_xm[nid][0] * scale, nodes_xm[nid][1] * scale, nodes_xm[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects2.append(hex_obj)

container2 = rad.ObjCnt(hex_objects2)
handle_ima = rad.BuildMatrix(container2, image='+x')
radia_mat_ima, _ = rad.GetInteractMatrix(handle_ima)

print(f"\nRadia normal matrix: {radia_mat_normal.shape}")
print(f"Radia IMA matrix: {radia_mat_ima.shape}")

# For element 0 (orig_idx=0, mirror_idx=39 in full model)
orig_idx = 0
mirror_idx = 39
xm_idx = 0

# ELF blocks from full model
K_AA = elf_mat_full[orig_idx*6:(orig_idx+1)*6, orig_idx*6:(orig_idx+1)*6]
K_AB = elf_mat_full[orig_idx*6:(orig_idx+1)*6, mirror_idx*6:(mirror_idx+1)*6]
K_BA = elf_mat_full[mirror_idx*6:(mirror_idx+1)*6, orig_idx*6:(orig_idx+1)*6]
K_BB = elf_mat_full[mirror_idx*6:(mirror_idx+1)*6, mirror_idx*6:(mirror_idx+1)*6]

# ELF IMA block
ELF_IMA_block = elf_mat_xm[xm_idx*6:(xm_idx+1)*6, xm_idx*6:(xm_idx+1)*6]

# Radia blocks
Radia_normal_block = radia_mat_normal[xm_idx*6:(xm_idx+1)*6, xm_idx*6:(xm_idx+1)*6]
Radia_IMA_block = radia_mat_ima[xm_idx*6:(xm_idx+1)*6, xm_idx*6:(xm_idx+1)*6]

print("\n" + "=" * 70)
print("Element 0 diagonal block")
print("=" * 70)

print("\n--- ELF Full Model Blocks ---")
print("K_AA (self at orig from orig):")
print(K_AA)

print("\nK_AB (self at orig from mirror):")
print(K_AB)

print("\nK_BA (self at mirror from orig):")
print(K_BA)

print("\n--- Radia Blocks ---")
print("Radia normal block (should match K_AA):")
print(Radia_normal_block)

diff_normal = np.abs(Radia_normal_block - K_AA).max()
print(f"\nMax |Radia_normal - K_AA|: {diff_normal:.2e}")

print("\n--- Mirror Contribution ---")
# ELF's IMA formula: IMA = K_AA + Q @ K_BA
Q = [3, 2, 1, 0, 4, 5]
ELF_mirror_contrib = K_BA[Q, :]  # Q @ K_BA
print("ELF mirror contribution (Q @ K_BA):")
print(ELF_mirror_contrib)

# Radia's computed mirror contribution
Radia_mirror_contrib = Radia_IMA_block - Radia_normal_block
print("\nRadia mirror contribution (IMA - normal):")
print(Radia_mirror_contrib)

print(f"\nMax |Radia_mirror - ELF_mirror|: {np.abs(Radia_mirror_contrib - ELF_mirror_contrib).max():.2e}")

print("\n--- Transformation analysis ---")
# The formula we're using in C++:
# result = K_AB[row_perm][:, col_perm_inv]
# where row_perm = [0,3,2,1,4,5], col_perm = [1,2,3,0,4,5]
row_perm = [0, 3, 2, 1, 4, 5]
col_perm_inv = [3, 0, 1, 2, 4, 5]  # This is what we pass to ApplyDOFPermutation

# ApplyDOFPermutation does: result[:, perm[j]] = input[:, j]
# With perm=[1,2,3,0,4,5]: result[:, 1] = input[:, 0], etc.
# So result columns are reordered as: [input[:,3], input[:,0], input[:,1], input[:,2], input[:,4], input[:,5]]
# which is input[:, [3,0,1,2,4,5]] in numpy terms

# Let me simulate exactly what C++ does:
# 1. K_AB_col[:, perm[j]] = K_AB[:, j]  where perm = [1,2,3,0,4,5]
K_AB_col = np.zeros((6, 6))
perm = [1, 2, 3, 0, 4, 5]
for j in range(6):
    K_AB_col[:, perm[j]] = K_AB[:, j]

print("After ApplyDOFPermutation with perm=[1,2,3,0,4,5]:")
print(K_AB_col)

# 2. K_AB_perm[i, :] = K_AB_col[rowPerm[i], :]  where rowPerm = [0,3,2,1,4,5]
K_AB_perm = np.zeros((6, 6))
rowPerm = [0, 3, 2, 1, 4, 5]
for i in range(6):
    K_AB_perm[i, :] = K_AB_col[rowPerm[i], :]

print("\nAfter ApplyRowPermutation with rowPerm=[0,3,2,1,4,5]:")
print(K_AB_perm)

print(f"\nMax |K_AB_perm - ELF_mirror|: {np.abs(K_AB_perm - ELF_mirror_contrib).max():.2e}")

# Also try the numpy-style permutation for reference
numpy_row_perm = [0, 3, 2, 1, 4, 5]
numpy_col_perm = [3, 0, 1, 2, 4, 5]
K_AB_numpy = K_AB[numpy_row_perm, :][:, numpy_col_perm]
print(f"\nNumpy-style K_AB[{numpy_row_perm}][:, {numpy_col_perm}]:")
print(K_AB_numpy)
print(f"Max |K_AB_numpy - ELF_mirror|: {np.abs(K_AB_numpy - ELF_mirror_contrib).max():.2e}")

# What exactly is Radia computing for the mirror contribution?
print("\n--- Detailed comparison ---")
print("ELF_mirror_contrib:")
print(ELF_mirror_contrib)
print("\nRadia_mirror_contrib:")
print(Radia_mirror_contrib)
print("\nDifference:")
print(Radia_mirror_contrib - ELF_mirror_contrib)
