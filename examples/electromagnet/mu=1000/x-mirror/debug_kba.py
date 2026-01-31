#!/usr/bin/env python
"""
Debug K_BA computation: compare Radia's Compute6x6BlockMirroredTarget with ELF's actual K_BA.
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
print("Debug K_BA computation")
print("=" * 70)

# Load ELF full model
nodes_full, elements_full = load_elf_geometry(ELF_FULL)
n_elem_full = len(elements_full)
elf_mat_full = load_elf_matrix(ELF_FULL, n_elem_full * 6)

# Load ELF x-mirror model
nodes_xm, elements_xm = load_elf_geometry(ELF_X_MIRROR)
n_elem_xm = len(elements_xm)
elf_mat_xm = load_elf_matrix(ELF_X_MIRROR, n_elem_xm * 6)

print(f"Full model: {n_elem_full} elements")
print(f"X-mirror model: {n_elem_xm} elements")

# Find element 0 in x-mirror and its correspondence in full model
xm_elem0_verts = [nodes_xm[nid] for nid in elements_xm[0][1]]
xm_elem0_centroid = np.mean(xm_elem0_verts, axis=0)
mirror_centroid = np.array([-xm_elem0_centroid[0], xm_elem0_centroid[1], xm_elem0_centroid[2]])

print(f"\nX-mirror element 0 centroid: {xm_elem0_centroid}")
print(f"Mirror centroid: {mirror_centroid}")

# Find corresponding elements in full model
orig_idx = -1
mirror_idx = -1
for f_idx, (f_id, f_nodes) in enumerate(elements_full):
    f_verts = [nodes_full[nid] for nid in f_nodes]
    f_centroid = np.mean(f_verts, axis=0)
    if np.linalg.norm(f_centroid - xm_elem0_centroid) < 1.0:
        orig_idx = f_idx
    if np.linalg.norm(f_centroid - mirror_centroid) < 1.0:
        mirror_idx = f_idx

print(f"Full model: orig_idx={orig_idx}, mirror_idx={mirror_idx}")

# Extract K_AA, K_AB, K_BA, K_BB from full model
K_AA = elf_mat_full[orig_idx*6:(orig_idx+1)*6, orig_idx*6:(orig_idx+1)*6]
K_AB = elf_mat_full[orig_idx*6:(orig_idx+1)*6, mirror_idx*6:(mirror_idx+1)*6]
K_BA = elf_mat_full[mirror_idx*6:(mirror_idx+1)*6, orig_idx*6:(orig_idx+1)*6]
K_BB = elf_mat_full[mirror_idx*6:(mirror_idx+1)*6, mirror_idx*6:(mirror_idx+1)*6]

ELF_IMA_block = elf_mat_xm[0:6, 0:6]

print("\n" + "=" * 70)
print("ELF Matrix Blocks (element 0)")
print("=" * 70)
print("\nK_AA (self-interaction of element 0):")
print(K_AA)
print("\nK_BA (field at mirror from original):")
print(K_BA)
print("\nELF IMA block (element 0,0):")
print(ELF_IMA_block)

# Compute Required = ELF_IMA - K_AA
Required = ELF_IMA_block - K_AA
print("\nRequired mirror contribution (ELF_IMA - K_AA):")
print(Required)

# Now create Radia model and get its K_BA equivalent
print("\n" + "=" * 70)
print("Radia Model")
print("=" * 70)

rad.UtiDelAll()
rad.FldUnits('m')

# Create hex elements for x-mirror model only
hex_objects = []
for elem_id, node_ids in elements_xm:
    verts = [[nodes_xm[nid][0] * scale, nodes_xm[nid][1] * scale, nodes_xm[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects.append(hex_obj)

container = rad.ObjCnt(hex_objects)

# Build Radia matrix with IMA
handle = rad.BuildMatrix(container, image='+x')
radia_mat, dof = rad.GetInteractMatrix(handle)

print(f"Radia IMA matrix shape: {radia_mat.shape}")

# Compare block 0,0
Radia_IMA_block = radia_mat[0:6, 0:6]
print("\nRadia IMA block (element 0,0):")
print(Radia_IMA_block)

# Compute Radia's implied mirror contribution
Radia_Required = Radia_IMA_block - K_AA
print("\nRadia's mirror contribution (Radia_IMA - K_AA):")
print(Radia_Required)

print("\n" + "=" * 70)
print("Comparison")
print("=" * 70)

diff = Radia_Required - Required
print(f"\nDifference (Radia - ELF) mirror contribution:")
print(diff)
print(f"\nMax |diff|: {np.abs(diff).max():.6e}")

# Check Q @ K_BA
Q = [3, 2, 1, 0, 4, 5]
Q_K_BA = K_BA[Q, :]
print("\nQ @ K_BA (expected mirror contribution):")
print(Q_K_BA)
print(f"\nQ @ K_BA matches Required: {np.allclose(Q_K_BA, Required)}")

# What does Radia's K_BA look like?
print("\n" + "=" * 70)
print("Analyzing Radia's K_BA")
print("=" * 70)

# The "K_BA" that Radia computes is different from ELF
# Let's see if we can find the pattern
print("\nK_BA (from ELF full model):")
print(K_BA)

print("\nRadia_Required (what Radia computes as mirror contribution):")
print(Radia_Required)

# Try to find if there's a row/column permutation that maps Radia to ELF
print("\n--- Trying to find relationship ---")

# Check each row of Radia_Required against K_BA
for i in range(6):
    min_dist = float('inf')
    best_match = -1
    for j in range(6):
        dist = np.linalg.norm(Radia_Required[i, :] - K_BA[j, :])
        if dist < min_dist:
            min_dist = dist
            best_match = j
    print(f"Radia_Required[{i}] closest to K_BA[{best_match}], dist={min_dist:.2e}")

# Check column-wise too
print("\nColumn-wise comparison:")
for i in range(6):
    min_dist = float('inf')
    best_match = -1
    for j in range(6):
        dist = np.linalg.norm(Radia_Required[:, i] - K_BA[:, j])
        if dist < min_dist:
            min_dist = dist
            best_match = j
    print(f"Radia_Required[:, {i}] closest to K_BA[:, {best_match}], dist={min_dist:.2e}")
