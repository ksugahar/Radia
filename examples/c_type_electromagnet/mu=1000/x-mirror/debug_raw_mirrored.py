#!/usr/bin/env python
"""
Debug: Compare raw Compute6x6BlockMirrored output with ELF's K_AB.
Check if element mappings are correct.
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
print("Debug: Element mapping and raw Compute6x6BlockMirrored")
print("=" * 70)

# Load ELF models
nodes_xm, elements_xm = load_elf_geometry(ELF_X_MIRROR)
nodes_full, elements_full = load_elf_geometry(ELF_FULL)

n_elem_xm = len(elements_xm)
n_elem_full = len(elements_full)

print(f"\nX-mirror model: {n_elem_xm} elements")
print(f"Full model: {n_elem_full} elements")

# Print first few element centroids from both models
print("\n--- X-mirror model element centroids (first 5) ---")
for idx in range(min(5, n_elem_xm)):
    verts = [nodes_xm[nid] for nid in elements_xm[idx][1]]
    centroid = np.mean(verts, axis=0)
    print(f"  XM elem {idx}: ({centroid[0]:8.2f}, {centroid[1]:8.2f}, {centroid[2]:8.2f})")

print("\n--- Full model element centroids (first 10) ---")
for idx in range(min(10, n_elem_full)):
    verts = [nodes_full[nid] for nid in elements_full[idx][1]]
    centroid = np.mean(verts, axis=0)
    print(f"  Full elem {idx}: ({centroid[0]:8.2f}, {centroid[1]:8.2f}, {centroid[2]:8.2f})")

# Find element mapping
print("\n--- Element mapping (x-mirror -> full) ---")
xm_to_full_orig = {}
xm_to_full_mirror = {}

for xm_idx in range(n_elem_xm):
    xm_verts = [nodes_xm[nid] for nid in elements_xm[xm_idx][1]]
    xm_centroid = np.mean(xm_verts, axis=0)
    mirror_centroid = np.array([-xm_centroid[0], xm_centroid[1], xm_centroid[2]])

    for f_idx, (f_id, f_nodes) in enumerate(elements_full):
        f_verts = [nodes_full[nid] for nid in f_nodes]
        f_centroid = np.mean(f_verts, axis=0)
        if np.linalg.norm(f_centroid - xm_centroid) < 1.0:
            xm_to_full_orig[xm_idx] = f_idx
        if np.linalg.norm(f_centroid - mirror_centroid) < 1.0:
            xm_to_full_mirror[xm_idx] = f_idx

for xm_idx in [0, 1, 2, 5, 10, 25]:
    if xm_idx < n_elem_xm:
        orig = xm_to_full_orig.get(xm_idx, -1)
        mirror = xm_to_full_mirror.get(xm_idx, -1)
        xm_verts = [nodes_xm[nid] for nid in elements_xm[xm_idx][1]]
        xm_centroid = np.mean(xm_verts, axis=0)
        print(f"  XM {xm_idx} at ({xm_centroid[0]:6.1f}, {xm_centroid[1]:6.1f}, {xm_centroid[2]:6.1f}) -> full orig={orig}, mirror={mirror}")

# Load matrices
elf_mat_full = load_elf_matrix(ELF_FULL, n_elem_full * 6)
elf_mat_xm = load_elf_matrix(ELF_X_MIRROR, n_elem_xm * 6)

# For element 0 of x-mirror model
xm_idx = 0
orig_idx = xm_to_full_orig.get(xm_idx, 0)
mirror_idx = xm_to_full_mirror.get(xm_idx, 39)

print(f"\n--- Analysis for x-mirror element {xm_idx} ---")
print(f"Full model: orig_idx={orig_idx}, mirror_idx={mirror_idx}")

# ELF blocks from full model
K_AA = elf_mat_full[orig_idx*6:(orig_idx+1)*6, orig_idx*6:(orig_idx+1)*6]
K_AB = elf_mat_full[orig_idx*6:(orig_idx+1)*6, mirror_idx*6:(mirror_idx+1)*6]
K_BA = elf_mat_full[mirror_idx*6:(mirror_idx+1)*6, orig_idx*6:(orig_idx+1)*6]

print("\nK_AA (self at orig):")
print(K_AA)

print("\nK_AB (at orig from mirror) - K[orig, mirror]:")
print(K_AB)

print("\nK_BA (at mirror from orig) - K[mirror, orig]:")
print(K_BA)

# Check if K_AB.T ≈ K_BA (reciprocity)
print(f"\nMax |K_AB.T - K_BA|: {np.abs(K_AB.T - K_BA).max():.2e}")

# Check the transformation
Q = [3, 2, 1, 0, 4, 5]
row_perm = [0, 3, 2, 1, 4, 5]
col_perm = [3, 0, 1, 2, 4, 5]

Required = K_BA[Q, :]  # Q @ K_BA
K_AB_transformed = K_AB[row_perm, :][:, col_perm]

print("\n--- Transformation verification ---")
print("Required = Q @ K_BA (ELF's expected mirror contribution):")
print(Required)

print("\nK_AB[row_perm][:, col_perm]:")
print(K_AB_transformed)

print(f"\nMax |K_AB_transformed - Required|: {np.abs(K_AB_transformed - Required).max():.2e}")

# Now check if K_BA is the same structure as the transpose
print("\n--- Structure comparison ---")
print("K_AB structure (first 4 elements of row 0):", K_AB[0, :4])
print("K_BA structure (first 4 elements of col 0):", K_BA[:4, 0])
print("K_AB.T structure (first 4 elements of row 0):", K_AB.T[0, :4])

# Try to find if there's a simpler relationship
print("\n--- Trying to find K_BA from K_AB ---")
# K_BA[i, j] should equal K_AB[j, i] if reciprocal
for i in range(3):
    for j in range(3):
        print(f"  K_BA[{i},{j}] = {K_BA[i,j]:.4f}, K_AB[{j},{i}] = {K_AB[j,i]:.4f}, diff = {K_BA[i,j] - K_AB[j,i]:.2e}")

# What if we use Compute6x6BlockMirroredTarget instead?
print("\n" + "=" * 70)
print("Checking if we should use Compute6x6BlockMirroredTarget instead")
print("=" * 70)

# Compute6x6BlockMirrored mirrors SOURCE -> should give K[orig, mirror] = K_AB
# Compute6x6BlockMirroredTarget mirrors TARGET -> should give K[mirror, orig] = K_BA

print("\nFor ELF IMA formula: K_IMA[i,i] = K[i,i] + Q @ K[mirror(i), i]")
print("                                = K_AA + Q @ K_BA")
print("\nSo we need K_BA = K[mirror, orig] from Compute6x6BlockMirroredTarget, NOT K_AB from Compute6x6BlockMirrored!")
