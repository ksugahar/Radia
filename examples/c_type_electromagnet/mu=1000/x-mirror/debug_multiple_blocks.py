#!/usr/bin/env python
"""
Debug IMA permutation by checking multiple blocks, not just element 0.
Compare ELF x-mirror matrix with Radia for several diagonal and off-diagonal blocks.
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
ELF_FULL = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\full"

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


def build_radia_model(nodes, elements):
    """Build Radia model from ELF geometry."""
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
    return container


print("=" * 70)
print("Debug IMA permutation for multiple blocks")
print("=" * 70)

# Load ELF models
nodes_xm, elements_xm = load_elf_geometry(ELF_X_MIRROR)
nodes_full, elements_full = load_elf_geometry(ELF_FULL)

n_elem_xm = len(elements_xm)
n_elem_full = len(elements_full)
print(f"X-mirror model: {n_elem_xm} elements")
print(f"Full model: {n_elem_full} elements")

elf_mat_xm = load_elf_matrix(ELF_X_MIRROR, n_elem_xm * 6)
elf_mat_full = load_elf_matrix(ELF_FULL, n_elem_full * 6)

# Build Radia model
container = build_radia_model(nodes_xm, elements_xm)
handle = rad.BuildMatrix(container, image='+x')
radia_mat, dof = rad.GetInteractMatrix(handle)

print(f"\nELF IMA matrix: {elf_mat_xm.shape}")
print(f"Radia IMA matrix: {radia_mat.shape}")

# Find element correspondences between x-mirror and full models
# For each x-mirror element, find original and mirror in full model
def find_element_in_full(xm_idx, nodes_xm, elements_xm, nodes_full, elements_full):
    """Find the corresponding element in full model."""
    xm_verts = [nodes_xm[nid] for nid in elements_xm[xm_idx][1]]
    xm_centroid = np.mean(xm_verts, axis=0)

    for f_idx, (f_id, f_nodes) in enumerate(elements_full):
        f_verts = [nodes_full[nid] for nid in f_nodes]
        f_centroid = np.mean(f_verts, axis=0)
        if np.linalg.norm(f_centroid - xm_centroid) < 1.0:
            return f_idx
    return -1


def find_mirror_in_full(xm_idx, nodes_xm, elements_xm, nodes_full, elements_full):
    """Find the X-mirrored element in full model."""
    xm_verts = [nodes_xm[nid] for nid in elements_xm[xm_idx][1]]
    xm_centroid = np.mean(xm_verts, axis=0)
    mirror_centroid = np.array([-xm_centroid[0], xm_centroid[1], xm_centroid[2]])

    for f_idx, (f_id, f_nodes) in enumerate(elements_full):
        f_verts = [nodes_full[nid] for nid in f_nodes]
        f_centroid = np.mean(f_verts, axis=0)
        if np.linalg.norm(f_centroid - mirror_centroid) < 1.0:
            return f_idx
    return -1


# Test several elements
test_elements = [0, 5, 10, 15, 20, 25]

print("\n" + "=" * 70)
print("Analyzing diagonal blocks (self-interaction with mirror)")
print("=" * 70)

Q = [3, 2, 1, 0, 4, 5]  # ELF reordering

for xm_idx in test_elements:
    if xm_idx >= n_elem_xm:
        continue

    orig_idx = find_element_in_full(xm_idx, nodes_xm, elements_xm, nodes_full, elements_full)
    mirror_idx = find_mirror_in_full(xm_idx, nodes_xm, elements_xm, nodes_full, elements_full)

    if orig_idx < 0 or mirror_idx < 0:
        print(f"Element {xm_idx}: Could not find correspondence")
        continue

    # ELF IMA block for this element
    ELF_IMA_block = elf_mat_xm[xm_idx*6:(xm_idx+1)*6, xm_idx*6:(xm_idx+1)*6]

    # Full model blocks
    K_AA = elf_mat_full[orig_idx*6:(orig_idx+1)*6, orig_idx*6:(orig_idx+1)*6]
    K_AB = elf_mat_full[orig_idx*6:(orig_idx+1)*6, mirror_idx*6:(mirror_idx+1)*6]
    K_BA = elf_mat_full[mirror_idx*6:(mirror_idx+1)*6, orig_idx*6:(orig_idx+1)*6]

    # Expected mirror contribution
    Required = K_BA[Q, :]  # Q @ K_BA

    # Radia IMA block
    Radia_IMA_block = radia_mat[xm_idx*6:(xm_idx+1)*6, xm_idx*6:(xm_idx+1)*6]

    # Differences
    elf_diff = np.abs(ELF_IMA_block - K_AA - Required).max()
    radia_diff = np.abs(Radia_IMA_block - K_AA - Required).max()

    print(f"\nElement {xm_idx} (full: orig={orig_idx}, mirror={mirror_idx}):")
    print(f"  ELF_IMA - K_AA - Q@K_BA: max diff = {elf_diff:.2e}")
    print(f"  Radia_IMA - K_AA - Q@K_BA: max diff = {radia_diff:.2e}")

    # Check if K_AB can be transformed to Required
    row_perm = [0, 3, 2, 1, 4, 5]
    col_perm = [3, 0, 1, 2, 4, 5]
    K_AB_transformed = K_AB[row_perm, :][:, col_perm]
    transform_diff = np.abs(K_AB_transformed - Required).max()
    print(f"  K_AB[row_perm][:, col_perm] - Required: max diff = {transform_diff:.2e}")

    # What did Radia compute as mirror contribution?
    Radia_mirror = Radia_IMA_block - K_AA
    radia_to_required_diff = np.abs(Radia_mirror - Required).max()
    print(f"  Radia_mirror - Required: max diff = {radia_to_required_diff:.2e}")


print("\n" + "=" * 70)
print("Analyzing off-diagonal blocks (element-to-element with mirrors)")
print("=" * 70)

# Test a few off-diagonal blocks
test_pairs = [(0, 1), (0, 5), (1, 10), (5, 15)]

for xm_i, xm_j in test_pairs:
    if xm_i >= n_elem_xm or xm_j >= n_elem_xm:
        continue

    orig_i = find_element_in_full(xm_i, nodes_xm, elements_xm, nodes_full, elements_full)
    orig_j = find_element_in_full(xm_j, nodes_xm, elements_xm, nodes_full, elements_full)
    mirror_j = find_mirror_in_full(xm_j, nodes_xm, elements_xm, nodes_full, elements_full)

    if orig_i < 0 or orig_j < 0 or mirror_j < 0:
        print(f"Pair ({xm_i}, {xm_j}): Could not find correspondence")
        continue

    # ELF IMA off-diagonal block
    ELF_IMA_block = elf_mat_xm[xm_i*6:(xm_i+1)*6, xm_j*6:(xm_j+1)*6]

    # Full model blocks
    K_ij = elf_mat_full[orig_i*6:(orig_i+1)*6, orig_j*6:(orig_j+1)*6]  # original i from original j
    K_i_mj = elf_mat_full[orig_i*6:(orig_i+1)*6, mirror_j*6:(mirror_j+1)*6]  # original i from mirror j

    # Radia IMA block
    Radia_IMA_block = radia_mat[xm_i*6:(xm_i+1)*6, xm_j*6:(xm_j+1)*6]

    # The IMA formula for off-diagonal should be: K_IMA[i,j] = K[i,j] + sign * perm(K[i, mirror(j)])
    row_perm = [0, 3, 2, 1, 4, 5]
    col_perm = [3, 0, 1, 2, 4, 5]
    K_i_mj_transformed = K_i_mj[row_perm, :][:, col_perm]

    # Expected IMA
    ELF_expected = K_ij + K_i_mj_transformed  # sign = +1 for symmetric

    elf_diff = np.abs(ELF_IMA_block - ELF_expected).max()
    radia_diff = np.abs(Radia_IMA_block - ELF_expected).max()

    print(f"\nPair ({xm_i}, {xm_j}) [full: i={orig_i}, j={orig_j}, mj={mirror_j}]:")
    print(f"  ELF_IMA - (K_ij + perm(K_i_mj)): max diff = {elf_diff:.2e}")
    print(f"  Radia_IMA - (K_ij + perm(K_i_mj)): max diff = {radia_diff:.2e}")


print("\n" + "=" * 70)
print("Summary")
print("=" * 70)

# Overall comparison
diff = np.abs(elf_mat_xm - radia_mat).max()
rel_err = diff / np.abs(elf_mat_xm).max() * 100
print(f"Max |ELF - Radia|: {diff:.2e}")
print(f"Relative error: {rel_err:.2f}%")
