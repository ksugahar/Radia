#!/usr/bin/env python
"""
Find the exact IMA formula by exhaustive search over permutations.
"""

import sys
import os
import numpy as np
import struct
from itertools import permutations

work_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(work_dir)
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(parent_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))

ELF_X_MIRROR = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\x-mirror"
ELF_FULL = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\full"

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


# Load data
nodes_full, elements_full = load_elf_geometry(ELF_FULL)
nodes_xmirror, elements_xmirror = load_elf_geometry(ELF_X_MIRROR)

elf_mat_full = load_elf_matrix(ELF_FULL, len(elements_full) * 6)
elf_mat_xmirror = load_elf_matrix(ELF_X_MIRROR, len(elements_xmirror) * 6)

# Find mirror element mapping
print("Finding mirror element correspondence...")
mirror_map = {}  # xmirror_elem -> (full_orig_elem, full_mirror_elem)

for xm_idx, (xm_id, xm_nodes) in enumerate(elements_xmirror):
    xm_verts = [nodes_xmirror[nid] for nid in xm_nodes]
    xm_centroid = np.mean(xm_verts, axis=0)
    mirror_centroid = np.array([-xm_centroid[0], xm_centroid[1], xm_centroid[2]])

    orig_idx = -1
    mirror_idx = -1
    for f_idx, (f_id, f_nodes) in enumerate(elements_full):
        f_verts = [nodes_full[nid] for nid in f_nodes]
        f_centroid = np.mean(f_verts, axis=0)
        if np.linalg.norm(f_centroid - xm_centroid) < 1.0:
            orig_idx = f_idx
        if np.linalg.norm(f_centroid - mirror_centroid) < 1.0:
            mirror_idx = f_idx

    mirror_map[xm_idx] = (orig_idx, mirror_idx)

print(f"Element 0: orig={mirror_map[0][0]}, mirror={mirror_map[0][1]}")

# Get key matrices for element 0
orig_idx, mirror_idx = mirror_map[0]
K_AA = elf_mat_full[orig_idx*6:(orig_idx+1)*6, orig_idx*6:(orig_idx+1)*6]
K_AB = elf_mat_full[orig_idx*6:(orig_idx+1)*6, mirror_idx*6:(mirror_idx+1)*6]
K_BA = elf_mat_full[mirror_idx*6:(mirror_idx+1)*6, orig_idx*6:(orig_idx+1)*6]
K_BB = elf_mat_full[mirror_idx*6:(mirror_idx+1)*6, mirror_idx*6:(mirror_idx+1)*6]

ELF_IMA_block = elf_mat_xmirror[0:6, 0:6]
Required = ELF_IMA_block - K_AA

print("\nRequired mirror contribution:")
print(Required)

# Standard permutations
P_X = np.array([0, 3, 2, 1, 4, 5])  # x-mirror permutation

def perm_to_matrix(p):
    """Convert permutation array to permutation matrix."""
    n = len(p)
    M = np.zeros((n, n))
    for i in range(n):
        M[p[i], i] = 1.0  # Column permutation
    return M

def apply_perm_rows(A, p):
    """Apply row permutation: result[i] = A[p[i]]"""
    return A[p, :]

def apply_perm_cols(A, p):
    """Apply column permutation: result[:, i] = A[:, p[i]]"""
    return A[:, p]

print("\n" + "=" * 70)
print("Testing known permutations and matrices")
print("=" * 70)

# Test: Required = K_BA with rows and columns permuted somehow
# Try to find the row permutation that makes Required rows match K_BA rows

print("\nSearching for row mapping Required -> K_BA...")
row_map = []
for req_row in range(6):
    best_match = -1
    best_dist = float('inf')
    for kba_row in range(6):
        dist = np.linalg.norm(Required[req_row, :] - K_BA[kba_row, :])
        if dist < best_dist:
            best_dist = dist
            best_match = kba_row
    row_map.append(best_match)
    print(f"  Required[{req_row}] matches K_BA[{best_match}], dist={best_dist:.2e}")

print(f"\nRow mapping: Required rows [0,1,2,3,4,5] come from K_BA rows {row_map}")

# Check if this works exactly
test_K = K_BA[row_map, :]
diff = Required - test_K
print(f"\nDifference after row permutation only: max={np.abs(diff).max():.2e}")

# Now try column permutations on K_BA[row_map]
print("\nSearching for column permutation on K_BA[row_map]...")
for col_perm in [[0,1,2,3,4,5], [0,3,2,1,4,5], list(P_X)]:
    test_K = K_BA[row_map, :][:, col_perm]
    diff = Required - test_K
    print(f"  Column perm {col_perm}: max diff = {np.abs(diff).max():.2e}")

# The row_map [3,2,1,0,4,5] is unusual. Let's check what permutation matrix this is.
row_perm_arr = np.array(row_map)
print(f"\nRow permutation array: {row_map}")

# Check if row_perm is related to P_X
# P_X swaps positions 1 and 3
# row_map swaps (0,3) and (1,2)

# Try formula: Required = Q @ K_BA where Q is the row permutation matrix
Q = np.zeros((6, 6))
for i in range(6):
    Q[i, row_map[i]] = 1.0

print("\nQ (row permutation matrix):")
print(Q.astype(int))

# Verify: Q @ K_BA should equal Required
test = Q @ K_BA
print(f"\nQ @ K_BA matches Required: {np.allclose(test, Required)}, max diff = {np.abs(test - Required).max():.2e}")

# Now, what is Q in terms of face operations?
# Q[i,j] = 1 means row j of K_BA goes to row i of Q@K_BA
# So Q moves K_BA rows according to row_map
print("\nQ analysis:")
print("  Q swaps row 0 <-> row 3")
print("  Q swaps row 1 <-> row 2")
print("  Q keeps rows 4, 5 in place")

# This is a 180-degree rotation around z-axis in face space!
# Let's call it R_z (rotation permutation around z)
R_z = [3, 2, 1, 0, 4, 5]  # 180-deg rotation around z

# Check: P_X @ R_z ?
P_X_arr = list(P_X)
P_X_R_z = [P_X_arr[R_z[i]] for i in range(6)]
print(f"\nP_X @ R_z = {P_X_R_z}")

R_z_P_X = [R_z[P_X_arr[i]] for i in range(6)]
print(f"R_z @ P_X = {R_z_P_X}")

# Hmm, neither gives the identity or a simple form.
# Let me check if Q = R_z = the transpose of something related to P_X

# Actually, looking at the face ordering:
# Face 0: y-, Face 1: x+, Face 2: y+, Face 3: x-
# Under 180-deg z-rotation:
# y- -> y+, y+ -> y-, x+ -> x-, x- -> x+
# So faces: 0 <-> 2, 1 <-> 3
# This gives permutation [2, 3, 0, 1, 4, 5]

R_z_180 = [2, 3, 0, 1, 4, 5]
print(f"\n180-deg rotation around z: {R_z_180}")
print(f"But we found Q = {row_map}")

# These don't match. Let me check the actual geometry transformation.
# For x-mirror: position (x,y,z) -> (-x,y,z)
# Face normals: (nx, ny, nz) -> (-nx, ny, nz)
# So:
# Face 0 (normal -y) -> (normal -y) stays face 0
# Face 1 (normal +x) -> (normal -x) becomes like face 3
# Face 2 (normal +y) -> (normal +y) stays face 2
# Face 3 (normal -x) -> (normal +x) becomes like face 1
# Face 4 (normal -z) -> (normal -z) stays face 4
# Face 5 (normal +z) -> (normal +z) stays face 5

print("\n" + "=" * 70)
print("Physical interpretation of Q")
print("=" * 70)
print("Q swaps: 0<->3, 1<->2")
print("This corresponds to swapping (y-, x-) and (x+, y+)")
print("This is NOT the same as x-mirror face permutation P_X")
print("")
print("Possible explanation:")
print("- ELF may use a different face ordering for the mirror element")
print("- OR the mirror element index 39 in full model has different face ordering")
print("")

# Let's check if the formula might be K_IMA = K_AA + Q @ K_BA
# where Q is this strange permutation
K_IMA_test = K_AA + Q @ K_BA
diff_to_elf = ELF_IMA_block - K_IMA_test
print(f"Testing K_IMA = K_AA + Q @ K_BA:")
print(f"  Max diff to ELF = {np.abs(diff_to_elf).max():.2e}")

# Perfect match? Let's verify
if np.abs(diff_to_elf).max() < 1e-10:
    print("  *** EXACT MATCH! ***")
    print(f"  Formula: K_IMA = K_AA + Q @ K_BA where Q = {row_map}")
