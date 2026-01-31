#!/usr/bin/env python
"""
Verify that Compute6x6BlockMirroredTarget already outputs in the correct order
and doesn't need Q permutation.

Hypothesis:
- Compute6x6BlockMirroredTarget uses original face indices to mirror coords/normals
- This produces output that's already in the "Q-permuted" order relative to ELF's K_BA
- So applying Q permutation is WRONG - it double-permutes!
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
print("Verify permutation hypothesis")
print("=" * 70)

# Load ELF data
nodes_full, elements_full = load_elf_geometry(ELF_FULL)
nodes_xm, elements_xm = load_elf_geometry(ELF_X_MIRROR)

n_elem_xm = len(elements_xm)
n_elem_full = len(elements_full)

elf_mat_full = load_elf_matrix(ELF_FULL, n_elem_full * 6)
elf_mat_xm = load_elf_matrix(ELF_X_MIRROR, n_elem_xm * 6)

# Build Radia models
print("\nBuilding Radia model (no IMA)...")
container_no_ima = build_radia_model(nodes_xm, elements_xm)
handle_no_ima = rad.BuildMatrix(container_no_ima)
radia_mat_no_ima, _ = rad.GetInteractMatrix(handle_no_ima)

print("Building Radia model (with IMA +x)...")
container_ima = build_radia_model(nodes_xm, elements_xm)
handle_ima = rad.BuildMatrix(container_ima, image='+x')
radia_mat_ima, _ = rad.GetInteractMatrix(handle_ima)

print(f"Radia no IMA: {radia_mat_no_ima.shape}")
print(f"Radia with IMA: {radia_mat_ima.shape}")

# For element 0
xm_idx = 0
orig_idx = 0  # From previous script output
mirror_idx = 39

# ELF blocks
K_AA = elf_mat_full[orig_idx*6:(orig_idx+1)*6, orig_idx*6:(orig_idx+1)*6]
K_BA = elf_mat_full[mirror_idx*6:(mirror_idx+1)*6, orig_idx*6:(orig_idx+1)*6]

# Radia's contribution
Radia_no_ima = radia_mat_no_ima[xm_idx*6:(xm_idx+1)*6, xm_idx*6:(xm_idx+1)*6]
Radia_ima = radia_mat_ima[xm_idx*6:(xm_idx+1)*6, xm_idx*6:(xm_idx+1)*6]
Radia_mirror = Radia_ima - Radia_no_ima  # What Radia computes as mirror contribution

print("\n" + "=" * 70)
print("Analyzing permutation")
print("=" * 70)

# ELF's expected mirror contribution: Q @ K_BA
Q = [3, 2, 1, 0, 4, 5]
ELF_mirror = K_BA[Q, :]

print("ELF expected (Q @ K_BA):")
print(ELF_mirror)

print("\nRadia computed mirror contribution:")
print(Radia_mirror)

print(f"\nMax |Radia - ELF_expected|: {np.abs(Radia_mirror - ELF_mirror).max():.2e}")

# Hypothesis: Radia already outputs Q @ K_BA, and then applies Q again
# So Radia actually outputs: Q @ (Q @ K_BA) = K_BA (since Q is self-inverse)
print("\n--- Testing hypothesis ---")

print("\nIf Radia outputs K_BA (without Q permutation):")
print(K_BA)
print(f"Max |Radia - K_BA|: {np.abs(Radia_mirror - K_BA).max():.2e}")

# Also try Q^(-1) @ Radia = Q @ Radia (since Q = Q^(-1))
Radia_Q = Radia_mirror[Q, :]
print("\nIf we apply Q to Radia's output:")
print(f"Max |Q @ Radia - ELF_expected|: {np.abs(Radia_Q - ELF_mirror).max():.2e}")

# What if Radia's output needs inverse permutation?
Q_inv = [3, 2, 1, 0, 4, 5]  # Q is its own inverse
Radia_Q_inv = Radia_mirror[Q_inv, :]
print(f"Max |Q^(-1) @ Radia - ELF_expected|: {np.abs(Radia_Q_inv - ELF_mirror).max():.2e}")

# Also check the column permutation
col_perm = [3, 0, 1, 2, 4, 5]
Radia_col = Radia_mirror[:, col_perm]
print(f"Max |Radia[:, col_perm] - ELF_expected|: {np.abs(Radia_col - ELF_mirror).max():.2e}")

# Try both row and column permutation
Radia_both = Radia_mirror[Q, :][:, col_perm]
print(f"Max |Q @ Radia[:, col_perm] - ELF_expected|: {np.abs(Radia_both - ELF_mirror).max():.2e}")

print("\n" + "=" * 70)
print("Detailed comparison: row by row")
print("=" * 70)

for i in range(6):
    print(f"\nRow {i}:")
    print(f"  ELF_expected: {ELF_mirror[i, :]}")
    print(f"  Radia:        {Radia_mirror[i, :]}")
    print(f"  K_BA[{i}]:     {K_BA[i, :]}")

    # Find which row of Radia matches which row of ELF
    for j in range(6):
        diff = np.abs(Radia_mirror[j, :] - ELF_mirror[i, :]).max()
        if diff < 0.01:
            print(f"  -> Radia row {j} matches ELF row {i} (diff={diff:.2e})")
