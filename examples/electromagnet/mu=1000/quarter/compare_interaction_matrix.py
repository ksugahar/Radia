#!/usr/bin/env python
"""
Compare Radia and ELF interaction matrices for quarter model.

ELF quarter model uses:
- MIMA X (symmetric mirror in X=0 plane)
- MIMA -Z (antisymmetric mirror in Z=0 plane)
- 13 hexahedral elements

The ELF approach computes:
  K_IMA[i,j] = K[i,j] + sign_x * K[i, mirror_x(j)] + sign_z * K[i, mirror_z(j)]
             + sign_x * sign_z * K[i, mirror_xz(j)]

where mirror_x(j) is the position of element j mirrored across x=0, etc.
"""

import sys
import os
import numpy as np
import struct

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(work_dir))))
sys.path.insert(0, os.path.join(repo_root, 'src'))

import radia as rad

# === Load ELF geometry ===
ELF_QUARTER_PATH = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\quater"
ELF_MEG_PATH = os.path.join(ELF_QUARTER_PATH, "ELF_magic.meg")
ELF_MAT_PATH = os.path.join(ELF_QUARTER_PATH, "ELF_magic.mat")

print("=" * 70)
print("Interaction Matrix Comparison: Radia vs ELF Quarter Model")
print("=" * 70)

# === Parse ELF geometry ===
print("\n--- Parsing ELF geometry ---")
nodes = {}
elements = []

with open(ELF_MEG_PATH, 'r') as f:
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

n_elem = len(elements)
n_dof = n_elem * 6
print(f"  Nodes: {len(nodes)}")
print(f"  Magnetic elements: {n_elem}")
print(f"  Total DOF: {n_dof}")

# Compute element centroids
centroids = []
for elem_id, node_ids in elements:
    verts = [nodes[nid] for nid in node_ids]
    centroid = np.mean(verts, axis=0)
    centroids.append(centroid)
    print(f"  Element {elem_id}: centroid = ({centroid[0]:.1f}, {centroid[1]:.1f}, {centroid[2]:.1f}) mm")

# === Load ELF interaction matrix ===
print("\n--- Loading ELF interaction matrix ---")
with open(ELF_MAT_PATH, 'rb') as f:
    data = f.read()

elf_mat = np.zeros((n_dof, n_dof))
offset = 0
for i in range(n_dof):
    rec_len = struct.unpack('<i', data[offset:offset+4])[0]
    row_data = np.frombuffer(data[offset+4:offset+4+rec_len], dtype=np.float64)
    elf_mat[i, :] = row_data
    offset += 4 + rec_len + 4

print(f"  ELF matrix shape: {elf_mat.shape}")
print(f"  ELF diagonal range: {elf_mat.diagonal().min():.4f} to {elf_mat.diagonal().max():.4f}")

# === Create Radia model and get interaction matrix ===
print("\n--- Creating Radia model ---")
MU_R = 1000
scale = 0.001  # mm to m

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
print(f"  Created {len(hex_objects)} hexahedral elements")

# Build interaction matrix with Image symmetry using new API
handle = rad.BuildMatrix(container, image='+x-z')
print(f"  BuildMatrix handle: {handle}")

# Get Radia interaction matrix for comparison
radia_mat, radia_dof = rad.GetInteractMatrix(handle)
print(f"  Radia matrix shape: {radia_mat.shape}, DOF: {radia_dof}")

print("\n--- Comparing diagonal blocks ---")

# Compare diagonal blocks between ELF and Radia
for i in range(min(3, n_elem)):
    block_start = i * 6
    block_end = (i + 1) * 6
    elf_block = elf_mat[block_start:block_end, block_start:block_end]
    radia_block = radia_mat[block_start:block_end, block_start:block_end]

    print(f"\nElement {i+1} diagonal block comparison:")
    print("  ELF:")
    for row in elf_block:
        print("    " + " ".join([f"{v:8.4f}" for v in row]))
    print("  Radia:")
    for row in radia_block:
        print("    " + " ".join([f"{v:8.4f}" for v in row]))

    diff = np.abs(elf_block - radia_block)
    print(f"  Max difference: {diff.max():.6f}")

# Check matrix structure
print("\n--- Matrix structure analysis ---")
print(f"\nELF matrix properties:")
print(f"  Diagonal sum: {np.trace(elf_mat):.4f}")
print(f"  Frobenius norm: {np.linalg.norm(elf_mat):.4f}")
print(f"  Max value: {np.max(np.abs(elf_mat)):.4f}")
print(f"  Min non-zero: {np.min(np.abs(elf_mat[elf_mat != 0])):.6f}")

# Check symmetry
sym_err_elf = np.abs(elf_mat - elf_mat.T).max()
print(f"  Symmetry error (max|A-A^T|): {sym_err_elf:.6f}")

print(f"\nRadia matrix properties:")
print(f"  Diagonal sum: {np.trace(radia_mat):.4f}")
print(f"  Frobenius norm: {np.linalg.norm(radia_mat):.4f}")
print(f"  Max value: {np.max(np.abs(radia_mat)):.4f}")
if np.any(radia_mat != 0):
    print(f"  Min non-zero: {np.min(np.abs(radia_mat[radia_mat != 0])):.6f}")

# Check symmetry
sym_err_radia = np.abs(radia_mat - radia_mat.T).max()
print(f"  Symmetry error (max|A-A^T|): {sym_err_radia:.6f}")

# Overall matrix comparison
print("\n--- Overall Matrix Comparison ---")
diff_mat = np.abs(elf_mat - radia_mat)
print(f"Max difference: {diff_mat.max():.6f}")
print(f"Mean difference: {diff_mat.mean():.6f}")
print(f"Relative error (Frobenius): {np.linalg.norm(diff_mat) / np.linalg.norm(elf_mat) * 100:.2f}%")

# Check block structure - ELF uses MIMA X and MIMA -Z
# This creates virtual elements through mirroring
print("\n--- ELF MIMA structure ---")
print("  MIMA X: symmetric mirror across X=0 plane")
print("  MIMA -Z: antisymmetric mirror across Z=0 plane")
print("  Combined effect: quarter model -> full model")

# DOF permutation for X and Z mirrors
# X-mirror: swaps DOF 1,3 (x- and x+ faces)
# Z-mirror: swaps DOF 4,5 (z- and z+ faces)
P_x = np.array([
    [1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, 1, 0, 0, 0],
    [0, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 0],
    [0, 0, 0, 0, 0, 1],
])

P_z = np.array([
    [1, 0, 0, 0, 0, 0],
    [0, 1, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 1, 0],
])

print("\nDOF permutation matrix P_x (X-mirror):")
print(P_x)
print("\nDOF permutation matrix P_z (Z-mirror):")
print(P_z)

# Save matrices for further analysis
np.save(os.path.join(work_dir, 'elf_quarter_matrix.npy'), elf_mat)
np.save(os.path.join(work_dir, 'radia_quarter_matrix.npy'), radia_mat)
print(f"\nSaved ELF matrix to elf_quarter_matrix.npy")
print(f"Saved Radia matrix to radia_quarter_matrix.npy")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
