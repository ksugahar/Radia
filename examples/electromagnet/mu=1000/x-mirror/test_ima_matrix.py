#!/usr/bin/env python
"""
Test IMA (Image) symmetry matrix against ELF x-mirror reference.
This validates the IMA implementation for MSC hexahedra.
"""

import sys
import os
import struct
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(work_dir))))
sys.path.insert(0, os.path.join(repo_root, 'src'))

import radia as rad

# === Load ELF x-mirror matrix (26 elements) ===
ELF_MAT_PATH = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\x-mirror\ELF_magic.mat"
ELF_MEG_PATH = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\x-mirror\ELF_magic.meg"

# Also load full model for IMA construction
ELF_FULL_MEG_PATH = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\full\ELF_magic.meg"

print("=" * 70)
print("IMA (Image) Symmetry Matrix Validation")
print("Comparing Radia IMA vs ELF x-mirror reference")
print("=" * 70)

# === Parse full model geometry (52 elements) ===
print("\nParsing full model geometry (52 elements)...")
nodes_full = {}
elements_full = []

with open(ELF_FULL_MEG_PATH, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('MGR1'):
            parts = line.split()
            node_id = int(parts[1])
            x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
            nodes_full[node_id] = np.array([x, y, z])
        elif line.startswith('MMB8T'):
            parts = line.split()
            elem_id = int(parts[1])
            node_ids = [int(parts[i]) for i in range(4, 12)]
            elements_full.append((elem_id, node_ids))

n_elem_full = len(elements_full)
print(f"  Full model: {n_elem_full} elements")

# === Parse x-mirror model geometry (26 elements) ===
print("\nParsing x-mirror model geometry (26 elements)...")
nodes_xm = {}
elements_xm = []

with open(ELF_MEG_PATH, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('MGR1'):
            parts = line.split()
            node_id = int(parts[1])
            x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
            nodes_xm[node_id] = np.array([x, y, z])
        elif line.startswith('MMB8T'):
            parts = line.split()
            elem_id = int(parts[1])
            node_ids = [int(parts[i]) for i in range(4, 12)]
            elements_xm.append((elem_id, node_ids))

n_elem_xm = len(elements_xm)
dof_xm = n_elem_xm * 6
print(f"  X-mirror model: {n_elem_xm} elements, {dof_xm} DOF")

# === Load ELF x-mirror matrix ===
print(f"\nLoading ELF x-mirror matrix ({dof_xm}x{dof_xm})...")
elf_matrix = np.zeros((dof_xm, dof_xm))
with open(ELF_MAT_PATH, 'rb') as f:
    for col in range(dof_xm):
        rec_len1 = struct.unpack('<I', f.read(4))[0]
        column_data = struct.unpack(f'<{dof_xm}d', f.read(dof_xm * 8))
        elf_matrix[:, col] = column_data
        rec_len2 = struct.unpack('<I', f.read(4))[0]

print(f"  ELF matrix loaded: {elf_matrix.shape}")

# === Build Radia model (full 52 elements) ===
print(f"\nBuilding Radia full model (52 elements)...")
rad.UtiDelAll()
rad.FldUnits('m')

MU_R = 1000
scale = 0.001  # mm to m

hex_objects = []
for elem_id, node_ids in elements_full:
    verts = [[nodes_full[nid][0] * scale, nodes_full[nid][1] * scale, nodes_full[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects.append(hex_obj)

print(f"  Created {len(hex_objects)} hexahedra")

# Create container
container = rad.ObjCnt(hex_objects)

# === Set IMA symmetry and build matrix ===
print("\nSetting IMA x-mirror symmetry...")
try:
    # PreRelax to setup interaction
    intrc_handle = rad.PreRelax(container, container)
    print(f"  Interaction handle: {intrc_handle}")

    # Set IMA symmetry (x-mirror)
    n_ima_elem = rad.SetIMASymmetry(intrc_handle, "x")  # X-mirror
    print(f"  IMA elements: {n_ima_elem}")

    # Build IMA matrix
    print("\nBuilding IMA matrix...")
    rad.BuildIMAMatrix(intrc_handle)

    # Get the matrix
    radia_ima_matrix, radia_dof = rad.GetInteractMatrix(intrc_handle)
    print(f"  Radia IMA matrix shape: {radia_ima_matrix.shape}, DOF: {radia_dof}")

except Exception as e:
    print(f"  Error: {e}")
    import traceback
    traceback.print_exc()
    radia_ima_matrix = None

if radia_ima_matrix is not None:
    # === Compare matrices ===
    print("\n" + "=" * 70)
    print("Matrix Comparison: Radia IMA vs ELF x-mirror")
    print("=" * 70)

    # Check dimensions
    if radia_ima_matrix.shape != elf_matrix.shape:
        print(f"  WARNING: Shape mismatch!")
        print(f"    Radia: {radia_ima_matrix.shape}")
        print(f"    ELF:   {elf_matrix.shape}")
    else:
        diff = radia_ima_matrix - elf_matrix
        max_diff = np.max(np.abs(diff))
        rms_diff = np.sqrt(np.mean(diff ** 2))
        max_elf = np.max(np.abs(elf_matrix))
        max_radia = np.max(np.abs(radia_ima_matrix))

        print(f"Max |ELF|:   {max_elf:.6e}")
        print(f"Max |Radia|: {max_radia:.6e}")
        print(f"Max |Radia - ELF|: {max_diff:.6e}")
        print(f"Relative max diff: {max_diff / max_elf * 100:.4f}%")
        print(f"RMS difference: {rms_diff:.6e}")

        # === Diagonal comparison ===
        print("\n--- Diagonal Block [1,1] Comparison ---")
        print("ELF:")
        for i in range(6):
            row_str = " ".join([f"{elf_matrix[i,j]:10.4f}" for j in range(6)])
            print(f"  [{row_str}]")

        print("\nRadia IMA:")
        for i in range(6):
            row_str = " ".join([f"{radia_ima_matrix[i,j]:10.4f}" for j in range(6)])
            print(f"  [{row_str}]")

        print("\nDifference:")
        for i in range(6):
            row_str = " ".join([f"{diff[i,j]:10.6f}" for j in range(6)])
            print(f"  [{row_str}]")

        # === Check if matrices match ===
        rel_diff = max_diff / max_elf
        if rel_diff < 0.0001:  # < 0.01%
            print(f"\n==> SUCCESS: Matrices match exactly! (rel diff < 0.01%)")
        elif rel_diff < 0.01:  # < 1%
            print(f"\n==> GOOD: Matrices match within 1%!")
        elif rel_diff < 0.1:  # < 10%
            print(f"\n==> WARNING: Matrices match within 10%")
        else:
            print(f"\n==> ERROR: Matrices differ significantly ({rel_diff*100:.1f}%)")

else:
    print("\nCould not build Radia IMA matrix")

print("\n" + "=" * 70)
