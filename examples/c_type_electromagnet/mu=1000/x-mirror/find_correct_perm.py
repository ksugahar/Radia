#!/usr/bin/env python
"""
Find the correct permutation formula by testing all combinations.
Compare K_AB @ P, P @ K_AB, P @ K_AB @ P, etc.
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

ELF_FULL = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\mu=1000\ELF_MMB8T_EIEM2_1x1x1\full"


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
print("Find correct permutation for IMA formula")
print("=" * 70)

# Load ELF full model
nodes_full, elements_full = load_elf_geometry(ELF_FULL)
n_elem_full = len(elements_full)
elf_mat_full = load_elf_matrix(ELF_FULL, n_elem_full * 6)

# Element 0 in x-mirror corresponds to element 0 in full model
# Its mirror is element 39 in full model
orig_idx = 0
mirror_idx = 39

# Extract blocks
K_AA = elf_mat_full[orig_idx*6:(orig_idx+1)*6, orig_idx*6:(orig_idx+1)*6]
K_AB = elf_mat_full[orig_idx*6:(orig_idx+1)*6, mirror_idx*6:(mirror_idx+1)*6]
K_BA = elf_mat_full[mirror_idx*6:(mirror_idx+1)*6, orig_idx*6:(orig_idx+1)*6]
K_BB = elf_mat_full[mirror_idx*6:(mirror_idx+1)*6, mirror_idx*6:(mirror_idx+1)*6]

# Required mirror contribution for IMA
# ELF_IMA[0,0] = K_AA + Required
# So Required = ELF_IMA[0,0] - K_AA
# We know from find_ima_formula.py that Required = Q @ K_BA where Q = [3,2,1,0,4,5]
Q = np.array([3, 2, 1, 0, 4, 5])
Required = K_BA[Q, :]  # This is Q @ K_BA

print("\nK_AA (self interaction):")
print(K_AA)

print("\nK_AB (field at orig from mirror source):")
print(K_AB)

print("\nK_BA (field at mirror from orig source):")
print(K_BA)

print("\nRequired = Q @ K_BA (expected mirror contribution):")
print(Required)

# Now find what formula using K_AB gives Required
# K_AB is what Compute6x6BlockMirrored computes (mirrors the source)
# Can we transform K_AB to get Required?

print("\n" + "=" * 70)
print("Testing formulas with K_AB")
print("=" * 70)

# Standard permutations to test
P_X = np.array([0, 3, 2, 1, 4, 5])  # Column swap 1 <-> 3
Q_X = np.array([3, 2, 1, 0, 4, 5])  # Row swap 0<->3, 1<->2

# Test 1: K_AB @ P (original Radia formula)
result1 = K_AB[:, P_X]
diff1 = np.abs(result1 - Required).max()
print(f"\n1. K_AB @ P_X (columns swapped): max diff = {diff1:.2e}")

# Test 2: Q @ K_AB (row swap only)
result2 = K_AB[Q_X, :]
diff2 = np.abs(result2 - Required).max()
print(f"2. Q @ K_AB (rows swapped): max diff = {diff2:.2e}")

# Test 3: Q @ K_AB @ P (both swaps)
result3 = K_AB[Q_X, :][:, P_X]
diff3 = np.abs(result3 - Required).max()
print(f"3. Q @ K_AB @ P (both swapped): max diff = {diff3:.2e}")

# Test 4: K_AB.T (transpose)
result4 = K_AB.T
diff4 = np.abs(result4 - Required).max()
print(f"4. K_AB.T: max diff = {diff4:.2e}")

# Test 5: Q @ K_AB.T
result5 = K_AB.T[Q_X, :]
diff5 = np.abs(result5 - Required).max()
print(f"5. Q @ K_AB.T: max diff = {diff5:.2e}")

# Test 6: K_AB.T @ P
result6 = K_AB.T[:, P_X]
diff6 = np.abs(result6 - Required).max()
print(f"6. K_AB.T @ P: max diff = {diff6:.2e}")

# Test 7: Q @ K_AB.T @ P
result7 = K_AB.T[Q_X, :][:, P_X]
diff7 = np.abs(result7 - Required).max()
print(f"7. Q @ K_AB.T @ P: max diff = {diff7:.2e}")

# Test 8: P @ K_AB
result8 = K_AB[P_X, :]
diff8 = np.abs(result8 - Required).max()
print(f"8. P @ K_AB (rows by P): max diff = {diff8:.2e}")

# Test 9: P @ K_AB @ P
result9 = K_AB[P_X, :][:, P_X]
diff9 = np.abs(result9 - Required).max()
print(f"9. P @ K_AB @ P: max diff = {diff9:.2e}")

# Test 10: P @ K_AB @ Q
result10 = K_AB[P_X, :][:, Q_X]
diff10 = np.abs(result10 - Required).max()
print(f"10. P @ K_AB @ Q: max diff = {diff10:.2e}")

# Find the best match
results = [
    ("K_AB @ P", diff1),
    ("Q @ K_AB", diff2),
    ("Q @ K_AB @ P", diff3),
    ("K_AB.T", diff4),
    ("Q @ K_AB.T", diff5),
    ("K_AB.T @ P", diff6),
    ("Q @ K_AB.T @ P", diff7),
    ("P @ K_AB", diff8),
    ("P @ K_AB @ P", diff9),
    ("P @ K_AB @ Q", diff10),
]

results.sort(key=lambda x: x[1])
print("\n" + "=" * 70)
print("Best matches:")
print("=" * 70)
for name, diff in results[:5]:
    print(f"  {name}: {diff:.2e}")

best_name, best_diff = results[0]
print(f"\n*** BEST: {best_name} with max diff = {best_diff:.2e} ***")

if best_diff < 1e-10:
    print("\n*** EXACT MATCH FOUND! ***")
    print(f"Formula: K_IMA = K_AA + {best_name}")
else:
    print("\nNo exact match found. Need different approach.")

    # Try exhaustive search over all possible row and column permutations
    print("\n" + "=" * 70)
    print("Exhaustive search (6! x 6! combinations)")
    print("=" * 70)

    best_diff = float('inf')
    best_row_perm = None
    best_col_perm = None

    # Only test permutations that keep z-faces (4,5) fixed
    base_perms = list(permutations([0,1,2,3]))

    for row_perm_base in base_perms:
        row_perm = list(row_perm_base) + [4, 5]
        for col_perm_base in base_perms:
            col_perm = list(col_perm_base) + [4, 5]
            result = K_AB[row_perm, :][:, col_perm]
            diff = np.abs(result - Required).max()
            if diff < best_diff:
                best_diff = diff
                best_row_perm = row_perm
                best_col_perm = col_perm

    print(f"Best row perm: {best_row_perm}")
    print(f"Best col perm: {best_col_perm}")
    print(f"Best diff: {best_diff:.2e}")

    if best_diff < 1e-10:
        print("\n*** EXACT MATCH with row/col permutation of K_AB ***")
        result = K_AB[best_row_perm, :][:, best_col_perm]
        print("\nResult:")
        print(result)
        print("\nRequired:")
        print(Required)
