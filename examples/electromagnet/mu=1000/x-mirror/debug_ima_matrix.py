#!/usr/bin/env python
"""
Debug IMA matrix differences between ELF and Radia.
Analyze specific blocks to identify the issue.
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
scale = 0.001  # mm to m


def load_elf_geometry(path):
    """Load ELF geometry from .meg file."""
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
    """Load ELF interaction matrix from .mat file (Fortran binary format)."""
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
print("IMA Matrix Debug Analysis")
print("=" * 70)

# Load ELF geometries
nodes_full, elements_full = load_elf_geometry(ELF_FULL)
nodes_xmirror, elements_xmirror = load_elf_geometry(ELF_X_MIRROR)

n_elem_full = len(elements_full)
n_elem_xmirror = len(elements_xmirror)
n_dof_full = n_elem_full * 6
n_dof_xmirror = n_elem_xmirror * 6

print(f"Full model: {n_elem_full} elements, {n_dof_full} DOF")
print(f"X-mirror model: {n_elem_xmirror} elements, {n_dof_xmirror} DOF")

# Load ELF matrices
elf_mat_full = load_elf_matrix(ELF_FULL, n_dof_full)
elf_mat_xmirror = load_elf_matrix(ELF_X_MIRROR, n_dof_xmirror)

# Create Radia models
rad.UtiDelAll()
rad.FldUnits('m')

# X-mirror model
hex_objects_xmirror = []
for elem_id, node_ids in elements_xmirror:
    verts = [[nodes_xmirror[nid][0] * scale, nodes_xmirror[nid][1] * scale, nodes_xmirror[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects_xmirror.append(hex_obj)

container_xmirror = rad.ObjCnt(hex_objects_xmirror)

# Build Radia x-mirror matrix
handle_xmirror = rad.BuildMatrix(container_xmirror, image='+x')
radia_mat_xmirror, _ = rad.GetInteractMatrix(handle_xmirror)

# Also build Radia full model matrix for comparison
rad.UtiDelAll()
rad.FldUnits('m')

hex_objects_full = []
for elem_id, node_ids in elements_full:
    verts = [[nodes_full[nid][0] * scale, nodes_full[nid][1] * scale, nodes_full[nid][2] * scale]
             for nid in node_ids]
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])
    mat = rad.MatLin(MU_R)
    rad.MatApl(hex_obj, mat)
    hex_objects_full.append(hex_obj)

container_full = rad.ObjCnt(hex_objects_full)

# Build Radia full model matrix
handle_full = rad.BuildMatrix(container_full)
radia_mat_full, _ = rad.GetInteractMatrix(handle_full)

print("\n" + "=" * 70)
print("Full Model Matrix Comparison")
print("=" * 70)

diff_full = np.abs(elf_mat_full - radia_mat_full)
print(f"Full model matrix: max diff = {diff_full.max():.6e}")
print(f"Full model matrix: mean diff = {diff_full.mean():.6e}")

print("\n" + "=" * 70)
print("X-Mirror Model Matrix Comparison")
print("=" * 70)

diff_xmirror = np.abs(elf_mat_xmirror - radia_mat_xmirror)
print(f"X-mirror model matrix: max diff = {diff_xmirror.max():.6e}")
print(f"X-mirror model matrix: mean diff = {diff_xmirror.mean():.6e}")

# Find worst differences
worst_idx = np.unravel_index(np.argmax(diff_xmirror), diff_xmirror.shape)
print(f"\nWorst difference at [{worst_idx[0]}, {worst_idx[1]}]:")
print(f"  ELF:   {elf_mat_xmirror[worst_idx[0], worst_idx[1]]:.10e}")
print(f"  Radia: {radia_mat_xmirror[worst_idx[0], worst_idx[1]]:.10e}")

# Analyze diagonal blocks
print("\n" + "=" * 70)
print("Diagonal Block Analysis (6x6 blocks)")
print("=" * 70)

for elem in range(min(5, n_elem_xmirror)):
    offset = elem * 6
    elf_block = elf_mat_xmirror[offset:offset+6, offset:offset+6]
    radia_block = radia_mat_xmirror[offset:offset+6, offset:offset+6]
    diff_block = np.abs(elf_block - radia_block)

    print(f"\nElement {elem} diagonal block:")
    print(f"  Max diff: {diff_block.max():.6e}")

    if diff_block.max() > 1e-6:
        worst_sub = np.unravel_index(np.argmax(diff_block), diff_block.shape)
        i, j = worst_sub
        print(f"  Worst at [{i}, {j}]: ELF={elf_block[i,j]:.6e}, Radia={radia_block[i,j]:.6e}")

# Analyze relationship between full model and x-mirror model
print("\n" + "=" * 70)
print("IMA Formula Verification")
print("=" * 70)
print("IMA formula: K_IMA[i,j] = K_full[i,j] + K_full[i, mirror(j)] @ P")
print("P is the DOF permutation matrix for x-mirror")

# For x-mirror, elements 0-25 in x-mirror correspond to elements 0-25 in full
# Their mirrors are elements 26-51 in full (for x>0 region)
# Check if this is correct by examining element centroids

print("\n--- Element Centroid Analysis ---")
for elem in range(min(3, n_elem_xmirror)):
    verts = [[nodes_xmirror[nid][0], nodes_xmirror[nid][1], nodes_xmirror[nid][2]]
             for nid in elements_xmirror[elem][1]]
    centroid = np.mean(verts, axis=0)
    print(f"X-mirror elem {elem}: centroid = ({centroid[0]:.1f}, {centroid[1]:.1f}, {centroid[2]:.1f}) mm")

print("\n--- Full Model Element Centroids ---")
for elem in [0, 1, 2, 26, 27, 28]:
    if elem < n_elem_full:
        verts = [[nodes_full[nid][0], nodes_full[nid][1], nodes_full[nid][2]]
                 for nid in elements_full[elem][1]]
        centroid = np.mean(verts, axis=0)
        print(f"Full elem {elem}: centroid = ({centroid[0]:.1f}, {centroid[1]:.1f}, {centroid[2]:.1f}) mm")

# Manually compute IMA matrix for element [0,0] block
print("\n" + "=" * 70)
print("Manual IMA Block Computation for Element 0 Self-Interaction")
print("=" * 70)

# Find mirror element index in full model
# For x-mirror: element at x should map to element at -x
elem_0_centroid = np.mean([[nodes_xmirror[nid][0] for nid in elements_xmirror[0][1]]]), \
                  np.mean([[nodes_xmirror[nid][1] for nid in elements_xmirror[0][1]]]), \
                  np.mean([[nodes_xmirror[nid][2] for nid in elements_xmirror[0][1]]])

verts_0 = np.array([[nodes_xmirror[nid][0], nodes_xmirror[nid][1], nodes_xmirror[nid][2]]
                    for nid in elements_xmirror[0][1]])
centroid_0 = np.mean(verts_0, axis=0)
mirror_centroid_0 = np.array([-centroid_0[0], centroid_0[1], centroid_0[2]])

print(f"Element 0 centroid: ({centroid_0[0]:.1f}, {centroid_0[1]:.1f}, {centroid_0[2]:.1f}) mm")
print(f"Mirror centroid:    ({mirror_centroid_0[0]:.1f}, {mirror_centroid_0[1]:.1f}, {mirror_centroid_0[2]:.1f}) mm")

# Find mirror element in full model
mirror_elem_idx = -1
for elem_idx, (elem_id, node_ids) in enumerate(elements_full):
    verts = [[nodes_full[nid][0], nodes_full[nid][1], nodes_full[nid][2]]
             for nid in node_ids]
    centroid = np.mean(verts, axis=0)
    dist = np.linalg.norm(centroid - mirror_centroid_0)
    if dist < 1.0:  # Within 1mm
        mirror_elem_idx = elem_idx
        print(f"Mirror element in full model: {mirror_elem_idx} (distance {dist:.3f} mm)")
        break

if mirror_elem_idx >= 0:
    # Get full model blocks
    K_ii = elf_mat_full[0:6, 0:6]  # Element 0 self-interaction
    K_i_mirror = elf_mat_full[0:6, mirror_elem_idx*6:(mirror_elem_idx+1)*6]  # Element 0 to mirror

    print("\nFull model K[0,0] block (self-interaction):")
    print(K_ii)

    print(f"\nFull model K[0,{mirror_elem_idx}] block (element 0 to mirror):")
    print(K_i_mirror)

    # DOF permutation for x-mirror: swap x+ (face 1) and x- (face 3)
    # Permutation: [0, 3, 2, 1, 4, 5] -> face index mapping
    perm = [0, 3, 2, 1, 4, 5]  # IMA_PERM_X

    # Apply permutation to columns
    K_i_mirror_perm = K_i_mirror[:, perm]

    print("\nK[0,mirror] after column permutation:")
    print(K_i_mirror_perm)

    # IMA formula: K_IMA = K_ii + K_i_mirror @ P
    K_IMA_computed = K_ii + K_i_mirror_perm

    print("\nComputed IMA block K_IMA[0,0] = K[0,0] + K[0,mirror] @ P:")
    print(K_IMA_computed)

    print("\nELF x-mirror K_IMA[0,0] block:")
    print(elf_mat_xmirror[0:6, 0:6])

    print("\nRadia x-mirror K_IMA[0,0] block:")
    print(radia_mat_xmirror[0:6, 0:6])

    print("\nDifference (ELF x-mirror - computed):")
    print(elf_mat_xmirror[0:6, 0:6] - K_IMA_computed)

    print("\nDifference (Radia x-mirror - computed):")
    print(radia_mat_xmirror[0:6, 0:6] - K_IMA_computed)

    # Try alternative IMA formulas
    print("\n" + "=" * 70)
    print("Alternative IMA Formulas")
    print("=" * 70)

    # Get K[mirror, 0] from full model
    K_mirror_i = elf_mat_full[mirror_elem_idx*6:(mirror_elem_idx+1)*6, 0:6]  # K[mirror(0), 0]
    K_mirror_mirror = elf_mat_full[mirror_elem_idx*6:(mirror_elem_idx+1)*6, mirror_elem_idx*6:(mirror_elem_idx+1)*6]

    print(f"\nFull model K[{mirror_elem_idx},0] block (mirror to element 0):")
    print(K_mirror_i)

    print(f"\nFull model K[{mirror_elem_idx},{mirror_elem_idx}] block (mirror self-interaction):")
    print(K_mirror_mirror)

    # Build permutation matrix P
    P = np.zeros((6, 6))
    for j in range(6):
        P[j, perm[j]] = 1.0

    print("\nPermutation matrix P:")
    print(P)

    # Try formula 2: K_IMA = K[i,j] + P^T @ K[mirror(i), j]
    K_IMA_formula2 = K_ii + P.T @ K_mirror_i
    print("\nFormula 2: K_IMA = K[i,j] + P^T @ K[mirror(i), j]:")
    print(K_IMA_formula2)
    print("Difference from ELF:", np.abs(elf_mat_xmirror[0:6, 0:6] - K_IMA_formula2).max())

    # Try formula 3: K_IMA = K[i,j] + P^T @ K[mirror(i), mirror(j)] @ P
    K_IMA_formula3 = K_ii + P.T @ K_mirror_mirror @ P
    print("\nFormula 3: K_IMA = K[i,j] + P^T @ K[mirror(i), mirror(j)] @ P:")
    print(K_IMA_formula3)
    print("Difference from ELF:", np.abs(elf_mat_xmirror[0:6, 0:6] - K_IMA_formula3).max())

    # Try formula 4: K_IMA = K[i,j] + K[i, mirror(j)] + K[mirror(i), j] + K[mirror(i), mirror(j)]
    # with permutations
    # K_IMA = K[i,j] + K[i,mirror(j)] @ P + P^T @ K[mirror(i),j] + P^T @ K[mirror(i),mirror(j)] @ P
    K_IMA_formula4 = K_ii + K_i_mirror_perm + P.T @ K_mirror_i + P.T @ K_mirror_mirror @ P
    print("\nFormula 4: K = K[i,j] + K[i,m(j)]@P + P^T@K[m(i),j] + P^T@K[m(i),m(j)]@P:")
    print(K_IMA_formula4)
    print("Difference from ELF:", np.abs(elf_mat_xmirror[0:6, 0:6] - K_IMA_formula4).max())

    # Try formula 5: The row DOF also needs permutation for the target face
    # K_IMA = K[i,j] + (P^T @ K[i, mirror(j)].T).T @ P = K[i,j] + K[i, mirror(j)] @ P @ P^T...?
    # Actually try: K_IMA = P^T @ (K[i,j] + K[i, mirror(j)] @ P)
    K_IMA_formula5 = P.T @ (K_ii + K_i_mirror_perm)
    print("\nFormula 5: K_IMA = P^T @ (K[i,j] + K[i, mirror(j)] @ P):")
    print(K_IMA_formula5)
    print("Difference from ELF:", np.abs(elf_mat_xmirror[0:6, 0:6] - K_IMA_formula5).max())

    # Try formula 6: What if we need to permute the output rows?
    # For IMA, the row index might need to be mapped through P as well
    # because the DOF ordering is different for the half model
    K_IMA_formula6 = K_ii + P.T @ K_i_mirror_perm
    print("\nFormula 6: K_IMA = K[i,j] + P^T @ K[i, mirror(j)] @ P:")
    print(K_IMA_formula6)
    print("Difference from ELF:", np.abs(elf_mat_xmirror[0:6, 0:6] - K_IMA_formula6).max())

    # Try a completely different approach: sum of two full model blocks
    # K_IMA[i,i] = K_full[i,i] + K_full[i,mirror(i)] + K_full[mirror(i),i] + K_full[mirror(i),mirror(i)]
    # with appropriate permutations
    K_IMA_formula7 = K_ii + K_i_mirror + K_mirror_i + K_mirror_mirror
    print("\nFormula 7: K = K[i,j] + K[i,m(j)] + K[m(i),j] + K[m(i),m(j)] (no permutation):")
    print(K_IMA_formula7)
    print("Difference from ELF:", np.abs(elf_mat_xmirror[0:6, 0:6] - K_IMA_formula7).max())

    # Formula 8: Try without column permutation but with row permutation
    K_IMA_formula8 = K_ii + P.T @ K_i_mirror
    print("\nFormula 8: K_IMA = K[i,j] + P^T @ K[i, mirror(j)] (row perm only):")
    print(K_IMA_formula8)
    print("Difference from ELF:", np.abs(elf_mat_xmirror[0:6, 0:6] - K_IMA_formula8).max())

    # Check specific element values to understand the pattern
    print("\n" + "=" * 70)
    print("Specific Element Analysis: K[1,1] and K[1,3]")
    print("=" * 70)

    print("\nELF x-mirror values:")
    print(f"  K_IMA[1,1] = {elf_mat_xmirror[1,1]:.6f}")
    print(f"  K_IMA[1,3] = {elf_mat_xmirror[1,3]:.6f}")
    print(f"  K_IMA[3,1] = {elf_mat_xmirror[3,1]:.6f}")
    print(f"  K_IMA[3,3] = {elf_mat_xmirror[3,3]:.6f}")

    print("\nRadia x-mirror values:")
    print(f"  K_IMA[1,1] = {radia_mat_xmirror[1,1]:.6f}")
    print(f"  K_IMA[1,3] = {radia_mat_xmirror[1,3]:.6f}")
    print(f"  K_IMA[3,1] = {radia_mat_xmirror[3,1]:.6f}")
    print(f"  K_IMA[3,3] = {radia_mat_xmirror[3,3]:.6f}")

    print("\nFull model components:")
    print(f"  K_AA[1,1] = {K_ii[1,1]:.6f}")
    print(f"  K_AA[1,3] = {K_ii[1,3]:.6f}")
    print(f"  K_AA[3,1] = {K_ii[3,1]:.6f}")
    print(f"  K_AA[3,3] = {K_ii[3,3]:.6f}")

    print(f"\n  K_AB[1,1] = K[0,39][1,1] = {K_i_mirror[1,1]:.6f}")
    print(f"  K_AB[1,3] = K[0,39][1,3] = {K_i_mirror[1,3]:.6f}")
    print(f"  K_BA[1,1] = K[39,0][1,1] = {K_mirror_i[1,1]:.6f}")
    print(f"  K_BA[1,3] = K[39,0][1,3] = {K_mirror_i[1,3]:.6f}")

    # Try to find the exact formula
    print("\n" + "=" * 70)
    print("Reverse Engineering the IMA Formula")
    print("=" * 70)

    # What is the required mirror contribution to match ELF?
    required = elf_mat_xmirror[0:6, 0:6] - K_ii
    print("\nRequired mirror contribution (ELF_IMA - K_AA):")
    print(required)

    # What does Radia compute?
    computed_mirror = K_i_mirror_perm
    print("\nRadia's mirror contribution (K_AB @ P):")
    print(computed_mirror)

    print("\nDifference (required - computed):")
    diff_mirror = required - computed_mirror
    print(diff_mirror)

    # Check if the difference equals P^T @ K_BA - K_AB @ P
    # This would mean: ELF uses K_AA + P^T @ K_BA instead of K_AA + K_AB @ P
    alt_contribution = P.T @ K_mirror_i
    print("\nAlternative contribution (P^T @ K_BA):")
    print(alt_contribution)

    print("\nDifference (required - P^T @ K_BA):")
    print(required - alt_contribution)

    # Maybe ELF uses average of both?
    avg_contribution = (K_i_mirror_perm + P.T @ K_mirror_i) / 2
    print("\nAverage contribution ((K_AB @ P + P^T @ K_BA) / 2):")
    print(avg_contribution)

    print("\nDifference (required - average):")
    print(required - avg_contribution)

    # Key insight: Try P^T @ K_BA @ P
    formula_new = P.T @ K_mirror_i @ P
    print("\n" + "=" * 70)
    print("NEW FORMULA: P^T @ K_BA @ P")
    print("=" * 70)
    print("\nComputed contribution (P^T @ K_BA @ P):")
    print(formula_new)

    print("\nDifference (required - P^T @ K_BA @ P):")
    diff_new = required - formula_new
    print(diff_new)
    print(f"\nMax difference: {np.abs(diff_new).max():.2e}")

    # Verify the full IMA block
    K_IMA_new_formula = K_ii + P.T @ K_mirror_i @ P
    print("\nFull K_IMA = K_AA + P^T @ K_BA @ P:")
    print(K_IMA_new_formula)

    print("\nDifference from ELF x-mirror:")
    diff_final = elf_mat_xmirror[0:6, 0:6] - K_IMA_new_formula
    print(diff_final)
    print(f"\n*** Max diff = {np.abs(diff_final).max():.2e} ***")

    # Compare all formulas
    print("\n" + "=" * 70)
    print("FORMULA COMPARISON SUMMARY")
    print("=" * 70)
    print(f"Original (K_AA + K_AB @ P):          max diff = {np.abs(elf_mat_xmirror[0:6, 0:6] - K_IMA_computed).max():.4e}")
    print(f"New (K_AA + P^T @ K_BA @ P):          max diff = {np.abs(diff_final).max():.4e}")

    # Try symmetrized formula: K_IMA = K_AA + (K_AB @ P + P^T @ K_BA @ P) / 2
    # No wait, this doesn't give right dimension. Let me try:
    # K_IMA = K_AA + (K_AB @ P + P^T @ K_BA) @ P
    # or just the average of K_AB and K_BA^T after permutation

    # Actually, by reciprocity in symmetric problems: K[i,j] = K[j,i]^T
    # So K_BA = K_AB^T, but with faces permuted
    # Let me check: K_BA ?= P @ K_AB^T @ P for mirror elements
    K_AB_T = K_i_mirror.T  # K[0,39]^T
    K_BA_from_reciprocity = P @ K_AB_T @ P

    print("\n" + "=" * 70)
    print("RECIPROCITY CHECK")
    print("=" * 70)
    print(f"\nK_BA (actual from full model):")
    print(K_mirror_i)
    print(f"\nK_BA from reciprocity (P @ K_AB^T @ P):")
    print(K_BA_from_reciprocity)
    print(f"\nDifference:")
    print(K_mirror_i - K_BA_from_reciprocity)
    print(f"Max diff: {np.abs(K_mirror_i - K_BA_from_reciprocity).max():.2e}")

    # If reciprocity holds, then:
    # K_IMA = K_AA + K_AB @ P should equal K_AA + P^T @ K_BA @ P
    # Let me verify...

    # Actually the key question is: what is the CORRECT IMA formula?
    # Let me try to match ELF exactly by finding what contribution is needed

    print("\n" + "=" * 70)
    print("EXACT ELF RECONSTRUCTION")
    print("=" * 70)

    # The required contribution to match ELF
    required_contribution = elf_mat_xmirror[0:6, 0:6] - K_ii

    # Try: K_AB @ P + correction
    correction = required_contribution - K_i_mirror_perm
    print("\nCorrection needed (Required - K_AB @ P):")
    print(correction)

    # Try: P^T @ K_BA @ P + correction
    correction2 = required_contribution - P.T @ K_mirror_i @ P
    print("\nCorrection needed (Required - P^T @ K_BA @ P):")
    print(correction2)

    # Check if correction is related to some other matrix
    # correction ≈ P^T @ K_BA @ P - K_AB @ P = P^T @ (K_BA - P @ K_AB) @ P ?
    diff_K = K_mirror_i - P @ K_i_mirror @ P.T  # K_BA - P @ K_AB @ P^T
    print("\nK_BA - P @ K_AB @ P^T:")
    print(diff_K)

    # Check pattern of correction
    print("\nPattern analysis of correction:")
    print(f"Row 0 sum: {correction[0,:].sum():.6f}")
    print(f"Row 1 sum: {correction[1,:].sum():.6f}")
    print(f"Row 2 sum: {correction[2,:].sum():.6f}")
    print(f"Row 3 sum: {correction[3,:].sum():.6f}")
    print(f"Row 4 sum: {correction[4,:].sum():.6f}")
    print(f"Row 5 sum: {correction[5,:].sum():.6f}")
