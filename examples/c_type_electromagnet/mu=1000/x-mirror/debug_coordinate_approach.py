#!/usr/bin/env python
"""
Debug IMA using coordinate-based approach.

The user suggested: "For each DOF, use node coordinates directly to create mirrors
and generate matrix elements, without relying on element/connectivity information."

This script tests computing the mirror contribution directly by:
1. For each target face, get its evaluation point and normal
2. Mirror the evaluation point coordinates: x' = -x
3. Mirror the normal vector component: nx' = -nx
4. Compute the field from the original source at the mirrored observation point
5. Dot product with mirrored normal

This should give us Q @ K_BA directly without needing permutation arrays.
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

# ELF face ordering for hexahedra: [y-, x+, y+, x-, z-, z+]
# Face 0: y- (nodes 0,1,5,4)
# Face 1: x+ (nodes 1,2,6,5)
# Face 2: y+ (nodes 2,3,7,6)
# Face 3: x- (nodes 3,0,4,7)
# Face 4: z- (nodes 0,3,2,1)
# Face 5: z+ (nodes 4,5,6,7)
ELF_HEX_FACES = [
    [0, 1, 5, 4],  # Face 0: y-
    [1, 2, 6, 5],  # Face 1: x+
    [2, 3, 7, 6],  # Face 2: y+
    [3, 0, 4, 7],  # Face 3: x-
    [0, 3, 2, 1],  # Face 4: z-
    [4, 5, 6, 7],  # Face 5: z+
]


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


def compute_face_centroid(vertices, face_indices):
    """Compute face centroid from vertex coordinates."""
    face_verts = [vertices[i] for i in face_indices]
    return np.mean(face_verts, axis=0)


def compute_face_normal(vertices, face_indices):
    """Compute outward face normal from vertex coordinates."""
    v0 = vertices[face_indices[0]]
    v1 = vertices[face_indices[1]]
    v2 = vertices[face_indices[2]]

    # Cross product of two edge vectors
    e1 = v1 - v0
    e2 = v2 - v0
    n = np.cross(e1, e2)
    norm = np.linalg.norm(n)
    if norm > 1e-12:
        n = n / norm
    return n


def get_element_faces(nodes, node_ids):
    """Get face data for a hexahedral element."""
    # Convert to 0-indexed
    vertices = [nodes[nid] * scale for nid in node_ids]

    faces = []
    for face_idx, face_nodes in enumerate(ELF_HEX_FACES):
        centroid = compute_face_centroid(vertices, face_nodes)
        normal = compute_face_normal(vertices, face_nodes)
        faces.append({
            'centroid': centroid,
            'normal': normal,
            'vertices': [vertices[i] for i in face_nodes]
        })
    return faces


def mirror_point_x(p):
    """Mirror a point about x=0 plane."""
    return np.array([-p[0], p[1], p[2]])


def mirror_normal_x(n):
    """Mirror a normal vector about x=0 plane (flip x component)."""
    return np.array([-n[0], n[1], n[2]])


def solid_angle_triangle(obs, v0, v1, v2):
    """Compute solid angle of a triangle as seen from observation point."""
    r0 = v0 - obs
    r1 = v1 - obs
    r2 = v2 - obs

    d0 = np.linalg.norm(r0)
    d1 = np.linalg.norm(r1)
    d2 = np.linalg.norm(r2)

    if d0 < 1e-12 or d1 < 1e-12 or d2 < 1e-12:
        return 0.0

    # Van Oosterom & Strackee formula
    triple = np.dot(r0, np.cross(r1, r2))
    denom = d0*d1*d2 + d0*np.dot(r1, r2) + d1*np.dot(r0, r2) + d2*np.dot(r0, r1)

    if abs(denom) < 1e-15:
        return 0.0

    omega = 2.0 * np.arctan2(triple, denom)
    return omega


def field_from_charged_face(obs, face_vertices, sigma=1.0):
    """
    Compute H field at obs from a uniformly charged quadrilateral face.
    Uses solid angle method: H = -sigma * grad(Omega) / (4*pi)
    For unit sigma: H = n * Omega / (4*pi) where Omega is solid angle
    """
    # Split quad into two triangles
    v0, v1, v2, v3 = face_vertices

    omega1 = solid_angle_triangle(obs, v0, v1, v2)
    omega2 = solid_angle_triangle(obs, v0, v2, v3)
    omega_total = omega1 + omega2

    # Face normal
    e1 = np.array(v1) - np.array(v0)
    e2 = np.array(v2) - np.array(v0)
    n = np.cross(e1, e2)
    norm = np.linalg.norm(n)
    if norm > 1e-12:
        n = n / norm

    # H = -sigma * n * Omega / (4*pi)
    # For target DOF: K[target, source] = n_target . H_source
    # H_source = -grad(phi) where phi = sigma_source * Omega / (4*pi)
    # So K[i,j] = n_i . (-grad_at_i(phi_from_j))
    # With solid angle: K[i,j] = n_i . (n_face_j * Omega_j / (4*pi))
    # Wait, this is getting confusing. Let me use the gradient approach.

    inv_4pi = 1.0 / (4.0 * np.pi)

    # Actually for MSC, K[target_face, source_face] = n_target . H_at_target_from_source
    # And H = (1/4pi) * sigma * n_source_face * Omega
    # This is an approximation - the actual field integral is more complex

    # Let's use a simpler approach: H ~= (1/4pi) * sigma * sum(solid angles from triangles) * normal
    # This is the "constant panel" approximation
    H = inv_4pi * sigma * omega_total * n

    return H


print("=" * 70)
print("Debug IMA using coordinate-based approach")
print("=" * 70)

# Load ELF data
nodes_full, elements_full = load_elf_geometry(ELF_FULL)
nodes_xm, elements_xm = load_elf_geometry(ELF_X_MIRROR)

n_elem_full = len(elements_full)
n_elem_xm = len(elements_xm)
print(f"Full model: {n_elem_full} elements")
print(f"X-mirror model: {n_elem_xm} elements")

elf_mat_full = load_elf_matrix(ELF_FULL, n_elem_full * 6)
elf_mat_xm = load_elf_matrix(ELF_X_MIRROR, n_elem_xm * 6)

# Find element 0's correspondence
elem_0_xm = elements_xm[0]
faces_0 = get_element_faces(nodes_xm, elem_0_xm[1])

print(f"\nElement 0 in x-mirror model:")
print(f"  Node IDs: {elem_0_xm[1]}")
for i, face in enumerate(faces_0):
    print(f"  Face {i}: centroid={face['centroid']}, normal={face['normal']}")

# Find the corresponding elements in full model
def find_element_by_centroid(target_centroid, nodes, elements):
    for idx, (elem_id, node_ids) in enumerate(elements):
        verts = [nodes[nid] * scale for nid in node_ids]
        centroid = np.mean(verts, axis=0)
        if np.linalg.norm(centroid - target_centroid) < 0.001:  # 1mm tolerance
            return idx
    return -1

elem_centroid = np.mean([nodes_xm[nid] * scale for nid in elem_0_xm[1]], axis=0)
mirrored_centroid = mirror_point_x(elem_centroid)

orig_idx = find_element_by_centroid(elem_centroid, nodes_full, elements_full)
mirror_idx = find_element_by_centroid(mirrored_centroid, nodes_full, elements_full)

print(f"\nElement correspondence:")
print(f"  X-mirror elem 0 centroid: {elem_centroid}")
print(f"  Full model orig idx: {orig_idx}")
print(f"  Full model mirror idx: {mirror_idx}")

# Extract ELF blocks
K_AA = elf_mat_full[orig_idx*6:(orig_idx+1)*6, orig_idx*6:(orig_idx+1)*6]
K_BA = elf_mat_full[mirror_idx*6:(mirror_idx+1)*6, orig_idx*6:(orig_idx+1)*6]

print("\n" + "=" * 70)
print("Computing mirror contribution using coordinate-based approach")
print("=" * 70)

# For each target face of element 0, compute the contribution from the mirrored source
# K_mirror[target_face, source_face] = n'_target . H(at mirrored obs point, from original source)
# where n'_target has flipped x component

# Get source element faces
source_faces = faces_0  # same element for diagonal block

K_mirror_computed = np.zeros((6, 6))

for target_face_idx in range(6):
    target_face = faces_0[target_face_idx]

    # Mirror the observation point
    obs_mirrored = mirror_point_x(target_face['centroid'])

    # Mirror the target normal
    n_mirrored = mirror_normal_x(target_face['normal'])

    for source_face_idx in range(6):
        source_face = source_faces[source_face_idx]

        # Compute field at mirrored observation point from original source face
        H = field_from_charged_face(obs_mirrored, source_face['vertices'], sigma=1.0)

        # Dot with mirrored target normal
        K_mirror_computed[target_face_idx, source_face_idx] = np.dot(n_mirrored, H)

print("\nComputed K_mirror (coordinate-based):")
print(K_mirror_computed)

# ELF's expected mirror contribution (Q @ K_BA)
Q = [3, 2, 1, 0, 4, 5]
ELF_mirror = K_BA[Q, :]

print("\nELF mirror contribution (Q @ K_BA):")
print(ELF_mirror)

print(f"\nMax |computed - ELF|: {np.abs(K_mirror_computed - ELF_mirror).max():.2e}")

# Also check the simpler relationship
print("\n" + "=" * 70)
print("Checking face-by-face what the Q permutation does")
print("=" * 70)

# For x-mirror, the face mapping should be:
# Original face 0 (y-) with normal [0, -1, 0] -> after x-mirror, normal is still [0, -1, 0]
# But the evaluation point moves from (x, y, z) to (-x, y, z)
# So which face of the mirrored element corresponds to this?

print("\nFace normal directions before/after x-mirror:")
for i in range(6):
    orig_n = faces_0[i]['normal']
    mirr_n = mirror_normal_x(orig_n)
    print(f"Face {i}: original={orig_n}, mirrored={mirr_n}")

print("\nELF Q permutation [3, 2, 1, 0, 4, 5] means:")
print("  Row 0 <- Row 3 (x- face)")
print("  Row 1 <- Row 2 (y+ face)")
print("  Row 2 <- Row 1 (x+ face)")
print("  Row 3 <- Row 0 (y- face)")
print("  Row 4 <- Row 4 (z- face, unchanged)")
print("  Row 5 <- Row 5 (z+ face, unchanged)")

# Wait, the permutation swaps 0<->3 and 1<->2, but the face ordering is:
# [y-, x+, y+, x-, z-, z+]
# For x-mirror, x+ becomes x- and x- becomes x+, so faces 1 and 3 should swap
# But Q = [3, 2, 1, 0, 4, 5] swaps 0<->3 and 1<->2, which is:
# Face 0 (y-) <-> Face 3 (x-)
# Face 1 (x+) <-> Face 2 (y+)
# This doesn't match the expected x-face swap!

# Let me reconsider. The Q permutation is applied to K_BA rows.
# K_BA[row, col] = field at mirrored_elem.face[row] from orig_elem.face[col]
# Q @ K_BA means reordering the rows: (Q @ K_BA)[i,:] = K_BA[Q[i], :]
# So row 0 of result = row 3 of K_BA = field at mirrored_elem.face[3] from orig_elem

print("\n" + "=" * 70)
print("Analyzing what ELF actually computes")
print("=" * 70)

# The key insight: in ELF's IMA, the formula is:
# K_IMA[i,i] = K_AA + Q @ K_BA
# where K_BA = K[mirror(i), i] = field at mirrored target from original source
# and Q reorders the rows to match the original DOF ordering

# Let me verify: for the mirrored element, face 3 (x-) should have the same
# physical location as the original element's face 1 (x+) after mirroring
# Because the original x+ face (at positive x side) becomes x- face (at negative x side)
# when the geometry is mirrored about x=0

# So Q[0] = 3 means: DOF 0 of the IMA result corresponds to mirrored element's face 3
# which is the original element's face 1 (x+) after geometric mirroring

# Let me verify this by checking the face centroids
elem_full_orig = elements_full[orig_idx]
elem_full_mirror = elements_full[mirror_idx]

faces_orig = get_element_faces(nodes_full, elem_full_orig[1])
faces_mirror = get_element_faces(nodes_full, elem_full_mirror[1])

print("\nOriginal element face centroids:")
for i, face in enumerate(faces_orig):
    print(f"  Face {i}: {face['centroid']}")

print("\nMirrored element face centroids:")
for i, face in enumerate(faces_mirror):
    print(f"  Face {i}: {face['centroid']}")

print("\nChecking which faces correspond after x-mirroring:")
for i in range(6):
    orig_c = faces_orig[i]['centroid']
    mirror_of_orig = mirror_point_x(orig_c)
    # Find which face of the mirrored element matches
    for j in range(6):
        if np.linalg.norm(faces_mirror[j]['centroid'] - mirror_of_orig) < 0.001:
            print(f"  Original face {i} ({orig_c}) -> Mirrored face {j} ({faces_mirror[j]['centroid']})")
            break
