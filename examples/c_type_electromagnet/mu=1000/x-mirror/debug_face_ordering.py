#!/usr/bin/env python
"""
Debug face ordering: Compare Radia vs ELF face conventions.

The issue might be that Radia and ELF use different face orderings for hexahedra.
"""

import sys
import os
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(work_dir)
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(parent_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))

import radia as rad

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


print("=" * 70)
print("Debug face ordering")
print("=" * 70)

# Load ELF geometry
nodes, elements = load_elf_geometry(ELF_X_MIRROR)

# Build Radia model and get face normals
rad.UtiDelAll()
rad.FldUnits('m')

elem_0 = elements[0]
print(f"\nElement 0 node IDs: {elem_0[1]}")
verts = [[nodes[nid][0] * scale, nodes[nid][1] * scale, nodes[nid][2] * scale]
         for nid in elem_0[1]]
print("Vertices (m):")
for i, v in enumerate(verts):
    print(f"  V{i}: {v}")

hex_obj = rad.ObjHexahedron(verts, [0, 0, 0])

# Get face info from Radia
face_info = rad.ObjGetFaces(hex_obj)
print(f"\nRadia face info: {len(face_info)} faces")

# For each face, compute normal direction
def compute_normal(v0, v1, v2):
    e1 = np.array(v1) - np.array(v0)
    e2 = np.array(v2) - np.array(v0)
    n = np.cross(e1, e2)
    norm = np.linalg.norm(n)
    if norm > 1e-12:
        n = n / norm
    return n


def classify_normal(n):
    """Classify normal direction."""
    tol = 0.9
    if n[0] > tol:
        return "x+"
    elif n[0] < -tol:
        return "x-"
    elif n[1] > tol:
        return "y+"
    elif n[1] < -tol:
        return "y-"
    elif n[2] > tol:
        return "z+"
    elif n[2] < -tol:
        return "z-"
    return "mixed"


print("\nRadia face analysis:")
for i, face in enumerate(face_info):
    # face is a list of vertex indices (1-indexed in Radia?)
    print(f"  Face {i}: vertex indices = {face}")
    if len(face) >= 3:
        # Get vertices (convert from 1-indexed if needed)
        face_verts = [verts[idx] if idx < len(verts) else verts[idx-1] for idx in face[:3]]
        n = compute_normal(*face_verts)
        direction = classify_normal(n)
        print(f"           normal = {n}, direction = {direction}")

# ELF face ordering (from CLAUDE.md): [y-, x+, y+, x-, z-, z+]
print("\nELF face ordering convention:")
print("  Face 0: y-")
print("  Face 1: x+")
print("  Face 2: y+")
print("  Face 3: x-")
print("  Face 4: z-")
print("  Face 5: z+")

# Compare orderings
print("\n" + "=" * 70)
print("Finding the permutation between Radia and ELF")
print("=" * 70)

# Map from direction to ELF face index
elf_direction_to_face = {
    'y-': 0, 'x+': 1, 'y+': 2, 'x-': 3, 'z-': 4, 'z+': 5
}

# Build Radia direction mapping
radia_face_directions = []
for i, face in enumerate(face_info):
    if len(face) >= 3:
        face_verts = [verts[idx] if idx < len(verts) else verts[idx-1] for idx in face[:3]]
        n = compute_normal(*face_verts)
        direction = classify_normal(n)
        radia_face_directions.append(direction)

print("\nRadia face directions:")
for i, d in enumerate(radia_face_directions):
    print(f"  Radia face {i}: {d}")

# Compute permutation
print("\nPermutation: Radia -> ELF")
radia_to_elf = []
for i, d in enumerate(radia_face_directions):
    elf_face = elf_direction_to_face.get(d, -1)
    radia_to_elf.append(elf_face)
    print(f"  Radia face {i} ({d}) = ELF face {elf_face}")

print(f"\nradia_to_elf permutation: {radia_to_elf}")

# Compute inverse permutation
elf_to_radia = [0] * 6
for radia_idx, elf_idx in enumerate(radia_to_elf):
    if 0 <= elf_idx < 6:
        elf_to_radia[elf_idx] = radia_idx

print(f"elf_to_radia permutation: {elf_to_radia}")

# What permutation should we use for x-mirror in Radia's ordering?
print("\n" + "=" * 70)
print("Deriving IMA permutation for Radia's face ordering")
print("=" * 70)

print("\nFor x-mirror, faces swap as follows (by direction):")
print("  x+ <-> x-")
print("  y-, y+, z-, z+ unchanged")

print("\nIn Radia's face ordering:")
# Find which Radia face is x+ and x-
x_plus_radia = radia_face_directions.index('x+') if 'x+' in radia_face_directions else -1
x_minus_radia = radia_face_directions.index('x-') if 'x-' in radia_face_directions else -1

print(f"  x+ is Radia face {x_plus_radia}")
print(f"  x- is Radia face {x_minus_radia}")

# Construct the Radia-specific Q permutation for x-mirror
# Q[i] tells which Radia face corresponds to face i after x-mirroring
radia_Q_x = list(range(6))
if x_plus_radia >= 0 and x_minus_radia >= 0:
    radia_Q_x[x_plus_radia] = x_minus_radia
    radia_Q_x[x_minus_radia] = x_plus_radia

print(f"\nRadia-specific Q for x-mirror: {radia_Q_x}")
print(f"(This swaps faces {x_plus_radia} and {x_minus_radia})")

# Compare with ELF's Q = [3, 2, 1, 0, 4, 5]
elf_Q_x = [3, 2, 1, 0, 4, 5]
print(f"\nELF Q for x-mirror: {elf_Q_x}")

# Convert ELF Q to Radia face ordering
# If K_elf is in ELF ordering and K_radia is in Radia ordering:
# K_elf = P @ K_radia @ P^T where P[i,j] = 1 if radia_to_elf[j] == i
# For row permutation: (P @ K_radia)[i,:] = K_radia[elf_to_radia[i], :]

# The ELF IMA formula: K_ima_elf = K_aa_elf + Q_elf @ K_ba_elf
# In Radia ordering: K_ima_radia = K_aa_radia + Q_radia @ K_ba_radia
# where Q_radia = P^(-1) @ Q_elf @ P in some sense

print("\n" + "=" * 70)
print("The correct permutation for Radia")
print("=" * 70)

# Actually, the simpler approach: what permutation does x-mirroring cause
# in Radia's face ordering?

# For x-mirror in Radia's ordering [based on face directions]:
# If Radia face i has direction D, and D' is the mirrored direction,
# then Q[i] = the Radia face with direction D'

mirror_direction = {
    'x+': 'x-', 'x-': 'x+',  # x-faces swap
    'y+': 'y+', 'y-': 'y-',  # y-faces unchanged
    'z+': 'z+', 'z-': 'z-',  # z-faces unchanged
}

radia_direction_to_face = {d: i for i, d in enumerate(radia_face_directions)}

correct_radia_Q = []
for i, d in enumerate(radia_face_directions):
    d_mirror = mirror_direction[d]
    j = radia_direction_to_face[d_mirror]
    correct_radia_Q.append(j)

print(f"Correct Radia Q for x-mirror (derived from face directions): {correct_radia_Q}")
