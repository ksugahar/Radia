#!/usr/bin/env python
"""
Debug script: Analyze A field computation per face
to understand why on-axis points give zero for ObjHexahedron.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'radia'))

import numpy as np

# Constants
MU_0 = 4 * np.pi * 1e-7
dx, dy, dz = 0.02, 0.02, 0.03  # half dimensions

# Magnet vertices (same as used in Radia)
vertices = np.array([
    [-dx, -dy, -dz],  # 0: bottom, front-left
    [ dx, -dy, -dz],  # 1: bottom, front-right
    [ dx,  dy, -dz],  # 2: bottom, back-right
    [-dx,  dy, -dz],  # 3: bottom, back-left
    [-dx, -dy,  dz],  # 4: top, front-left
    [ dx, -dy,  dz],  # 5: top, front-right
    [ dx,  dy,  dz],  # 6: top, back-right
    [-dx,  dy,  dz],  # 7: top, back-left
])

# Face definitions (counter-clockwise when viewed from outside)
# Each face has 4 vertices
HEX_FACES = [
    [0, 3, 2, 1],  # bottom face (z=-dz) - CCW from below
    [4, 5, 6, 7],  # top face (z=+dz) - CCW from above
    [0, 1, 5, 4],  # front face (y=-dy)
    [2, 3, 7, 6],  # back face (y=+dy)
    [0, 4, 7, 3],  # left face (x=-dx)
    [1, 2, 6, 5],  # right face (x=+dx)
]

def compute_face_normal(v0, v1, v2, v3, center):
    """Compute outward-pointing face normal."""
    e1 = v1 - v0
    e2 = v2 - v0
    n = np.cross(e1, e2)
    n = n / np.linalg.norm(n)

    # Ensure outward-pointing
    face_center = (v0 + v1 + v2 + v3) / 4
    if np.dot(n, face_center - center) < 0:
        n = -n
    return n

def compute_BufVect_triangle(V0, V1, V2, obs_point, elem_center):
    """
    Compute BufVect contribution from a single triangle face.
    BufVect = n * integral_face(1/|r-r'|) dS

    Uses the analytical formula for the integral:
    I = sum_edges[ L * ln((R_i + R_j + L) / (R_i + R_j - L)) ] + |z| * omega

    where omega is the solid angle subtended by the triangle.
    """
    EPS = 1e-15

    # Face edges
    e1 = V1 - V0
    e2 = V2 - V0

    # Face normal
    normal = np.cross(e1, e2)
    normal_len = np.linalg.norm(normal)
    if normal_len < EPS:
        return np.array([0., 0., 0.])
    normal = normal / normal_len

    # Check outward direction
    face_center = (V0 + V1 + V2) / 3
    if np.dot(normal, face_center - elem_center) < 0:
        normal = -normal

    # Build local coordinate system
    basis_a = e1 / np.linalg.norm(e1)
    basis_c = normal
    basis_b = np.cross(basis_c, basis_a)
    basis_b = basis_b / np.linalg.norm(basis_b)

    # Transform vertices to local 2D (V0 as origin)
    xy0 = np.array([0., 0.])
    xy1 = np.array([np.dot(e1, basis_a), np.dot(e1, basis_b)])
    e2_vec = V2 - V0
    xy2 = np.array([np.dot(e2_vec, basis_a), np.dot(e2_vec, basis_b)])

    # Transform observation point to local coordinates
    DD = obs_point - V0
    X = np.dot(DD, basis_a)
    Y = np.dot(DD, basis_b)
    Z = np.dot(DD, basis_c)  # Signed distance to face plane

    # Distance to each vertex in local coords
    x0, y0 = X - xy0[0], Y - xy0[1]
    x1, y1 = X - xy1[0], Y - xy1[1]
    x2, y2 = X - xy2[0], Y - xy2[1]

    R0 = np.sqrt(x0**2 + y0**2 + Z**2)
    R1 = np.sqrt(x1**2 + y1**2 + Z**2)
    R2 = np.sqrt(x2**2 + y2**2 + Z**2)

    # Edge lengths
    L01 = np.linalg.norm(xy1 - xy0)
    L12 = np.linalg.norm(xy2 - xy1)
    L20 = np.linalg.norm(xy0 - xy2)

    I_scalar = 0.0

    # Edge contributions (log terms)
    if L01 > EPS:
        Rplus = R0 + R1 + L01
        Rminus = R0 + R1 - L01
        if Rminus < EPS:
            Rminus = EPS
        I_scalar += L01 * np.log(Rplus / Rminus)

    if L12 > EPS:
        Rplus = R1 + R2 + L12
        Rminus = R1 + R2 - L12
        if Rminus < EPS:
            Rminus = EPS
        I_scalar += L12 * np.log(Rplus / Rminus)

    if L20 > EPS:
        Rplus = R2 + R0 + L20
        Rminus = R2 + R0 - L20
        if Rminus < EPS:
            Rminus = EPS
        I_scalar += L20 * np.log(Rplus / Rminus)

    # Solid angle contribution
    absZ = abs(Z)
    if absZ > EPS:
        # van Oosterom & Strackee solid angle formula
        r0 = np.array([x0, y0, Z])
        r1 = np.array([x1, y1, Z])
        r2 = np.array([x2, y2, Z])

        num = np.dot(r0, np.cross(r1, r2))
        denom = (R0 * R1 * R2
                + R0 * np.dot(r1, r2)
                + R1 * np.dot(r0, r2)
                + R2 * np.dot(r0, r1))

        if abs(denom) > EPS:
            omega = 2.0 * np.arctan2(num, denom)
            I_scalar += absZ * omega

    # BufVect = n * I_scalar
    BufVect = normal * I_scalar
    return BufVect

def compute_A_from_hexahedron(obs_point, magnetization):
    """Compute vector potential A from hexahedron using face integration."""
    INV_FOUR_PI = 1.0 / (4.0 * np.pi)

    center = np.array([0., 0., 0.])  # Element center

    BufVect_total = np.array([0., 0., 0.])

    print("\nPer-face BufVect contributions:")
    print("-" * 70)

    for face_idx, face in enumerate(HEX_FACES):
        V = [vertices[i] for i in face]

        # Split quad into 2 triangles: [0,1,2] and [0,2,3]
        BufVect_t1 = compute_BufVect_triangle(V[0], V[1], V[2], obs_point, center)
        BufVect_t2 = compute_BufVect_triangle(V[0], V[2], V[3], obs_point, center)
        BufVect_face = BufVect_t1 + BufVect_t2

        BufVect_total += BufVect_face

        print(f"Face {face_idx}: BufVect = [{BufVect_face[0]:12.6e}, {BufVect_face[1]:12.6e}, {BufVect_face[2]:12.6e}]")

    print("-" * 70)
    print(f"Total BufVect = [{BufVect_total[0]:12.6e}, {BufVect_total[1]:12.6e}, {BufVect_total[2]:12.6e}]")

    # A = (1/4pi) * M x BufVect
    M = np.array(magnetization)
    MxBuf = np.cross(M, BufVect_total)
    A = INV_FOUR_PI * MxBuf

    return A

# Test
print("=" * 70)
print("Debug: Per-face BufVect analysis for A field")
print("=" * 70)

Br = 1.2  # T
Mr = Br / MU_0  # A/m
magnetization = [0, 0, Mr]

test_points = [
    [0.05, 0.0, 0.0],   # Off-axis (should work)
    [0.0, 0.0, 0.05],   # On z-axis (problem case)
]

for pt in test_points:
    print(f"\n{'='*70}")
    print(f"Observation point: {pt}")
    print(f"{'='*70}")

    A = compute_A_from_hexahedron(np.array(pt), magnetization)
    print(f"\nComputed A = [{A[0]:12.6e}, {A[1]:12.6e}, {A[2]:12.6e}]")
    print(f"|A| = {np.linalg.norm(A):12.6e}")
