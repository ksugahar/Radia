#!/usr/bin/env python
"""
Diagnose V304: Check if ObjHexahedron face normals point outward.

Tests:
1. A standard unit cube - all normals should point outward
2. V304 elements - check all face normals
3. Compare Radia's field from single element vs analytical expectation
"""

import sys
import os
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(work_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, work_dir)

import radia as rad

scale = 0.001  # mm -> m

# ObjHexahedron face ordering (from radia_pybind.cpp, 1-indexed)
# Face 0: z- = {1,2,3,4}
# Face 1: x+ = {2,6,7,3}
# Face 2: y- = {1,5,6,2}
# Face 3: x- = {1,4,8,5}
# Face 4: y+ = {3,7,8,4}
# Face 5: z+ = {5,8,7,6}
FACES_0idx = [
    [0, 1, 2, 3],  # z-
    [1, 5, 6, 2],  # x+
    [0, 4, 5, 1],  # y-
    [0, 3, 7, 4],  # x-
    [2, 6, 7, 3],  # y+
    [4, 7, 6, 5],  # z+
]

FACE_NAMES = ['z-', 'x+', 'y-', 'x-', 'y+', 'z+']


def cross_product_normal(v0, v1, v2):
    """Get normal from triangle winding: (v1-v0) x (v2-v0)."""
    e1 = np.array(v1) - np.array(v0)
    e2 = np.array(v2) - np.array(v0)
    n = np.cross(e1, e2)
    norm = np.linalg.norm(n)
    if norm > 1e-20:
        n /= norm
    return n


def check_normal_outward(face_verts, elem_center):
    """Check if face normal (from cross product of first 3 vertices) points outward."""
    face_center = np.mean(face_verts, axis=0)
    normal = cross_product_normal(face_verts[0], face_verts[1], face_verts[2])
    to_face = face_center - elem_center
    dot = np.dot(normal, to_face)
    return dot > 0, normal, dot


def test_unit_cube():
    """Test with standard unit cube."""
    print("=" * 60)
    print("Test 1: Standard unit cube")
    print("=" * 60)

    # Standard hex8 vertices
    verts = np.array([
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],  # bottom
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],   # top
    ], dtype=float) * 0.01  # 1cm cube in meters

    elem_center = np.mean(verts, axis=0)
    print(f"  Element center: [{elem_center[0]*100:.1f}, {elem_center[1]*100:.1f}, {elem_center[2]*100:.1f}] cm")

    all_outward = True
    for fi, face_ids in enumerate(FACES_0idx):
        fv = np.array([verts[i] for i in face_ids])
        fc = np.mean(fv, axis=0)
        outward, normal, dot = check_normal_outward(fv, elem_center)
        status = "OUTWARD" if outward else "**INWARD**"
        print(f"  Face {fi} ({FACE_NAMES[fi]:>2s}): normal=[{normal[0]:+.3f},{normal[1]:+.3f},{normal[2]:+.3f}]  dot={dot:+.4f}  {status}")
        if not outward:
            all_outward = False

    print(f"\n  Result: {'ALL outward' if all_outward else 'SOME INWARD - PROBLEM!'}")

    # Now create in Radia and check field
    rad.UtiDelAll()
    rad.FldUnits('m')
    obj = rad.ObjHexahedron(verts.tolist(), [0, 0, 954930])  # PM in z
    B = rad.Fld(obj, 'b', [0, 0, 0.03])  # 3cm above
    print(f"  Bz at (0,0,3cm) = {B[2]*1000:.4f} mT")

    return all_outward


def test_v304_elements():
    """Test V304 mesh elements."""
    print("\n" + "=" * 60)
    print("Test 2: V304 mesh elements")
    print("=" * 60)

    ELF_V304 = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT\ELF_MMB8T_EIEM2_V304\MMB8T_only"
    nodes = {}
    elements = []
    with open(os.path.join(ELF_V304, "ELF_magic.meg"), 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('MGR1'):
                parts = line.split()
                nid = int(parts[1])
                x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                nodes[nid] = np.array([x, y, z]) * scale
            elif line.startswith('MMB8T'):
                parts = line.split()
                eid = int(parts[1])
                nids = [int(parts[i]) for i in range(4, 12)]
                elements.append((eid, nids))

    total_faces = 0
    inward_faces = 0
    worst_elems = []

    for eid, nids in elements:
        verts = [nodes[nid] for nid in nids]
        elem_center = np.mean(verts, axis=0)

        n_inward = 0
        for fi, face_ids in enumerate(FACES_0idx):
            fv = np.array([verts[i] for i in face_ids])
            outward, normal, dot = check_normal_outward(fv, elem_center)
            total_faces += 1
            if not outward:
                inward_faces += 1
                n_inward += 1

        if n_inward > 0:
            worst_elems.append((eid, n_inward))

    print(f"  Total elements: {len(elements)}")
    print(f"  Total faces: {total_faces}")
    print(f"  Inward-pointing normals: {inward_faces} / {total_faces}")
    print(f"  Elements with inward normals: {len(worst_elems)} / {len(elements)}")

    if worst_elems:
        print(f"\n  Elements with inward normals:")
        for eid, n_inv in sorted(worst_elems, key=lambda x: -x[1])[:10]:
            print(f"    Element #{eid}: {n_inv} inward faces")

        # Show details for worst element
        worst_eid = worst_elems[0][0]
        for eid, nids in elements:
            if eid == worst_eid:
                verts = [nodes[nid] for nid in nids]
                elem_center = np.mean(verts, axis=0)
                print(f"\n  Details for element #{worst_eid}:")
                print(f"  Center: [{elem_center[0]*1000:.2f}, {elem_center[1]*1000:.2f}, {elem_center[2]*1000:.2f}] mm")
                for fi, face_ids in enumerate(FACES_0idx):
                    fv = np.array([verts[i] for i in face_ids])
                    fc = np.mean(fv, axis=0)
                    outward, normal, dot = check_normal_outward(fv, elem_center)
                    status = "OUTWARD" if outward else "**INWARD**"
                    print(f"    Face {fi} ({FACE_NAMES[fi]:>2s}): "
                          f"center=[{fc[0]*1000:.2f},{fc[1]*1000:.2f},{fc[2]*1000:.2f}] mm "
                          f"normal=[{normal[0]:+.3f},{normal[1]:+.3f},{normal[2]:+.3f}] "
                          f"dot={dot:+.6f} {status}")
                break


def test_1x1x1_elements():
    """Test 1x1x1 mesh elements (should all be correct)."""
    print("\n" + "=" * 60)
    print("Test 3: 1x1x1 mesh elements")
    print("=" * 60)

    ELF_1x1x1 = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT\ELF_MMB8T_EIEM2_1x1x1\full"
    nodes = {}
    elements = []
    with open(os.path.join(ELF_1x1x1, "ELF_magic.meg"), 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('MGR1'):
                parts = line.split()
                nid = int(parts[1])
                x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                nodes[nid] = np.array([x, y, z]) * scale
            elif line.startswith('MMB8T'):
                parts = line.split()
                eid = int(parts[1])
                nids = [int(parts[i]) for i in range(4, 12)]
                elements.append((eid, nids))

    total_faces = 0
    inward_faces = 0

    for eid, nids in elements:
        verts = [nodes[nid] for nid in nids]
        elem_center = np.mean(verts, axis=0)

        for fi, face_ids in enumerate(FACES_0idx):
            fv = np.array([verts[i] for i in face_ids])
            outward, normal, dot = check_normal_outward(fv, elem_center)
            total_faces += 1
            if not outward:
                inward_faces += 1

    print(f"  Total elements: {len(elements)}")
    print(f"  Total faces: {total_faces}")
    print(f"  Inward-pointing normals: {inward_faces} / {total_faces}")
    print(f"  Result: {'ALL outward' if inward_faces == 0 else f'{inward_faces} INWARD FACES!'}")


def test_field_sign():
    """Test if the sign of field from a single hex element is correct."""
    print("\n" + "=" * 60)
    print("Test 4: Field sign consistency")
    print("=" * 60)

    rad.UtiDelAll()
    rad.FldUnits('m')

    # Create a simple cube centered at origin
    s = 0.01  # 1cm half-size
    verts = [
        [-s, -s, -s], [s, -s, -s], [s, s, -s], [-s, s, -s],
        [-s, -s, s], [s, -s, s], [s, s, s], [-s, s, s],
    ]

    # Test 1: Hex with z-magnetization (should match ObjRecMag)
    hex_obj = rad.ObjHexahedron(verts, [0, 0, 954930])
    rec_obj = rad.ObjRecMag([0, 0, 0], [2*s, 2*s, 2*s], [0, 0, 954930])

    # Compare fields at several points
    test_pts = [
        [0, 0, 0.03],   # Along z-axis
        [0.03, 0, 0],   # Along x-axis
        [0, 0.03, 0],   # Along y-axis
        [0.02, 0.02, 0.02],  # Diagonal
    ]

    print(f"\n  Comparison: ObjHexahedron vs ObjRecMag (z-magnetized cube)")
    print(f"  {'Point':>20s}  {'Hex Bz (mT)':>12s}  {'Rec Bz (mT)':>12s}  {'Diff %':>8s}")

    for pt in test_pts:
        B_hex = rad.Fld(hex_obj, 'b', pt)
        B_rec = rad.Fld(rec_obj, 'b', pt)
        diff = np.linalg.norm(np.array(B_hex) - np.array(B_rec))
        rel_diff = diff / np.linalg.norm(B_rec) * 100 if np.linalg.norm(B_rec) > 1e-20 else 0
        pt_str = f"[{pt[0]*100:.0f},{pt[1]*100:.0f},{pt[2]*100:.0f}]cm"
        print(f"  {pt_str:>20s}  {B_hex[2]*1000:12.4f}  {B_rec[2]*1000:12.4f}  {rel_diff:8.4f}")


if __name__ == '__main__':
    test_unit_cube()
    test_v304_elements()
    test_1x1x1_elements()
    test_field_sign()
