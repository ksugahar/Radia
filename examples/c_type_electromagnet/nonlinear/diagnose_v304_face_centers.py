#!/usr/bin/env python
"""
Diagnose V304 5% error: Compare face center computation methods.

Hypothesis: Radia uses 2D polygon centroid (area centroid) for FaceCenter,
while ELF uses simple geometric average of 4 face vertices.
For distorted V304 elements, these differ, affecting evaluation points
and thus the interaction matrix.

This script:
1. Loads V304 mesh elements
2. For each quad face, computes both face center methods
3. Quantifies the difference
4. Compares with a structured mesh (1x1x1) where they should be identical
"""

import sys
import os
import numpy as np

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(os.path.dirname(work_dir)))
sys.path.insert(0, os.path.join(repo_root, 'src'))
sys.path.insert(0, work_dir)

import radia as rad

# Mesh paths
ELF_V304 = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT\ELF_MMB8T_EIEM2_V304\MMB8T_only"
ELF_1x1x1 = r"S:\ELF_MAGIC\2020_03_07_CEFC_2020\model_C-Type\nonlinear_20000AT\ELF_MMB8T_EIEM2_1x1x1\full"

scale = 0.001  # mm -> m

# Radia/Netgen face ordering (0-indexed)
HEX_FACES_0idx = [
    [0, 3, 2, 1],  # Bottom (z-)
    [4, 5, 6, 7],  # Top (z+)
    [0, 1, 5, 4],  # Front (y-)
    [2, 3, 7, 6],  # Back (y+)
    [0, 4, 7, 3],  # Left (x-)
    [1, 2, 6, 5],  # Right (x+)
]


def load_mesh(path):
    """Load ELF mesh, return nodes and hex elements."""
    nodes = {}
    elements = []
    with open(os.path.join(path, "ELF_magic.meg"), 'r') as f:
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
    return nodes, elements


def geometric_face_center(verts, face_ids):
    """Simple geometric average of face vertices (ELF convention)."""
    pts = np.array([verts[i] for i in face_ids])
    return np.mean(pts, axis=0)


def area_centroid_2d(pts_3d):
    """Compute 2D area centroid of a quad face (Radia convention).

    Projects 3D quad to its best-fit 2D plane, computes area centroid
    using the shoelace formula, then transforms back to 3D.
    """
    # Compute face normal (cross product of diagonals)
    d1 = pts_3d[2] - pts_3d[0]
    d2 = pts_3d[3] - pts_3d[1]
    normal = np.cross(d1, d2)
    norm_len = np.linalg.norm(normal)
    if norm_len < 1e-20:
        return np.mean(pts_3d, axis=0)
    normal = normal / norm_len

    # Build local coordinate frame
    e1 = pts_3d[1] - pts_3d[0]
    e1_len = np.linalg.norm(e1)
    if e1_len < 1e-20:
        return np.mean(pts_3d, axis=0)
    e1 = e1 / e1_len
    e2 = np.cross(normal, e1)
    e2 = e2 / np.linalg.norm(e2)

    # Project to 2D
    origin = pts_3d[0]
    pts_2d = np.array([
        [np.dot(p - origin, e1), np.dot(p - origin, e2)]
        for p in pts_3d
    ])

    # Shoelace formula for centroid
    n = len(pts_2d)
    A = 0.0
    cx = 0.0
    cy = 0.0
    for i in range(n):
        j = (i + 1) % n
        xi, yi = pts_2d[i]
        xj, yj = pts_2d[j]
        cross = xi * yj - xj * yi
        A += cross
        cx += (xi + xj) * cross
        cy += (yi + yj) * cross
    A *= 0.5
    if abs(A) < 1e-30:
        return np.mean(pts_3d, axis=0)
    cx /= (6.0 * A)
    cy /= (6.0 * A)

    # Transform back to 3D
    center_3d = origin + cx * e1 + cy * e2
    return center_3d


def quad_planarity(pts_3d):
    """Compute planarity error of quad face.

    Returns max distance from best-fit plane, normalized by face size.
    """
    # Best-fit plane through centroid
    centroid = np.mean(pts_3d, axis=0)
    d1 = pts_3d[2] - pts_3d[0]
    d2 = pts_3d[3] - pts_3d[1]
    normal = np.cross(d1, d2)
    norm_len = np.linalg.norm(normal)
    if norm_len < 1e-20:
        return 0.0
    normal = normal / norm_len

    # Max distance of vertices from plane
    max_dist = 0.0
    for p in pts_3d:
        dist = abs(np.dot(p - centroid, normal))
        max_dist = max(max_dist, dist)

    # Normalize by face diagonal
    diag = max(np.linalg.norm(d1), np.linalg.norm(d2))
    if diag < 1e-20:
        return 0.0
    return max_dist / diag


def quad_aspect_ratio(pts_3d):
    """Compute aspect ratio of quad face."""
    edges = [
        np.linalg.norm(pts_3d[1] - pts_3d[0]),
        np.linalg.norm(pts_3d[2] - pts_3d[1]),
        np.linalg.norm(pts_3d[3] - pts_3d[2]),
        np.linalg.norm(pts_3d[0] - pts_3d[3]),
    ]
    if min(edges) < 1e-20:
        return float('inf')
    return max(edges) / min(edges)


def analyze_mesh(name, path):
    """Analyze face center differences for a mesh."""
    nodes, elements = load_mesh(path)

    print(f"\n{'='*70}")
    print(f"  {name}: {len(elements)} elements")
    print(f"{'='*70}")

    max_center_diff = 0.0
    max_eval_diff = 0.0
    max_planarity = 0.0
    max_aspect = 0.0
    worst_elem = -1

    all_center_diffs = []
    all_eval_diffs = []
    all_planarities = []

    for eid, nids in elements:
        verts = [nodes[nid] for nid in nids]

        # Element center = geometric average of 8 vertices
        elem_center = np.mean(verts, axis=0)

        elem_max_cdiff = 0.0
        elem_max_ediff = 0.0
        elem_max_plan = 0.0

        for fi, face_ids in enumerate(HEX_FACES_0idx):
            face_verts = np.array([verts[i] for i in face_ids])

            # Method 1: Simple geometric average (ELF)
            fc_geom = np.mean(face_verts, axis=0)

            # Method 2: Area centroid of projected 2D polygon (Radia)
            fc_area = area_centroid_2d(face_verts)

            # Difference
            cdiff = np.linalg.norm(fc_geom - fc_area)
            all_center_diffs.append(cdiff)

            # Evaluation point difference
            ep_geom = 0.5 * (fc_geom + elem_center)
            ep_area = 0.5 * (fc_area + elem_center)
            ediff = np.linalg.norm(ep_geom - ep_area)
            all_eval_diffs.append(ediff)

            # Planarity
            plan = quad_planarity(face_verts)
            all_planarities.append(plan)

            # Aspect ratio
            ar = quad_aspect_ratio(face_verts)

            elem_max_cdiff = max(elem_max_cdiff, cdiff)
            elem_max_ediff = max(elem_max_ediff, ediff)
            elem_max_plan = max(elem_max_plan, plan)
            max_aspect = max(max_aspect, ar)

        if elem_max_cdiff > max_center_diff:
            max_center_diff = elem_max_cdiff
            worst_elem = eid
        max_eval_diff = max(max_eval_diff, elem_max_ediff)
        max_planarity = max(max_planarity, elem_max_plan)

    all_center_diffs = np.array(all_center_diffs)
    all_eval_diffs = np.array(all_eval_diffs)
    all_planarities = np.array(all_planarities)

    # Find characteristic element size for normalization
    sizes = []
    for eid, nids in elements:
        verts = [nodes[nid] for nid in nids]
        max_dim = 0.0
        for i in range(8):
            for j in range(i+1, 8):
                d = np.linalg.norm(np.array(verts[i]) - np.array(verts[j]))
                max_dim = max(max_dim, d)
        sizes.append(max_dim)
    avg_size = np.mean(sizes)

    print(f"\n  Element size (avg diagonal): {avg_size*1000:.2f} mm")
    print(f"  Max face aspect ratio: {max_aspect:.2f}")
    print(f"\n  Face center difference (geometric avg vs area centroid):")
    print(f"    Max:  {max_center_diff*1000:.6f} mm ({max_center_diff/avg_size*100:.4f}% of elem size)")
    print(f"    Mean: {np.mean(all_center_diffs)*1000:.6f} mm")
    print(f"    Std:  {np.std(all_center_diffs)*1000:.6f} mm")
    n_nonzero = np.sum(all_center_diffs > 1e-12)
    print(f"    Faces with difference > 0: {n_nonzero}/{len(all_center_diffs)}")

    print(f"\n  Evaluation point difference:")
    print(f"    Max:  {max_eval_diff*1000:.6f} mm ({max_eval_diff/avg_size*100:.4f}% of elem size)")
    print(f"    Mean: {np.mean(all_eval_diffs)*1000:.6f} mm")

    print(f"\n  Face planarity (max vertex dist / diagonal):")
    print(f"    Max:  {max_planarity:.6f} ({max_planarity*100:.4f}%)")
    print(f"    Mean: {np.mean(all_planarities):.6f}")
    n_nonplanar = np.sum(all_planarities > 0.001)
    print(f"    Faces with >0.1% non-planarity: {n_nonplanar}/{len(all_planarities)}")

    print(f"\n  Worst element: #{worst_elem}")

    return all_center_diffs, all_eval_diffs, all_planarities


def test_single_element_field():
    """Test: compare field from single element with different face center methods."""
    print(f"\n{'='*70}")
    print(f"  Single Element Field Test")
    print(f"{'='*70}")

    nodes, elements = load_mesh(ELF_V304)

    # Find the element with worst face center difference
    from coil_model import create_racetrack_coil
    bh_data = np.loadtxt(os.path.join(work_dir, 'BH.txt'), comments='#').tolist()

    # Pick a few elements: first, last, and one near the gap
    test_elems = [elements[0], elements[len(elements)//2], elements[-1]]

    for eid, nids in test_elems:
        verts = [nodes[nid] for nid in nids]
        elem_center = np.mean(verts, axis=0)

        # Compute max face center difference
        max_diff = 0.0
        for fi, face_ids in enumerate(HEX_FACES_0idx):
            face_verts = np.array([verts[i] for i in face_ids])
            fc_geom = np.mean(face_verts, axis=0)
            fc_area = area_centroid_2d(face_verts)
            diff = np.linalg.norm(fc_geom - fc_area)
            max_diff = max(max_diff, diff)

        print(f"\n  Element #{eid}: center = [{elem_center[0]*1000:.2f}, "
              f"{elem_center[1]*1000:.2f}, {elem_center[2]*1000:.2f}] mm")
        print(f"  Max face center diff: {max_diff*1000:.6f} mm")


def compare_field_at_multiple_points():
    """Compare Radia V304 field at multiple observation points vs ELF."""
    print(f"\n{'='*70}")
    print(f"  Field Comparison at Multiple Points")
    print(f"{'='*70}")

    from coil_model import create_racetrack_coil
    bh_data = np.loadtxt(os.path.join(work_dir, 'BH.txt'), comments='#').tolist()
    nodes, elements, _ = load_mesh_full(ELF_V304)

    rad.UtiDelAll()
    rad.FldUnits('m')

    mat = rad.MatSatIsoTab(bh_data)
    all_objs = []
    for eid, nids in elements:
        v = [[nodes[nid][0], nodes[nid][1], nodes[nid][2]] for nid in nids]
        obj = rad.ObjHexahedron(v, [0, 0, 0])
        rad.MatApl(obj, mat)
        all_objs.append(obj)

    yoke = rad.ObjCnt(all_objs)
    coil = create_racetrack_coil(20000.0)
    model = rad.ObjCnt([yoke, coil])

    rad.Solve(model, 0.001, 100, 0, image='+x-z')

    # Test at multiple points along y-axis at origin x,z
    y_pts = [0.0, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20]
    print(f"\n  {'y (mm)':>8s} {'Bz (mT)':>10s}")
    print(f"  {'-'*20}")
    for y in y_pts:
        B = rad.Fld(model, 'b', [0, y, 0])
        print(f"  {y*1000:8.1f} {B[2]*1000:10.3f}")


def load_mesh_full(path):
    """Load mesh (same as load_mesh but compatible with original load_elf_geometry)."""
    nodes = {}
    elements = []
    wedges = []
    with open(os.path.join(path, "ELF_magic.meg"), 'r') as f:
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
    return nodes, elements, wedges


if __name__ == '__main__':
    print("V304 Face Center Diagnostic")
    print("Hypothesis: Radia uses area centroid, ELF uses geometric average")
    print("For distorted elements, these differ -> evaluation point mismatch")

    # Analyze both meshes
    cd_v304, ed_v304, pl_v304 = analyze_mesh("V304 (unstructured, 74 hex)", ELF_V304)
    cd_1x1x1, ed_1x1x1, pl_1x1x1 = analyze_mesh("1x1x1 (structured, 52 hex)", ELF_1x1x1)

    # Summary comparison
    print(f"\n{'='*70}")
    print(f"  Summary Comparison")
    print(f"{'='*70}")
    print(f"\n  {'Metric':>35s} {'V304':>15s} {'1x1x1':>15s}")
    print(f"  {'-'*65}")
    print(f"  {'Max face center diff (mm)':>35s} {np.max(cd_v304)*1000:15.6f} {np.max(cd_1x1x1)*1000:15.6f}")
    print(f"  {'Mean face center diff (mm)':>35s} {np.mean(cd_v304)*1000:15.6f} {np.mean(cd_1x1x1)*1000:15.6f}")
    print(f"  {'Max eval point diff (mm)':>35s} {np.max(ed_v304)*1000:15.6f} {np.max(ed_1x1x1)*1000:15.6f}")
    print(f"  {'Max non-planarity':>35s} {np.max(pl_v304):15.6f} {np.max(pl_1x1x1):15.6f}")

    print(f"\n  Conclusion:")
    if np.max(cd_v304) > 10 * np.max(cd_1x1x1):
        print(f"  -> V304 has significantly larger face center differences")
        print(f"     This could explain the 5% error vs ELF")
    else:
        print(f"  -> Face center differences are similar")
        print(f"     Root cause is likely elsewhere")
