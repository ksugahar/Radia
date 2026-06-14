"""cad_topology_edges.py -- IGA-style CAD-direct edge extraction (mesh-free).

Given an OCC shape (the CAD primitive used to seed the mesh, e.g.
`Box(...)`, `Compound([...])`, the loaded STEP file), extract the
TRUE topological edges + their dihedrals + lengths -- bypassing
the BBND mesh-segment artifact.

For a cuboid the CAD topology has 12 edges (one per geometric edge),
whereas a BBND-based extractor reports N_per_edge x 12 segments
(typically 72-200 depending on mesh density), some of which can
land on the FLAT face interior with dihedral approximated as
something other than pi (numerical artifact). The CAD-direct
extraction gives the exact 12 edges with dihedral pi/2 each.

For an L-section the CAD has 18 edges (12 outer + 6 from the
notch), including 2 re-entrant edges at dihedral 3pi/2.

This is the IGA-flavored part of Mixed Galerkin SIBC: bulk
Foster modes still need a volumetric FEM mesh (no way around
that for arbitrary geometry), but the SURFACE / EDGE / VERTEX
contributions can be IGA-direct from the CAD topology and stay
mesh-independent (Marussig 2014, Dolz 2019 isogeometric BEM
spirit).
"""
from __future__ import annotations

import math
from typing import Iterable

import numpy as np


def _vec(p):
    """Convert OCC Pnt or Vec to numpy array."""
    return np.array([p[0], p[1], p[2]], dtype=float)


def _edge_length(edge) -> float:
    """Edge length from OCC RefEdge.mass."""
    return float(edge.mass)


def _edge_midpoint(edge) -> np.ndarray:
    """Approximate midpoint by averaging start and end vertices.
    For a straight edge this is exact; for a curved edge it is a
    proxy good enough to bin edges by location.
    """
    verts = list(edge.vertices)
    if len(verts) >= 2:
        p0 = _vec(verts[0].p)
        p1 = _vec(verts[-1].p)
        return 0.5 * (p0 + p1)
    return _vec(verts[0].p) if verts else np.zeros(3)


def _face_normal_at(face, point) -> np.ndarray:
    """Approximate outward face normal at a given 3D point.

    Strategy:
      1. Use face center and 2 vertex offsets to build a tangent basis.
      2. Cross product gives the unit normal (sign chosen so it points
         away from the geometric center of the parent shape, set later).
    """
    verts = list(face.vertices)
    if len(verts) < 3:
        # Degenerate -- fall back to face.center based estimate
        c = _vec(face.center)
        return (c - point) / (np.linalg.norm(c - point) + 1e-30)
    p0 = _vec(verts[0].p)
    p1 = _vec(verts[1].p)
    p2 = _vec(verts[2].p)
    n = np.cross(p1 - p0, p2 - p0)
    nn = np.linalg.norm(n)
    if nn < 1e-30:
        # try other vertex pairs
        for i in range(1, len(verts) - 1):
            p2 = _vec(verts[i + 1].p)
            n = np.cross(p1 - p0, p2 - p0)
            nn = np.linalg.norm(n)
            if nn > 1e-30:
                break
    return n / (nn + 1e-30)


def _orient_outward(normals, face_centers, body_center) -> np.ndarray:
    """Flip each normal so it points AWAY from body_center."""
    out = []
    for n, c in zip(normals, face_centers):
        v = c - body_center
        if np.dot(n, v) < 0:
            out.append(-n)
        else:
            out.append(n)
    return np.array(out)


def cad_topology_edges(shape) -> list[dict]:
    """Return one dict per CAD edge with exact length + interior dihedral.

    Each dict: {'edge_idx', 'length', 'dihedral', 'midpoint', 'adjacent_faces'}

    Dihedral convention:
      - pi for flat (no actual edge)
      - pi/2 for right-angle (cuboid)
      - 3*pi/2 for re-entrant 270 deg corner (L-shape notch)
      - between 0 and 2*pi
    """
    # OCC reports each edge once per adjacent face -- deduplicate by midpoint.
    raw_edges = list(shape.edges)
    faces = list(shape.faces)
    if not raw_edges or not faces:
        return []

    seen_keys = {}
    edges = []
    for e in raw_edges:
        m = _edge_midpoint(e)
        key = (round(m[0]*1e9), round(m[1]*1e9), round(m[2]*1e9))
        if key in seen_keys:
            continue
        seen_keys[key] = True
        edges.append(e)

    # Body center (used to orient face normals outward)
    body_center = np.mean([_vec(f.center) for f in faces], axis=0)

    # For each face, store outward-oriented normal at face center
    face_centers = [_vec(f.center) for f in faces]
    face_normals_raw = [_face_normal_at(f, c) for f, c in zip(faces, face_centers)]
    face_normals = _orient_outward(face_normals_raw, face_centers, body_center)

    # For each edge, find adjacent faces by checking shared vertices.
    # OCC face.edges enumerates edges of each face -- safer to use that.
    # But the OCC python binding may not provide it directly. Use vertex containment.
    edge_vert_keys = []
    for e in edges:
        v_keys = []
        for v in e.vertices:
            p = v.p
            v_keys.append((round(p[0]*1e9), round(p[1]*1e9), round(p[2]*1e9)))
        edge_vert_keys.append(set(v_keys))

    face_vert_keys = []
    for f in faces:
        v_keys = []
        for v in f.vertices:
            p = v.p
            v_keys.append((round(p[0]*1e9), round(p[1]*1e9), round(p[2]*1e9)))
        face_vert_keys.append(set(v_keys))

    info = []
    for ei, e in enumerate(edges):
        L_e = _edge_length(e)
        mid = _edge_midpoint(e)

        # Adjacent faces: faces that contain BOTH endpoints of this edge
        e_keys = edge_vert_keys[ei]
        adj = []
        for fi, f in enumerate(faces):
            if e_keys.issubset(face_vert_keys[fi]):
                adj.append(fi)
        if len(adj) < 2:
            # Degenerate (open edge or single-face) -- skip
            info.append({
                "edge_idx": ei,
                "length": L_e,
                "dihedral": math.pi,
                "midpoint": mid.tolist(),
                "adjacent_faces": adj,
                "degenerate": True,
            })
            continue

        # Take the 2 face normals
        n1 = face_normals[adj[0]]
        n2 = face_normals[adj[1]]

        # Dihedral (interior of body):
        #   For convex right-angle corner: n1 . n2 = 0 -> dihedral = pi/2
        #   For flat (no edge): n1 = n2 -> dot = 1 -> dihedral = pi
        #   For re-entrant 270 deg: n1 . n2 = 0 but signed differently
        dot = float(np.clip(np.dot(n1, n2), -1.0, 1.0))
        # Distinguish convex (interior < pi) from re-entrant (interior > pi)
        # via the sign of (n1 + n2) . (mid - body_center).
        ext_dir = mid - body_center
        side = np.dot(n1 + n2, ext_dir)
        if side >= 0:
            dihedral = math.pi - math.acos(dot)         # convex
        else:
            dihedral = math.pi + math.acos(dot)         # re-entrant
        info.append({
            "edge_idx": ei,
            "length": L_e,
            "dihedral": dihedral,
            "midpoint": mid.tolist(),
            "adjacent_faces": adj,
            "degenerate": False,
        })
    return info


def cad_topology_total_area(shape) -> float:
    """Sum of face areas (CAD-exact)."""
    return float(sum(f.mass for f in shape.faces))


def cad_topology_c1(shape, mu) -> tuple[float, float, int]:
    """Polyhedral c_1 from CAD topology:
        c_1 = -(1/mu) sum_e L_e * W(alpha_e),  W(alpha) = (4/pi) cot(alpha/2)

    Returns (c_1, total_edge_length, n_edges).
    """
    edges_info = cad_topology_edges(shape)
    s = 0.0
    L_total = 0.0
    n = 0
    for e in edges_info:
        if e.get("degenerate", False):
            continue
        alpha = e["dihedral"]
        # Avoid singularity at alpha = 0 (knife edge, not physical for solid)
        if alpha < 1e-3 or alpha > 2 * math.pi - 1e-3:
            continue
        W = (4.0 / math.pi) * (1.0 / math.tan(alpha / 2.0))
        s += e["length"] * W
        L_total += e["length"]
        n += 1
    c1 = -s / mu
    return c1, L_total, n


def _self_test():
    """Compare mesh-derived vs CAD-derived c_1 for cube, cuboid, L-shape."""
    from netgen.occ import Box, Pnt, Glue, OCCGeometry, Vec, Prism, WorkPlane
    from ngsolve import Mesh

    MU_0 = 4 * math.pi * 1e-7

    print("=== CAD-topology edge extractor self-test ===\n")

    # Case 1: cube
    L = 5e-3
    shape = Box(Pnt(0, 0, 0), Pnt(L, L, L))
    edges_info = cad_topology_edges(shape)
    print(f"--- Cube L = {L*1e3:.1f} mm ---")
    print(f"  CAD edges: {len(edges_info)}")
    diheds = [e['dihedral'] for e in edges_info if not e.get('degenerate', False)]
    print(f"  dihedrals (deg): {[f'{math.degrees(d):.1f}' for d in diheds[:6]]}{'...' if len(diheds) > 6 else ''}")
    c1, Ltot, n = cad_topology_c1(shape, MU_0)
    expected_c1 = -16.0 * (3 * L) / (math.pi * MU_0)
    print(f"  CAD c_1 = {c1:.4e},  expected (cube formula) = {expected_c1:.4e},  ratio = {c1/expected_c1:.6f}")
    print()

    # Case 2: non-cubic cuboid 5x7x3 mm
    a, b, c = 5e-3, 7e-3, 3e-3
    shape = Box(Pnt(0, 0, 0), Pnt(a, b, c))
    edges_info = cad_topology_edges(shape)
    print(f"--- Cuboid {a*1e3:.0f}x{b*1e3:.0f}x{c*1e3:.0f} mm ---")
    print(f"  CAD edges: {len(edges_info)}")
    c1, Ltot, n = cad_topology_c1(shape, MU_0)
    expected_c1 = -16.0 * (a + b + c) / (math.pi * MU_0)
    print(f"  CAD c_1 = {c1:.4e},  expected = {expected_c1:.4e},  ratio = {c1/expected_c1:.6f}")
    print()

    # Case 3: sheared cuboid (parallelepiped) -- non-pi/2 dihedrals
    alpha_shr = 0.3
    bot = WorkPlane().Rectangle(L, L).Face()
    shape = Prism(bot, Vec(alpha_shr * L, 0, L))
    edges_info = cad_topology_edges(shape)
    print(f"--- Parallelepiped (shear={alpha_shr}) ---")
    print(f"  CAD edges: {len(edges_info)}")
    diheds = sorted([math.degrees(e['dihedral']) for e in edges_info if not e.get('degenerate', False)])
    print(f"  unique dihedrals (deg): {sorted(set(round(d, 1) for d in diheds))}")
    c1, Ltot, n = cad_topology_c1(shape, MU_0)
    print(f"  CAD c_1 = {c1:.4e},  total edge length = {Ltot*1e3:.4f} mm,  n = {n}")
    print()


if __name__ == "__main__":
    _self_test()
