"""cad_edges.py -- IGA-style CAD-direct edge extraction (mesh-free).

Given an OCC shape (the CAD primitive used to seed the mesh, e.g.
`Box(...)`, `Compound([...])`, the loaded STEP file), extract the
TRUE topological edges + their dihedrals + lengths -- bypassing the
BBND mesh-segment artifact.

For a cuboid the CAD topology has 12 edges (one per geometric edge),
whereas a BBND-based extractor reports N_per_edge x 12 segments
(typically 72-200 depending on mesh density), some of which can
land on the FLAT face interior with dihedral approximated as
something other than pi (numerical artifact).  The CAD-direct
extraction gives the exact 12 edges with dihedral pi/2 each.

For an L-section the CAD has 18 edges (12 outer + 6 from the notch),
including 2 re-entrant edges at dihedral 3*pi/2.

This is the IGA-flavored part of Mixed Galerkin SIBC: bulk Foster
modes still need a volumetric FEM mesh (no way around that for
arbitrary geometry), but the SURFACE / EDGE / VERTEX contributions
can be IGA-direct from the CAD topology and stay mesh-independent
(Marussig 2014, Dolz 2019 isogeometric BEM spirit).
"""
from __future__ import annotations

import math
import numpy as np


def _vec(p):
    """Convert OCC Pnt or Vec to numpy array."""
    return np.array([p[0], p[1], p[2]], dtype=float)


def _edge_length(edge) -> float:
    return float(edge.mass)


def _edge_midpoint(edge) -> np.ndarray:
    verts = list(edge.vertices)
    if len(verts) >= 2:
        p0 = _vec(verts[0].p)
        p1 = _vec(verts[-1].p)
        return 0.5 * (p0 + p1)
    return _vec(verts[0].p) if verts else np.zeros(3)


def _edge_endpoints(edge):
    """Return (p0, p1) endpoint coordinates of a (straight) edge."""
    verts = list(edge.vertices)
    if len(verts) >= 2:
        return _vec(verts[0].p), _vec(verts[-1].p)
    p = _vec(verts[0].p) if verts else np.zeros(3)
    return p, p


def _face_normal_at(face, point) -> np.ndarray:
    verts = list(face.vertices)
    if len(verts) < 3:
        c = _vec(face.center)
        return (c - point) / (np.linalg.norm(c - point) + 1e-30)
    p0 = _vec(verts[0].p)
    p1 = _vec(verts[1].p)
    p2 = _vec(verts[2].p)
    n = np.cross(p1 - p0, p2 - p0)
    nn = np.linalg.norm(n)
    if nn < 1e-30:
        for i in range(1, len(verts) - 1):
            p2 = _vec(verts[i + 1].p)
            n = np.cross(p1 - p0, p2 - p0)
            nn = np.linalg.norm(n)
            if nn > 1e-30:
                break
    return n / (nn + 1e-30)


def _orient_outward(normals, face_centers, body_center) -> np.ndarray:
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

    Each dict: `{'edge_idx', 'length', 'dihedral', 'midpoint', 'adjacent_faces'}`

    Dihedral convention:
        pi for flat (no actual edge)
        pi/2 for right-angle (cuboid)
        3*pi/2 for re-entrant 270 deg corner (L-shape notch)
        between 0 and 2*pi
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

    body_center = np.mean([_vec(f.center) for f in faces], axis=0)
    face_centers = [_vec(f.center) for f in faces]
    face_normals_raw = [_face_normal_at(f, c) for f, c in zip(faces, face_centers)]
    face_normals = _orient_outward(face_normals_raw, face_centers, body_center)

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
        p0, p1 = _edge_endpoints(e)

        e_keys = edge_vert_keys[ei]
        adj = []
        for fi, f in enumerate(faces):
            if e_keys.issubset(face_vert_keys[fi]):
                adj.append(fi)
        if len(adj) < 2:
            info.append({
                "edge_idx": ei,
                "length": L_e,
                "dihedral": math.pi,
                "midpoint": mid.tolist(),
                "p0": p0.tolist(),
                "p1": p1.tolist(),
                "adjacent_faces": adj,
                "degenerate": True,
            })
            continue

        n1 = face_normals[adj[0]]
        n2 = face_normals[adj[1]]
        dot = float(np.clip(np.dot(n1, n2), -1.0, 1.0))
        ext_dir = mid - body_center
        side = np.dot(n1 + n2, ext_dir)
        if side >= 0:
            dihedral = math.pi - math.acos(dot)
        else:
            dihedral = math.pi + math.acos(dot)
        info.append({
            "edge_idx": ei,
            "length": L_e,
            "dihedral": dihedral,
            "midpoint": mid.tolist(),
            "p0": p0.tolist(),
            "p1": p1.tolist(),
            "adjacent_faces": adj,
            "degenerate": False,
        })
    return info


def cad_topology_faces(shape) -> list[dict]:
    """Return one dict per CAD face with exact area + outward normal + center.

    Each dict: `{'face_idx', 'area', 'normal', 'center'}` (normal oriented
    outward from the body centroid).  This is the per-face partition of the
    wetted surface used by the multi-port SIBC envelope: for a polyhedron the
    flat faces are disjoint, so the partition-of-unity weights theta_f are
    face indicators and the surface integral

        K_geom[p,q] = integral_{dOmega} f_p f_q dS

    decomposes as the plain per-face sum sum_F integral_F f_p f_q dS (the
    `alpha.surface_moment_matrix` boundary integral).  Curved bodies tiled
    with OVERLAPPING patches need genuine theta_f weights (Paper 1 SIII.B
    partition-of-unity SIBC); that smooth-surface extension is not yet
    verified here -- this function gives the verified flat-face partition.
    """
    faces = list(shape.faces)
    if not faces:
        return []
    centers = [_vec(f.center) for f in faces]
    body_center = np.mean(centers, axis=0)
    normals_raw = [_face_normal_at(f, c) for f, c in zip(faces, centers)]
    normals = _orient_outward(normals_raw, centers, body_center)
    out = []
    for fi, f in enumerate(faces):
        out.append({
            "face_idx": fi,
            "area": float(f.mass),
            "normal": normals[fi].tolist(),
            "center": centers[fi].tolist(),
        })
    return out


def edge_moment_matrix(shape, drive_fns, mu) -> np.ndarray:
    """Matrix (multi-port) edge coefficient C1[p,q] (CAD-direct, mesh-free).

        C1[p,q] = -(1/mu) sum_e W(alpha_e) * integral_e f_p(x) f_q(x) dl

    the matrix generalization of `cad_topology_c1` (the scalar c_1 Mellin
    edge coefficient).  `drive_fns` is the list of P port functions as
    callables R^3 -> R (numpy 3-vector in); they must be the SAME ports used
    for the bulk residue (`alpha.bulk_foster_matrix_via_eigen`) and the
    surface moment (`alpha.surface_moment_matrix`), expressed in the global
    frame.

    For the single monopole drive f=1 this returns [[c1]] identical to
    `cad_topology_c1(shape, mu)[0]`.  The straight-edge integral of the
    product of two affine port functions is integrated exactly with a 2-point
    Gauss rule.

    Returns a symmetric (P, P) ndarray.
    """
    edges_info = cad_topology_edges(shape)
    P = len(drive_fns)
    C1 = np.zeros((P, P))
    g = 1.0 / math.sqrt(3.0)            # 2-point Gauss node offset on [-1, 1]
    for e in edges_info:
        if e.get("degenerate", False):
            continue
        alpha = e["dihedral"]
        if alpha < 1e-3 or alpha > 2 * math.pi - 1e-3:
            continue
        W = (4.0 / math.pi) * (1.0 / math.tan(alpha / 2.0))
        L_e = e["length"]
        p0 = np.asarray(e["p0"], dtype=float)
        p1 = np.asarray(e["p1"], dtype=float)
        mid = 0.5 * (p0 + p1)
        half = 0.5 * (p1 - p0)
        nodes = [mid - g * half, mid + g * half]   # weights L_e/2 each
        fvals = np.array([[fn(xg) for fn in drive_fns] for xg in nodes])
        # integral_e f_p f_q dl = (L_e/2) sum_g f_p(x_g) f_q(x_g)
        moment = (L_e / 2.0) * (np.outer(fvals[0], fvals[0])
                                + np.outer(fvals[1], fvals[1]))
        C1 += -W * moment / mu
    return C1


def cad_topology_total_area(shape) -> float:
    """Sum of face areas (CAD-exact)."""
    return float(sum(f.mass for f in shape.faces))


def cad_topology_c1(shape, mu) -> tuple[float, float, int]:
    """Polyhedral c_1 from CAD topology.

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
        if alpha < 1e-3 or alpha > 2 * math.pi - 1e-3:
            continue
        W = (4.0 / math.pi) * (1.0 / math.tan(alpha / 2.0))
        s += e["length"] * W
        L_total += e["length"]
        n += 1
    c1 = -s / mu
    return c1, L_total, n
