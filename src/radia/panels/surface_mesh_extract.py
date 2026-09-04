"""Surface-mesh extraction utilities shared by BEM application backends.

Moved out of the legacy BEM reference calculator on 2026-04-24 because the
production BEM backend was importing from a validation copy, which fails on a
wheel-only install and allows the two implementations to drift.

Consumers:
  - src/radia/panels/calc_inductance.py         (production: PEEC|BEM-A
                                                  coil + scalar BEM-SIBC
                                                  workpiece, weak coupling)
"""
from __future__ import annotations


def orient_surface_triangles(points, tris):
    """Make a triangulated surface's winding globally consistent + outward.

    Two passes per connected component:
      1. face-BFS flip propagation: every interior edge must be listed in
         opposite directions by its two triangles (manifold orientation);
      2. outward normalisation: if the component's signed volume
         (divergence theorem, sum of p0.(p1 x p2)/6) is negative, flip
         the whole component.

    Surface extractors can return triangles whose winding is not globally
    consistent.  The
    old per-triangle "centroid-outward" heuristic in the hole extractor
    actively CREATES the inconsistency on a genus-1 tube (the bore-wall
    outward normal points toward the centroid, so the whole inner wall is
    flipped).  The inconsistent winding corrupts the double-layer
    operator of the scalar BIE + SIBC.  On already-consistent meshes this
    is a no-op.

    Args:
        points: (nv, 3) float array of vertex coordinates.
        tris: (nt, 3) int array of vertex indices (any winding).

    Returns:
        (tris_oriented, stats) -- a NEW (nt, 3) int64 array and a dict
        with n_flips / n_components / components_flipped /
        conflicts_before / conflicts_after.  For a closed manifold
        surface conflicts_after is always 0.
    """
    import numpy as np
    from collections import defaultdict, deque

    pts = np.asarray(points, dtype=float)
    tris = np.array(tris, dtype=np.int64, copy=True)
    nt = len(tris)

    def _directed_conflicts(tt):
        seen_dir = {}
        n = 0
        for t in tt:
            for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
                key = (int(a), int(b))
                if key in seen_dir:
                    n += 1
                seen_dir[key] = True
        return n

    conflicts_before = _directed_conflicts(tris)

    e2t = defaultdict(list)
    for ti, t in enumerate(tris):
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            e2t[(min(a, b), max(a, b))].append(ti)

    def _has_dir(t, a, b):
        return any(t[k] == a and t[(k + 1) % 3] == b for k in range(3))

    seen = set()
    n_flips = 0
    components = []
    for seed in range(nt):
        if seed in seen:
            continue
        comp = [seed]
        seen.add(seed)
        dq = deque([seed])
        while dq:
            ti = dq.popleft()
            t = tris[ti]
            for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
                for tj in e2t[(min(a, b), max(a, b))]:
                    if tj == ti or tj in seen:
                        continue
                    if _has_dir(tris[tj], a, b):
                        tris[tj][1], tris[tj][2] = tris[tj][2], tris[tj][1]
                        n_flips += 1
                    seen.add(tj)
                    comp.append(tj)
                    dq.append(tj)
        components.append(comp)

    comps_flipped = 0
    for comp in components:
        v6 = 0.0
        for ti in comp:
            t = tris[ti]
            v6 += float(np.dot(pts[t[0]], np.cross(pts[t[1]], pts[t[2]])))
        if v6 < 0.0:
            for ti in comp:
                tris[ti][1], tris[ti][2] = tris[ti][2], tris[ti][1]
            comps_flipped += 1

    stats = {
        "n_flips": int(n_flips),
        "n_components": len(components),
        "components_flipped": int(comps_flipped),
        "conflicts_before": int(conflicts_before),
        "conflicts_after": int(_directed_conflicts(tris)),
    }
    return tris, stats


def _extract_surface_mesh_filtered(vol_mesh, keep_label="",
                                    return_vertex_map=False):
    """Extract a clean 2D surface mesh containing only boundary elements
    that touch the specified material.

    Args:
        vol_mesh: NGSolve Mesh (3D volume mesh)
        keep_label: material name to keep (e.g. "coil").  Empty string
            means "no filter" (keep all boundary labels).
        return_vertex_map: if True, also return the ``new_to_old``
            mapping (extracted-mesh vertex nr -> volume-mesh vertex nr)
            so the caller can re-evaluate per-vertex fields on the
            parent (curved) volume mesh.

    Returns:
        NGSolve Mesh (surface only, no orphan vertices, boundary labels
        renumbered consecutively).  If ``return_vertex_map`` is True a
        ``(mesh, new_to_old)`` tuple is returned.

    Raises:
        ValueError if ``keep_label`` is non-empty and does not match any
        material in the mesh, OR if the matched material has no
        adjacent boundary elements.
    """
    from ngsolve import Mesh, BND
    import netgen.meshing as ngm

    keep_dom = 0  # 0 = no filter
    materials = vol_mesh.GetMaterials()
    if keep_label:
        for i, m in enumerate(materials, 1):
            if m == keep_label:
                keep_dom = i
                break
        if keep_dom == 0:
            raise ValueError(
                f"Surface mesh extractor: requested keep_label "
                f"{keep_label!r} is not in the .vol's materials list "
                f"({sorted(set(materials))}). Fix the .jou (block name) "
                f"or pass an existing material name.")

    ngmesh_new = ngm.Mesh(dim=3)
    bnd_labels = list(vol_mesh.GetBoundaries())

    # Pre-scan: which boundary labels have elements adjacent to keep_dom?
    used_labels = []
    label_to_fd = {}
    for el in vol_mesh.Elements(BND):
        if keep_dom > 0:
            fd = vol_mesh.ngmesh.FaceDescriptor(el.index + 1)
            if fd.domin != keep_dom and fd.domout != keep_dom:
                continue
        lbl = bnd_labels[el.index]
        if lbl not in label_to_fd:
            new_idx = len(used_labels) + 1
            fd_new = ngm.FaceDescriptor(bc=new_idx)
            fd_idx = ngmesh_new.Add(fd_new)
            ngmesh_new.SetBCName(new_idx - 1, lbl)
            label_to_fd[lbl] = fd_idx
            used_labels.append(lbl)

    if not used_labels:
        raise ValueError(
            f"Surface mesh extractor: material {keep_label!r} has no "
            f"adjacent boundary elements. The .vol export is "
            f"inconsistent — re-run the Cubit .jou and re-export.")

    # Pass 1: collect the filtered triangles (old vertex ids) + labels.
    kept = []                        # (fd_idx, [old vertex nrs])
    for el in vol_mesh.Elements(BND):
        if keep_dom > 0:
            fd = vol_mesh.ngmesh.FaceDescriptor(el.index + 1)
            if fd.domin != keep_dom and fd.domout != keep_dom:
                continue
        lbl = bnd_labels[el.index]
        fd_idx = label_to_fd.get(lbl)
        if fd_idx is None:
            continue
        kept.append((fd_idx, [v.nr for v in el.vertices]))

    # Pass 2: make the winding globally consistent + outward BEFORE
    # emitting (the BND winding of a .vol is not globally consistent, and
    # an inconsistent winding corrupts the double-layer operator -- see
    # orient_surface_triangles).
    import numpy as np
    tri_old = np.array([verts for _fd, verts in kept], dtype=np.int64)
    used = sorted({int(v) for t in tri_old for v in t})
    old_to_compact = {v: i for i, v in enumerate(used)}
    pts_used = np.array([vol_mesh.vertices[v].point for v in used])
    tri_compact = np.vectorize(old_to_compact.get)(tri_old)
    tri_oriented, _stats = orient_surface_triangles(pts_used, tri_compact)

    old_to_new = {}
    compact_to_old = {i: v for v, i in old_to_compact.items()}
    for (fd_idx, _verts), t_or in zip(kept, tri_oriented):
        new_verts = []
        for ci in t_or:
            old_nr = compact_to_old[int(ci)]
            if old_nr not in old_to_new:
                pt = vol_mesh.vertices[old_nr].point
                old_to_new[old_nr] = ngmesh_new.Add(
                    ngm.MeshPoint(ngm.Pnt(pt[0], pt[1], pt[2])))
            new_verts.append(old_to_new[old_nr])
        ngmesh_new.Add(ngm.Element2D(fd_idx, new_verts))

    surf_mesh = Mesh(ngmesh_new)
    if return_vertex_map:
        # ngmesh.Add returns a netgen PointId (1-indexed).  NGSolve Mesh
        # exposes vertices via 0-indexed ``v.nr``.  PointId is NOT
        # directly castable to int — use its ``.nr`` attribute (learned
        # the hard way 2026-04-12 when ``int(new_id)`` crashed BEM
        # with TypeError).
        new_to_old = {new_id.nr - 1: int(old_nr)
                       for old_nr, new_id in old_to_new.items()}
        return surf_mesh, new_to_old
    return surf_mesh
