"""Surface-mesh extraction utilities shared between BEM-based panels.

Moved out of the legacy BEM reference calc_heating_bem.py
on 2026-04-24 because the production BEM panel was importing from the
demoted examples path, which fails silently on any wheel-only install
(no repo checkout on the user's box).

Consumers:
  - src/radia/panels/calc_inductance.py         (production: PEEC|BEM-A
                                                  coil + scalar BEM-SIBC
                                                  workpiece, weak coupling)
  - validation_test/induction_heating/bem_reference/calc_heating_bem.py
    (research validation script; keeps its own local copy for back-compat)
"""
from __future__ import annotations


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

    # Add boundary elements (and their vertices) for the filtered set.
    old_to_new = {}
    for el in vol_mesh.Elements(BND):
        if keep_dom > 0:
            fd = vol_mesh.ngmesh.FaceDescriptor(el.index + 1)
            if fd.domin != keep_dom and fd.domout != keep_dom:
                continue
        lbl = bnd_labels[el.index]
        fd_idx = label_to_fd.get(lbl)
        if fd_idx is None:
            continue
        new_verts = []
        for v in el.vertices:
            if v.nr not in old_to_new:
                pt = vol_mesh.vertices[v.nr].point
                old_to_new[v.nr] = ngmesh_new.Add(
                    ngm.MeshPoint(ngm.Pnt(pt[0], pt[1], pt[2])))
            new_verts.append(old_to_new[v.nr])
        se = ngm.Element2D(fd_idx, new_verts)
        ngmesh_new.Add(se)

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
