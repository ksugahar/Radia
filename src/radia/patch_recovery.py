"""Region-separated linear patch recovery for affine P1 tetrahedra.

Zienkiewicz and Zhu (1992), DOI 10.1002/nme.1620330702. This is a
centroid-sampled SPR-style least-squares recovery, not a guarantee of
superconvergence on arbitrary unstructured meshes. Recovered fields are
postprocessing: they do not replace the discrete solution or its residual.
"""

import operator
import time

import numpy as np


def recover_vertex_patches(vertices, tetrahedra, samples, *, max_patch_rings=3):
    """Fit a linear vector polynomial per vertex within ONE material region.

    Samples are (ne, 3) gradients at the affine tetrahedron centroids. Patches
    with insufficient rank expand through vertex adjacency, never beyond the
    supplied region. Unresolved rank deficiency raises instead of silently
    substituting nodal averaging. Returns active vertex IDs, recovered values
    and diagnostics. The original arrays are not modified.
    """
    points = np.asarray(vertices, dtype=float)
    cells = np.asarray(tetrahedra)
    values = np.asarray(samples)
    if (points.ndim != 2 or points.shape[1] != 3
            or not np.isfinite(points).all()):
        raise ValueError("vertices must be finite (nv, 3) coordinates")
    if (cells.ndim != 2 or cells.shape[1] != 4 or not len(cells)
            or not np.issubdtype(cells.dtype, np.integer)
            or np.any(cells < 0) or np.any(cells >= len(points))):
        raise ValueError("tetrahedra must be a nonempty (ne, 4) vertex-index array")
    if (values.shape != (len(cells), 3) or not np.isfinite(values).all()):
        raise ValueError("samples must be finite (ne, 3) vectors")
    try:
        if isinstance(max_patch_rings, (bool, np.bool_)):
            raise TypeError
        max_patch_rings = operator.index(max_patch_rings)
    except TypeError as exc:
        raise ValueError("max_patch_rings must be a positive integer") from exc
    if max_patch_rings < 1:
        raise ValueError("max_patch_rings must be a positive integer")
    edges = points[cells[:, 1:]] - points[cells[:, :1]]
    volume6 = np.abs(np.linalg.det(edges))
    scale = np.max(np.linalg.norm(edges, axis=2), axis=1)
    if np.any(volume6 <= 1e-14 * scale**3):
        raise ValueError("degenerate tetrahedron in recovery region")
    centers = points[cells].mean(axis=1)
    active = np.unique(cells)
    incidence = {int(vertex): [] for vertex in active}
    for element, cell in enumerate(cells):
        for vertex in cell:
            incidence[int(vertex)].append(element)
    recovered = np.empty((len(active), 3), dtype=np.result_type(values, float))
    expanded = 0
    max_condition = 0.0
    for row, vertex in enumerate(active):
        patch = set(incidence[int(vertex)])
        for ring in range(int(max_patch_rings)):
            indices = np.array(sorted(patch), dtype=int)
            offsets = centers[indices] - points[vertex]
            length = np.max(np.linalg.norm(offsets, axis=1))
            design = np.column_stack((np.ones(len(indices)), offsets / length))
            coefficients, _, rank, singular = np.linalg.lstsq(
                design, values[indices], rcond=1e-12)
            if rank == 4:
                recovered[row] = coefficients[0]
                expanded += int(ring > 0)
                max_condition = max(max_condition, float(singular[0] / singular[-1]))
                break
            if ring + 1 == int(max_patch_rings):
                raise ValueError(
                    f"rank-deficient SPR patch at vertex {vertex}: rank={rank}, "
                    f"elements={len(indices)}, rings={max_patch_rings}")
            patch.update(element for node in np.unique(cells[indices])
                         for element in incidence[int(node)])
        else:  # pragma: no cover
            raise RuntimeError("patch recovery did not run")
    return {
        "vertex_ids": active, "values": recovered,
        "expanded_patches": expanded, "max_design_condition": max_condition,
        "sample_count": len(cells), "vertex_count": len(active),
    }


def recover_mixed_omega_p1(mesh, result, *, reduced_materials, total_materials,
                           max_patch_rings=3):
    """Recover potential gradients separately for every physical material.

    Requires a linear mixed-Omega result with phi_reduced and phi_total
    potentials on an affine, undeformed tetrahedral P1 mesh. On reduced regions H_s is added back *after* recovery, retaining
    the analytic source. Harmonic total-region sources are likewise retained.
    Interface vertices have distinct values on each material side. Kelvin
    exteriors and nonlinear material states are deliberately unsupported.
    Returned raw fields alias the original solution; recovered fields are
    distinct and are not asserted to satisfy flux continuity or div(B)=0.
    """
    import ngsolve as ng

    started = time.perf_counter()
    if (mesh.dim != 3 or mesh.GetCurveOrder() > 1
            or getattr(mesh, "deformation", None) is not None):
        raise ValueError("SPR requires an affine, undeformed 3D mesh")
    if result.get("nonlinear_stats") is not None or result.get("kelvin_materials"):
        raise ValueError("SPR currently supports linear physical regions only")
    if result["H_cf"].is_complex:
        raise ValueError("Mixed Omega SPR currently requires real fields")
    space = result["fes"]
    spaces = (result.get("fes_reduced", space), result.get("fes_total", space))
    if any(int(item.globalorder) != 1 for item in spaces):
        raise ValueError("SPR requires P1 potentials")
    reduced, total = set(reduced_materials), set(total_materials)
    if not reduced or not total or reduced & total or reduced | total != set(mesh.GetMaterials()):
        raise ValueError("reduced/total materials must partition the mesh")
    if not all(name in result for name in ("phi_reduced", "phi_total")):
        raise ValueError("SPR requires the solved phi_reduced and phi_total potentials")
    elements = list(mesh.Elements(ng.VOL))
    if any(element.type != ng.ET.TET for element in elements):
        raise ValueError("SPR currently supports tetrahedra only")
    vertices = np.array([vertex.point for vertex in mesh.vertices], dtype=float)
    gradients = {"total": ng.grad(result["phi_total"])}
    if "source_lift" in result:
        gradients["reduced"] = gradients["total"] - ng.grad(result["source_lift"])
    else:
        gradients["reduced"] = ng.grad(result["phi_reduced"])
    recovered_fields, patches, summaries = {}, {}, {}
    for material in sorted(reduced | total):
        cells = np.array([[v.nr for v in element.vertices]
                          for element in elements if element.mat == material], dtype=int)
        gradient = gradients["reduced" if material in reduced else "total"]
        # Preserve every analytic source contribution, including total-region
        # harmonic fields, without requiring optional solver metadata.
        source = result["H_cf"] + gradient
        centers = vertices[cells].mean(axis=1)
        samples = np.asarray([gradient(mesh(*point)) for point in centers])
        if np.iscomplexobj(samples):
            raise ValueError("Mixed Omega SPR currently requires real fields")
        patch = recover_vertex_patches(vertices, cells, samples,
                                       max_patch_rings=max_patch_rings)
        fes = ng.H1(mesh, order=1, definedon=mesh.Materials(material))
        gradient_gf = ng.GridFunction(fes ** 3, name=f"spr_gradient_{material}")
        for vertex, value in zip(patch["vertex_ids"], patch["values"]):
            dofs = fes.GetDofNrs(ng.NodeId(ng.VERTEX, int(vertex)))
            if len(dofs) != 1 or dofs[0] < 0:
                raise ValueError("SPR requires one active H1 DOF per region vertex")
            for component in range(3):
                gradient_gf.components[component].vec[dofs[0]] = float(value[component])
        recovered_fields[material] = source - gradient_gf
        patches[material] = gradient_gf
        summaries[material] = {key: value for key, value in patch.items()
                               if key not in ("vertex_ids", "values")}
    recovered_h = ng.CF(tuple(mesh.MaterialCF({
        name: field[component] for name, field in recovered_fields.items()
    }) for component in range(3)))
    return {
        "raw_H_cf": result["H_cf"], "raw_B_cf": result["B_cf"],
        "H_cf": recovered_h, "B_cf": result["mu_cf"] * recovered_h,
        "gradient_gridfunctions": patches,
        "diagnostics": {"method": "region-separated-centroid-linear-patch-recovery",
                        "regions": summaries, "seconds": time.perf_counter() - started,
                        "accuracy_accepted": False, "equilibrium_preserved": False},
    }
