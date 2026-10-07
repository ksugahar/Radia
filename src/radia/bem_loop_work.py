"""Distributed loop work on undeformed, flat P1 surface triangles.

The caller owns TaskManager. Electric work is bilinear in phasor
coefficients; dissipated power uses the real impedance and conjugation.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import spsolve

MU_0 = 4e-7 * np.pi


def _p0_single_layer(bem_solver, areas, *, quadrature_bonus=4):
    from ngsolve import BND, BilinearForm, SurfaceL2, ds
    from ngsolve.bem import LaplaceSL
    from .surface_impedance import panel_mesh_identity

    mesh = bem_solver.mesh
    if mesh.GetCurveOrder() > 1 or mesh.deformation is not None:
        raise ValueError("Loop work requires undeformed flat surface triangles")
    if not isinstance(quadrature_bonus, int) or quadrature_bonus < 0:
        raise ValueError("Loop work quadrature bonus must be a nonnegative integer")
    identity, _ = panel_mesh_identity(mesh)
    frozen = getattr(bem_solver, "_loop_work_geometry_identity", identity)
    if frozen != identity:
        raise ValueError("Loop work geometry changed; rebuild the BIE solver")
    bem_solver._loop_work_geometry_identity = identity
    key = (identity, quadrature_bonus)
    cached = getattr(bem_solver, "_loop_work_single_layer", None)
    if cached is not None and cached[0] == key:
        return cached[1]

    if len(areas) > 14000:
        raise ValueError("Dense loop work exceeds 14000 surface faces "
                         "(about 1.5 GiB for this matrix alone); use a coarser validated mesh")
    import time
    started = time.perf_counter()
    space = SurfaceL2(mesh, order=0, dual_mapping=False)
    elements = list(mesh.Elements(BND))
    dofs = [space.GetDofNrs(el)[0] for el in elements]
    if (len(dofs) != len(areas) or space.ndof != len(dofs)
            or len(set(dofs)) != len(dofs) or min(dofs) < 0):
        raise ValueError("Loop work requires one constant P0 DOF per triangle")
    u, v = space.TnT()
    mass = BilinearForm(space)
    mass += u*v*ds
    mass.Assemble()
    rows, cols, values = mass.mat.COO()
    diagonal = coo_matrix((values, (rows, cols)),
                          shape=(space.ndof, space.ndof)).tocsr().diagonal()[dofs]
    if not np.allclose(diagonal, areas, rtol=1e-10, atol=1e-16):
        raise ValueError("P0 basis scaling does not match the flat triangle areas")

    measure = ds(bonus_intorder=quadrature_bonus)
    operator = LaplaceSL(u*measure, use_fmm=False)*v*measure
    rows, cols, values = operator.mat.COO()
    matrix = coo_matrix((values, (rows, cols)),
                        shape=(space.ndof, space.ndof)).toarray()
    matrix = np.ascontiguousarray(matrix[np.ix_(dofs, dofs)])
    if not np.all(np.isfinite(matrix)) or np.any(matrix.diagonal() <= 0):
        raise RuntimeError("Invalid singular Galerkin P0 single-layer matrix")
    matrix.flags.writeable = False
    bem_solver._loop_work_single_layer = (key, matrix)
    bem_solver._loop_work_assembly_seconds = time.perf_counter()-started
    return matrix


def _loop_work_row(bem_solver, points, closed_triangles, open_triangles,
                   vertex_map, gradients, normals, areas, theta, impedance,
                   incident_vector_potential, neumann_trace, incident_potential,
                   *, quadrature_bonus=4):
    """Return electric/magnetic loop rows and their measured work diagnostics.

    Magnetic work includes the distributed Cartesian surface current and
    normal magnetic flux with the same exterior Calderon convention as the BIE.
    Seam/conormal diagnostics rebuild edge pairs and solve a surface dual
    problem per call; only the singular P0 operator is cached.
    Its source pairs with the scattered field; electric work uses the total
    field. An independent reverse contraction checks work reciprocity.
    """
    n = len(points)
    nt = len(areas)
    impedance = np.asarray(impedance, dtype=complex)
    if impedance.ndim == 0:
        impedance = np.full(nt, impedance, dtype=complex)
    if (impedance.shape != (nt,) or not np.all(np.isfinite(impedance))
            or np.any(impedance.real < 0)):
        raise ValueError("Loop work requires finite passive face-ordered impedance")
    currents = np.cross(normals[:, None, :], -gradients)
    loop_current = np.cross(normals, -np.einsum("tik,ti->tk", gradients,
                                               theta[open_triangles]))
    if np.iscomplexobj(loop_current) and np.any(loop_current.imag != 0):
        raise ValueError("Loop lift must be real")
    loop_current = np.asarray(loop_current.real)
    row_ids = np.repeat(np.arange(nt), 3)
    column_ids = vertex_map[open_triangles].ravel()
    maps = [coo_matrix((currents[:, :, c].ravel(), (row_ids, column_ids)),
                       shape=(nt, n)).tocsr() for c in range(3)]

    electric = sum(mapping.T @ (areas*impedance*loop_current[:, c])
                   for c, mapping in enumerate(maps))
    electric_diagonal = np.sum(areas*impedance*np.sum(loop_current**2, axis=1))
    single_layer = _p0_single_layer(bem_solver, areas,
                                    quadrature_bonus=quadrature_bonus)
    magnetic = MU_0*sum(mapping.T @ (single_layer.T @ loop_current[:, c])
                       for c, mapping in enumerate(maps))
    reverse = MU_0*sum(mapping.T @ (single_layer @ loop_current[:, c])
                      for c, mapping in enumerate(maps))
    q = np.asarray(neumann_trace)
    if np.iscomplexobj(q) and np.any(q.imag != 0):
        raise ValueError("Unit carrier normal trace must be real")
    q = np.asarray(q.real, dtype=float)
    incident_phi = np.asarray(incident_potential, dtype=complex)
    if (q.shape != (n,) or incident_phi.shape != (n,) or not np.all(np.isfinite(q))
            or not np.all(np.isfinite(incident_phi))):
        raise ValueError("Loop work requires finite matching normal and incident traces")
    B = .5*bem_solver.M-bem_solver.DL
    # Exterior work: Q_g,u = j_g^T S0 C - q_g^T (.5 M-DL),
    # Q_g,g = j_g^T S0 j_g + q_g^T SL q_g. Current-only work
    # omits finite-SIBC normal flux and fails a generic lift change.
    magnetic -= MU_0*(q @ B)
    reverse -= MU_0*(B.T @ q)
    magnetic_diagonal = MU_0*sum(loop_current[:, c] @ single_layer
                                 @ loop_current[:, c] for c in range(3))
    magnetic_diagonal += MU_0*(q @ bem_solver.SL @ q)
    if not np.isfinite(magnetic_diagonal) or magnetic_diagonal <= 0:
        raise RuntimeError("Unit carrier exterior magnetic work must be positive")
    scale = max(np.max(np.abs(magnetic)), abs(magnetic_diagonal),
                np.finfo(float).tiny)
    symmetry = float(np.max(np.abs(magnetic-reverse))/scale)
    if not np.isfinite(symmetry) or symmetry > 1e-6:
        raise RuntimeError(f"Magnetic loop work reciprocity {symmetry:.3e} exceeds 1e-6; "
                           "refine quadrature")

    centroids = points[closed_triangles].mean(axis=1)
    incident = np.asarray(incident_vector_potential(centroids), dtype=complex)
    if incident.shape != (nt, 3) or not np.all(np.isfinite(incident)):
        raise ValueError("Incident vector potential must be finite on every face")
    # Magnetic work acts on the scattered potential; electric work acts
    # on total potential. The caller multiplies this source by -i omega.
    forcing = (np.sum(areas*np.einsum("ij,ij->i", incident, loop_current))
               - magnetic @ incident_phi)

    gradient_load = np.zeros(n)
    for k in range(3):
        np.add.at(gradient_load, vertex_map[open_triangles[:, k]],
                  areas*np.einsum("ij,ij->i", gradients[:, k], loop_current))
    dual = spsolve(csc_matrix(bem_solver.K[1:, 1:]), gradient_load[1:])
    norm = float(np.sum(areas*np.sum(loop_current**2, axis=1)))
    if not np.isfinite(norm) or norm <= 0 or not np.all(np.isfinite(dual)):
        raise RuntimeError("Invalid loop lift or surface-gradient dual solve")
    gradient_work = float(np.sqrt(max(0., gradient_load[1:] @ dual)/norm))
    edge_flux = {}
    for ti, triangle in enumerate(closed_triangles):
        for a, b in zip(triangle, np.roll(triangle, -1)):
            edge = points[b]-points[a]
            flux = float(loop_current[ti] @ np.cross(edge, normals[ti]))
            edge_flux.setdefault(tuple(sorted((int(a), int(b)))), []).append(flux)
    if any(len(values) != 2 for values in edge_flux.values()):
        raise ValueError("Loop work requires a closed manifold triangular surface")
    flux_scale = max(max(sum(abs(v) for v in values)
                         for values in edge_flux.values()), np.finfo(float).tiny)
    flux_error = max(abs(sum(values)) for values in edge_flux.values())/flux_scale
    diagnostics = dict(magnetic_work_symmetry=symmetry,
                       gradient_work_residual=gradient_work,
                       conormal_flux_relative_mismatch=float(flux_error),
                       quadrature_bonus=quadrature_bonus,
                       single_layer_matrix_bytes=int(single_layer.nbytes),
                       single_layer_assembly_seconds=bem_solver._loop_work_assembly_seconds)
    return electric, electric_diagonal, magnetic, magnetic_diagonal, forcing, diagnostics
