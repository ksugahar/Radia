"""Distributed loop work on undeformed, flat P1 surface triangles.

The caller owns TaskManager. Electric work is bilinear in phasor
coefficients; dissipated power uses the real impedance and conjugation.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix, csc_matrix, csr_matrix
from scipy.sparse.linalg import spsolve, LinearOperator

MU_0 = 4e-7 * np.pi
# Native Laplace Galerkin FMM compresses separated interactions; singular and
# near interactions use the same direct quadrature. These controls are fixed.
_FMM_OPTIONS = dict(fmm_minorder=20, fmm_maxdirect=100, fmm_separation=2,
                    fmm_eval_separation=3, fmm_maxlevel=20,
                    fmm_split_kr=5, fmm_order_factor=2)
_FMM_FACE_THRESHOLD = 512


def _loop_work_route(bem_solver, face_count):
    requested = getattr(bem_solver, "loop_work_backend", "dense")
    if requested not in ("auto", "dense", "fmm"):
        raise ValueError("loop_work_backend must be auto, dense, or fmm")
    return ("fmm" if face_count >= _FMM_FACE_THRESHOLD else "dense") if requested == "auto" else requested


class _P0SingleLayer(LinearOperator):
    """Face-ordered products retaining the native operator, without dense copies."""
    def __init__(self, operator, space, dofs, native_bytes, *, backend="dense"):
        self.operator, self.space = operator, space
        self.matrix, self.dofs = operator.mat, np.asarray(dofs)
        self.dense_array_bytes = 0  # no NumPy dense matrix is retained
        self.native_bytes = self.nbytes = None if native_bytes is None else int(native_bytes)
        self.backend = backend
        super().__init__(dtype=np.dtype(float), shape=(len(dofs), len(dofs)))

    def _component(self, values, transpose):
        x = self.matrix.CreateColVector()
        y = self.matrix.CreateRowVector()
        x.FV().NumPy()[self.dofs] = values
        try:
            (self.matrix.T if transpose else self.matrix).Mult(x, y)
        except Exception as error:
            if self.backend == "fmm":
                raise RuntimeError("Loop-work FMM apply failed; select loop_work_backend='dense' "
                                   "or --wp-loop-work-backend dense") from error
            raise
        result = y.FV().NumPy()[self.dofs].copy()
        if not np.all(np.isfinite(result)):
            raise RuntimeError("Nonfinite loop-work single-layer product")
        return result

    def _apply(self, values, transpose):
        values = np.asarray(values).reshape(-1)
        if np.iscomplexobj(values):
            return self._component(values.real, transpose)+1j*self._component(values.imag, transpose)
        return self._component(values, transpose)

    def _matvec(self, values):
        return self._apply(values, False)

    def _rmatvec(self, values):
        return self._apply(values, True)


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
    backend = _loop_work_route(bem_solver, len(areas))
    key = (identity, quadrature_bonus, backend)
    cached = getattr(bem_solver, "_loop_work_single_layer", None)
    if cached is not None and cached[0] == key:
        return cached[1]

    if len(areas) > 14000:
        raise ValueError("Loop work exceeds 14000 surface faces; use a coarser validated mesh")
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
    if backend == "fmm":
        try:
            operator = LaplaceSL(u*measure, use_fmm=True, **_FMM_OPTIONS)*v*measure
        except Exception as error:
            raise RuntimeError("Loop-work Galerkin FMM construction failed; select "
                               "loop_work_backend='dense' or --wp-loop-work-backend dense") from error
        # The native near+FMM operator exposes no complete storage inventory.
        native_bytes = None
    else:
        operator = LaplaceSL(u*measure, use_fmm=False)*v*measure
        values, columns, offsets = (np.asarray(a) for a in operator.mat.CSR())
        native_bytes = values.nbytes+columns.nbytes+offsets.nbytes
        native = csr_matrix((values, columns, offsets), shape=(space.ndof, space.ndof), copy=False)
        if not np.all(np.isfinite(values)) or np.any(native.diagonal() <= 0):
            raise RuntimeError("Invalid singular Galerkin P0 single-layer matrix")
    matrix = _P0SingleLayer(operator, space, dofs, native_bytes, backend=backend)
    bem_solver._loop_work_single_layer = (key, matrix)
    bem_solver._loop_work_assembly_seconds = time.perf_counter()-started
    return matrix


def _prepare_loop_magnetic_work(bem_solver, points, closed_triangles, open_triangles,
                   vertex_map, gradients, normals, areas, theta,
                   neumann_trace, *, quadrature_bonus=4):
    """Return electric/magnetic loop rows and their measured work diagnostics.

    Magnetic work includes the distributed Cartesian surface current and
    normal magnetic flux with the same exterior Calderon convention as the BIE.
    Seam/conormal diagnostics rebuild edge pairs and solve a surface dual
    problem per unprepared call. Strong prepared solves also retain geometry-only
    diagnostics and magnetic work; impedance-dependent work has a separate key.
    Its source pairs with the scattered field; electric work uses the total
    field. An independent reverse contraction checks work reciprocity.
    """
    n = len(points)
    nt = len(areas)
    currents = np.cross(normals[:, None, :], -gradients)
    loop_current = np.cross(normals, -np.einsum("tik,ti->tk", gradients,
                                               theta[open_triangles]))
    if np.iscomplexobj(loop_current) and np.any(loop_current.imag != 0):
        raise ValueError("Loop lift must be real")
    loop_current = np.asarray(loop_current.real)
    a, b = closed_triangles, np.roll(closed_triangles, -1, axis=1)
    edges = np.sort(np.stack((a, b), axis=-1).reshape(-1, 2), axis=1)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    ordered_edges = edges[order]
    starts = np.flatnonzero(np.r_[True, np.any(np.diff(ordered_edges, axis=0), axis=1), True])
    if np.any(np.diff(starts) != 2):
        raise ValueError("Loop work requires a closed manifold triangular surface")
    row_ids = np.repeat(np.arange(nt), 3)
    column_ids = vertex_map[open_triangles].ravel()
    maps = [coo_matrix((currents[:, :, c].ravel(), (row_ids, column_ids)),
                       shape=(nt, n)).tocsr() for c in range(3)]

    single_layer = _p0_single_layer(bem_solver, areas,
                                    quadrature_bonus=quadrature_bonus)
    applied = single_layer @ loop_current
    transposed = single_layer.T @ loop_current
    magnetic = MU_0*sum(mapping.T @ transposed[:, c] for c, mapping in enumerate(maps))
    reverse = MU_0*sum(mapping.T @ applied[:, c] for c, mapping in enumerate(maps))
    q = np.asarray(neumann_trace)
    if np.iscomplexobj(q) and np.any(q.imag != 0):
        raise ValueError("Unit carrier normal trace must be real")
    q = np.asarray(q.real, dtype=float)
    if q.shape != (n,) or not np.all(np.isfinite(q)):
        raise ValueError("Loop work requires finite matching normal and incident traces")
    B = .5*bem_solver.M-bem_solver.DL
    # Exterior work: Q_g,u = j_g^T S0 C - q_g^T (.5 M-DL),
    # Q_g,g = j_g^T S0 j_g + q_g^T SL q_g. Current-only work
    # omits finite-SIBC normal flux and fails a generic lift change.
    magnetic -= MU_0*(q @ B)
    reverse -= MU_0*(B.T @ q)
    magnetic_diagonal = MU_0*sum(transposed[:, c] @ loop_current[:, c] for c in range(3))
    magnetic_diagonal += MU_0*(q @ bem_solver.SL @ q)
    if not np.isfinite(magnetic_diagonal) or magnetic_diagonal <= 0:
        raise RuntimeError("Unit carrier exterior magnetic work must be positive")
    scale = max(np.max(np.abs(magnetic)), abs(magnetic_diagonal),
                np.finfo(float).tiny)
    symmetry = float(np.max(np.abs(magnetic-reverse))/scale)
    if not np.isfinite(symmetry) or symmetry > 1e-6:
        raise RuntimeError(f"Magnetic loop work reciprocity {symmetry:.3e} exceeds 1e-6; "
                           "refine quadrature")

    gradient_load = np.zeros(n)
    for k in range(3):
        np.add.at(gradient_load, vertex_map[open_triangles[:, k]],
                  areas*np.einsum("ij,ij->i", gradients[:, k], loop_current))
    dual = spsolve(csc_matrix(bem_solver.K[1:, 1:]), gradient_load[1:])
    norm = float(np.sum(areas*np.sum(loop_current**2, axis=1)))
    if not np.isfinite(norm) or norm <= 0 or not np.all(np.isfinite(dual)):
        raise RuntimeError("Invalid loop lift or surface-gradient dual solve")
    gradient_work = float(np.sqrt(max(0., gradient_load[1:] @ dual)/norm))
    flux = np.einsum("ti,tki->tk", loop_current,
                     np.cross(points[b]-points[a], normals[:, None, :])).ravel()[order]
    pairs = flux.reshape(-1, 2)
    flux_scale = max(np.max(np.sum(np.abs(pairs), axis=1)), np.finfo(float).tiny)
    flux_error = np.max(np.abs(np.sum(pairs, axis=1)))/flux_scale
    diagnostics = dict(magnetic_work_symmetry=symmetry,
                       gradient_work_residual=gradient_work,
                       conormal_flux_relative_mismatch=float(flux_error),
                       quadrature_bonus=quadrature_bonus,
                       single_layer_matrix_bytes=single_layer.nbytes,
                       single_layer_native_storage_bytes=single_layer.native_bytes,
                       single_layer_dense_array_bytes=single_layer.dense_array_bytes,
                       single_layer_dense_equivalent_bytes=8*nt**2,
                       single_layer_backend="ngsolve-galerkin-fmm" if single_layer.backend == "fmm" else "native-exact",
                       single_layer_requested_backend=getattr(bem_solver, "loop_work_backend", "dense"),
                       single_layer_auto_face_threshold=_FMM_FACE_THRESHOLD,
                       single_layer_fmm_parameters=_FMM_OPTIONS.copy() if single_layer.backend == "fmm" else None,
                       ngsolve_version=__import__("ngsolve").__version__,
                       single_layer_assembly_seconds=bem_solver._loop_work_assembly_seconds)
    return maps, loop_current, magnetic, magnetic_diagonal, diagnostics

def _prepare_loop_work(bem_solver, points, closed_triangles, open_triangles,
                   vertex_map, gradients, normals, areas, theta, impedance,
                   neumann_trace, *, quadrature_bonus=4, _geometry_key=None):
    """Static magnetic work and value-dependent electric work are cached separately."""
    nt = len(areas)
    impedance = np.asarray(impedance, dtype=complex)
    if impedance.ndim == 0:
        impedance = np.full(nt, impedance, dtype=complex)
    if (impedance.shape != (nt,) or not np.all(np.isfinite(impedance))
            or np.any(impedance.real < 0)):
        raise ValueError("Loop work requires finite passive face-ordered impedance")
    key = (_geometry_key, quadrature_bonus, _loop_work_route(bem_solver, len(areas)))
    cached = getattr(bem_solver, '_prepared_loop_magnetic_work', None) if _geometry_key is not None else None
    if cached is not None and cached[0] == key:
        magnetic_work = cached[1]
    else:
        magnetic_work = _prepare_loop_magnetic_work(bem_solver, points, closed_triangles,
            open_triangles, vertex_map, gradients, normals, areas, theta, neumann_trace,
            quadrature_bonus=quadrature_bonus)
        if _geometry_key is not None:
            bem_solver._prepared_loop_magnetic_work = (key, magnetic_work)
    maps, loop_current, magnetic, magnetic_diagonal, diagnostics = magnetic_work
    electric = sum(mapping.T @ (areas*impedance*loop_current[:, c])
                   for c, mapping in enumerate(maps))
    electric_diagonal = np.sum(areas*impedance*np.sum(loop_current**2, axis=1))
    return electric, electric_diagonal, magnetic, magnetic_diagonal, loop_current, diagnostics



def _loop_work_row(bem_solver, points, closed_triangles, open_triangles,
                   vertex_map, gradients, normals, areas, theta, impedance,
                   incident_vector_potential, neumann_trace, incident_potential,
                   *, quadrature_bonus=4, _prepared_key=None):
    """Reuse static pairings only; evaluate the current source on every call."""
    phi = np.asarray(incident_potential, dtype=complex)
    if phi.shape != (len(points),) or not np.all(np.isfinite(phi)):
        raise ValueError("Loop work requires a finite matching incident trace")
    cached = getattr(bem_solver, '_prepared_loop_work', None) if _prepared_key is not None else None
    key = (_prepared_key, quadrature_bonus, _loop_work_route(bem_solver, len(areas)))
    if cached is not None and cached[0] == key:
        work = cached[1]
    else:
        work = _prepare_loop_work(bem_solver, points, closed_triangles, open_triangles,
            vertex_map, gradients, normals, areas, theta, impedance, neumann_trace,
            quadrature_bonus=quadrature_bonus,
            _geometry_key=_prepared_key[0] if _prepared_key is not None else None)
        if _prepared_key is not None:
            bem_solver._prepared_loop_work = (key, work)
    electric, diagonal, magnetic, magnetic_diagonal, loop_current, diagnostics = work
    centers = points[closed_triangles].mean(axis=1)
    incident = np.asarray(incident_vector_potential(centers), dtype=complex)
    if incident.shape != (len(areas), 3) or not np.all(np.isfinite(incident)):
        raise ValueError("Incident vector potential must be finite on every face")
    forcing = np.sum(areas*np.einsum('td,td->t', incident, loop_current))-magnetic @ phi
    return electric, diagonal, magnetic, magnetic_diagonal, forcing, diagnostics.copy()
