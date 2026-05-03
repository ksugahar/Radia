"""Coil self-inductance via ngsolve.bem LaplaceSL on HDivSurface RT₀.

Production replacement for the intree (Python) BEM-A path retired
2026-05-03 after benchmarking showed ngsolve.bem was 50-60x faster
than our pure-Python double-loop assembler at every measured N (302
to 1014 triangles).  ngsolve.bem stores the dense Galerkin matrix
internally and exposes it via COO for O(1) extraction; intree's
``radia.bem.efie_rwg.solve_inductance_source_sink_intree`` was deleted.

Solves the constrained EFIE on the coil surface:
  [SL  D^T] [J]   [0]
  [D   0  ] [p] = [g]

where SL = LaplaceSL (single-layer BEM operator on surface currents),
D = divergence matrix (HDivSurface -> SurfaceL2), g = source/sink
current injection (+1/A_src at source, -1/A_snk at sink), J = surface
current coefficients, p = Lagrange multiplier enforcing current
conservation.

Inductance: ``L = mu_0 * J^T @ SL @ J``.
AC SIBC resistance: ``R = Re(Z_s) * J^T @ M @ J`` where M is the
HDivSurface mass matrix and Re(Z_s) = 1/(σ δ) is the real part of
the surface impedance.

Per-triangle J sampling: ``compute_centroids_areas_J(...)`` returns
the centroid + area + averaged J vector per BND triangle, ready to
feed into ``radia.bem_sibc_solver.compute_phi_inc_from_surface_J`` for
the workpiece weak-coupling bridge.

Requires: NGSolve >= 6.2.2603 with ngsolve.bem, scipy.
"""
from __future__ import annotations

import time

import numpy as np
from scipy.linalg import solve as scipy_solve
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import minres as scipy_minres
from scipy.sparse.linalg import gmres as scipy_gmres


MU_0 = 4e-7 * np.pi


def _solve_saddle(K, rhs, method, tol=1e-8, maxiter=2000):
    """Solve the saddle-point system K x = rhs.

    method:
      "lu"     dense LAPACK direct solve (scipy.linalg.solve, partial-
               pivot LU). Best for ndof <~ 5000. Memory O(N^2),
               work O(N^3).
      "minres" Krylov solver for symmetric indefinite systems
               (scipy.sparse.linalg.minres). Best when ndof is large
               and the matrix is symmetric (the saddle point IS
               symmetric: K^T = [[SL, D'],[D, 0]] = K).
      "gmres"  General iterative Krylov (scipy.sparse.linalg.gmres).
               Use only if MINRES stalls.
    """
    if method == "lu":
        x = scipy_solve(K, rhs)
        info = {"method": "lu", "iterations": 1,
                "residual": float(np.linalg.norm(K @ x - rhs))}
        return x, info
    if method == "minres":
        x, code = scipy_minres(K, rhs, rtol=tol, maxiter=maxiter)
        if code != 0:
            raise RuntimeError(
                f"MINRES did not converge (info={code}, "
                f"maxiter={maxiter}, tol={tol})")
        info = {"method": "minres", "iterations": maxiter,
                "residual": float(np.linalg.norm(K @ x - rhs))}
        return x, info
    if method == "gmres":
        x, code = scipy_gmres(K, rhs, rtol=tol, maxiter=maxiter,
                              restart=50)
        if code != 0:
            raise RuntimeError(
                f"GMRES did not converge (info={code}, "
                f"maxiter={maxiter}, tol={tol})")
        info = {"method": "gmres", "iterations": maxiter,
                "residual": float(np.linalg.norm(K @ x - rhs))}
        return x, info
    raise ValueError(
        f"Unknown saddle-point solver: {method!r}.  "
        f"Choices: lu, minres, gmres.")


def _to_dense(mat):
    """Extract dense NumPy array from NGSolve BaseMatrix via COO.

    NGSolve BEM operators (LaplaceSL etc.) store the full dense matrix
    internally as a SparseMatrix with 100% fill.  The built-in
    ToDense() is O(N) slower than necessary (~144 s vs 0.06 s at
    N=5085) because it performs N column-by-column MatVecs instead
    of a direct memory copy.

    This function extracts the COO triplets and converts via scipy,
    which is ~2500x faster.  Note: requires the operator to have
    been constructed with ``use_fmm=False`` so the underlying
    storage is dense (FMM = lazy operator, no COO).
    """
    rows, cols, vals = mat.COO()
    return coo_matrix((vals, (rows, cols)),
                      shape=(mat.height, mat.width)).toarray()


def compute_inductance_source_sink(
        mesh, source_label="source", sink_label="sink",
        fes_order=0, solver="lu", Z_s_re=0.0):
    """Compute self-inductance + (optional) AC SIBC R via ngsolve.bem saddle EFIE.

    Args:
        mesh: NGSolve Mesh (volume mesh with boundary, or surface-only).
        source_label, sink_label: BND labels for current injection /
            extraction faces.  Set in the OCC face.name = "source"/"sink"
            BEFORE meshing, or rely on the panel's smallest-2-PLANE
            auto-detection.
        fes_order: HDivSurface polynomial order.  0 = RT₀ (= RWG, the
            production setting); 1+ for higher-order is supported by
            NGSolve but not validated against radia goldens here.
        solver: saddle-point solver — "lu" (dense direct, default),
            "minres" (symmetric Krylov), or "gmres".
        Z_s_re: real part of surface impedance Z_s = 1/(σ δ) for
            AC SIBC closure on the conductor.  Pass 0.0 to skip
            (vacuum L only).  Used to compute R = Z_s_re * J^T M J.

    Returns:
        dict with keys
        ``L``         : float [H], self-inductance for unit terminal current
        ``R``         : float [Ω], AC SIBC dissipation (0.0 if Z_s_re=0)
        ``J``         : (n_J,) HDivSurface coefficients (unit terminal current)
        ``gf_J``      : NGSolve GridFunction(HDivSurface) wrapping J
        ``SL``, ``D`` : dense Galerkin matrices (post-processing)
        ``A_source``, ``A_sink`` : port areas [m²]
        ``residual``  : max|D J - g|, machine precision for LU solve
        ``n_J``, ``n_f`` : DOF counts
        ``t_assembly``, ``t_solve``, ``t_total`` : timings [s]
    """
    from ngsolve import (HDivSurface, SurfaceL2, TaskManager, ds, BND,
                         BilinearForm, LinearForm, div, GridFunction)
    from ngsolve.bem import LaplaceSL

    t_start = time.perf_counter()

    fes_J = HDivSurface(mesh, order=fes_order)
    fes_L2 = SurfaceL2(mesh, order=max(0, fes_order - 1))
    n_J = fes_J.ndof
    n_f = fes_L2.ndof

    # --- Divergence matrix D: n_f x n_J ---
    u_J = fes_J.TrialFunction()
    q = fes_L2.TestFunction()
    bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
    bf_D += div(u_J.Trace()) * q * ds
    bf_D.Assemble()
    D = _to_dense(bf_D.mat)

    # --- LaplaceSL matrix: n_J x n_J ---
    t0 = time.perf_counter()
    jt, jv = fes_J.TnT()
    with TaskManager():
        V_op = LaplaceSL(jt.Trace() * ds, use_fmm=False) * jv.Trace() * ds
    SL = _to_dense(V_op.mat)
    t_assembly = time.perf_counter() - t0

    # --- Source/sink RHS ---
    f_src = LinearForm(fes_L2)
    f_src += q * ds(source_label)
    f_src.Assemble()
    g_src = f_src.vec.FV().NumPy().copy()
    A_src = np.sum(g_src)

    f_snk = LinearForm(fes_L2)
    f_snk += q * ds(sink_label)
    f_snk.Assemble()
    g_snk = f_snk.vec.FV().NumPy().copy()
    A_snk = np.sum(g_snk)

    if A_src < 1e-30 or A_snk < 1e-30:
        return {"error": (f"Source/sink faces not found "
                          f"(A_src={A_src}, A_snk={A_snk}).  "
                          f"Check BND labels {source_label!r} and "
                          f"{sink_label!r}.  Available: "
                          f"{list(set(mesh.GetBoundaries()))}")}

    g = g_src / A_src - g_snk / A_snk

    # --- Saddle point system (drop one constraint to deflate the
    #     1-D null space of the constant function on closed surfaces) ---
    D_red = D[:-1, :]
    g_red = g[:-1]
    n_constraint = n_f - 1

    t0 = time.perf_counter()
    K = np.block([
        [SL,              D_red.T],
        [D_red, np.zeros((n_constraint, n_constraint))]
    ])
    rhs = np.zeros(n_J + n_constraint)
    rhs[n_J:] = g_red

    x, solve_info = _solve_saddle(K, rhs, method=solver)
    J = x[:n_J]
    t_lu = time.perf_counter() - t0

    # --- Inductance: L = mu_0 * J^T @ SL @ J ---
    L = MU_0 * J @ SL @ J
    residual = np.max(np.abs(D @ J - g))

    # --- (Optional) AC SIBC R = Z_s_re * J^T M J ---
    R_coil = 0.0
    if Z_s_re != 0.0:
        # Mass matrix on HDivSurface: ∫_S (jt · jv) dS
        bf_M = BilinearForm(fes_J)
        bf_M += jt.Trace() * jv.Trace() * ds
        bf_M.Assemble()
        M = _to_dense(bf_M.mat)
        R_coil = float(Z_s_re * (J @ M @ J))

    t_total = time.perf_counter() - t_start

    # --- GridFunction for post-processing ---
    gf_J = GridFunction(fes_J)
    gf_J.vec.FV().NumPy()[:] = J

    return {
        'L': float(L),
        'R': R_coil,
        'n_J': n_J,
        'n_f': n_f,
        'A_source': float(A_src),
        'A_sink': float(A_snk),
        't_assembly': round(t_assembly, 2),
        't_solve': round(t_lu, 2),
        't_total': round(t_total, 2),
        'residual': float(residual),
        'solver': solve_info,
        'J': J,
        'SL': SL,
        'D': D,
        'gf_J': gf_J,
    }


def compute_centroids_areas_J(mesh, gf_J):
    """Sample the HDivSurface coil current at each BND triangle's centroid.

    Returns (centroids, areas, J_per_tri) -- (n_t, *)-shaped float64
    arrays ready to feed into
    ``radia.bem_sibc_solver.compute_phi_inc_from_surface_J`` for the
    workpiece weak-coupling bridge.

    Per-element averaging via NGSolve ``Integrate(..., element_wise=True)``
    over each component of the HDivSurface field.  This avoids the
    HDivSurface-specific reference-coord evaluation API and Just Works
    for any fes_order (RT₀ or higher).
    """
    from ngsolve import Integrate, BND, CF
    n_bnd = mesh.GetNE(BND)
    # Per-element area (denominator for the average).
    elem_A = Integrate(CF(1), mesh, VOL_or_BND=BND, element_wise=True)
    # Per-element integral of each J component.
    elem_J = [
        Integrate(gf_J[i], mesh, VOL_or_BND=BND, element_wise=True)
        for i in range(3)
    ]
    centroids = np.zeros((n_bnd, 3), dtype=np.float64)
    areas = np.zeros(n_bnd, dtype=np.float64)
    J_per_tri = np.zeros((n_bnd, 3), dtype=np.float64)
    for i, el in enumerate(mesh.Elements(BND)):
        a = abs(elem_A[el.nr])
        if a < 1e-30:
            continue
        pts = [mesh.vertices[v.nr].point for v in el.vertices]
        centroids[i] = np.mean(pts, axis=0)
        areas[i] = a
        J_per_tri[i] = [elem_J[k][el.nr] / a for k in range(3)]
    return centroids, areas, J_per_tri
