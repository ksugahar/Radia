"""
BEM inductance extraction via source/sink saddle point EFIE.

Core solver for self-inductance extraction using ngsolve.bem LaplaceSL.

Solves the constrained EFIE on the conductor surface:
  [SL  D^T] [J] = [0]
  [D   0  ] [p] = [g]

where:
  SL = LaplaceSL (single layer BEM operator on surface currents)
  D  = divergence matrix (HDivSurface -> SurfaceL2)
  g  = source/sink current injection (+1 at source, -1 at sink)
  J  = surface current (unknowns)
  p  = Lagrange multiplier enforcing current conservation

Then inductance: L = mu_0 * J^T @ SL @ J

Requires: NGSolve with ngsolve.bem, scipy
"""

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
    Returns ``(x, info)`` where info is a dict with keys
    ``method, iterations, residual``.
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
        f"Unknown saddle-point solver: {method!r}. "
        f"Choices: lu, minres, gmres.")


def _to_dense(mat):
    """Extract dense NumPy array from NGSolve BaseMatrix via COO.

    NGSolve BEM operators (LaplaceSL etc.) store the full dense matrix
    internally as a SparseMatrix with 100% fill.  The built-in ToDense()
    is O(N) slower than necessary (~144s vs 0.06s at N=5085) because it
    performs N column-by-column MatVecs instead of a direct memory copy.

    This function extracts the COO triplets and converts via scipy,
    which is ~2500x faster.
    """
    rows, cols, vals = mat.COO()
    return coo_matrix((vals, (rows, cols)),
                      shape=(mat.height, mat.width)).toarray()


def compute_inductance_source_sink(mesh, source_label="source", sink_label="sink",
                                    fes_order=0, solver="lu"):
    """Compute self-inductance via saddle point EFIE with source/sink ports.

    Args:
        mesh: NGSolve Mesh (volume mesh with boundary, or surface-only mesh)
        source_label: Boundary label for current injection face
        sink_label: Boundary label for current extraction face
        fes_order: HDivSurface polynomial order (0=RWG, 1=higher-order)
        solver: saddle-point solver — "lu" (dense direct, default),
            "minres" (Krylov, symmetric indefinite), or "gmres".

    Returns:
        dict with keys:
            L: inductance [H]
            n_J: number of HDivSurface DOFs (edges)
            n_f: number of SurfaceL2 DOFs (faces)
            A_source, A_sink: port areas [m^2]
            t_assembly: BEM integral assembly time [s]
            t_solve: saddle point LU solve time [s]
            t_total: total solver time [s]
            residual: max|D*J - g|
            J: surface current coefficients (ndarray)
            SL: LaplaceSL dense matrix (ndarray, for post-processing)
            gf_J: GridFunction(HDivSurface) with solved J
        or dict with 'error' key on failure.
    """
    from ngsolve import (HDivSurface, SurfaceL2, TaskManager, ds, BND,
                         BilinearForm, LinearForm, div, GridFunction)
    from ngsolve.bem import LaplaceSL

    t_start = time.perf_counter()

    nse = mesh.GetNE(BND)
    nv = mesh.nv

    fes_J = HDivSurface(mesh, order=fes_order)
    # Lagrange multiplier space must be one polynomial order LOWER than the
    # current space so that div(HDivSurface) maps onto it (de Rham complex):
    #     div: HDivSurface(order=p) -> SurfaceL2(order=p-1)
    # For order=0 (RWG) the divergence is piecewise constant, so order=0
    # works. For order>=1 we need order=p-1 to capture the full image.
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
        return {"error": f"Source/sink faces not found (A_src={A_src}, A_snk={A_snk}). "
                f"Check boundary labels: '{source_label}', '{sink_label}'. "
                f"Available: {list(set(mesh.GetBoundaries()))}"}

    g = g_src / A_src - g_snk / A_snk

    # --- Saddle point system ---
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

    t_total = time.perf_counter() - t_start

    # --- GridFunction for post-processing ---
    gf_J = GridFunction(fes_J)
    gf_J.vec.FV().NumPy()[:] = J

    return {
        'L': float(L),
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
        'gf_J': gf_J,
    }
