"""Small Krylov solvers for optimization and inverse-design teaching examples.

The routines here are intentionally compact and NumPy-only so they can be read
side-by-side with MATLAB scripts.  The first building block is linear conjugate
gradient for symmetric positive-definite systems, the matrix-free workhorse
behind quadratic minimization and many PDE-constrained optimization inner loops.
"""
import numpy as np


def _as_matvec(operator):
    if callable(operator):
        return operator
    A = np.asarray(operator, dtype=float)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("operator matrix must be square")
    return lambda x: A @ x


def linear_conjugate_gradient(operator, b, x0=None, tol=1e-10, max_iter=None,
                              preconditioner=None):
    """Solve ``A x = b`` for symmetric positive-definite ``A`` by linear CG.

    Args:
        operator       : callable ``v -> A v`` or a square dense matrix.
        b              : right-hand side vector.
        x0             : optional initial guess, default zero.
        tol            : relative residual tolerance ``||r||/max(||b||,1)``.
        max_iter       : iteration cap, default ``len(b)``.
        preconditioner : optional callable ``r -> M^{-1} r``.

    Returns ``{x, n_iter, converged, residual_norm, relative_residual,
    residual_history}``.  For exact arithmetic and an n-dimensional SPD matrix,
    unpreconditioned CG terminates in at most n steps; in floating point the same
    property is a useful regression check on small well-conditioned examples.
    """
    matvec = _as_matvec(operator)
    b = np.asarray(b, dtype=float)
    if b.ndim != 1:
        raise ValueError("b must be a vector")
    if not np.all(np.isfinite(b)):
        raise ValueError("b must be finite")
    n = b.size
    if n == 0:
        raise ValueError("b must be nonempty")
    if not (tol > 0.0):
        raise ValueError("tol must be positive")
    if max_iter is None:
        max_iter = n
    if max_iter < 1:
        raise ValueError("max_iter must be at least 1")

    x = np.zeros_like(b) if x0 is None else np.asarray(x0, dtype=float).copy()
    if x.shape != b.shape:
        raise ValueError("x0 must have the same shape as b")
    if not np.all(np.isfinite(x)):
        raise ValueError("x0 must be finite")

    def apply_preconditioner(r):
        if preconditioner is None:
            return r.copy()
        z = np.asarray(preconditioner(r), dtype=float)
        if z.shape != r.shape:
            raise ValueError("preconditioner output must have the same shape as b")
        if not np.all(np.isfinite(z)):
            raise ValueError("preconditioner output must be finite")
        return z

    r = b - np.asarray(matvec(x), dtype=float)
    if r.shape != b.shape:
        raise ValueError("operator output must have the same shape as b")
    if not np.all(np.isfinite(r)):
        raise ValueError("operator output must be finite")
    norm_scale = max(float(np.linalg.norm(b)), 1.0)
    residual_norm = float(np.linalg.norm(r))
    history = [residual_norm]
    if residual_norm / norm_scale <= tol:
        return {"x": x, "n_iter": 0, "converged": True,
                "residual_norm": residual_norm,
                "relative_residual": residual_norm / norm_scale,
                "residual_history": history}

    z = apply_preconditioner(r)
    p = z.copy()
    rz_old = float(r @ z)
    if rz_old <= 0.0:
        raise ValueError("preconditioner must be positive definite")

    for it in range(1, int(max_iter) + 1):
        Ap = np.asarray(matvec(p), dtype=float)
        if Ap.shape != b.shape:
            raise ValueError("operator output must have the same shape as b")
        if not np.all(np.isfinite(Ap)):
            raise ValueError("operator output must be finite")
        denom = float(p @ Ap)
        if denom <= 0.0:
            raise ValueError("operator must be symmetric positive definite")
        alpha = rz_old / denom
        x = x + alpha * p
        r = r - alpha * Ap
        residual_norm = float(np.linalg.norm(r))
        history.append(residual_norm)
        rel = residual_norm / norm_scale
        if rel <= tol:
            return {"x": x, "n_iter": it, "converged": True,
                    "residual_norm": residual_norm,
                    "relative_residual": rel,
                    "residual_history": history}
        z = apply_preconditioner(r)
        rz_new = float(r @ z)
        if rz_new <= 0.0:
            raise ValueError("preconditioner must be positive definite")
        beta = rz_new / rz_old
        p = z + beta * p
        rz_old = rz_new

    return {"x": x, "n_iter": int(max_iter), "converged": False,
            "residual_norm": residual_norm,
            "relative_residual": residual_norm / norm_scale,
            "residual_history": history}
