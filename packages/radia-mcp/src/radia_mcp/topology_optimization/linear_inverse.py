"""Regularized linear inversion for magnetostatic FIELD SYNTHESIS (TSVD + Tikhonov).

The field-synthesis / source-design inverse problem is LINEAR: a forward matrix ``A``
maps a source distribution ``x`` (element magnetizations, coil stream-function weights,
shim moments) to the field samples ``b = A x`` at the target points.  ``A`` is almost
always ILL-CONDITIONED (the field is smooth, high source modes radiate weakly), so the
naive least-squares solution amplifies noise.  The classical cures are spectral filters
of the SVD ``A = U diag(s) V^T``:

    minimum-norm LS  x = sum_n (u_n^T b / s_n) v_n                 (= pinv(A) b)
    TSVD (rank k)    x = sum_{n<=k} (u_n^T b / s_n) v_n           (hard cut at mode k)
    Tikhonov (lam)   x = sum_n phi_n (u_n^T b / s_n) v_n,  phi_n = s_n^2/(s_n^2+lam^2)

TSVD/Tikhonov differ only in the FILTER FACTORS phi_n (1/0 step vs smooth roll-off); both
trade residual ||A x - b|| against solution norm ||x|| -- the L-curve.  This is the linear
core under :func:`topology_opt_applications` topic ``field_synthesis`` and the standard
tool for planar gradient / shim-coil design.

PUBLIC: textbook regularization (Hansen).  All quantities verified analytically in
tests/mcp_server/test_topology_linear_inverse.py.
"""
import numpy as np


def _svd(A):
    U, s, Vt = np.linalg.svd(np.asarray(A, dtype=float), full_matrices=False)
    return U, s, Vt


def tsvd_solve(A, b, k):
    """Truncated-SVD regularized solution of ``A x ~= b`` keeping the ``k`` largest
    singular values:  x = sum_{n<=k} (u_n^T b / s_n) v_n.  k = rank(A) recovers the
    minimum-norm least-squares (pseudo-inverse) solution."""
    U, s, Vt = _svd(A)
    b = np.asarray(b, dtype=float)
    k = int(max(0, min(k, len(s))))
    if k == 0:
        return np.zeros(Vt.shape[1])
    coeffs = (U[:, :k].T @ b) / s[:k]
    return Vt[:k, :].T @ coeffs


def filter_factors(s, lam):
    """Tikhonov filter factors  phi_n = s_n^2 / (s_n^2 + lam^2)  (lam = regularization
    parameter).  phi -> 1 for s >> lam (kept), phi -> 0 for s << lam (damped)."""
    s = np.asarray(s, dtype=float)
    return s * s / (s * s + lam * lam)


def tikhonov_solve(A, b, lam):
    """Tikhonov-regularized solution  x = argmin ||A x - b||^2 + lam^2 ||x||^2
    = sum_n phi_n (u_n^T b / s_n) v_n  with phi_n = s_n^2/(s_n^2+lam^2).  lam -> 0
    recovers the minimum-norm least-squares solution."""
    U, s, Vt = _svd(A)
    b = np.asarray(b, dtype=float)
    phi = filter_factors(s, lam)
    coeffs = phi * (U.T @ b) / s
    return Vt.T @ coeffs


def lcurve(A, b, lams):
    """L-curve trade-off points ``[(residual_norm, solution_norm), ...]`` over a list of
    Tikhonov ``lams`` -- residual ||A x - b|| vs solution norm ||x||.  As lam decreases the
    residual decreases (better fit) and the solution norm grows (more noise amplification);
    the corner balances the two."""
    A = np.asarray(A, dtype=float); b = np.asarray(b, dtype=float)
    out = []
    for lam in lams:
        x = tikhonov_solve(A, b, lam)
        out.append((float(np.linalg.norm(A @ x - b)), float(np.linalg.norm(x))))
    return out


def lcurve_corner(A, b, lams, eps=1.0e-300):
    """Select the L-curve corner by maximum discrete curvature in log-log space.

    The return value is a small dictionary so notebook/MATLAB teaching gates can
    record the chosen lambda, its index, and the full curvature vector without a
    plotting dependency.  Endpoints have zero curvature by definition; at least
    three positive lambda samples are required.
    """
    lams = np.asarray(lams, dtype=float)
    if lams.ndim != 1 or lams.size < 3:
        raise ValueError("lams must be a one-dimensional array with at least three entries")
    if np.any(lams <= 0.0):
        raise ValueError("lams must be positive for log-log L-curve curvature")

    pts = lcurve(A, b, lams)
    residuals = np.asarray([p[0] for p in pts], dtype=float)
    solution_norms = np.asarray([p[1] for p in pts], dtype=float)
    xy = np.column_stack((np.log(np.maximum(residuals, eps)),
                          np.log(np.maximum(solution_norms, eps))))

    curvature = np.zeros(lams.shape, dtype=float)
    for i in range(1, lams.size - 1):
        a = xy[i] - xy[i - 1]
        c = xy[i + 1] - xy[i]
        chord = xy[i + 1] - xy[i - 1]
        denom = np.linalg.norm(a) * np.linalg.norm(c) * np.linalg.norm(chord)
        if denom > 0.0:
            curvature[i] = 2.0 * abs(a[0] * c[1] - a[1] * c[0]) / denom

    index = int(np.argmax(curvature))
    return {
        "lambda": float(lams[index]),
        "index": index,
        "curvature": curvature.tolist(),
        "residual_norm": float(residuals[index]),
        "solution_norm": float(solution_norms[index]),
        "points": pts,
    }
