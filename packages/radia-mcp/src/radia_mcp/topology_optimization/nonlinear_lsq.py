"""Levenberg-Marquardt nonlinear least squares -- the ``lsqnonlin`` member of the
optimisation toolkit, the NONLINEAR complement to the linear inverse solvers in
:mod:`linear_inverse` (TSVD / Tikhonov) and the derivative-free outer loop
(Nelder-Mead, ``outer_loop`` knowledge topic).

Minimises ``0.5 * ||r(x)||^2`` for a vector residual ``r`` by interpolating between
Gauss-Newton (small damping) and gradient descent (large damping):

    (J^T J + lam I) delta = -J^T r ,   x <- x + delta ,

adapting ``lam`` down on a successful step and up on a rejected one. The Jacobian is
supplied or formed by forward differences. Typical uses in magnet/coil work: recover
geometry or excitation from measured field samples (inverse design / fitting), B-H curve
fits, and any over-determined nonlinear calibration.

SCALING CAVEAT (learned): LM's gradient stop ``||J^T r||_inf < gtol`` is an ABSOLUTE test, so
a residual carried in tiny physical units (e.g. teslas ~ 1e-5) can stop in a flat valley far
from the true minimum. NORMALISE the residual to O(1) (divide by a characteristic scale) for
well-conditioned recovery -- then noiseless data is recovered to machine precision.
"""
import numpy as np


def levenberg_marquardt(residual, x0, jac=None, max_iter=300, gtol=1e-12,
                        xtol=1e-14, lam0=1e-3, eps=1e-7):
    """Levenberg-Marquardt minimiser of ``0.5*||residual(x)||^2``.

    Args:
        residual : callable x -> 1-D residual vector r(x).
        x0       : initial parameter vector.
        jac      : callable x -> Jacobian dr/dx (m x n); forward differences if None.
        max_iter : outer iteration cap.
        gtol     : stop when ``||J^T r||_inf < gtol`` (ABSOLUTE -- scale residuals to O(1)).
        xtol     : stop when a step is smaller than ``xtol*(1+||x||)``.
        lam0     : initial damping.
        eps      : forward-difference step scale (numerical Jacobian).

    Returns ``{x, cost, n_iter, converged, grad_norm}`` (cost = 0.5||r||^2)."""
    x = np.asarray(x0, dtype=float).copy()

    def numjac(xx):
        r0 = np.asarray(residual(xx), float)
        J = np.empty((r0.size, xx.size))
        for j in range(xx.size):
            dx = eps * max(1.0, abs(xx[j]))
            xp = xx.copy(); xp[j] += dx
            J[:, j] = (np.asarray(residual(xp), float) - r0) / dx
        return J

    jfun = jac if jac is not None else numjac
    r = np.asarray(residual(x), float)
    cost = 0.5 * float(r @ r)
    lam = lam0
    converged = False
    g = np.zeros_like(x)
    it = 0
    for it in range(max_iter):
        J = np.asarray(jfun(x), float)
        g = J.T @ r
        if np.linalg.norm(g, np.inf) < gtol:
            converged = True
            break
        H = J.T @ J
        improved = False
        for _ in range(40):
            try:
                delta = np.linalg.solve(H + lam * np.eye(H.shape[0]), -g)
            except np.linalg.LinAlgError:
                lam *= 5.0
                continue
            x_new = x + delta
            r_new = np.asarray(residual(x_new), float)
            cost_new = 0.5 * float(r_new @ r_new)
            if cost_new < cost:
                x, r, cost = x_new, r_new, cost_new
                lam = max(lam / 3.0, 1e-14)
                improved = True
                if np.linalg.norm(delta) < xtol * (1.0 + np.linalg.norm(x)):
                    converged = True
                break
            lam *= 3.0
            if lam > 1e14:
                break
        if converged or not improved:
            converged = converged or not improved
            break
    return {"x": x, "cost": float(cost), "n_iter": it + 1,
            "converged": bool(converged), "grad_norm": float(np.linalg.norm(g, np.inf))}
