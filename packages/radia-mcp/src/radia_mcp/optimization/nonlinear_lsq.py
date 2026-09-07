"""Solver-neutral Levenberg-Marquardt teaching helper -- the ``lsqnonlin`` member of the
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

    Returns ``{x, cost, n_iter, converged, grad_norm, termination_reason}``.
    Cost and gradient describe the RETURNED x. ``converged`` means only that
    the absolute gradient tolerance holds, not that x is a minimum. A small
    step or failed improvement alone no longer reports convergence. Callbacks
    must be deterministic; this compact helper does not support constraints.
    """
    x = np.asarray(x0, dtype=float).copy()
    if x.ndim != 1 or x.size == 0 or not np.all(np.isfinite(x)):
        raise ValueError("x0 must be a finite nonempty vector")
    if isinstance(max_iter, bool) or not isinstance(max_iter, (int, np.integer)) or max_iter < 1:
        raise ValueError("max_iter must be a positive integer")
    if any(not np.isfinite(v) or v <= 0 for v in (gtol, xtol, lam0, eps)):
        raise ValueError("gtol, xtol, lam0 and eps must be finite and positive")

    def evaluate(xx):
        rr = np.asarray(residual(xx), dtype=float)
        if rr.ndim != 1 or rr.size == 0 or not np.all(np.isfinite(rr)):
            raise ValueError("residual must return a finite nonempty vector")
        return rr

    def numjac(xx):
        r0 = evaluate(xx)
        J = np.empty((r0.size, xx.size))
        for j in range(xx.size):
            dx = eps * max(1.0, abs(xx[j]))
            xp = xx.copy(); xp[j] += dx
            rp = evaluate(xp)
            if rp.shape != r0.shape:
                raise ValueError("residual shape must remain constant")
            J[:, j] = (rp - r0) / dx
        return J

    jfun = jac if jac is not None else numjac
    r = evaluate(x)
    def checked_jac(xx):
        jj = np.asarray(jfun(xx), dtype=float)
        if jj.shape != (r.size, x.size) or not np.all(np.isfinite(jj)):
            raise ValueError("Jacobian must be finite with shape (residuals, variables)")
        return jj

    cost = 0.5 * float(r @ r)
    if not np.isfinite(cost):
        raise ValueError("residual cost must be finite; rescale the problem")
    lam = lam0
    reason = "iteration_limit"
    g = np.zeros_like(x)
    it = 0
    for it in range(max_iter):
        J = checked_jac(x)
        g = J.T @ r
        if not np.all(np.isfinite(g)):
            raise ValueError("gradient must be finite; rescale the problem")
        if np.linalg.norm(g, np.inf) < gtol:
            reason = "gradient_tolerance"
            break
        H = J.T @ J
        if not np.all(np.isfinite(H)):
            raise ValueError("normal matrix must be finite; rescale the problem")
        improved = False
        small_step = False
        for _ in range(40):
            try:
                delta = np.linalg.solve(H + lam * np.eye(H.shape[0]), -g)
            except np.linalg.LinAlgError:
                lam *= 5.0
                continue
            x_new = x + delta
            if not np.all(np.isfinite(x_new)):
                raise ValueError("trial point must be finite; rescale the problem")
            r_new = evaluate(x_new)
            if r_new.shape != r.shape:
                raise ValueError("residual shape must remain constant")
            cost_new = 0.5 * float(r_new @ r_new)
            if not np.isfinite(cost_new):
                raise ValueError("trial cost must be finite; rescale the problem")
            if cost_new < cost:
                x, r, cost = x_new, r_new, cost_new
                lam = max(lam / 3.0, 1e-14)
                improved = True
                if np.linalg.norm(delta) < xtol * (1.0 + np.linalg.norm(x)):
                    small_step = True
                break
            lam *= 3.0
            if lam > 1e14:
                break
        if small_step or not improved:
            reason = "step_tolerance" if small_step else "no_improvement"
            break
    # An accepted last step changes x: never return its predecessor's gradient.
    g = checked_jac(x).T @ r
    if not np.all(np.isfinite(g)):
        raise ValueError("final gradient must be finite; rescale the problem")
    grad_norm = float(np.linalg.norm(g, np.inf))
    converged = grad_norm < gtol
    if converged:
        reason = "gradient_tolerance"
    return {"x": x, "cost": float(cost), "n_iter": it + 1,
            "converged": bool(converged), "grad_norm": grad_norm,
            "termination_reason": reason}
