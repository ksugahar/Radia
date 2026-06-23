"""Readable line-search helpers for local optimization.

This module intentionally keeps the algorithm small enough to compare with
MATLAB teaching scripts: given a descent direction ``p``, shrink ``alpha`` until

    f(x + alpha p) <= f(x) + c1 alpha grad(f)(x)^T p.

That Armijo sufficient-decrease gate is the usual stabilizer behind steepest
descent, nonlinear CG, and quasi-Newton methods.
"""
import numpy as np


def armijo_backtracking(func, x, direction, grad=None, f0=None, g0=None,
                        alpha0=1.0, rho=0.5, c1=1e-4, min_alpha=1e-12,
                        max_iter=50):
    """Backtracking line search for a descent direction.

    Args:
        func      : callable ``x -> scalar objective``.
        x         : current point.
        direction : search direction ``p``.
        grad      : optional callable ``x -> gradient``.
        f0        : optional cached ``func(x)``.
        g0        : optional cached gradient at ``x``.
        alpha0    : first trial step length.
        rho       : shrink factor in ``(0, 1)``.
        c1        : Armijo constant in ``(0, 1)``.
        min_alpha : smallest useful trial step.
        max_iter  : maximum number of trial steps.

    Returns ``{alpha, x_new, f_new, n_iter, nfev, accepted, slope0}``.
    ``g0`` or ``grad`` must be supplied; numerical gradients are intentionally
    left to callers so the scaling and perturbation policy stays explicit.
    """
    if not (alpha0 > 0.0):
        raise ValueError("alpha0 must be positive")
    if not (0.0 < rho < 1.0):
        raise ValueError("rho must be in (0, 1)")
    if not (0.0 < c1 < 1.0):
        raise ValueError("c1 must be in (0, 1)")
    if not (min_alpha > 0.0):
        raise ValueError("min_alpha must be positive")
    if max_iter < 1:
        raise ValueError("max_iter must be at least 1")

    x = np.asarray(x, dtype=float)
    p = np.asarray(direction, dtype=float)
    if x.shape != p.shape:
        raise ValueError("x and direction must have the same shape")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(p)):
        raise ValueError("x and direction must be finite")

    nfev = 0
    f_base = float(func(x)) if f0 is None else float(f0)
    if f0 is None:
        nfev += 1
    if not np.isfinite(f_base):
        raise ValueError("f0 must be finite")

    if g0 is None:
        if grad is None:
            raise ValueError("g0 or grad must be supplied")
        g = np.asarray(grad(x), dtype=float)
    else:
        g = np.asarray(g0, dtype=float)
    if g.shape != x.shape:
        raise ValueError("gradient must have the same shape as x")
    if not np.all(np.isfinite(g)):
        raise ValueError("gradient must be finite")

    slope0 = float(g @ p)
    if slope0 >= 0.0:
        raise ValueError("direction must be a strict descent direction")

    alpha = float(alpha0)
    last_x = x.copy()
    last_f = f_base
    last_alpha = 0.0

    for it in range(1, max_iter + 1):
        x_trial = x + alpha * p
        f_trial = float(func(x_trial))
        nfev += 1
        last_x = x_trial
        last_f = f_trial
        last_alpha = alpha
        if np.isfinite(f_trial) and f_trial <= f_base + c1 * alpha * slope0:
            return {"alpha": alpha, "x_new": x_trial, "f_new": f_trial,
                    "n_iter": it, "nfev": nfev, "accepted": True,
                    "slope0": slope0}
        alpha *= rho
        if alpha < min_alpha:
            break

    return {"alpha": last_alpha, "x_new": last_x, "f_new": float(last_f),
            "n_iter": it, "nfev": nfev, "accepted": False,
            "slope0": slope0}
