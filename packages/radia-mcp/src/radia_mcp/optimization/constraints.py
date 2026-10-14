"""Scaled first-order constraint diagnostics on caller-supplied evidence."""
from __future__ import annotations

import math
from typing import Any

from .diagnostics import _number, _scaled_gradient, _vector


def _norm(values):
    return max((abs(_number(v, "residual")) for v in values), default=0.0)


def _constraint_data(values, jacobian, multipliers, scales, n, name):
    if not all(isinstance(a, list) for a in (values, jacobian, multipliers, scales)):
        raise ValueError(f"{name}: values, row Jacobian, multipliers and scales must be lists")
    if not len(values) == len(jacobian) == len(multipliers) == len(scales):
        raise ValueError(f"{name}: one Jacobian row, multiplier and scale per constraint")
    rows = [_vector(row, f"{name} Jacobian row") for row in jacobian]
    if any(len(row) != n for row in rows):
        raise ValueError(f"{name}: Jacobian rows must match the variable dimension")
    return ([_number(v, name) for v in values], rows,
            [_number(v, f"{name} multiplier") for v in multipliers],
            [_number(v, f"{name} scale", True) for v in scales])


def optimization_kkt_audit(
    gradient: list[float], variable_scales: list[float], objective_scale: float,
    inequality_values: list[float], inequality_jacobian: list[list[float]],
    inequality_multipliers: list[float], inequality_scales: list[float],
    equality_values: list[float], equality_jacobian: list[list[float]],
    equality_multipliers: list[float], equality_scales: list[float],
    problem_class: str, tolerance: float = 1e-6,
) -> dict[str, Any]:
    """Audit smooth minimization with g(x)<=0, h(x)=0 and L=f+mu*g+nu*h.

    Supply physical objective gradient and constraint Jacobian ROWS at the same
    point; mu and nu are physical (unscaled) Lagrangian multipliers. mu>=0;
    nu is unrestricted. Include ALL constraints, including bounds as rows.
    Empty lists denote an absent constraint family. This does not solve for
    multipliers or verify constraint qualifications, smoothness or sufficiency.
    """
    if problem_class != "smooth_constrained":
        raise ValueError("Only smooth_constrained is supported")
    scales = _vector(variable_scales, "variable_scales", True)
    q = _number(objective_scale, "objective_scale", True)
    tol = _number(tolerance, "tolerance", True)
    grad = _scaled_gradient(_vector(gradient, "gradient"), scales, q)
    gi, ji, mu, ci = _constraint_data(inequality_values, inequality_jacobian,
        inequality_multipliers, inequality_scales, len(scales), "inequality")
    he, je, nu, ce = _constraint_data(equality_values, equality_jacobian,
        equality_multipliers, equality_scales, len(scales), "equality")
    if not gi and not he:
        raise ValueError("Supply at least one constraint; use stopping_audit for unconstrained problems")
    g = [_number(v/c, "scaled inequality") for v, c in zip(gi, ci)]
    h = [_number(v/c, "scaled equality") for v, c in zip(he, ce)]
    lam = [_number(m*c/q, "scaled multiplier") for m, c in zip(mu, ci)]
    eta = [_number(m*c/q, "scaled multiplier") for m, c in zip(nu, ce)]
    rows = [[_number(v*s/c, "scaled Jacobian") for v, s in zip(row, scales)]
            for row, c in zip(ji+je, ci+ce)]
    stationarity = [_number(math.fsum([v] + [
        _number(m*row[i], "Lagrangian contribution")
        for m, row in zip(lam+eta, rows)]), "stationarity")
        for i, v in enumerate(grad)]
    metrics = {
        "inequality_violation_inf": _norm(max(v, 0.0) for v in g),
        "equality_violation_inf": _norm(h),
        "dual_violation_inf": _norm(min(m, 0.0) for m in lam),
        "complementarity_inf": _norm(m*v for m, v in zip(lam, g)),
        "stationarity_inf": _norm(stationarity),
    }
    checks = {key: value <= tol for key, value in metrics.items()}
    feasible = checks["inequality_violation_inf"] and checks["equality_violation_inf"]
    return {
        "schema": "radia.optimization.kkt-audit.v1",
        "status": "infeasible" if not feasible else (
            "first_order_candidate" if all(checks.values()) else "needs_attention"),
        "metrics": metrics, "checks": checks, "tolerance": tol,
        "scaled_lagrangian_gradient": stationarity,
        "constraint_qualification": "not_verified",
        "optimality_certified": False,
        "evidence_origin": "caller-supplied; completeness and derivatives not verified",
        "caution": "KKT necessity requires appropriate constraint qualifications. Residuals alone do not certify a minimum or convexity; a failed audit is not proof that no minimum exists.",
    }


def optimization_projected_gradient_audit(
    x: list[float], gradient: list[float], lower_bounds: list[float | None],
    upper_bounds: list[float | None], variable_scales: list[float],
    objective_scale: float, problem_class: str, tolerance: float = 1e-6,
) -> dict[str, Any]:
    """Audit smooth BOX-ONLY minimization; null means an absent bound.

    Report G=z-project_box(z-grad_z(F)) at fixed unit dimensionless step,
    with x=S*z and F=f/objective_scale. Fixed variables are allowed. Supply
    all bounds; nonlinear/general constraints require kkt_audit instead.
    This is a stationarity residual, not an optimizer or minimum certificate.
    """
    if problem_class != "smooth_box_only":
        raise ValueError("Only smooth_box_only is supported; general constraints require kkt_audit")
    current = _vector(x, "x")
    scales = _vector(variable_scales, "variable_scales", True)
    q = _number(objective_scale, "objective_scale", True)
    tol = _number(tolerance, "tolerance", True)
    if not isinstance(lower_bounds, list) or not isinstance(upper_bounds, list):
        raise ValueError("bounds must be lists; use null for an absent bound")
    if not len(current) == len(scales) == len(lower_bounds) == len(upper_bounds):
        raise ValueError("x, scales and bounds must have matching lengths")
    g = _scaled_gradient(_vector(gradient, "gradient"), scales, q)
    residual, violations = [], []
    for value, s, v, lower, upper in zip(current, scales, g, lower_bounds, upper_bounds):
        lo = None if lower is None else _number(lower, "lower bound")
        hi = None if upper is None else _number(upper, "upper bound")
        if lo is not None and hi is not None and lo > hi:
            raise ValueError("lower bound must not exceed upper bound")
        # Evaluate in displacement coordinates: avoids subtracting two large z's.
        a = None if hi is None else _number((value-hi)/s, "upper displacement")
        b = None if lo is None else _number((value-lo)/s, "lower displacement")
        projected = v if a is None else max(v, a)
        projected = projected if b is None else min(projected, b)
        residual.append(projected)
        violations.append(max(0.0, a if a is not None else 0.0, -b if b is not None else 0.0))
    metrics = {"bound_violation_inf": _norm(violations),
               "projected_gradient_inf": _norm(residual)}
    feasible = metrics["bound_violation_inf"] <= tol
    return {
        "schema": "radia.optimization.projected-gradient-audit.v1",
        "status": "infeasible" if not feasible else (
            "first_order_candidate" if metrics["projected_gradient_inf"] <= tol else "needs_attention"),
        "metrics": metrics, "scaled_projected_gradient": residual,
        "tolerance": tol, "dimensionless_step": 1.0,
        "optimality_certified": False,
        "evidence_origin": "caller-supplied; box completeness and derivatives not verified",
        "caution": "A small residual is scale- and tolerance-dependent first-order evidence, not a minimum certificate. Check derivatives and curvature; no general constraints are handled.",
    }
