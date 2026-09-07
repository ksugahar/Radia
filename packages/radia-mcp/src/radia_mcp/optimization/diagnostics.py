"""Read-only diagnostics on supplied evidence, not a numerical solver."""
from __future__ import annotations

import math
from typing import Any


def _number(value, name, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite real number")
    result = float(value)
    if not math.isfinite(result) or (positive and result <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return result


def _vector(values, name, positive=False):
    if not isinstance(values, list) or not values:
        raise ValueError(f"{name} must be a nonempty list")
    return [_number(v, name, positive) for v in values]


def _scaled_gradient(gradient, scales, objective_scale):
    if len(gradient) != len(scales):
        raise ValueError("gradient and variable_scales lengths must agree")
    return [_number(g * s / objective_scale, "scaled gradient") for g, s in zip(gradient, scales)]


def optimization_gradient_check(
    gradient: list[float], reference_gradients: list[list[float]],
    variable_scales: list[float], objective_scale: float,
    reference_method: str, steps: list[float], atol: float = 1e-7,
    rtol: float = 1e-4,
) -> dict[str, Any]:
    """Compare supplied gradients in dimensionless coordinates over >=2 steps.

    References are computed independently by the caller at the SAME point and
    objective. Steps are distinct positive dimensionless perturbations.
    This checks agreement of evidence, not independence of its origin.
    """
    if reference_method not in {"central_difference", "complex_step"}:
        raise ValueError("reference_method must be central_difference or complex_step")
    scales = _vector(variable_scales, "variable_scales", True)
    fscale = _number(objective_scale, "objective_scale", True)
    g = _scaled_gradient(_vector(gradient, "gradient"), scales, fscale)
    hs = _vector(steps, "steps", True)
    if len(hs) < 2 or len(set(hs)) != len(hs):
        raise ValueError("provide at least two distinct positive steps")
    if len(reference_gradients) != len(hs):
        raise ValueError("one reference gradient is required per step")
    atol = _number(atol, "atol", True)
    rtol = _number(rtol, "rtol", True)
    rows = []
    for h, reference in zip(hs, reference_gradients):
        ref = _scaled_gradient(_vector(reference, "reference gradient"), scales, fscale)
        errors = [_number(abs(a-b), "gradient error") for a, b in zip(g, ref)]
        limits = [_number(atol + rtol * max(abs(a), abs(b)), "error limit") for a, b in zip(g, ref)]
        rows.append({"step": h, "max_abs_error": max(errors),
                     "agrees": all(e <= limit for e, limit in zip(errors, limits))})
    return {"schema": "radia.optimization.gradient-check.v1",
            "status": "consistent" if all(r["agrees"] for r in rows) else "needs_attention",
            "method": reference_method, "samples": rows,
            "evidence_origin": "caller-supplied; independence not verified",
            "caution": "Complex-step requires a holomorphic evaluation path; clipping, abs and branch-dependent solves may invalidate it."}


def optimization_stopping_audit(
    x: list[float], previous_x: list[float], objective: float,
    previous_objective: float, gradient: list[float],
    variable_scales: list[float], objective_scale: float,
    problem_class: str, gradient_tolerance: float = 1e-6,
    step_tolerance: float = 1e-8, objective_tolerance: float = 1e-10,
) -> dict[str, Any]:
    """Audit smooth unconstrained minimization with explicit fixed scales.

    x=S*z, F=f/f_scale, grad_z(F)=S*grad_x(f)/f_scale. First-order evidence
    only: stationary points may be saddles/maxima. Constraints need KKT or
    projected-gradient evidence; nonsmooth objectives need other criteria.
    """
    if problem_class != "smooth_unconstrained":
        raise ValueError("Only smooth_unconstrained is supported; bounds/constraints require projected-gradient or KKT evidence")
    current = _vector(x, "x")
    previous = _vector(previous_x, "previous_x")
    scales = _vector(variable_scales, "variable_scales", True)
    if len(current) != len(previous) or len(current) != len(scales):
        raise ValueError("x, previous_x and variable_scales lengths must agree")
    fscale = _number(objective_scale, "objective_scale", True)
    f = _number(objective, "objective")
    fprevious = _number(previous_objective, "previous_objective")
    g = _scaled_gradient(_vector(gradient, "gradient"), scales, fscale)
    step = max(_number(abs(a-b)/s, "scaled step") for a, b, s in zip(current, previous, scales))
    change = _number(abs(f-fprevious)/fscale, "scaled objective change")
    checks = {
        "gradient_small": max(map(abs, g)) <= _number(gradient_tolerance, "gradient_tolerance", True),
        "step_small": step <= _number(step_tolerance, "step_tolerance", True),
        "objective_change_small": change <= _number(objective_tolerance, "objective_tolerance", True),
        "objective_not_increased": f <= fprevious,
    }
    status = "continue"
    if all(checks.values()):
        status = "first_order_candidate"
    elif checks["step_small"] and checks["objective_change_small"] and not checks["gradient_small"]:
        status = "stalled_nonstationary"
    elif checks["gradient_small"]:
        status = "stationarity_only"
    return {"schema": "radia.optimization.stopping-audit.v1", "status": status,
            "checks": checks, "metrics": {"scaled_gradient_inf": max(map(abs, g)),
                "scaled_step_inf": step, "scaled_objective_change": change},
            "optimality_certified": False,
            "next_check": "Independently check the gradient; inspect curvature and feasibility before claiming a minimum.",
            "evidence_origin": "caller-supplied"}
