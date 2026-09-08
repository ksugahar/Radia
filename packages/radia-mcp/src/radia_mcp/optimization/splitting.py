"""Read-only normalized-coordinate proximal and consensus-ADMM diagnostics."""
from __future__ import annotations

import math
from typing import Any

from .diagnostics import _number, _vector


def _norm(values):
    return _number(math.hypot(*values), "residual norm")


def _coordinates(coordinate_system):
    if coordinate_system != "dimensionless":
        raise ValueError("Use dimensionless coordinates; transform objective, gradient and penalty consistently before the run")


def optimization_proximal_gradient_audit(
    x: list[float], smooth_gradient: list[float], proximal_point: list[float],
    step_size: float, penalty: str, weights: list[float],
    coordinate_system: str, groups: list[list[int]] | None = None,
    mapping_tolerance: float = 1e-6, proximal_tolerance: float = 1e-8,
) -> dict[str, Any]:
    """Check supplied p=prox_(alpha R)(x-alpha grad f) for L1/disjoint group L2.

    All inputs belong to a pre-normalized objective f+R. L1 uses one weight per
    coordinate and groups=None. group_l2 uses a disjoint complete partition of
    zero-based coordinate indices, with one weight per group. Weights are
    nonnegative. Only f is differentiated; do NOT include R in smooth_gradient.
    This verifies the proximal subproblem condition, not the supplied gradient,
    a descent step, Lipschitz constant, or a minimum of f+R. No solver executes.
    """
    _coordinates(coordinate_system)
    current = _vector(x, "x")
    gradient = _vector(smooth_gradient, "smooth_gradient")
    point = _vector(proximal_point, "proximal_point")
    if not len(current) == len(gradient) == len(point):
        raise ValueError("x, smooth_gradient and proximal_point lengths must agree")
    alpha = _number(step_size, "step_size", True)
    mtol = _number(mapping_tolerance, "mapping_tolerance", True)
    ptol = _number(proximal_tolerance, "proximal_tolerance", True)
    if penalty == "l1":
        if groups is not None:
            raise ValueError("l1 requires groups=None")
        partition = [[i] for i in range(len(current))]
    elif penalty == "group_l2":
        if not isinstance(groups, list) or not groups or any(
            not isinstance(group, list) or not group for group in groups
        ):
            raise ValueError("group_l2 requires nonempty groups")
        indices = [i for group in groups for i in group]
        if any(isinstance(i, bool) or not isinstance(i, int) for i in indices):
            raise ValueError("group indices must be integers")
        if sorted(indices) != list(range(len(current))):
            raise ValueError("groups must be a disjoint complete partition of zero-based indices")
        partition = groups
    else:
        raise ValueError("Only l1 and disjoint group_l2 are supported; TV/overlap need different proximal operators")
    coefficients = _vector(weights, "weights")
    if len(coefficients) != len(partition) or any(w < 0 for w in coefficients):
        raise ValueError("provide one nonnegative weight per coordinate/group")
    mapping = [_number((v-p)/alpha, "gradient mapping") for v, p in zip(current, point)]
    # Proximal optimality: 0 in (p-x)/alpha + grad f(x) + partial R(p).
    # Displacement form avoids cancellation in p - (x-alpha*gradient).
    direction = [_number(g-m, "proximal stationarity") for g, m in zip(gradient, mapping)]
    residuals = []
    for group, weight in zip(partition, coefficients):
        norm_p = _norm([point[i] for i in group])
        if norm_p == 0:
            residuals.append(max(0.0, _norm([direction[i] for i in group])-weight))
        else:
            residuals.append(_norm([
                _number(direction[i] + weight*(point[i]/norm_p), "proximal residual")
                for i in group]))
    error = max(residuals)
    mapping_norm = max(map(abs, mapping))
    valid = error <= ptol
    return {
        "schema": "radia.optimization.proximal-gradient-audit.v1",
        "status": "invalid_proximal_evidence" if not valid else (
            "first_order_candidate" if mapping_norm <= mtol else "continue"),
        "metrics": {"mapping_inf": mapping_norm, "proximal_stationarity_max_group_l2": error},
        "checks": {"proximal_condition": valid, "mapping_small": mapping_norm <= mtol},
        "group_residuals": residuals,
        "tolerances": {"mapping": mtol, "proximal": ptol},
        "coordinate_system": coordinate_system, "optimality_certified": False,
        "evidence_origin": "caller-supplied; smooth gradient and physical scaling not verified",
        "caution": "A small mapping requires a valid proximal result. Step size, descent, smoothness and curvature are not certified. A tiny nonzero group is not treated as zero.",
    }


def optimization_admm_consensus_audit(
    x: list[float], z: list[float], previous_z: list[float],
    scaled_dual: list[float], rho: float, coordinate_system: str,
    variant: str, absolute_tolerance: float = 1e-6,
    relative_tolerance: float = 1e-4,
) -> dict[str, Any]:
    """Audit standard two-block ADMM for min f(x)+g(z), x=z (fixed rho).

    Supply post-update x^k,z^k,u^k, previous z^(k-1), with scaled dual u=y/rho,
    in one normalized coordinate/objective system. variant must be
    two_block_consensus_fixed_rho: no relaxation, general A/B, or multi-block
    formulas. Only residuals are checked; subproblem solves and dual updates
    are caller-owned and unverified. No optimum or convergence is certified.
    """
    _coordinates(coordinate_system)
    if variant != "two_block_consensus_fixed_rho":
        raise ValueError("Only two_block_consensus_fixed_rho is supported")
    current = _vector(x, "x")
    auxiliary = _vector(z, "z")
    previous = _vector(previous_z, "previous_z")
    dual = _vector(scaled_dual, "scaled_dual")
    if not len(current) == len(auxiliary) == len(previous) == len(dual):
        raise ValueError("x, z, previous_z and scaled_dual lengths must agree")
    penalty = _number(rho, "rho", True)
    atol = _number(absolute_tolerance, "absolute_tolerance", True)
    rtol = _number(relative_tolerance, "relative_tolerance", True)
    primal = [_number(a-b, "primal residual") for a, b in zip(current, auxiliary)]
    # General convention s=rho A^T B(z^k-z^(k-1)), here A=I, B=-I.
    dual_residual = [_number(-penalty*(b-c), "dual residual") for b, c in zip(auxiliary, previous)]
    multiplier = [_number(penalty*u, "unscaled dual") for u in dual]
    base = math.sqrt(len(current))*atol
    epri = _number(base + rtol*max(_norm(current), _norm(auxiliary)), "primal tolerance")
    edual = _number(base + rtol*_norm(multiplier), "dual tolerance")
    rnorm, snorm = _norm(primal), _norm(dual_residual)
    checks = {"primal_small": rnorm <= epri, "dual_small": snorm <= edual}
    return {
        "schema": "radia.optimization.admm-consensus-audit.v1",
        "status": "residuals_small" if all(checks.values()) else "continue",
        "metrics": {"primal_l2": rnorm, "dual_l2": snorm},
        "tolerances": {"primal": epri, "dual": edual}, "checks": checks,
        "coordinate_system": coordinate_system, "optimality_certified": False,
        "subproblem_optimality": "not_verified", "dual_update": "not_verified",
        "evidence_origin": "caller-supplied; iteration history and fixed-rho assumption not verified",
        "caution": "Both residuals are required. Small residuals alone do not prove correct subproblem solves, convexity, a saddle point, or convergence. Do not use this formula for adaptive-rho or relaxed/general/multiblock ADMM.",
    }
