"""Solver-neutral reverse-mode adjoint scaling validation."""

from __future__ import annotations

import math
from typing import Any


def adjoint_gradient_scaling_gate(
    rows: list[dict[str, Any]],
    *,
    max_gradient_relative_error: float = 1.0e-6,
    max_forward_affine_residual: float = 1.0e-10,
    min_final_objective_ratio: float = 1.0,
) -> dict[str, Any]:
    """Gate gradient correctness, direction and solve-count independence."""

    if not isinstance(rows, list) or len(rows) < 2:
        raise ValueError("rows must contain at least two design-variable counts")
    tolerances = (
        float(max_gradient_relative_error),
        float(max_forward_affine_residual),
        float(min_final_objective_ratio),
    )
    if any(not math.isfinite(value) or value < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")

    def finite(row: dict[str, Any], key: str) -> float:
        try:
            value = float(row[key])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"row field {key!r} must be numeric") from exc
        if not math.isfinite(value):
            raise ValueError(f"row field {key!r} must be finite")
        return value

    counts = [int(finite(row, "designVariableCount")) for row in rows]
    adjoint_solves = [int(finite(row, "adjointSolves")) for row in rows]
    gradient_errors = [finite(row, "gradientCheckRelativeError") for row in rows]
    affine_residuals = [finite(row, "forwardAffineResidual") for row in rows]
    plus_ratios = [finite(row, "plusAscentObjectiveRatio") for row in rows]
    minus_ratios = [finite(row, "minusAscentObjectiveRatio") for row in rows]
    final_ratios = [finite(row, "fiftyStepObjectiveRatio") for row in rows]
    monotone = [row.get("fiftyStepMonotone") is True for row in rows]

    checks = {
        "design_variable_counts_strictly_increase": all(a < b for a, b in zip(counts, counts[1:])),
        "one_adjoint_solve_for_every_size": all(value == 1 for value in adjoint_solves),
        "gradients_match_finite_difference": all(value <= max_gradient_relative_error for value in gradient_errors),
        "forward_map_affine_residual_small": all(value <= max_forward_affine_residual for value in affine_residuals),
        "positive_direction_raises_objective": all(value > 1.0 for value in plus_ratios),
        "negative_direction_lowers_objective": all(value < 1.0 for value in minus_ratios),
        "iterative_ascent_is_monotone": all(monotone),
        "iterative_ascent_improves_objective": all(value > min_final_objective_ratio for value in final_ratios),
    }
    return {
        "policy": "adjoint_gradient_scaling_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "design_variable_counts": counts,
            "adjoint_solves": adjoint_solves,
            "max_gradient_relative_error": max(gradient_errors),
            "max_forward_affine_residual": max(affine_residuals),
            "minimum_plus_ascent_ratio": min(plus_ratios),
            "maximum_minus_ascent_ratio": max(minus_ratios),
            "minimum_final_objective_ratio": min(final_ratios),
        },
        "notes": [
            "finite-difference agreement alone does not establish the sign of a complex ascent direction",
            "reverse-mode scaling requires one adjoint solve as the design-variable count grows",
        ],
    }
