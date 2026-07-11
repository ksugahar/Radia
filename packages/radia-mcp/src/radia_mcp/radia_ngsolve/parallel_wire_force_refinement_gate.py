"""Mesh-refinement validation for a reciprocal two-wire force pair."""
from __future__ import annotations

import math


def parallel_wire_force_refinement_gate(
    refinement_levels,
    force_wire1_rows,
    force_wire2_rows,
    *,
    expected_force_magnitude: float,
    separation_direction=(1.0, 0.0),
    expected_wire2_radial_sign: int | None = None,
    min_sample_count: int = 3,
    max_final_relative_error: float = 0.01,
    max_final_pair_relative_residual: float = 0.01,
    max_final_transverse_relative_force: float = 0.01,
    min_initial_to_final_error_ratio: float = 1.2,
):
    """Gate a non-monotone force-convergence sweep.

    Refinement errors often oscillate, so the contract compares the first and
    final samples instead of requiring every intermediate error to decrease.
    """

    levels = [float(value) for value in refinement_levels]
    f1 = [[float(component) for component in row] for row in force_wire1_rows]
    f2 = [[float(component) for component in row] for row in force_wire2_rows]
    if len(levels) != len(f1) or len(levels) != len(f2):
        raise ValueError("refinement_levels and both force row sets must have the same length")
    if any(len(row) != 2 for row in f1 + f2):
        raise ValueError("force rows must contain exactly two Cartesian components")
    if min_sample_count < 2:
        raise ValueError("min_sample_count must be >= 2")
    expected = float(expected_force_magnitude)
    if not math.isfinite(expected) or expected <= 0.0:
        raise ValueError("expected_force_magnitude must be finite and positive")
    if expected_wire2_radial_sign not in (None, -1, 1):
        raise ValueError("expected_wire2_radial_sign must be -1, 1, or None")
    thresholds = (
        max_final_relative_error,
        max_final_pair_relative_residual,
        max_final_transverse_relative_force,
        min_initial_to_final_error_ratio,
    )
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in thresholds):
        raise ValueError("all tolerances and ratios must be finite and nonnegative")

    direction = [float(value) for value in separation_direction]
    if len(direction) != 2 or not all(math.isfinite(value) for value in direction):
        raise ValueError("separation_direction must contain two finite components")
    norm = math.hypot(*direction)
    if norm <= 0.0:
        raise ValueError("separation_direction must be nonzero")
    radial = [direction[0] / norm, direction[1] / norm]
    transverse = [-radial[1], radial[0]]

    finite = all(math.isfinite(value) for value in levels)
    finite = finite and all(math.isfinite(value) for row in f1 + f2 for value in row)
    increasing = finite and all(right > left for left, right in zip(levels, levels[1:]))

    def project(rows, axis):
        return [row[0] * axis[0] + row[1] * axis[1] for row in rows]

    radial1 = project(f1, radial)
    radial2 = project(f2, radial)
    transverse1 = project(f1, transverse)
    transverse2 = project(f2, transverse)
    relative_errors = [abs(abs(value) - expected) / expected for value in radial2]
    pair_residuals = [
        math.hypot(left[0] + right[0], left[1] + right[1]) / expected
        for left, right in zip(f1, f2)
    ]
    transverse_ratios = [
        max(abs(left), abs(right)) / expected
        for left, right in zip(transverse1, transverse2)
    ]
    initial_error = relative_errors[0] if relative_errors else math.inf
    final_error = relative_errors[-1] if relative_errors else math.inf
    improvement_ratio = (
        initial_error / max(final_error, 1.0e-300) if relative_errors else 0.0
    )
    opposite_pair_sign = all(left * right < 0.0 for left, right in zip(radial1, radial2))
    expected_sign_ok = expected_wire2_radial_sign is None or all(
        math.copysign(1.0, value) == float(expected_wire2_radial_sign)
        for value in radial2
    )
    checks = {
        "sample_count_sufficient": len(levels) >= int(min_sample_count),
        "all_finite": finite,
        "refinement_levels_strictly_increase": increasing,
        "opposite_pair_force_sign": opposite_pair_sign,
        "wire2_radial_sign_matches_expectation": expected_sign_ok,
        "final_analytic_error_ok": final_error <= float(max_final_relative_error),
        "final_action_reaction_residual_ok": bool(pair_residuals) and pair_residuals[-1] <= float(max_final_pair_relative_residual),
        "final_transverse_force_ok": bool(transverse_ratios) and transverse_ratios[-1] <= float(max_final_transverse_relative_force),
        "first_to_final_error_improves": improvement_ratio >= float(min_initial_to_final_error_ratio),
    }
    monotone_error = all(right <= left for left, right in zip(relative_errors, relative_errors[1:]))
    return {
        "policy": "parallel_wire_force_refinement_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "sample_count": len(levels),
        "refinement_levels": levels,
        "expected_force_magnitude": expected,
        "wire1_radial_force": radial1,
        "wire2_radial_force": radial2,
        "relative_errors": relative_errors,
        "pair_relative_residuals": pair_residuals,
        "transverse_relative_forces": transverse_ratios,
        "initial_relative_error": initial_error,
        "final_relative_error": final_error,
        "initial_to_final_error_ratio": improvement_ratio,
        "error_is_monotone": monotone_error,
        "checks": checks,
        "lesson": (
            "A force-refinement sweep may converge non-monotonically. Gate the final analytic error, "
            "first-to-final improvement, reciprocal pair balance, transverse leakage, and sign together."
        ),
    }
