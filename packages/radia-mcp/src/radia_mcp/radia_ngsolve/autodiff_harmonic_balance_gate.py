"""Convergence gate for automatic-differentiation harmonic balance."""
from __future__ import annotations

import math


def _number(mapping: dict, key: str) -> float:
    value = float(mapping.get(key, math.nan))
    return value


def _four_stage_timing(value: object) -> bool:
    if not isinstance(value, dict) or len(value) != 4:
        return False
    try:
        return all(
            math.isfinite(float(seconds)) and float(seconds) >= 0.0
            for seconds in value.values()
        )
    except (TypeError, ValueError):
        return False


def autodiff_harmonic_balance_convergence_gate(
    summary: dict,
    *,
    jacobian_rtol: float = 1.0e-11,
    coefficient_rtol: float = 1.0e-8,
    exact_residual_rms_max: float = 1.0e-10,
    replay_rtol: float = 1.0e-14,
    target_improvement_ratio: float = 0.15,
    gradient_inf_max: float = 1.0e-6,
) -> dict:
    """Reject mean-only false convergence and gate a nonlinear AD solve.

    A zero-mean periodic residual can have a large norm, so its arithmetic mean
    is not a convergence metric.  This gate requires residual RMS and gradient
    stationarity, validates the AD Jacobian by two independent routes, and
    closes an exactly generated coefficient-recovery problem against a second
    nonlinear least-squares implementation and a deterministic replay.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be an object")
    tolerances = (
        jacobian_rtol,
        coefficient_rtol,
        exact_residual_rms_max,
        replay_rtol,
        target_improvement_ratio,
        gradient_inf_max,
    )
    if any(not math.isfinite(value) or value < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")

    problem = summary.get("problem")
    false_control = summary.get("false_convergence_control")
    jacobian = summary.get("jacobian")
    recovery = summary.get("exact_recovery")
    target = summary.get("target_fit")
    stopping = summary.get("stopping")
    sections = (problem, false_control, jacobian, recovery, target, stopping)
    if any(not isinstance(section, dict) for section in sections):
        raise ValueError("problem, convergence, Jacobian, recovery, target, and stopping objects are required")

    point_count = int(problem.get("periodic_point_count", 0))
    variable_count = int(problem.get("variable_count", 0))
    harmonic_max = int(problem.get("odd_harmonic_maximum", 0))
    endpoint_duplicated = bool(problem.get("periodic_endpoint_duplicated", True))

    false_mean = _number(false_control, "mean_residual_abs")
    false_rms = _number(false_control, "rms_residual")
    mean_threshold = _number(false_control, "mean_criterion_threshold")

    ad_analytic = _number(jacobian, "ad_analytic_relative_error")
    ad_complex = _number(jacobian, "ad_complex_step_relative_error")
    jacobian_rows = int(jacobian.get("rows", 0))
    jacobian_columns = int(jacobian.get("columns", 0))

    ad_coefficient = _number(recovery, "ad_coefficient_relative_error")
    reference_coefficient = _number(recovery, "reference_coefficient_relative_error")
    cross_coefficient = _number(recovery, "ad_reference_coefficient_relative_error")
    replay_coefficient = _number(recovery, "ad_replay_coefficient_relative_error")
    ad_rms = _number(recovery, "ad_final_rms_residual")
    reference_rms = _number(recovery, "reference_final_rms_residual")
    reference_exitflag = int(recovery.get("reference_exitflag", 0))

    initial_rms = _number(target, "initial_rms_residual")
    final_rms = _number(target, "final_rms_residual")
    final_gradient = _number(target, "final_gradient_inf_norm")

    finite_values = (
        false_mean,
        false_rms,
        mean_threshold,
        ad_analytic,
        ad_complex,
        ad_coefficient,
        reference_coefficient,
        cross_coefficient,
        replay_coefficient,
        ad_rms,
        reference_rms,
        initial_rms,
        final_rms,
        final_gradient,
    )
    checks = {
        "all_metrics_are_finite": all(math.isfinite(value) for value in finite_values),
        "periodic_grid_has_no_duplicate_endpoint": point_count >= 32
        and variable_count >= 3
        and harmonic_max > 0
        and not endpoint_duplicated,
        "mean_only_rule_is_a_demonstrated_false_convergence": false_mean <= mean_threshold
        and false_rms >= 1.0e-3
        and false_rms >= 1.0e6 * max(false_mean, 1.0e-300),
        "jacobian_dimensions_match_problem": jacobian_rows == point_count
        and jacobian_columns == variable_count,
        "ad_matches_analytic_jacobian": ad_analytic <= jacobian_rtol,
        "ad_matches_complex_step_jacobian": ad_complex <= jacobian_rtol,
        "exact_coefficients_are_recovered": ad_coefficient <= coefficient_rtol
        and reference_coefficient <= coefficient_rtol,
        "independent_nonlinear_solvers_agree": cross_coefficient <= coefficient_rtol
        and reference_exitflag > 0,
        "exact_residuals_close": ad_rms <= exact_residual_rms_max
        and reference_rms <= exact_residual_rms_max,
        "fresh_ad_replay_is_deterministic": replay_coefficient <= replay_rtol,
        "nonexact_target_improves_and_is_stationary": initial_rms > 0.0
        and final_rms <= target_improvement_ratio * initial_rms
        and final_gradient <= gradient_inf_max,
        "rms_and_gradient_stopping_replace_mean_only_stopping": bool(
            stopping.get("residual_rms_used")
        )
        and bool(stopping.get("gradient_inf_used"))
        and not bool(stopping.get("residual_mean_only_used", True)),
        "exactly_four_timing_stages": _four_stage_timing(
            summary.get("timing_breakdown_s")
        ),
    }
    return {
        "policy": "autodiff_harmonic_balance_convergence_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "false_mean_residual_abs": false_mean,
            "false_rms_residual": false_rms,
            "maximum_jacobian_relative_error": max(ad_analytic, ad_complex),
            "maximum_exact_coefficient_relative_error": max(
                ad_coefficient, reference_coefficient, cross_coefficient
            ),
            "ad_replay_coefficient_relative_error": replay_coefficient,
            "target_rms_improvement_ratio": final_rms / max(initial_rms, 1.0e-300),
            "target_gradient_inf_norm": final_gradient,
        },
        "lesson": (
            "A periodic residual may average to zero while remaining large. "
            "Stop on residual RMS and gradient stationarity, triple-check the "
            "AD Jacobian, and require exact recovery, an independent nonlinear "
            "solver, and a fresh deterministic replay."
        ),
    }
