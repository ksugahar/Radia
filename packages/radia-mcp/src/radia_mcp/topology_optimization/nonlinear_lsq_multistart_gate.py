"""Validate nonlinear least-squares multistart and Jacobian evidence."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence


def _number(value: object, name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def _vector(value: object, name: str) -> list[float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise ValueError(f"{name} must be a non-empty array")
    return [_number(item, name) for item in value]


def _relative_error(value: float, reference: float) -> float:
    return abs(value - reference) / max(abs(value), abs(reference), 1.0e-300)


def evaluate_nonlinear_lsq_multistart(summary: Mapping[str, object]) -> dict:
    """Gate corrected multistart diversity, stationarity, and Jacobians."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    contract = summary.get("contract")
    tolerances = summary.get("tolerances")
    legacy_starts = summary.get("legacy_starts")
    runs = summary.get("runs")
    if not isinstance(contract, Mapping) or not isinstance(tolerances, Mapping):
        raise ValueError("contract and tolerances must be objects")
    if not isinstance(legacy_starts, Sequence) or isinstance(legacy_starts, (str, bytes)):
        raise ValueError("legacy_starts must be an array")
    if not isinstance(runs, Sequence) or isinstance(runs, (str, bytes)):
        raise ValueError("runs must be an array")

    residual_tolerance = _number(
        tolerances.get("residual_norm"), "tolerances.residual_norm"
    )
    gradient_tolerance = _number(
        tolerances.get("projected_gradient_inf_norm"),
        "tolerances.projected_gradient_inf_norm",
    )
    solver_jacobian_tolerance = _number(
        tolerances.get("solver_jacobian_relative_error"),
        "tolerances.solver_jacobian_relative_error",
    )
    finite_difference_tolerance = _number(
        tolerances.get("finite_difference_jacobian_relative_error"),
        "tolerances.finite_difference_jacobian_relative_error",
    )
    resnorm_tolerance = _number(
        tolerances.get("resnorm_identity_relative_error"),
        "tolerances.resnorm_identity_relative_error",
    )
    if min(
        residual_tolerance,
        gradient_tolerance,
        solver_jacobian_tolerance,
        finite_difference_tolerance,
        resnorm_tolerance,
    ) < 0.0:
        raise ValueError("tolerances must be nonnegative")

    legacy = [
        _vector(row, f"legacy_starts[{index}]")
        for index, row in enumerate(legacy_starts)
    ]
    dimension = len(legacy[0]) if legacy else 0
    if any(len(row) != dimension for row in legacy):
        raise ValueError("legacy start dimensions must match")

    parsed_runs = []
    parse_errors: list[str] = []
    for index, run in enumerate(runs):
        if not isinstance(run, Mapping):
            parse_errors.append(f"runs[{index}] must be an object")
            continue
        try:
            start = _vector(run.get("start"), f"runs[{index}].start")
            solution = _vector(run.get("solution"), f"runs[{index}].solution")
            residual = _vector(run.get("residual"), f"runs[{index}].residual")
            if len(start) != dimension or len(solution) != dimension:
                raise ValueError("start/solution dimensions must match legacy starts")
            residual_norm = _number(
                run.get("residual_norm"), f"runs[{index}].residual_norm"
            )
            resnorm = _number(run.get("resnorm"), f"runs[{index}].resnorm")
            projected_gradient = _number(
                run.get("projected_gradient_inf_norm"),
                f"runs[{index}].projected_gradient_inf_norm",
            )
            solver_jacobian_error = _number(
                run.get("solver_jacobian_relative_error"),
                f"runs[{index}].solver_jacobian_relative_error",
            )
            finite_difference_error = _number(
                run.get("finite_difference_jacobian_relative_error"),
                f"runs[{index}].finite_difference_jacobian_relative_error",
            )
            exitflag = int(run["exitflag"])
            iterations = int(run["iterations"])
            func_count = int(run["func_count"])
            residual_vector_norm = math.sqrt(sum(value * value for value in residual))
            parsed_runs.append(
                {
                    "start": start,
                    "solution": solution,
                    "residual_norm": residual_norm,
                    "residual_vector_norm": residual_vector_norm,
                    "resnorm": resnorm,
                    "projected_gradient_inf_norm": projected_gradient,
                    "solver_jacobian_relative_error": solver_jacobian_error,
                    "finite_difference_jacobian_relative_error": finite_difference_error,
                    "exitflag": exitflag,
                    "iterations": iterations,
                    "func_count": func_count,
                    "resnorm_identity_relative_error": _relative_error(
                        resnorm, residual_vector_norm**2
                    ),
                }
            )
        except (KeyError, TypeError, ValueError) as exc:
            parse_errors.append(f"runs[{index}]: {exc}")

    starts = [tuple(row["start"]) for row in parsed_runs]
    checks = {
        "units_explicit": bool(str(summary.get("parameter_unit") or "").strip())
        and bool(str(summary.get("residual_unit") or "").strip()),
        "least_squares_contract_recorded": contract.get("objective")
        == "0.5*||r||^2"
        and contract.get("gradient") == "J^T*r"
        and contract.get("resnorm_identity") == "resnorm=||r||^2"
        and contract.get("solver_jacobian_semantics") == "residual_jacobian",
        "legacy_zero_multiplicative_start_collapse_detected": len(legacy) >= 2
        and dimension > 0
        and all(all(value == 0.0 for value in row) for row in legacy),
        "runs_parsed_and_finite": not parse_errors and len(parsed_runs) == len(runs),
        "corrected_multistart_count_sufficient": len(parsed_runs) >= max(4, dimension + 1),
        "corrected_starts_are_distinct": bool(starts) and len(set(starts)) == len(starts),
        "all_solvers_report_success": bool(parsed_runs)
        and all(row["exitflag"] > 0 for row in parsed_runs),
        "iteration_and_evaluation_counts_positive": bool(parsed_runs)
        and all(row["iterations"] > 0 and row["func_count"] > 0 for row in parsed_runs),
        "residual_vectors_match_reported_norms": bool(parsed_runs)
        and all(
            _relative_error(row["residual_norm"], row["residual_vector_norm"])
            <= resnorm_tolerance
            for row in parsed_runs
        ),
        "resnorm_square_identity": bool(parsed_runs)
        and all(
            row["resnorm_identity_relative_error"] <= resnorm_tolerance
            for row in parsed_runs
        ),
        "residuals_within_tolerance": bool(parsed_runs)
        and all(row["residual_norm"] <= residual_tolerance for row in parsed_runs),
        "projected_stationarity_within_tolerance": bool(parsed_runs)
        and all(
            row["projected_gradient_inf_norm"] <= gradient_tolerance
            for row in parsed_runs
        ),
        "solver_jacobians_match_analytic": bool(parsed_runs)
        and all(
            row["solver_jacobian_relative_error"] <= solver_jacobian_tolerance
            for row in parsed_runs
        ),
        "finite_difference_jacobians_match_analytic": bool(parsed_runs)
        and all(
            row["finite_difference_jacobian_relative_error"]
            <= finite_difference_tolerance
            for row in parsed_runs
        ),
    }
    issues = parse_errors + [name for name, ok in checks.items() if not ok]
    return {
        "schema": "radia-topology-nonlinear-lsq-multistart/v1",
        "policy": "independent_multistart_jacobian_and_projected_stationarity_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "dimension": dimension,
            "legacy_start_count": len(legacy),
            "corrected_run_count": len(parsed_runs),
            "distinct_corrected_start_count": len(set(starts)),
            "max_residual_norm": max(
                (row["residual_norm"] for row in parsed_runs), default=math.inf
            ),
            "max_projected_gradient_inf_norm": max(
                (row["projected_gradient_inf_norm"] for row in parsed_runs),
                default=math.inf,
            ),
            "max_solver_jacobian_relative_error": max(
                (row["solver_jacobian_relative_error"] for row in parsed_runs),
                default=math.inf,
            ),
            "max_finite_difference_jacobian_relative_error": max(
                (
                    row["finite_difference_jacobian_relative_error"]
                    for row in parsed_runs
                ),
                default=math.inf,
            ),
            "max_resnorm_identity_relative_error": max(
                (row["resnorm_identity_relative_error"] for row in parsed_runs),
                default=math.inf,
            ),
        },
        "tolerances": {
            "residual_norm": residual_tolerance,
            "projected_gradient_inf_norm": gradient_tolerance,
            "solver_jacobian_relative_error": solver_jacobian_tolerance,
            "finite_difference_jacobian_relative_error": finite_difference_tolerance,
            "resnorm_identity_relative_error": resnorm_tolerance,
        },
        "notes": [
            "Multiplying an all-zero start by random factors does not create a multistart ensemble.",
            "For 0.5*||r||^2, the first-order gradient is J^T*r; resnorm itself is ||r||^2.",
            "Validate the solver-returned residual Jacobian against both an analytic Jacobian and an independent finite difference.",
            "Use a projected gradient for bound-constrained stationarity rather than an unconstrained gradient norm alone.",
        ],
    }


def nonlinear_lsq_multistart_gate(summary_json: str) -> str:
    try:
        summary = json.loads(summary_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"summary_json must be valid JSON: {exc.msg}") from exc
    return json.dumps(evaluate_nonlinear_lsq_multistart(summary), indent=2, sort_keys=True)
