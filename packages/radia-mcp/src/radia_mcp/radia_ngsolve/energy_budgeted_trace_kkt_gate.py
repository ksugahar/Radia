"""Solver-neutral KKT gate for an energy-budgeted FEM/BEM trace fit."""

from __future__ import annotations

import math
from typing import Any


def energy_budgeted_trace_kkt_gate(
    payload: dict[str, Any],
    *,
    max_gradient_relative_error: float = 1.0e-6,
    max_solution_relative_error: float = 1.0e-5,
    max_stationarity_inf: float = 1.0e-6,
    max_complementarity_abs: float = 1.0e-7,
    max_constraint_relative: float = 1.0e-8,
) -> dict[str, Any]:
    """Gate primal, dual, stationarity, complementarity, and method agreement."""

    if not isinstance(payload, dict):
        raise ValueError("payload must be an object")
    tolerances = [
        max_gradient_relative_error,
        max_solution_relative_error,
        max_stationarity_inf,
        max_complementarity_abs,
        max_constraint_relative,
    ]
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")

    mesh = payload.get("mesh") or {}
    problem = payload.get("problem") or {}
    analytic = payload.get("analytic") or {}
    metrics = payload.get("metrics") or {}
    solvers = payload.get("solvers")
    if not isinstance(solvers, list) or len(solvers) != 2:
        raise ValueError("solvers must contain exactly two records")
    by_algorithm = {
        str(row.get("algorithm") or "").strip().lower(): row
        for row in solvers
        if isinstance(row, dict)
    }
    if set(by_algorithm) != {"sqp", "interior-point"}:
        raise ValueError("solver algorithms must be sqp and interior-point")

    try:
        budget = float(problem["energy_budget"])
        budget_fraction = float(problem["budget_fraction"])
        analytic_energy = float(analytic["energy"])
        analytic_constraint = float(analytic["constraint_value"])
        analytic_stationarity = float(analytic["stationarity_inf"])
        analytic_complementarity = float(analytic["complementarity_abs"])
        objective_gradient_error = float(metrics["objective_gradient_relative_error"])
        constraint_gradient_error = float(metrics["constraint_gradient_relative_error"])
        pair_difference = float(metrics["solver_pair_solution_relative_difference"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("problem, analytic, or metric fields are invalid") from exc
    numeric = [
        budget,
        budget_fraction,
        analytic_energy,
        analytic_constraint,
        analytic_stationarity,
        analytic_complementarity,
        objective_gradient_error,
        constraint_gradient_error,
        pair_difference,
    ]
    if not all(math.isfinite(value) for value in numeric):
        raise ValueError("KKT metrics must be finite")

    solver_metrics = {}
    solver_rows_ok = True
    for algorithm, row in by_algorithm.items():
        try:
            parsed = {
                "exitflag": int(row["exitflag"]),
                "objective": float(row["objective"]),
                "energy": float(row["energy"]),
                "constraint_value": float(row["constraint_value"]),
                "dual_lambda": float(row["dual_lambda"]),
                "stationarity_inf": float(row["stationarity_inf"]),
                "complementarity_abs": float(row["complementarity_abs"]),
                "solution_relative_error": float(row["solution_relative_error"]),
                "objective_relative_error": float(row["objective_relative_error"]),
            }
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"invalid solver record: {algorithm}") from exc
        solver_rows_ok = solver_rows_ok and all(
            math.isfinite(float(value)) for value in parsed.values()
        )
        solver_metrics[algorithm] = parsed

    energy_scale = max(abs(budget), 1.0e-30)
    checks = {
        "vol_tri_tet_p1_trace_contract": (
            mesh.get("volume_element") == "tet4"
            and mesh.get("boundary_element") == "tri3"
            and mesh.get("volume_basis") == "H1_P1"
            and mesh.get("trace_basis") == "scalar_P1"
            and int(mesh.get("points", 0)) > int(mesh.get("trace_dofs", 0)) > 0
            and int(mesh.get("tets", 0)) > 0
            and int(mesh.get("triangles", 0)) > 0
        ),
        "strict_positive_active_budget": budget > 0.0
        and 0.0 < budget_fraction < 1.0
        and abs(analytic_energy - budget) / energy_scale <= max_constraint_relative,
        "analytic_kkt_stationary": analytic_stationarity <= max_stationarity_inf,
        "analytic_complementarity": analytic_complementarity
        <= max_complementarity_abs
        and abs(analytic_constraint) / energy_scale <= max_constraint_relative,
        "objective_gradient_matches_finite_difference": objective_gradient_error
        <= max_gradient_relative_error,
        "constraint_gradient_matches_finite_difference": constraint_gradient_error
        <= max_gradient_relative_error,
        "solver_records_are_finite": solver_rows_ok,
        "both_solvers_converged": all(
            row["exitflag"] > 0 for row in solver_metrics.values()
        ),
        "both_solvers_primal_feasible": all(
            row["constraint_value"] / energy_scale <= max_constraint_relative
            for row in solver_metrics.values()
        ),
        "both_solvers_positive_dual": all(
            row["dual_lambda"] > 0.0 for row in solver_metrics.values()
        ),
        "both_solvers_stationary": all(
            row["stationarity_inf"] <= max_stationarity_inf
            for row in solver_metrics.values()
        ),
        "both_solvers_complementary": all(
            row["complementarity_abs"] <= max_complementarity_abs
            for row in solver_metrics.values()
        ),
        "both_solvers_match_analytic_solution": all(
            row["solution_relative_error"] <= max_solution_relative_error
            and row["objective_relative_error"] <= max_solution_relative_error
            for row in solver_metrics.values()
        ),
        "independent_solver_pair_agrees": pair_difference
        <= max_solution_relative_error,
    }
    return {
        "policy": "energy_budgeted_trace_kkt_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "energy_budget": budget,
            "budget_fraction": budget_fraction,
            "objective_gradient_relative_error": objective_gradient_error,
            "constraint_gradient_relative_error": constraint_gradient_error,
            "solver_pair_solution_relative_difference": pair_difference,
            "solvers": solver_metrics,
        },
        "lesson": (
            "An energy-budgeted trace fit is not validated by a small objective alone. "
            "Check primal feasibility, a nonnegative active dual, stationarity, "
            "complementarity, finite-difference gradients, and agreement with an "
            "independent analytic KKT route."
        ),
    }
