"""Audit derivative-free convergence reports against stationarity evidence."""
from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence


def _number(value: object, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _vector(value: object, name: str) -> list[float]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        raise ValueError(f"{name} must be a non-empty array")
    return [_number(item, name) for item in value]


def _norm(values: Sequence[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def evaluate_simplex_stationarity_audit(summary: Mapping[str, object]) -> dict:
    """Classify reported convergence using independent first-order checks.

    A small simplex or objective spread is not a stationarity certificate.
    Every reported result is compared with an independently evaluated
    gradient and, when available, a trusted reference solution.
    """
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    reference = summary.get("reference")
    methods = summary.get("methods")
    tolerances = summary.get("tolerances")
    if not isinstance(reference, Mapping) or not isinstance(tolerances, Mapping):
        raise ValueError("reference and tolerances must be objects")
    if not isinstance(methods, Sequence) or isinstance(methods, (str, bytes)):
        raise ValueError("methods must be an array")
    reference_x = _vector(reference.get("x"), "reference.x")
    reference_f = _number(reference.get("objective"), "reference.objective")
    stationarity_tolerance = _number(tolerances.get("gradient_norm"), "tolerances.gradient_norm")
    objective_tolerance = _number(tolerances.get("objective_gap"), "tolerances.objective_gap")
    parameter_tolerance = _number(tolerances.get("parameter_distance"), "tolerances.parameter_distance")
    if min(stationarity_tolerance, objective_tolerance, parameter_tolerance) < 0.0:
        raise ValueError("tolerances must be nonnegative")

    rows = []
    ids: list[str] = []
    roles: list[str] = []
    for index, method in enumerate(methods):
        if not isinstance(method, Mapping):
            raise ValueError(f"methods[{index}] must be an object")
        method_id = str(method.get("method_id") or "").strip()
        if not method_id:
            raise ValueError(f"methods[{index}].method_id is required")
        ids.append(method_id)
        role = str(method.get("role") or "").strip().lower()
        if role not in {"candidate", "control"}:
            raise ValueError(f"methods[{index}].role must be candidate or control")
        roles.append(role)
        x = _vector(method.get("x"), f"methods[{index}].x")
        gradient = _vector(method.get("gradient"), f"methods[{index}].gradient")
        if len(x) != len(reference_x) or len(gradient) != len(reference_x):
            raise ValueError("method x/gradient dimensions must match the reference")
        objective = _number(method.get("objective"), f"methods[{index}].objective")
        try:
            evaluations = int(method["function_evaluations"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"methods[{index}].function_evaluations must be an integer") from exc
        gradient_norm = _norm(gradient)
        objective_gap = abs(objective - reference_f)
        parameter_distance = _norm([value - target for value, target in zip(x, reference_x)])
        reported = method.get("reported_converged") is True
        accepted = (
            reported
            and evaluations > 0
            and gradient_norm <= stationarity_tolerance
            and objective_gap <= objective_tolerance
            and parameter_distance <= parameter_tolerance
        )
        reasons = []
        if not reported:
            reasons.append("solver_did_not_report_convergence")
        if evaluations <= 0:
            reasons.append("function_evaluations_not_positive")
        if gradient_norm > stationarity_tolerance:
            reasons.append("gradient_norm_exceeds_tolerance")
        if objective_gap > objective_tolerance:
            reasons.append("objective_gap_exceeds_tolerance")
        if parameter_distance > parameter_tolerance:
            reasons.append("parameter_distance_exceeds_tolerance")
        rows.append({
            "method_id": method_id,
            "role": role,
            "reported_converged": reported,
            "accepted": accepted,
            "function_evaluations": evaluations,
            "gradient_norm": gradient_norm,
            "objective_gap": objective_gap,
            "parameter_distance": parameter_distance,
            "rejection_reasons": reasons,
        })

    accepted_ids = [row["method_id"] for row in rows if row["accepted"]]
    accepted_control_ids = [
        row["method_id"] for row in rows if row["role"] == "control" and row["accepted"]
    ]
    false_ids = [row["method_id"] for row in rows if row["reported_converged"] and not row["accepted"]]
    false_candidate_ids = [
        row["method_id"]
        for row in rows
        if row["role"] == "candidate" and row["reported_converged"] and not row["accepted"]
    ]
    checks = {
        "units_explicit": bool(str(summary.get("parameter_unit") or "").strip()) and bool(str(summary.get("objective_unit") or "").strip()),
        "at_least_two_methods": len(rows) >= 2,
        "method_ids_unique": len(set(ids)) == len(ids),
        "candidate_and_control_roles_present": {"candidate", "control"}.issubset(set(roles)),
        "accepted_independent_control_present": bool(accepted_control_ids),
        "candidate_false_convergence_detected": bool(false_candidate_ids),
        "every_reported_result_has_secondary_checks": all(row["function_evaluations"] > 0 for row in rows),
    }
    return {
        "schema": "radia-topology-simplex-stationarity-audit/v1",
        "policy": "simplex_or_objective_spread_is_not_a_stationarity_certificate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "accepted_method_ids": accepted_ids,
        "accepted_control_method_ids": accepted_control_ids,
        "false_convergence_method_ids": false_ids,
        "false_convergence_candidate_ids": false_candidate_ids,
        "methods": rows,
        "reference": {"x": reference_x, "objective": reference_f},
        "tolerances": {
            "gradient_norm": stationarity_tolerance,
            "objective_gap": objective_tolerance,
            "parameter_distance": parameter_tolerance,
        },
    }


def simplex_stationarity_audit_gate(summary_json: str) -> str:
    try:
        summary = json.loads(summary_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"summary_json must be valid JSON: {exc.msg}") from exc
    return json.dumps(evaluate_simplex_stationarity_audit(summary), indent=2, sort_keys=True)
