"""Solver-neutral gate for a regularized FEM/BEM trace inverse path."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

from .slot_gates import lcurve_corner_choice, morozov_discrepancy_choice


def regularized_trace_inverse_path_gate(
    summary: Mapping[str, Any],
    *,
    max_solution_relative_error: float = 1.0e-10,
    max_trace_relative_error: float = 1.0e-10,
    max_regularized_objective_relative_error: float = 1.0e-10,
    max_zero_alpha_objective_absolute_error: float = 1.0e-20,
    max_normal_equation_residual: float = 1.0e-9,
    max_gradient_check_absolute_error: float = 2.0e-7,
    max_replay_relative_error: float = 1.0e-12,
    monotonicity_tolerance: float = 1.0e-12,
) -> dict[str, Any]:
    """Validate P1 trace regularization, parameter choices, and replay.

    The input is deliberately solver-neutral.  It records a first-order
    tetrahedron/triangle trace path, two regularization choices, independent
    linear-solver checks, and deterministic replay.  The gate recomputes both
    L-curve and Morozov choices instead of trusting reported indices.
    """

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    tolerances = {
        "max_solution_relative_error": _nonnegative(
            max_solution_relative_error, "max_solution_relative_error"
        ),
        "max_trace_relative_error": _nonnegative(
            max_trace_relative_error, "max_trace_relative_error"
        ),
        "max_regularized_objective_relative_error": _nonnegative(
            max_regularized_objective_relative_error,
            "max_regularized_objective_relative_error",
        ),
        "max_zero_alpha_objective_absolute_error": _nonnegative(
            max_zero_alpha_objective_absolute_error,
            "max_zero_alpha_objective_absolute_error",
        ),
        "max_normal_equation_residual": _nonnegative(
            max_normal_equation_residual, "max_normal_equation_residual"
        ),
        "max_gradient_check_absolute_error": _nonnegative(
            max_gradient_check_absolute_error, "max_gradient_check_absolute_error"
        ),
        "max_replay_relative_error": _nonnegative(
            max_replay_relative_error, "max_replay_relative_error"
        ),
        "monotonicity_tolerance": _nonnegative(
            monotonicity_tolerance, "monotonicity_tolerance"
        ),
    }

    mesh = _mapping(summary, "mesh")
    path = _mapping(summary, "path")
    problem = _mapping(summary, "problem")
    reported_lcurve = _mapping(summary, "lcurve")
    reported_morozov = _mapping(summary, "morozov")
    crosscheck = _mapping(summary, "crosscheck")
    replay = _mapping(summary, "replay")

    alphas = _float_list(path, "alphas")
    solution_norms = _float_list(path, "solution_norms")
    trace_residuals = _float_list(path, "trace_residual_norms")
    weighted_residuals = _float_list(path, "weighted_trace_residuals")
    normal_residuals = _float_list(path, "normal_equation_residuals")
    gradient_errors = _float_list(path, "gradient_check_max_abs_errors")
    lengths = {
        len(alphas),
        len(solution_norms),
        len(trace_residuals),
        len(weighted_residuals),
        len(normal_residuals),
        len(gradient_errors),
    }
    if len(lengths) != 1 or not alphas:
        raise ValueError("all regularization path arrays must have one nonzero common length")
    if len(alphas) < 5:
        raise ValueError("at least five regularization path rows are required")
    if any(value < 0.0 for value in alphas + solution_norms + trace_residuals + weighted_residuals):
        raise ValueError("path weights and norms must be non-negative")

    lcurve = lcurve_corner_choice(
        alphas[1:],
        trace_residuals[1:],
        solution_norms[1:],
        tol=tolerances["monotonicity_tolerance"],
    )
    morozov = morozov_discrepancy_choice(
        alphas,
        weighted_residuals,
        _finite_float(problem, "noise_norm"),
        tol=tolerances["monotonicity_tolerance"],
    )

    alpha_path_ok = alphas[0] == 0.0 and all(
        right > left for left, right in zip(alphas[1:], alphas[2:])
    )
    solution_monotone = all(
        right <= left + tolerances["monotonicity_tolerance"]
        for left, right in zip(solution_norms, solution_norms[1:])
    )
    trace_monotone = all(
        right + tolerances["monotonicity_tolerance"] >= left
        for left, right in zip(trace_residuals, trace_residuals[1:])
    )
    weighted_monotone = all(
        right + tolerances["monotonicity_tolerance"] >= left
        for left, right in zip(weighted_residuals, weighted_residuals[1:])
    )
    row_identity_ok = _optional_row_identity_is_aligned(path, len(alphas))
    gradient_resolution_ok = _optional_gradient_check_is_resolved(
        path, len(alphas)
    )
    gradient_generation_ok, solution_run_generation_ok = (
        _optional_generation_ids_are_aligned(path, len(alphas))
    )

    checks = {
        "schema_is_regularized_trace_inverse_v1": (
            str(summary.get("schema", "")) == "regularized_trace_inverse_path/v1"
        ),
        "mesh_is_first_order_tri_tet_trace": (
            str(mesh.get("volume_element", "")) == "tetrahedron"
            and str(mesh.get("boundary_element", "")) == "triangle"
            and _integer(mesh, "polynomial_order") == 1
            and _positive_integer(mesh, "tetrahedra") > 0
            and _positive_integer(mesh, "triangles") > 0
            and _positive_integer(mesh, "surface_nodes")
            < _positive_integer(mesh, "volume_nodes")
            and _positive_integer(mesh, "trace_rows")
            == _positive_integer(mesh, "surface_nodes")
            and _positive_integer(mesh, "fem_unknowns")
            == _positive_integer(mesh, "volume_nodes")
            and _positive_integer(mesh, "trace_nnz")
            == _positive_integer(mesh, "surface_nodes")
        ),
        "zero_then_strictly_increasing_alpha_path": alpha_path_ok,
        "solution_norm_decreases_along_path": solution_monotone,
        "trace_residual_increases_along_path": trace_monotone,
        "weighted_residual_increases_along_path": weighted_monotone,
        "normal_equations_close": (
            max(normal_residuals) <= tolerances["max_normal_equation_residual"]
        ),
        "finite_difference_gradients_close": (
            max(gradient_errors) <= tolerances["max_gradient_check_absolute_error"]
        ),
        "alpha_path_row_identity_is_aligned": row_identity_ok,
        "finite_difference_steps_are_numerically_resolved": gradient_resolution_ok,
        "gradient_uses_current_parameter_generation": gradient_generation_ok,
        "solution_rows_share_regularization_run_generation": (
            solution_run_generation_ok
        ),
        "lcurve_recomputation_passes": lcurve["status"] == "ok",
        "reported_lcurve_choice_matches": (
            _integer(reported_lcurve, "selected_index") == lcurve["selected_index"]
            and _close(
                _finite_float(reported_lcurve, "selected_alpha"),
                lcurve["selected_alpha"],
            )
        ),
        "morozov_recomputation_passes": morozov["status"] == "ok",
        "reported_morozov_choice_matches": (
            _integer(reported_morozov, "selected_index") == morozov["selected_index"]
            and _close(
                _finite_float(reported_morozov, "selected_alpha"),
                morozov["selected_alpha"],
            )
        ),
        "two_independent_linear_references_close": (
            _positive_integer(crosscheck, "reference_solver_count") >= 2
            and _finite_float(crosscheck, "max_solution_relative_error")
            <= tolerances["max_solution_relative_error"]
            and _finite_float(crosscheck, "max_trace_relative_error")
            <= tolerances["max_trace_relative_error"]
            and _finite_float(
                crosscheck, "max_regularized_objective_relative_error"
            )
            <= tolerances["max_regularized_objective_relative_error"]
            and _finite_float(crosscheck, "zero_alpha_objective_absolute_error")
            <= tolerances["max_zero_alpha_objective_absolute_error"]
        ),
        "deterministic_replay_closes": (
            _positive_integer(replay, "count") >= 2
            and bool(replay.get("selectors_identical", False))
            and _finite_float(replay, "max_relative_error")
            <= tolerances["max_replay_relative_error"]
        ),
    }
    issues = [name for name, passed in checks.items() if not passed]
    return {
        "policy": "regularized_trace_inverse_path_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "lcurve": lcurve,
        "morozov": morozov,
        "metrics": {
            "path_count": len(alphas),
            "max_normal_equation_residual": max(normal_residuals),
            "max_gradient_check_absolute_error": max(gradient_errors),
            "max_solution_relative_error": _finite_float(
                crosscheck, "max_solution_relative_error"
            ),
            "max_replay_relative_error": _finite_float(
                replay, "max_relative_error"
            ),
        },
        "tolerances": tolerances,
        "lesson": (
            "Promote a trace inverse path only after the zero-weight minimum-norm "
            "limit, monotone Tikhonov trade-off, independently recomputed L-curve "
            "and Morozov choices, two linear references, and replay all close."
        ),
    }


def _mapping(parent: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    value = parent.get(key)
    if not isinstance(value, Mapping):
        raise ValueError(f"{key} must be a mapping")
    return value


def _float_list(parent: Mapping[str, Any], key: str) -> list[float]:
    value = parent.get(key)
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{key} must be an array")
    out = [float(item) for item in value]
    if not all(math.isfinite(item) for item in out):
        raise ValueError(f"{key} must contain only finite values")
    return out


def _finite_float(parent: Mapping[str, Any], key: str) -> float:
    value = float(parent[key])
    if not math.isfinite(value):
        raise ValueError(f"{key} must be finite")
    return value


def _optional_row_identity_is_aligned(
    path: Mapping[str, Any], expected_length: int
) -> bool:
    names = ("alpha_row_ids", "solution_row_ids", "residual_row_ids")
    if not any(name in path for name in names):
        return True
    rows = []
    for name in names:
        value = path.get(name)
        if (
            not isinstance(value, Sequence)
            or isinstance(value, (str, bytes))
            or len(value) != expected_length
        ):
            return False
        row = [str(item).strip() for item in value]
        if not all(row) or len(set(row)) != expected_length:
            return False
        rows.append(row)
    return rows[0] == rows[1] == rows[2]


def _optional_gradient_check_is_resolved(
    path: Mapping[str, Any], expected_length: int
) -> bool:
    names = (
        "gradient_check_step_sizes",
        "gradient_check_parameter_scales",
        "gradient_check_objective_pair_deltas",
    )
    if not any(name in path for name in names):
        return True
    try:
        steps, scales, objective_deltas = (
            _float_list(path, name) for name in names
        )
    except (KeyError, TypeError, ValueError):
        return False
    if any(
        len(row) != expected_length
        for row in (steps, scales, objective_deltas)
    ):
        return False
    resolution = math.sqrt(math.ulp(1.0))
    return all(
        scale > 0.0
        and step >= resolution * max(scale, 1.0)
        and abs(objective_delta) > 0.0
        for step, scale, objective_delta in zip(steps, scales, objective_deltas)
    )


def _optional_generation_ids_are_aligned(
    path: Mapping[str, Any], expected_length: int
) -> tuple[bool, bool]:
    pairs = (
        ("parameter_generation_ids", "gradient_parameter_generation_ids"),
        ("path_run_generation_ids", "solution_run_generation_ids"),
    )
    results = []
    for left_name, right_name in pairs:
        if left_name not in path and right_name not in path:
            results.append(True)
            continue
        rows = []
        for name in (left_name, right_name):
            value = path.get(name)
            if (
                not isinstance(value, Sequence)
                or isinstance(value, (str, bytes))
                or len(value) != expected_length
            ):
                rows = []
                break
            row = [str(item).strip() for item in value]
            if not all(row):
                rows = []
                break
            rows.append(row)
        results.append(len(rows) == 2 and rows[0] == rows[1])
    return results[0], results[1]


def _integer(parent: Mapping[str, Any], key: str) -> int:
    value = parent[key]
    if isinstance(value, bool) or int(value) != float(value):
        raise ValueError(f"{key} must be an integer")
    return int(value)


def _positive_integer(parent: Mapping[str, Any], key: str) -> int:
    value = _integer(parent, key)
    if value <= 0:
        raise ValueError(f"{key} must be positive")
    return value


def _nonnegative(value: float, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result < 0.0:
        raise ValueError(f"{name} must be finite and non-negative")
    return result


def _close(left: float, right: float, tol: float = 1.0e-12) -> bool:
    return abs(left - right) <= tol * max(abs(left), abs(right), 1.0)
