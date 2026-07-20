"""Transform and solver-normalization lineage checks for v48 summaries."""

from __future__ import annotations

import math
from collections.abc import Mapping


_ALE = "ale_reference_current_force_quadrature_owner_identity"
_SEGREGATED = "segregated_variable_scaling_residual_iteration_solution_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation_closed(row: Mapping[str, object], *names: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _result_identity_ok(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _ale_ok(row: Mapping[str, object]) -> bool:
    return (
        _generation_closed(
            row,
            "reference_mesh_generation",
            "current_mesh_generation",
            "quadrature_generation",
            "normal_generation",
            "result_generation",
        )
        and str(row.get("reference_configuration_id") or "").startswith("ale/reference-")
        and row.get("result_reference_configuration_id") == row.get("reference_configuration_id")
        and str(row.get("current_configuration_id") or "").startswith("ale/current-")
        and row.get("result_current_configuration_id") == row.get("current_configuration_id")
        and str(row.get("quadrature_rule") or "").startswith("gauss-surface-")
        and row.get("result_quadrature_rule") == row.get("quadrature_rule")
        and _digest(row.get("normal_orientation_sha256"))
        and row.get("result_normal_orientation_sha256") == row.get("normal_orientation_sha256")
        and str(row.get("body_owner") or "").startswith("body:")
        and row.get("result_body_owner") == row.get("body_owner")
        and _result_identity_ok(row)
    )


def _unique_strings(value: object) -> bool:
    return isinstance(value, list) and bool(value) and all(isinstance(item, str) and item for item in value) and len(value) == len(set(value))


def _segregated_ok(row: Mapping[str, object]) -> bool:
    groups = row.get("variable_groups")
    scaling = row.get("variable_scaling")
    iterations = row.get("iteration_rows")
    return (
        _generation_closed(
            row,
            "variable_group_generation",
            "scaling_generation",
            "residual_generation",
            "iteration_generation",
            "solution_generation",
            "result_generation",
        )
        and _unique_strings(groups)
        and row.get("result_variable_groups") == groups
        and isinstance(scaling, Mapping)
        and set(scaling) == set(groups)
        and all(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) and float(value) > 0.0 for value in scaling.values())
        and row.get("result_variable_scaling") == scaling
        and row.get("residual_norm") == row.get("result_residual_norm") == "scaled_l2"
        and _unique_strings(iterations)
        and len(iterations) >= len(groups)
        and row.get("result_iteration_rows") == iterations
        and str(row.get("solution_owner") or "").startswith("solution:")
        and row.get("result_solution_owner") == row.get("solution_owner")
        and _result_identity_ok(row)
    )


def validate_public_v48_identity(payload: object) -> dict[str, object]:
    """Validate optional v48 ALE and segregated-solver records."""
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    ale = payload.get(_ALE)
    segregated = payload.get(_SEGREGATED)
    if ale is not None:
        checks["v48_ale_configuration_quadrature_normal_owner"] = isinstance(ale, Mapping) and _ale_ok(ale)
    if segregated is not None:
        checks["v48_segregated_scaling_residual_iteration_owner"] = isinstance(segregated, Mapping) and _segregated_ok(segregated)
    if not checks:
        return {}
    return {
        "policy": "transform_normalization_v48",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
    }
