"""Solver-state and material lineage checks for COMSOL-derived v49 summaries."""

from __future__ import annotations

import math
from collections.abc import Mapping


_MATERIAL = "nonlinear_material_interpolation_branch_unit_temperature_extrapolation_owner_identity"
_SLIDING = "moving_mesh_sliding_interface_frame_time_remesh_solution_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation_closed(row: Mapping[str, object], *names: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _result_identity_ok(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _nonempty_string_mapping(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and bool(value)
        and all(isinstance(key, str) and key and isinstance(item, str) and item for key, item in value.items())
    )


def _material_ok(row: Mapping[str, object]) -> bool:
    units = row.get("input_units")
    temperature = row.get("temperature_value")
    return (
        _generation_closed(
            row,
            "material_generation",
            "branch_generation",
            "temperature_generation",
            "interpolation_generation",
            "result_generation",
        )
        and bool(str(row.get("interpolation_branch") or ""))
        and row.get("result_interpolation_branch") == row.get("interpolation_branch")
        and _nonempty_string_mapping(units)
        and row.get("result_input_units") == units
        and isinstance(temperature, (int, float))
        and math.isfinite(float(temperature))
        and float(temperature) > 0.0
        and row.get("result_temperature_value") == temperature
        and row.get("extrapolation_policy") in {"reject", "constant", "linear"}
        and row.get("result_extrapolation_policy") == row.get("extrapolation_policy")
        and str(row.get("material_owner") or "").startswith("material:")
        and row.get("result_material_owner") == row.get("material_owner")
        and _result_identity_ok(row)
    )


def _sliding_ok(row: Mapping[str, object]) -> bool:
    time_value = row.get("time_value_s")
    return (
        _generation_closed(
            row,
            "mesh_generation",
            "interface_generation",
            "frame_generation",
            "time_generation",
            "remesh_generation",
            "solution_generation",
            "result_generation",
        )
        and _digest(row.get("sliding_interface_map_sha256"))
        and row.get("result_sliding_interface_map_sha256") == row.get("sliding_interface_map_sha256")
        and row.get("coordinate_frame") in {"spatial", "material"}
        and row.get("result_coordinate_frame") == row.get("coordinate_frame")
        and isinstance(time_value, (int, float))
        and math.isfinite(float(time_value))
        and float(time_value) >= 0.0
        and row.get("result_time_value_s") == time_value
        and bool(str(row.get("remesh_revision") or ""))
        and row.get("result_remesh_revision") == row.get("remesh_revision")
        and str(row.get("solution_owner") or "").startswith("solution:")
        and row.get("result_solution_owner") == row.get("solution_owner")
        and _result_identity_ok(row)
    )


def validate_public_v49_identity(payload: object) -> dict[str, object]:
    """Validate optional nonlinear-material and moving-mesh identity records."""
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    material = payload.get(_MATERIAL)
    sliding = payload.get(_SLIDING)
    if material is not None:
        checks["v49_nonlinear_material_branch_unit_temperature_extrapolation_owner"] = (
            isinstance(material, Mapping) and _material_ok(material)
        )
    if sliding is not None:
        checks["v49_sliding_interface_frame_time_remesh_solution_owner"] = (
            isinstance(sliding, Mapping) and _sliding_ok(sliding)
        )
    if not checks:
        return {}
    return {
        "policy": "solver_state_identity_v49",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
    }

