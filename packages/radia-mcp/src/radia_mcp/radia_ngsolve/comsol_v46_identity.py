"""Neutral COMSOL-derived replay identity checks for v46 public artifacts.

The records are optional so older solver summaries retain their v1-v45
behavior. When present, a record must close the accepted partial/restarted
state or the field unit/frame convention all the way to the result digest.
"""

from __future__ import annotations

import math
from collections.abc import Mapping


_TIME = "time_adaptive_partial_solution_nan_inf_restart_window_identity"
_FIELD = "unit_scale_coordinate_frame_complex_field_vector_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _same(row: Mapping[str, object], *names: str) -> bool:
    return all(row.get(f"result_{name}") == row.get(name) for name in names)


def _time_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation", "")).strip()
    restart_window = row.get("restart_window_s")
    accepted_steps = row.get("accepted_step_count")
    result_steps = row.get("result_accepted_step_count")
    return (
        bool(generation)
        and row.get("adaptive_step_generation") == generation
        and row.get("restart_window_generation") == generation
        and row.get("finite_value_generation") == generation
        and row.get("result_generation") == generation
        and isinstance(restart_window, list)
        and len(restart_window) == 2
        and all(isinstance(value, (int, float)) and math.isfinite(float(value)) for value in restart_window)
        and 0.0 <= float(restart_window[0]) < float(restart_window[1])
        and row.get("result_restart_window_s") == restart_window
        and isinstance(accepted_steps, int)
        and accepted_steps > 0
        and result_steps == accepted_steps
        and row.get("adaptive_step_policy") in {"accepted_only", "bdf_adaptive_accepted_only"}
        and row.get("result_adaptive_step_policy") == row.get("adaptive_step_policy")
        and row.get("finite_value_status") == row.get("result_finite_value_status") == "finite"
        and row.get("restart_state") == row.get("result_restart_state") == "partial_restartable"
        and bool(str(row.get("owner") or ""))
        and row.get("accepted_owner") == row.get("owner")
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def _field_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation", "")).strip()
    scale = row.get("unit_scale_to_si")
    return (
        bool(generation)
        and row.get("unit_scale_generation") == generation
        and row.get("coordinate_frame_generation") == generation
        and row.get("complex_vector_generation") == generation
        and row.get("result_generation") == generation
        and isinstance(scale, (int, float))
        and math.isfinite(float(scale))
        and float(scale) > 0.0
        and row.get("result_unit_scale_to_si") == scale
        and row.get("coordinate_frame") in {"global_cartesian", "global_cylindrical", "local_right_handed"}
        and row.get("result_coordinate_frame") == row.get("coordinate_frame")
        and row.get("complex_vector_convention") in {"real_imag_components", "phasor_real_imag"}
        and row.get("result_complex_vector_convention") == row.get("complex_vector_convention")
        and _same(row, "unit_name", "coordinate_frame", "complex_vector_convention")
        and bool(str(row.get("owner") or ""))
        and row.get("accepted_owner") == row.get("owner")
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def validate_public_identity(payload: object) -> dict[str, object]:
    """Return v46 checks; an absent v46 record is a compatibility pass."""

    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    time_record = payload.get(_TIME)
    field_record = payload.get(_FIELD)
    if time_record is not None:
        checks["v46_time_partial_restart_identity"] = isinstance(time_record, Mapping) and _time_ok(time_record)
    if field_record is not None:
        checks["v46_field_unit_frame_identity"] = isinstance(field_record, Mapping) and _field_ok(field_record)
    if not checks:
        return {}
    return {
        "policy": "comsol_v46_public_identity_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
    }
