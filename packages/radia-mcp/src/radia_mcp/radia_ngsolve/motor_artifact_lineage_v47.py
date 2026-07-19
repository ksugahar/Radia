"""Solver-neutral motor transform and aggregate lineage checks."""

from __future__ import annotations

import math
from collections.abc import Mapping


DQ = "v47_public_dq_phase_order_electrical_angle_pole_pair_mapping_mismatch"
WINDOW = "v47_public_torque_loss_integration_window_parameter_row_key_mismatch"


def _sha(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation(row: Mapping[str, object], names: tuple[str, ...]) -> bool:
    value = str(row.get("generation") or "")
    return bool(value) and all(row.get(name) == value for name in names)


def _digest(row: Mapping[str, object]) -> bool:
    return _sha(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _dq_ok(row: Mapping[str, object]) -> bool:
    phases = row.get("phase_order")
    pole_pairs = row.get("pole_pairs")
    angle = row.get("electrical_angle_origin_deg")
    return (
        _generation(
            row,
            (
                "dq_generation",
                "phase_order_generation",
                "electrical_angle_generation",
                "pole_pair_generation",
                "transform_generation",
                "result_generation",
            ),
        )
        and isinstance(phases, list)
        and phases == ["A", "B", "C"]
        and row.get("result_phase_order") == phases
        and isinstance(pole_pairs, int)
        and pole_pairs > 0
        and row.get("result_pole_pairs") == pole_pairs
        and isinstance(angle, (int, float))
        and math.isfinite(float(angle))
        and row.get("result_electrical_angle_origin_deg") == angle
        and row.get("transform_identity") == row.get("result_transform_identity") == "power_invariant_park"
        and row.get("angle_direction") == row.get("result_angle_direction") == "electrical_ccw"
        and _digest(row)
    )


def _window_ok(row: Mapping[str, object]) -> bool:
    window = row.get("integration_window_s")
    torque = row.get("torque_mean_nm")
    loss = row.get("loss_total_w")
    return (
        _generation(
            row,
            (
                "integration_window_generation",
                "parameter_row_generation",
                "torque_generation",
                "loss_generation",
                "result_generation",
            ),
        )
        and isinstance(window, list)
        and len(window) == 2
        and all(isinstance(value, (int, float)) and math.isfinite(float(value)) for value in window)
        and float(window[0]) < float(window[1])
        and row.get("result_integration_window_s") == window
        and isinstance(row.get("parameter_row_key"), str)
        and bool(row.get("parameter_row_key"))
        and row.get("result_parameter_row_key") == row.get("parameter_row_key")
        and isinstance(torque, (int, float))
        and math.isfinite(float(torque))
        and row.get("result_torque_mean_nm") == torque
        and isinstance(loss, (int, float))
        and math.isfinite(float(loss))
        and row.get("result_loss_total_w") == loss
        and str(row.get("study_owner") or "").startswith("study:")
        and row.get("result_study_owner") == row.get("study_owner")
        and _digest(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    dq = identity.get(DQ)
    window = identity.get(WINDOW)
    if dq is not None:
        checks["motor_v47_dq_phase_angle_pole_pair_transform"] = isinstance(dq, Mapping) and _dq_ok(dq)
    if window is not None:
        checks["motor_v47_torque_loss_window_parameter_row"] = isinstance(window, Mapping) and _window_ok(window)
    return checks
