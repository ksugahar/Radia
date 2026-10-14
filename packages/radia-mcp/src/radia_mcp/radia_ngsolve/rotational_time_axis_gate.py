"""Rotational kinematics gate for result-table time axes."""
from __future__ import annotations

import math


_DISPLAY_TIME_SCALES = {
    "s": 1.0,
    "ms": 1.0e-3,
    "us": 1.0e-6,
    "ns": 1.0e-9,
}


def rotational_kinematics_time_axis_gate(
    time_values,
    angles_deg,
    speeds_rpm,
    *,
    reported_time_unit: str,
    time_value_basis: str = "si_seconds",
    max_central_relative_error: float = 1.0e-8,
    min_sample_count: int = 5,
):
    """Check ``d(angle)/dt = 6 * speed_rpm`` on interior rows.

    Some result APIs return time values in SI seconds while separately exposing
    the graph's display unit.  ``time_value_basis`` makes that contract explicit
    instead of silently scaling values from the unit label.
    """
    time = [float(value) for value in time_values]
    angle = [float(value) for value in angles_deg]
    speed = [float(value) for value in speeds_rpm]
    if not (len(time) == len(angle) == len(speed)):
        raise ValueError("time, angle, and speed must have the same length")
    if min_sample_count < 5:
        raise ValueError("min_sample_count must be >= 5")
    if max_central_relative_error < 0.0:
        raise ValueError("max_central_relative_error must be >= 0")
    if time_value_basis not in {"si_seconds", "display_unit"}:
        raise ValueError("time_value_basis must be 'si_seconds' or 'display_unit'")
    unit = str(reported_time_unit).strip().lower()
    if unit not in _DISPLAY_TIME_SCALES:
        raise ValueError(f"unsupported reported_time_unit: {reported_time_unit!r}")

    finite = all(math.isfinite(value) for value in time + angle + speed)
    increasing = finite and all(right > left for left, right in zip(time, time[1:]))
    scale = 1.0 if time_value_basis == "si_seconds" else _DISPLAY_TIME_SCALES[unit]

    rows = []
    errors = []
    alternate_display_errors = []
    if finite and increasing:
        for index in range(1, len(time) - 1):
            delta_angle = angle[index + 1] - angle[index - 1]
            delta_time = (time[index + 1] - time[index - 1]) * scale
            implied_rpm = delta_angle / delta_time / 6.0
            denom = max(abs(implied_rpm), abs(speed[index]), 1.0e-30)
            relative_error = abs(implied_rpm - speed[index]) / denom
            errors.append(relative_error)

            display_delta_time = (
                time[index + 1] - time[index - 1]
            ) * _DISPLAY_TIME_SCALES[unit]
            display_implied_rpm = delta_angle / display_delta_time / 6.0
            display_denom = max(abs(display_implied_rpm), abs(speed[index]), 1.0e-30)
            alternate_display_errors.append(
                abs(display_implied_rpm - speed[index]) / display_denom
            )
            rows.append(
                {
                    "index": index,
                    "time_value": time[index],
                    "angle_deg": angle[index],
                    "reported_speed_rpm": speed[index],
                    "implied_speed_rpm": implied_rpm,
                    "relative_error": relative_error,
                }
            )

    max_error = max(errors) if errors else math.inf
    checks = {
        "sample_count_sufficient": len(time) >= min_sample_count,
        "all_finite": finite,
        "time_values_strictly_increase": increasing,
        "angle_nontrivial": finite and bool(angle) and max(angle) > min(angle),
        "speed_nontrivial": finite and bool(speed) and max(abs(value) for value in speed) > 0.0,
        "interior_rows_available": len(errors) >= 3,
        "rotational_kinematics_match": max_error <= max_central_relative_error,
        "time_value_basis_recorded": time_value_basis in {"si_seconds", "display_unit"},
        "reported_time_unit_recorded": bool(unit),
    }
    return {
        "policy": "rotational_kinematics_time_axis_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "sample_count": len(time),
        "central_sample_count": len(errors),
        "time_value_basis": time_value_basis,
        "reported_time_unit": reported_time_unit,
        "applied_time_scale_to_seconds": scale,
        "max_central_relative_error": max_error,
        "mean_central_relative_error": sum(errors) / len(errors) if errors else None,
        "display_unit_interpretation_relative_error": (
            max(alternate_display_errors) if alternate_display_errors else None
        ),
        "display_unit_label_is_metadata_only": (
            time_value_basis == "si_seconds" and unit != "s"
        ),
        "initial_and_final_rows_are_diagnostic_only": True,
        "checks": checks,
        "rows": rows,
        "lesson": (
            "Validate result-table time values by kinematics before applying a display-unit label. "
            "Keep the numeric value basis and the graph display unit as separate metadata."
        ),
    }
