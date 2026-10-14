"""Solver-neutral contracts for prescribed rotating-conductor transients."""

from __future__ import annotations

import math
from typing import Any


def _rows(value: Any, name: str) -> list[tuple[float, float]]:
    if not isinstance(value, list):
        raise ValueError(f"{name} must be a list")
    result = []
    for index, row in enumerate(value):
        if not isinstance(row, (list, tuple)) or len(row) != 2:
            raise ValueError(f"{name}[{index}] must contain time and value")
        pair = (float(row[0]), float(row[1]))
        if not all(math.isfinite(item) for item in pair):
            raise ValueError(f"{name}[{index}] must be finite")
        result.append(pair)
    return result


def _loss_rows(value: Any) -> list[dict[str, Any]]:
    if not isinstance(value, list):
        raise ValueError("loss_rows must be a list")
    result = []
    for index, row in enumerate(value):
        if not isinstance(row, dict):
            raise ValueError(f"loss_rows[{index}] must be a mapping")
        time_s = float(row["time_s"])
        total_w = float(row["total_w"])
        parts_w = [float(item) for item in row["parts_w"]]
        if not parts_w or not all(
            math.isfinite(item) for item in [time_s, total_w, *parts_w]
        ):
            raise ValueError(f"loss_rows[{index}] must be finite with loss parts")
        result.append({"time_s": time_s, "total_w": total_w, "parts_w": parts_w})
    return result


def _relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(abs(actual), abs(expected), 1.0e-300)


def rotating_conductor_transient_gate(
    summary: dict[str, Any],
    *,
    min_sample_count: int = 20,
    max_time_alignment_error_s: float = 1.0e-12,
    max_kinematic_relative_error: float = 1.0e-8,
    max_loss_partition_relative_error: float = 1.0e-10,
) -> dict[str, Any]:
    """Gate a migrated moving-mesh run without asserting a false power law.

    Angle and torque include the initial row.  Speed, loss and current-flux
    tables use right-endpoint rows, so their time axes must equal
    ``angle_time[1:]``.  The prescribed rotation identity is
    ``delta(theta) = 2*pi*speed_rpm*delta(t)/60``.
    """
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    boundary = summary.get("moving_axis_boundary")
    units = summary.get("units")
    diagnostic = summary.get("energy_balance_contract")
    if not all(isinstance(value, dict) for value in (boundary, units, diagnostic)):
        raise ValueError("moving_axis_boundary, units and energy_balance_contract must be mappings")

    tolerances = {
        "min_sample_count": int(min_sample_count),
        "max_time_alignment_error_s": float(max_time_alignment_error_s),
        "max_kinematic_relative_error": float(max_kinematic_relative_error),
        "max_loss_partition_relative_error": float(max_loss_partition_relative_error),
    }
    if tolerances["min_sample_count"] < 3 or any(
        not math.isfinite(value) or value < 0.0
        for key, value in tolerances.items()
        if key != "min_sample_count"
    ):
        raise ValueError("tolerances must be finite and nonnegative")

    angle = _rows(summary.get("angle_rows"), "angle_rows")
    speed = _rows(summary.get("speed_rows"), "speed_rows")
    torque = _rows(summary.get("torque_rows"), "torque_rows")
    current_flux = _rows(summary.get("current_flux_rows"), "current_flux_rows")
    losses = _loss_rows(summary.get("loss_rows"))
    angle_time = [row[0] for row in angle]
    endpoint_time = [row[0] for row in speed]
    loss_time = [row["time_s"] for row in losses]

    increasing = all(right > left for left, right in zip(angle_time, angle_time[1:]))
    steps = [right - left for left, right in zip(angle_time, angle_time[1:])]
    mean_step = sum(steps) / len(steps) if steps else math.inf
    step_spread = (
        (max(steps) - min(steps)) / mean_step if steps and mean_step > 0.0 else math.inf
    )
    alignment_errors = [
        abs(left - right) for left, right in zip(angle_time[1:], endpoint_time)
    ]
    loss_alignment_errors = [
        abs(left - right) for left, right in zip(endpoint_time, loss_time)
    ]
    flux_alignment_errors = [
        abs(left - right)
        for left, right in zip(endpoint_time, [row[0] for row in current_flux])
    ]
    torque_alignment_errors = [
        abs(left - right) for left, right in zip(angle_time, [row[0] for row in torque])
    ]

    kinematic_errors = []
    if len(angle) == len(speed) + 1:
        for index, (_, speed_rpm) in enumerate(speed, start=1):
            dt = angle[index][0] - angle[index - 1][0]
            actual = angle[index][1] - angle[index - 1][1]
            expected = 2.0 * math.pi * speed_rpm * dt / 60.0
            kinematic_errors.append(_relative_error(actual, expected))
    loss_errors = [
        _relative_error(row["total_w"], sum(row["parts_w"])) for row in losses
    ]
    mechanical_power = []
    if len(torque) == len(angle):
        for index, (_, speed_rpm) in enumerate(speed, start=1):
            mechanical_power.append(torque[index][1] * 2.0 * math.pi * speed_rpm / 60.0)

    boundary_before = str(boundary.get("before") or "").strip().lower()
    boundary_after = str(boundary.get("after") or "").strip().lower()
    checks = {
        "moving_axis_recorded": str(boundary.get("axis") or "").strip().lower()
        in {"x", "y", "z"},
        "legacy_open_axis_boundary_recorded": boundary_before == "open",
        "moving_axis_boundary_migrated_closed": boundary_after in {"magnetic", "electric"},
        "source_model_left_immutable": boundary.get("source_modified") is False,
        "si_units_explicit": units
        == {
            "time": "s",
            "angle": "rad",
            "speed": "r/min",
            "torque": "N*m",
            "loss": "W",
            "current_flux": "A",
        },
        "sample_count_sufficient": len(speed) >= tolerances["min_sample_count"],
        "table_cardinality_contract": len(angle) == len(torque) == len(speed) + 1
        and len(speed) == len(losses) == len(current_flux),
        "angle_time_strictly_increases": increasing,
        "time_step_uniform": step_spread <= 1.0e-10,
        "right_endpoint_tables_align": bool(alignment_errors)
        and max(alignment_errors + loss_alignment_errors + flux_alignment_errors)
        <= tolerances["max_time_alignment_error_s"],
        "initial_row_tables_align": bool(torque_alignment_errors)
        and max(torque_alignment_errors) <= tolerances["max_time_alignment_error_s"],
        "angle_matches_right_endpoint_speed": bool(kinematic_errors)
        and max(kinematic_errors) <= tolerances["max_kinematic_relative_error"],
        "loss_partition_closes": bool(loss_errors)
        and max(loss_errors) <= tolerances["max_loss_partition_relative_error"],
        "losses_nonnegative": bool(losses)
        and all(row["total_w"] >= 0.0 and min(row["parts_w"]) >= 0.0 for row in losses),
        "torque_and_current_flux_finite_nontrivial": bool(torque and current_flux)
        and max(abs(row[1]) for row in torque) > 0.0
        and max(abs(row[1]) for row in current_flux) > 0.0,
        "mechanical_power_vs_joule_loss_is_diagnostic": diagnostic.get(
            "mechanical_power_vs_joule_loss"
        )
        == "diagnostic_only"
        and diagnostic.get("external_drive_power_available") is False,
    }
    return {
        "policy": "rotating_conductor_transient_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "angle_row_count": len(angle),
            "right_endpoint_row_count": len(speed),
            "time_step_s": mean_step,
            "time_step_relative_spread": step_spread,
            "maximum_time_alignment_error_s": max(
                alignment_errors + loss_alignment_errors + flux_alignment_errors,
                default=math.inf,
            ),
            "maximum_kinematic_relative_error": max(kinematic_errors, default=math.inf),
            "maximum_loss_partition_relative_error": max(loss_errors, default=math.inf),
            "final_angle_rad": angle[-1][1] if angle else None,
            "mean_speed_rpm": sum(row[1] for row in speed) / len(speed) if speed else None,
            "peak_torque_abs_nm": max((abs(row[1]) for row in torque), default=0.0),
            "peak_total_loss_w": max((row["total_w"] for row in losses), default=0.0),
            "peak_current_flux_abs_a": max((abs(row[1]) for row in current_flux), default=0.0),
            "mechanical_power_abs_w_diagnostic": {
                "minimum": min((abs(value) for value in mechanical_power), default=0.0),
                "mean": (
                    sum(abs(value) for value in mechanical_power) / len(mechanical_power)
                    if mechanical_power
                    else 0.0
                ),
                "maximum": max((abs(value) for value in mechanical_power), default=0.0),
            },
        },
        "tolerances": tolerances,
        "notes": [
            "Moving-axis boundary compatibility is a solver contract separate from field accuracy.",
            "Right-endpoint speed rows explain the one-row offset from angle and torque tables.",
            "Do not equate torque times speed with Joule loss unless source power and magnetic-energy rate are both available.",
        ],
    }
