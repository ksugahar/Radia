"""Solver-neutral range-angle localization gate for wideband array data."""
from __future__ import annotations

import math
from statistics import median
from typing import Any


SPEED_OF_LIGHT_M_S = 299_792_458.0


def radar_range_angle_localization_gate(
    frequency_hz: list[float],
    targets: list[dict[str, float]],
    *,
    max_range_resolution_multiples: float = 1.0,
    max_angle_error_deg: float = 2.0,
    max_frequency_step_relative_drift: float = 1.0e-7,
) -> dict[str, Any]:
    """Gate target localization against physical range and angular tolerances."""

    frequency = [float(value) for value in frequency_hz]
    if len(frequency) < 3 or not all(math.isfinite(value) and value > 0.0 for value in frequency):
        raise ValueError("frequency_hz must contain at least three finite positive samples")
    if not isinstance(targets, list) or len(targets) < 2:
        raise ValueError("targets must contain at least two localization rows")
    if max_range_resolution_multiples <= 0.0 or max_angle_error_deg <= 0.0:
        raise ValueError("localization tolerances must be positive")

    steps = [b - a for a, b in zip(frequency, frequency[1:])]
    increasing = all(step > 0.0 for step in steps)
    representative_step_hz = median(steps)
    step_drift = (
        max(abs(step - representative_step_hz) for step in steps) / representative_step_hz
        if increasing else math.inf
    )
    bandwidth_hz = frequency[-1] - frequency[0]
    resolution_m = SPEED_OF_LIGHT_M_S / (2.0 * bandwidth_hz) if bandwidth_hz > 0.0 else math.inf
    unambiguous_range_m = SPEED_OF_LIGHT_M_S / (2.0 * representative_step_hz) if representative_step_hz > 0.0 else 0.0

    target_metrics = []
    ids = []
    rows_valid = True
    for index, row in enumerate(targets):
        if not isinstance(row, dict):
            raise ValueError("each target row must be an object")
        target_id = str(row.get("target_id", "")).strip()
        values = {
            key: float(row[key])
            for key in ("expected_range_m", "estimated_range_m", "expected_angle_deg", "estimated_angle_deg")
        }
        if not target_id or not all(math.isfinite(value) for value in values.values()):
            rows_valid = False
        range_error = abs(values["estimated_range_m"] - values["expected_range_m"])
        angle_error = abs(values["estimated_angle_deg"] - values["expected_angle_deg"])
        row_ok = (
            0.0 <= values["expected_range_m"] < unambiguous_range_m
            and range_error <= max_range_resolution_multiples * resolution_m
            and angle_error <= max_angle_error_deg
        )
        rows_valid = rows_valid and row_ok
        ids.append(target_id)
        target_metrics.append({
            "target_id": target_id or f"target_{index}",
            "range_error_m": range_error,
            "range_error_resolution_multiples": range_error / resolution_m,
            "angle_error_deg": angle_error,
            "ok": row_ok,
        })

    checks = {
        "frequency_grid_strictly_increasing": increasing,
        "frequency_grid_equispaced": step_drift <= max_frequency_step_relative_drift,
        "positive_bandwidth": bandwidth_hz > 0.0,
        "target_ids_unique": len(ids) == len(set(ids)) and all(ids),
        "all_targets_localized": rows_valid,
    }
    return {
        "policy": "radar_range_angle_localization_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "frequency_count": len(frequency),
            "frequency_step_hz": representative_step_hz,
            "frequency_step_relative_drift": step_drift,
            "bandwidth_hz": bandwidth_hz,
            "physical_range_resolution_m": resolution_m,
            "unambiguous_range_m": unambiguous_range_m,
            "targets": target_metrics,
        },
        "notes": [
            "range tolerance is tied to c/(2*bandwidth), not display-grid spacing",
            "angle tolerance is declared separately because array aperture and estimator regularization set angular resolution",
        ],
    }
