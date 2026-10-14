"""Agreement gate for two independently evaluated static-torque curves."""

from __future__ import annotations

import json
import math
from typing import Any


def _finite_vector(values: Any, name: str) -> list[float]:
    if not isinstance(values, list):
        raise ValueError(f"{name} must be a list")
    parsed = [float(value) for value in values]
    if not parsed or not all(math.isfinite(value) for value in parsed):
        raise ValueError(f"{name} must contain finite values")
    return parsed


def dual_torque_method_curve_gate(
    summary: dict[str, Any],
    *,
    max_rms_relative_difference: float = 1.0e-2,
    max_point_relative_difference: float = 2.0e-2,
    max_endpoint_relative: float = 1.0e-2,
    max_peak_angle_difference_deg: float = 5.0,
    min_sample_count: int = 9,
) -> dict[str, Any]:
    """Gate air-gap-volume and stress-contour torque on one angle grid."""
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    angles = _finite_vector(summary.get("angles_deg"), "angles_deg")
    primary = _finite_vector(summary.get("primary_torque_nm"), "primary_torque_nm")
    secondary = _finite_vector(summary.get("secondary_torque_nm"), "secondary_torque_nm")
    same_length = len(angles) == len(primary) == len(secondary)
    if not same_length:
        raise ValueError("angle and torque arrays must have equal length")

    scale = max(max(map(abs, primary)), max(map(abs, secondary)), 1.0e-30)
    differences = [left - right for left, right in zip(primary, secondary)]
    rms_relative = math.sqrt(sum(value * value for value in differences) / len(differences)) / scale
    maximum_relative = max(map(abs, differences)) / scale
    primary_peak_index = max(range(len(primary)), key=primary.__getitem__)
    secondary_peak_index = max(range(len(secondary)), key=secondary.__getitem__)
    peak_angle_difference = abs(angles[primary_peak_index] - angles[secondary_peak_index])
    endpoint_relative = max(
        abs(primary[0]), abs(primary[-1]), abs(secondary[0]), abs(secondary[-1])
    ) / scale
    increasing = all(right > left for left, right in zip(angles, angles[1:]))

    checks = {
        "sample_count_sufficient": len(angles) >= int(min_sample_count),
        "common_finite_grid": same_length,
        "angle_strictly_increases": increasing,
        "angle_unit_degrees": summary.get("angle_unit") == "deg",
        "torque_unit_newton_metre": summary.get("torque_unit") == "N*m",
        "nontrivial_torque_scale": scale > 1.0e-12,
        "rms_method_agreement": rms_relative <= float(max_rms_relative_difference),
        "pointwise_method_agreement": maximum_relative <= float(max_point_relative_difference),
        "low_torque_endpoints": endpoint_relative <= float(max_endpoint_relative),
        "primary_peak_is_interior": 0 < primary_peak_index < len(primary) - 1,
        "secondary_peak_is_interior": 0 < secondary_peak_index < len(secondary) - 1,
        "peak_angles_agree": peak_angle_difference <= float(max_peak_angle_difference_deg),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "dual_torque_method_curve_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "sample_count": len(angles),
            "torque_scale_nm": scale,
            "rms_method_difference_relative": rms_relative,
            "maximum_method_difference_relative": maximum_relative,
            "endpoint_relative": endpoint_relative,
            "primary_peak_angle_deg": angles[primary_peak_index],
            "secondary_peak_angle_deg": angles[secondary_peak_index],
            "peak_angle_difference_deg": peak_angle_difference,
        },
        "tolerances": {
            "max_rms_relative_difference": float(max_rms_relative_difference),
            "max_point_relative_difference": float(max_point_relative_difference),
            "max_endpoint_relative": float(max_endpoint_relative),
            "max_peak_angle_difference_deg": float(max_peak_angle_difference_deg),
            "min_sample_count": int(min_sample_count),
        },
        "notes": [
            "Agreement between two postprocessors is a consistency check, not an analytic torque reference.",
            "Keep both method identifiers, the common angle grid, and units with the result artifact.",
        ],
    }


def dual_torque_method_curve_gate_json(summary_json: str) -> dict[str, Any]:
    return dual_torque_method_curve_gate(json.loads(summary_json))
