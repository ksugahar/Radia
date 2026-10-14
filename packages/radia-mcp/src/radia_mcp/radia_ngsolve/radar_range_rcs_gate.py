"""Solver-neutral validation for wideband range-resolved RCS profiles."""

from __future__ import annotations

import math
from statistics import median
from typing import Any


SPEED_OF_LIGHT_M_S = 299_792_458.0


def radar_range_rcs_profile_gate(
    frequency_hz: list[float],
    *,
    target_range_m: float,
    radar_peak_range_m: float,
    radar_peak_rcs_m2: float,
    generalized_peak_range_m: float,
    generalized_peak_rcs_m2: float,
    analytic_peak_rcs_m2: float,
    profile_relative_l2: float,
    max_frequency_step_relative_drift: float = 1.0e-7,
    max_peak_range_resolution_multiples: float = 1.0,
    max_profile_relative_l2: float = 1.0e-4,
    max_method_peak_relative_error: float = 1.0e-4,
    max_analytic_peak_relative_error: float = 0.05,
) -> dict[str, Any]:
    """Gate range localization and RCS amplitude without confusing display sampling with resolution."""

    frequency = [float(value) for value in frequency_hz]
    scalars = {
        "target_range_m": float(target_range_m),
        "radar_peak_range_m": float(radar_peak_range_m),
        "radar_peak_rcs_m2": float(radar_peak_rcs_m2),
        "generalized_peak_range_m": float(generalized_peak_range_m),
        "generalized_peak_rcs_m2": float(generalized_peak_rcs_m2),
        "analytic_peak_rcs_m2": float(analytic_peak_rcs_m2),
        "profile_relative_l2": float(profile_relative_l2),
    }
    tolerances = {
        "max_frequency_step_relative_drift": float(max_frequency_step_relative_drift),
        "max_peak_range_resolution_multiples": float(max_peak_range_resolution_multiples),
        "max_profile_relative_l2": float(max_profile_relative_l2),
        "max_method_peak_relative_error": float(max_method_peak_relative_error),
        "max_analytic_peak_relative_error": float(max_analytic_peak_relative_error),
    }
    if len(frequency) < 3:
        raise ValueError("frequency_hz must contain at least three samples")
    if not all(math.isfinite(value) and value > 0.0 for value in frequency):
        raise ValueError("frequency_hz must be finite and positive")
    if not all(math.isfinite(value) for value in scalars.values()):
        raise ValueError("RCS and range metrics must be finite")
    if any(not math.isfinite(value) or value < 0.0 for value in tolerances.values()):
        raise ValueError("tolerances must be finite and nonnegative")

    steps = [b - a for a, b in zip(frequency, frequency[1:])]
    strictly_increasing = all(step > 0.0 for step in steps)
    representative_step_hz = median(steps)
    step_relative_drift = (
        max(abs(step - representative_step_hz) for step in steps) / representative_step_hz
        if strictly_increasing
        else math.inf
    )
    bandwidth_hz = frequency[-1] - frequency[0]
    physical_range_resolution_m = (
        SPEED_OF_LIGHT_M_S / (2.0 * bandwidth_hz) if bandwidth_hz > 0.0 else math.inf
    )
    unambiguous_range_m = (
        SPEED_OF_LIGHT_M_S / (2.0 * representative_step_hz)
        if representative_step_hz > 0.0
        else 0.0
    )

    radar_range_error_m = abs(scalars["radar_peak_range_m"] - scalars["target_range_m"])
    generalized_range_error_m = abs(
        scalars["generalized_peak_range_m"] - scalars["target_range_m"]
    )
    method_peak_relative_error = abs(
        scalars["radar_peak_rcs_m2"] - scalars["generalized_peak_rcs_m2"]
    ) / max(
        abs(scalars["radar_peak_rcs_m2"]),
        abs(scalars["generalized_peak_rcs_m2"]),
        math.ulp(1.0),
    )
    analytic_errors = {
        "radar": abs(scalars["radar_peak_rcs_m2"] - scalars["analytic_peak_rcs_m2"])
        / max(abs(scalars["analytic_peak_rcs_m2"]), math.ulp(1.0)),
        "generalized": abs(
            scalars["generalized_peak_rcs_m2"] - scalars["analytic_peak_rcs_m2"]
        )
        / max(abs(scalars["analytic_peak_rcs_m2"]), math.ulp(1.0)),
    }
    max_range_error_m = (
        tolerances["max_peak_range_resolution_multiples"] * physical_range_resolution_m
    )

    checks = {
        "frequency_grid_strictly_increasing": strictly_increasing,
        "frequency_grid_equispaced": (
            step_relative_drift <= tolerances["max_frequency_step_relative_drift"]
        ),
        "positive_bandwidth": bandwidth_hz > 0.0,
        "target_inside_unambiguous_range": 0.0 <= scalars["target_range_m"] < unambiguous_range_m,
        "radar_peak_localizes_target": radar_range_error_m <= max_range_error_m,
        "generalized_peak_localizes_target": generalized_range_error_m <= max_range_error_m,
        "peak_rcs_positive": (
            scalars["radar_peak_rcs_m2"] > 0.0
            and scalars["generalized_peak_rcs_m2"] > 0.0
            and scalars["analytic_peak_rcs_m2"] > 0.0
        ),
        "two_reconstruction_profiles_agree": (
            0.0 <= scalars["profile_relative_l2"] <= tolerances["max_profile_relative_l2"]
        ),
        "two_reconstruction_peaks_agree": (
            method_peak_relative_error <= tolerances["max_method_peak_relative_error"]
        ),
        "radar_peak_matches_analytic_reference": (
            analytic_errors["radar"] <= tolerances["max_analytic_peak_relative_error"]
        ),
        "generalized_peak_matches_analytic_reference": (
            analytic_errors["generalized"] <= tolerances["max_analytic_peak_relative_error"]
        ),
    }
    return {
        "policy": "radar_range_rcs_profile_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "frequency_count": len(frequency),
            "frequency_step_hz": representative_step_hz,
            "frequency_step_relative_drift": step_relative_drift,
            "bandwidth_hz": bandwidth_hz,
            "physical_range_resolution_m": physical_range_resolution_m,
            "unambiguous_range_m": unambiguous_range_m,
            "radar_peak_range_error_m": radar_range_error_m,
            "generalized_peak_range_error_m": generalized_range_error_m,
            "profile_relative_l2": scalars["profile_relative_l2"],
            "method_peak_relative_error": method_peak_relative_error,
            "analytic_peak_relative_errors": analytic_errors,
        },
        "notes": [
            "Physical range resolution is c/(2*bandwidth); an interpolated display-grid step is not resolution.",
            "Unambiguous range is c/(2*frequency_step) for an equispaced frequency sweep.",
        ],
    }
