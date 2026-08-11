"""Solver-neutral gate for a second-order analog all-pass response."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


def _finite(value: object, name: str, *, positive: bool = False) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed) or (positive and parsed <= 0.0):
        requirement = "positive and finite" if positive else "finite"
        raise ValueError(f"{name} must be {requirement}")
    return parsed


def _complex_pair(value: object, name: str) -> complex:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 2:
        raise ValueError(f"{name} must contain [real, imag]")
    parsed = complex(_finite(value[0], name), _finite(value[1], name))
    return parsed


def second_order_allpass_gate(summary: Mapping[str, object]) -> dict[str, Any]:
    """Gate flat magnitude, phase winding, group delay, and pole-zero mirroring."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    model = summary.get("model_contract")
    runs = summary.get("runs")
    pole_zero = summary.get("pole_zero")
    metrics = summary.get("metrics")
    timing = summary.get("timing_breakdown_s")
    if not isinstance(model, Mapping) or not isinstance(pole_zero, Mapping) or not isinstance(metrics, Mapping):
        raise ValueError("model_contract, pole_zero, and metrics must be objects")
    if not isinstance(runs, Sequence) or isinstance(runs, (str, bytes)) or len(runs) != 2:
        raise ValueError("runs must contain exactly two objects")
    if any(not isinstance(run, Mapping) for run in runs):
        raise ValueError("run entries must be objects")

    f0 = _finite(model.get("center_frequency_hz"), "center_frequency_hz", positive=True)
    quality = _finite(model.get("quality_factor"), "quality_factor", positive=True)
    gain = _finite(model.get("nominal_gain"), "nominal_gain", positive=True)
    expected_delay = 4.0 * quality / (2.0 * math.pi * f0)
    normalized_runs = []
    for index, run in enumerate(runs):
        normalized_runs.append({
            "point_count": int(_finite(run.get("point_count"), f"runs[{index}].point_count", positive=True)),
            "frequency_min_hz": _finite(run.get("frequency_min_hz"), f"runs[{index}].frequency_min_hz", positive=True),
            "frequency_max_hz": _finite(run.get("frequency_max_hz"), f"runs[{index}].frequency_max_hz", positive=True),
            "minimum_magnitude": _finite(run.get("minimum_magnitude"), f"runs[{index}].minimum_magnitude", positive=True),
            "maximum_magnitude": _finite(run.get("maximum_magnitude"), f"runs[{index}].maximum_magnitude", positive=True),
            "center_sample_hz": _finite(run.get("center_frequency_sample_hz"), f"runs[{index}].center_frequency_sample_hz", positive=True),
            "center_phase_error_deg": _finite(run.get("center_phase_error_deg"), f"runs[{index}].center_phase_error_deg"),
            "low_high_phase_sum_error_deg": _finite(run.get("low_high_phase_sum_error_deg"), f"runs[{index}].low_high_phase_sum_error_deg"),
            "group_delay_at_center_s": _finite(run.get("group_delay_at_center_s"), f"runs[{index}].group_delay_at_center_s", positive=True),
            "phase_monotonic_violation_rad": _finite(run.get("phase_monotonic_violation_rad"), f"runs[{index}].phase_monotonic_violation_rad"),
        })
    poles = pole_zero.get("poles")
    zeros = pole_zero.get("zeros")
    if not isinstance(poles, Sequence) or not isinstance(zeros, Sequence) or len(poles) != 2 or len(zeros) != 2:
        raise ValueError("pole_zero must contain two poles and two zeros")
    parsed_poles = [_complex_pair(value, "pole") for value in poles]
    parsed_zeros = [_complex_pair(value, "zero") for value in zeros]
    recorded_mirror_error = _finite(pole_zero.get("mirror_relative_error"), "mirror_relative_error")
    omega0 = 2.0 * math.pi * f0
    mirror_error = min(
        max(
            abs(parsed_zeros[0] + parsed_poles[order].conjugate()),
            abs(parsed_zeros[1] + parsed_poles[1 - order].conjugate()),
        )
        for order in (0, 1)
    ) / omega0
    analytic_l2 = _finite(metrics.get("maximum_analytic_complex_relative_l2"), "maximum_analytic_complex_relative_l2")
    replay_error = _finite(metrics.get("maximum_complex_replay_relative_error"), "maximum_complex_replay_relative_error")

    timing_ok = False
    if isinstance(timing, Mapping) and len(timing) == 4:
        try:
            timing_ok = all(_finite(value, "timing") >= 0.0 for value in timing.values())
        except ValueError:
            timing_ok = False
    checks = {
        "second_order_allpass_model_contract": model.get("topology") == "second_order_allpass",
        "two_dense_ac_replays_cover_two_decades": all(
            run["point_count"] >= 200
            and run["frequency_min_hz"] <= f0 / 10.0
            and run["frequency_max_hz"] >= f0 * 10.0
            for run in normalized_runs
        ),
        "magnitude_is_flat_at_nominal_gain": all(
            max(abs(run["minimum_magnitude"] / gain - 1.0), abs(run["maximum_magnitude"] / gain - 1.0)) <= 5.0e-6
            for run in normalized_runs
        ),
        "center_sample_and_minus_180_phase_are_resolved": all(
            abs(run["center_sample_hz"] / f0 - 1.0) <= 1.0e-9
            and abs(run["center_phase_error_deg"]) <= 0.1
            for run in normalized_runs
        ),
        "reciprocal_frequency_phase_sum_is_minus_360_degrees": all(
            abs(run["low_high_phase_sum_error_deg"]) <= 0.1 for run in normalized_runs
        ),
        "group_delay_matches_four_q_over_omega0": all(
            abs(run["group_delay_at_center_s"] - expected_delay) / expected_delay <= 0.01
            for run in normalized_runs
        ),
        "unwrapped_phase_is_monotone": all(
            run["phase_monotonic_violation_rad"] <= 1.0e-12 for run in normalized_runs
        ),
        "complex_transfer_matches_analytic_allpass": analytic_l2 <= 5.0e-6,
        "stable_poles_and_right_half_plane_zeros_are_mirrored": (
            all(value.real < 0.0 for value in parsed_poles)
            and all(value.real > 0.0 for value in parsed_zeros)
            and mirror_error <= 1.0e-12
            and abs(recorded_mirror_error - mirror_error) <= 1.0e-12
        ),
        "complex_observable_replay_is_deterministic": replay_error <= 1.0e-12,
        "exactly_four_timing_stages": timing_ok,
    }
    return {
        "policy": "second_order_allpass_phase_group_delay_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "expected_group_delay_at_f0_s": expected_delay,
            "maximum_analytic_complex_relative_l2": analytic_l2,
            "maximum_complex_replay_relative_error": replay_error,
            "pole_zero_mirror_relative_error": mirror_error,
        },
        "lesson": (
            "An all-pass section is not validated by flat magnitude alone. Check the full complex "
            "transfer, monotone phase winding, phase symmetry about f0, group delay 4Q/omega0, "
            "and the right-half-plane zero mirror of stable poles."
        ),
    }
