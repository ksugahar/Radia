"""Periodic torque-table sampling and duplicate-endpoint gate."""
from __future__ import annotations

import math


def periodic_torque_sampling_gate(
    *,
    period_deg: float,
    sample_count: int,
    endpoint_included: bool,
    spectrum_excludes_duplicate_endpoint: bool,
    torque_min_Nm: float,
    torque_max_Nm: float,
    speed_rps: float,
    expected_step_deg: float | None = None,
    step_tolerance_deg: float = 1.0e-9,
) -> dict:
    period = float(period_deg)
    count = int(sample_count)
    low = float(torque_min_Nm)
    high = float(torque_max_Nm)
    speed = float(speed_rps)
    tolerance = float(step_tolerance_deg)
    if period <= 0.0:
        raise ValueError("period_deg must be > 0")
    if count < 3:
        raise ValueError("sample_count must be >= 3")
    if not isinstance(endpoint_included, bool) or not isinstance(spectrum_excludes_duplicate_endpoint, bool):
        raise ValueError("endpoint policy values must be booleans")
    if tolerance < 0.0:
        raise ValueError("step_tolerance_deg must be >= 0")

    interval_count = count - 1 if endpoint_included else count
    step = period / interval_count
    unique_count = count - 1 if endpoint_included else count
    expected_step = None if expected_step_deg is None else float(expected_step_deg)
    finite = all(math.isfinite(value) for value in (period, low, high, speed, step))
    checks = {
        "all_finite": finite,
        "torque_range_ordered": finite and high >= low,
        "torque_profile_nontrivial": finite and high > low,
        "speed_positive": finite and speed > 0.0,
        "expected_step_matches": expected_step is None or abs(step - expected_step) <= tolerance,
        "duplicate_endpoint_excluded_from_spectrum": (
            not endpoint_included or spectrum_excludes_duplicate_endpoint
        ),
    }
    return {
        "policy": "periodic_torque_sampling_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "period_deg": period,
        "sample_count": count,
        "interval_count": interval_count,
        "unique_spectrum_sample_count": unique_count,
        "endpoint_included": endpoint_included,
        "spectrum_excludes_duplicate_endpoint": spectrum_excludes_duplicate_endpoint,
        "step_deg": step,
        "expected_step_deg": expected_step,
        "torque_min_Nm": low,
        "torque_max_Nm": high,
        "torque_ripple_peak_to_peak_Nm": high - low,
        "speed_rps": speed,
        "checks": checks,
        "lesson": (
            "A periodic torque table may include both endpoints for plotting, but "
            "the duplicate endpoint must be removed before FFT or ripple harmonics."
        ),
    }
