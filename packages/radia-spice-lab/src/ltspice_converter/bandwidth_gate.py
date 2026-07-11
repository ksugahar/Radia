"""Solver-neutral -3 dB bandwidth and measure-log consistency gate."""

from __future__ import annotations

import math
from typing import Any


def _finite_vector(values: Any, name: str) -> list[float]:
    if not isinstance(values, list):
        raise ValueError(f"{name} must be a list")
    parsed = [float(value) for value in values]
    if not parsed or not all(math.isfinite(value) for value in parsed):
        raise ValueError(f"{name} must contain finite values")
    return parsed


def _crossing(
    frequencies: list[float], response_db: list[float], threshold_db: float, *, rising: bool
) -> float:
    inside = [index for index, value in enumerate(response_db) if value >= threshold_db]
    if not inside:
        raise ValueError("response never reaches the -3 dB threshold")
    if rising:
        right = inside[0]
        left = right - 1
    else:
        left = inside[-1]
        right = left + 1
    if left < 0 or right >= len(frequencies):
        raise ValueError("threshold crossing is outside the sampled frequency band")
    denominator = response_db[right] - response_db[left]
    if denominator == 0.0:
        raise ValueError("flat samples cannot locate a threshold crossing")
    fraction = (threshold_db - response_db[left]) / denominator
    return frequencies[left] + fraction * (frequencies[right] - frequencies[left])


def measure_bandwidth_crossing_gate(
    summary: dict[str, Any],
    *,
    max_peak_db_absolute_error: float = 1.0e-6,
    max_crossing_relative_error: float = 1.0e-4,
    max_bandwidth_relative_error: float = 1.0e-4,
    min_sample_count: int = 20,
) -> dict[str, Any]:
    """Recompute a rise/fall -3 dB bandwidth from sampled AC magnitude."""
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    frequencies = _finite_vector(summary.get("frequency_hz"), "frequency_hz")
    magnitudes = _finite_vector(summary.get("magnitude"), "magnitude")
    if len(frequencies) != len(magnitudes):
        raise ValueError("frequency_hz and magnitude must have equal length")
    if any(value <= 0.0 for value in frequencies + magnitudes):
        raise ValueError("frequency and magnitude values must be positive")

    measured_peak_db = float(summary["measured_peak_db"])
    measured_lower = float(summary["measured_lower_3db_hz"])
    measured_upper = float(summary["measured_upper_3db_hz"])
    measured_bandwidth = float(summary["measured_bandwidth_hz"])
    if not all(
        math.isfinite(value)
        for value in (measured_peak_db, measured_lower, measured_upper, measured_bandwidth)
    ):
        raise ValueError("measured bandwidth values must be finite")

    response_db = [20.0 * math.log10(value) for value in magnitudes]
    raw_peak_db = max(response_db)
    peak_index = response_db.index(raw_peak_db)
    threshold_db = raw_peak_db - 20.0 * math.log10(math.sqrt(2.0))
    raw_lower = _crossing(frequencies, response_db, threshold_db, rising=True)
    raw_upper = _crossing(frequencies, response_db, threshold_db, rising=False)
    raw_bandwidth = raw_upper - raw_lower

    def relative_error(actual: float, expected: float) -> float:
        return abs(actual - expected) / max(abs(actual), abs(expected), 1.0e-30)

    lower_error = relative_error(raw_lower, measured_lower)
    upper_error = relative_error(raw_upper, measured_upper)
    bandwidth_error = relative_error(raw_bandwidth, measured_bandwidth)
    measured_closure_error = relative_error(measured_upper - measured_lower, measured_bandwidth)
    increasing = all(right > left for left, right in zip(frequencies, frequencies[1:]))
    checks = {
        "sample_count_sufficient": len(frequencies) >= int(min_sample_count),
        "frequency_strictly_increases": increasing,
        "peak_is_interior": 0 < peak_index < len(frequencies) - 1,
        "measured_crossings_ordered": 0.0 < measured_lower < measured_upper,
        "measured_bandwidth_closes": measured_closure_error <= 1.0e-12,
        "peak_measure_matches_samples": abs(raw_peak_db - measured_peak_db)
        <= float(max_peak_db_absolute_error),
        "lower_crossing_matches_samples": lower_error <= float(max_crossing_relative_error),
        "upper_crossing_matches_samples": upper_error <= float(max_crossing_relative_error),
        "bandwidth_matches_samples": bandwidth_error <= float(max_bandwidth_relative_error),
        "peak_lies_between_crossings": measured_lower < frequencies[peak_index] < measured_upper,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "measure_bandwidth_crossing_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "sample_count": len(frequencies),
            "raw_peak_db": raw_peak_db,
            "threshold_db": threshold_db,
            "peak_frequency_hz": frequencies[peak_index],
            "raw_lower_3db_hz": raw_lower,
            "raw_upper_3db_hz": raw_upper,
            "raw_bandwidth_hz": raw_bandwidth,
            "peak_db_absolute_error": abs(raw_peak_db - measured_peak_db),
            "lower_crossing_relative_error": lower_error,
            "upper_crossing_relative_error": upper_error,
            "bandwidth_relative_error": bandwidth_error,
        },
        "tolerances": {
            "max_peak_db_absolute_error": float(max_peak_db_absolute_error),
            "max_crossing_relative_error": float(max_crossing_relative_error),
            "max_bandwidth_relative_error": float(max_bandwidth_relative_error),
            "min_sample_count": int(min_sample_count),
        },
        "notes": [
            "Interpolate the response in dB on the sampled linear-frequency interval to replay rise/fall -3 dB measures.",
            "A bandwidth result needs both threshold crossings; a single low-pass cutoff is a different contract.",
        ],
    }
