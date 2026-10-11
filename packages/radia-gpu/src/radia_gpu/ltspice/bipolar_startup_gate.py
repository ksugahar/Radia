"""Solver-neutral bipolar supply startup and power-good gate."""

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


def _quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    left = int(math.floor(position))
    right = int(math.ceil(position))
    if left == right:
        return ordered[left]
    weight = position - left
    return ordered[left] * (1.0 - weight) + ordered[right] * weight


def _first_crossing(
    time_s: list[float], values: list[float], level: float, *, rising: bool
) -> float | None:
    for index in range(1, len(values)):
        left, right = values[index - 1], values[index]
        crossed = left < level <= right if rising else left > level >= right
        if not crossed:
            continue
        if right == left:
            return time_s[index]
        fraction = (level - left) / (right - left)
        return time_s[index - 1] + fraction * (time_s[index] - time_s[index - 1])
    return None


def bipolar_supply_startup_gate(
    summary: dict[str, Any],
    *,
    max_final_relative_error: float = 0.03,
    max_rail_balance_relative: float = 0.03,
    max_tail_ripple_fraction: float = 0.03,
    max_overshoot_fraction: float = 0.10,
    max_t90_skew_s: float = 1.0e-3,
    max_power_good_delay_s: float = 5.0e-4,
    min_sample_count: int = 100,
) -> dict[str, Any]:
    """Gate dual-rail startup from sampled voltages and power-good traces."""
    if not isinstance(summary, dict):
        raise ValueError("summary must be a mapping")
    names = (
        "time_s",
        "positive_v",
        "negative_v",
        "power_good_positive_v",
        "power_good_negative_v",
    )
    arrays = {name: _finite_vector(summary.get(name), name) for name in names}
    lengths = {len(values) for values in arrays.values()}
    if len(lengths) != 1:
        raise ValueError("all sampled traces must have equal length")
    time_s = arrays["time_s"]
    positive = arrays["positive_v"]
    negative = arrays["negative_v"]
    pg_positive = arrays["power_good_positive_v"]
    pg_negative = arrays["power_good_negative_v"]
    target = float(summary["target_magnitude_v"])
    if not math.isfinite(target) or target <= 0.0:
        raise ValueError("target_magnitude_v must be positive and finite")

    strictly_increasing = all(right > left for left, right in zip(time_s, time_s[1:]))
    if not strictly_increasing:
        tail_start = time_s[0]
    else:
        tail_start = time_s[0] + 0.9 * (time_s[-1] - time_s[0])
    tail_indices = [index for index, value in enumerate(time_s) if value >= tail_start]
    if not tail_indices:
        raise ValueError("tail window is empty")
    positive_tail = [positive[index] for index in tail_indices]
    negative_tail = [negative[index] for index in tail_indices]
    pg_positive_tail = [pg_positive[index] for index in tail_indices]
    pg_negative_tail = [pg_negative[index] for index in tail_indices]
    positive_final = _quantile(positive_tail, 0.5)
    negative_final = _quantile(negative_tail, 0.5)
    pg_positive_final = _quantile(pg_positive_tail, 0.5)
    pg_negative_final = _quantile(pg_negative_tail, 0.5)

    positive_error = abs(positive_final - target) / target
    negative_error = abs(abs(negative_final) - target) / target
    rail_balance = abs(abs(positive_final) - abs(negative_final)) / target
    positive_ripple = (_quantile(positive_tail, 0.99) - _quantile(positive_tail, 0.01)) / target
    negative_ripple = (_quantile(negative_tail, 0.99) - _quantile(negative_tail, 0.01)) / target
    positive_overshoot = max(max(positive) / target - 1.0, 0.0)
    negative_overshoot = max(abs(min(negative)) / target - 1.0, 0.0)

    positive_t10 = _first_crossing(time_s, positive, 0.1 * target, rising=True)
    positive_t90 = _first_crossing(time_s, positive, 0.9 * target, rising=True)
    negative_t10 = _first_crossing(time_s, negative, -0.1 * target, rising=False)
    negative_t90 = _first_crossing(time_s, negative, -0.9 * target, rising=False)
    pg_positive_t = (
        _first_crossing(time_s, pg_positive, 0.5 * pg_positive_final, rising=True)
        if pg_positive_final > 0.0
        else None
    )
    pg_negative_t = (
        _first_crossing(time_s, pg_negative, 0.5 * pg_negative_final, rising=True)
        if pg_negative_final > 0.0
        else None
    )
    crossings = (positive_t10, positive_t90, negative_t10, negative_t90, pg_positive_t, pg_negative_t)
    all_crossings_present = all(value is not None for value in crossings)
    t90_skew = abs(positive_t90 - negative_t90) if positive_t90 is not None and negative_t90 is not None else math.inf
    pg_positive_delay = pg_positive_t - positive_t90 if pg_positive_t is not None and positive_t90 is not None else math.inf
    pg_negative_delay = pg_negative_t - negative_t90 if pg_negative_t is not None and negative_t90 is not None else math.inf

    checks = {
        "units_explicit": summary.get("time_unit") == "s" and summary.get("voltage_unit") == "V",
        "sample_count_sufficient": len(time_s) >= int(min_sample_count),
        "time_strictly_increases": strictly_increasing,
        "initial_rails_near_zero": abs(positive[0]) <= 0.1 * target and abs(negative[0]) <= 0.1 * target,
        "final_rail_polarities_correct": positive_final > 0.0 and negative_final < 0.0,
        "final_positive_regulated": positive_error <= float(max_final_relative_error),
        "final_negative_regulated": negative_error <= float(max_final_relative_error),
        "rail_magnitudes_balanced": rail_balance <= float(max_rail_balance_relative),
        "positive_tail_ripple_bounded": positive_ripple <= float(max_tail_ripple_fraction),
        "negative_tail_ripple_bounded": negative_ripple <= float(max_tail_ripple_fraction),
        "positive_overshoot_bounded": positive_overshoot <= float(max_overshoot_fraction),
        "negative_overshoot_bounded": negative_overshoot <= float(max_overshoot_fraction),
        "startup_crossings_present": all_crossings_present,
        "ten_before_ninety_percent": all_crossings_present and positive_t10 < positive_t90 and negative_t10 < negative_t90,
        "t90_skew_bounded": t90_skew <= float(max_t90_skew_s),
        "power_good_asserts_after_rail_t90": pg_positive_delay >= 0.0 and pg_negative_delay >= 0.0,
        "power_good_delay_bounded": pg_positive_delay <= float(max_power_good_delay_s)
        and pg_negative_delay <= float(max_power_good_delay_s),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "bipolar_supply_startup_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "sample_count": len(time_s),
            "positive_final_v": positive_final,
            "negative_final_v": negative_final,
            "positive_final_relative_error": positive_error,
            "negative_final_relative_error": negative_error,
            "rail_balance_relative": rail_balance,
            "positive_tail_ripple_fraction": positive_ripple,
            "negative_tail_ripple_fraction": negative_ripple,
            "positive_overshoot_fraction": positive_overshoot,
            "negative_overshoot_fraction": negative_overshoot,
            "positive_t10_s": positive_t10,
            "positive_t90_s": positive_t90,
            "negative_t10_s": negative_t10,
            "negative_t90_s": negative_t90,
            "t90_skew_s": t90_skew,
            "power_good_positive_assert_s": pg_positive_t,
            "power_good_negative_assert_s": pg_negative_t,
            "power_good_positive_delay_s": pg_positive_delay,
            "power_good_negative_delay_s": pg_negative_delay,
        },
        "tolerances": {
            "max_final_relative_error": float(max_final_relative_error),
            "max_rail_balance_relative": float(max_rail_balance_relative),
            "max_tail_ripple_fraction": float(max_tail_ripple_fraction),
            "max_overshoot_fraction": float(max_overshoot_fraction),
            "max_t90_skew_s": float(max_t90_skew_s),
            "max_power_good_delay_s": float(max_power_good_delay_s),
            "min_sample_count": int(min_sample_count),
        },
        "notes": [
            "Evaluate signed rails separately; taking absolute values too early hides polarity faults.",
            "Power-good must follow each rail's 90 percent crossing and remain within a bounded delay.",
        ],
    }
