"""Solver-neutral propagation and loss gate for a distributed RLC line pair."""

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


def _relative_gap(first: float, second: float) -> float:
    return abs(first - second) / max(abs(first), abs(second), 1.0e-300)


def distributed_line_delay_loss_gate(summary: Mapping[str, object]) -> dict[str, Any]:
    """Separate LC propagation delay from series-resistance attenuation."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    cases = summary.get("cases")
    if (
        not isinstance(cases, Sequence)
        or isinstance(cases, (str, bytes))
        or len(cases) != 2
        or any(not isinstance(case, Mapping) for case in cases)
    ):
        raise ValueError("cases must contain exactly two objects")
    replay = summary.get("replay")
    if not isinstance(replay, Mapping):
        raise ValueError("replay must be an object")

    parsed: dict[str, dict[str, float]] = {}
    for index, case in enumerate(cases):
        label = str(case.get("label") or "").strip()
        if label in parsed:
            raise ValueError(f"duplicate cases[{index}].label")
        parsed[label] = {
            "length": _finite(case.get("length"), f"cases[{index}].length", positive=True),
            "resistance": _finite(
                case.get("series_resistance_ohm_per_length"),
                f"cases[{index}].series_resistance_ohm_per_length",
                positive=True,
            ),
            "inductance": _finite(
                case.get("inductance_h_per_length"),
                f"cases[{index}].inductance_h_per_length",
                positive=True,
            ),
            "capacitance": _finite(
                case.get("capacitance_f_per_length"),
                f"cases[{index}].capacitance_f_per_length",
                positive=True,
            ),
            "delay": _finite(
                case.get("measured_one_way_delay_s"),
                f"cases[{index}].measured_one_way_delay_s",
                positive=True,
            ),
            "peak": _finite(
                case.get("first_pulse_peak_output_v"),
                f"cases[{index}].first_pulse_peak_output_v",
                positive=True,
            ),
        }
    required_labels = {"low_resistance", "high_resistance"}
    if set(parsed) != required_labels:
        raise ValueError(f"case labels must be {sorted(required_labels)}")

    low = parsed["low_resistance"]
    high = parsed["high_resistance"]
    expected_delay = low["length"] * math.sqrt(low["inductance"] * low["capacitance"])
    characteristic_impedance = math.sqrt(low["inductance"] / low["capacitance"])
    resistance_ratio = high["resistance"] / low["resistance"]
    peak_ratio = high["peak"] / low["peak"]
    delay_gap = _relative_gap(low["delay"], high["delay"])
    replay_errors = []
    for label in required_labels:
        row = replay.get(label)
        if not isinstance(row, Mapping):
            raise ValueError(f"replay.{label} must be an object")
        replay_errors.extend(
            [
                _finite(row.get("delay_relative_gap"), f"replay.{label}.delay_relative_gap"),
                _finite(row.get("peak_relative_gap"), f"replay.{label}.peak_relative_gap"),
            ]
        )

    same_lc = (
        low["length"] == high["length"]
        and low["inductance"] == high["inductance"]
        and low["capacitance"] == high["capacitance"]
    )
    checks = {
        "same_length_inductance_and_capacitance": same_lc,
        "spice_m_suffix_resistance_ratio_is_one_thousand": math.isclose(
            resistance_ratio, 1000.0, rel_tol=1.0e-12
        ),
        "both_arrivals_match_lc_delay": max(
            abs(row["delay"] - expected_delay) / expected_delay for row in (low, high)
        )
        <= 0.15,
        "resistance_change_preserves_delay_scale": delay_gap <= 0.10,
        "high_resistance_attenuates_first_pulse": 0.0 < peak_ratio <= 0.30,
        "observable_replay_is_deterministic": max(replay_errors) <= 1.0e-12,
    }
    return {
        "policy": "distributed_line_delay_loss_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "expected_lossless_one_way_delay_s": expected_delay,
            "lossless_characteristic_impedance_ohm": characteristic_impedance,
            "series_resistance_ratio_high_to_low": resistance_ratio,
            "measured_delay_cross_case_relative_gap": delay_gap,
            "first_pulse_peak_ratio_high_to_low": peak_ratio,
            "maximum_replay_relative_gap": max(replay_errors),
        },
        "lesson": (
            "In SPICE notation, suffix m means milli. For two lines with the same length, "
            "L, and C, the first-arrival scale follows length*sqrt(LC), while a large change "
            "in series R should primarily appear as attenuation rather than a new delay law."
        ),
    }
