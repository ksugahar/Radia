"""Solver-neutral complex-impedance gate for a coil self-resonance sweep."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any



def _numpy():
    import numpy as np

    return np


def _array(value: object, name: str):
    np = _numpy()
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    result = np.asarray(value, dtype=float).ravel()
    if result.size < 3 or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain at least three finite values")
    return result


def _finite(value: object, name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def coil_self_resonance_sweep_gate(summary: Mapping[str, object]) -> dict[str, Any]:
    """Gate passivity, reactance sign change, equivalent LC, and replay."""
    np = _numpy()
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    frequency = _array(summary.get("frequency_hz"), "frequency_hz")
    resistance = _array(summary.get("resistance_ohm"), "resistance_ohm")
    reactance = _array(summary.get("reactance_ohm"), "reactance_ohm")
    if not (frequency.size == resistance.size == reactance.size):
        raise ValueError("frequency, resistance, and reactance arrays must have equal length")
    if not np.all(np.diff(frequency) > 0.0):
        raise ValueError("frequency_hz must be strictly increasing")

    crossing_indices = np.flatnonzero(reactance[:-1] * reactance[1:] < 0.0)
    if len(crossing_indices) == 1:
        index = int(crossing_indices[0])
        lower_frequency = float(frequency[index])
        upper_frequency = float(frequency[index + 1])
        lower_reactance = float(reactance[index])
        upper_reactance = float(reactance[index + 1])
        resonance = lower_frequency - lower_reactance * (
            upper_frequency - lower_frequency
        ) / (upper_reactance - lower_reactance)
    else:
        index = -1
        lower_frequency = math.nan
        upper_frequency = math.nan
        resonance = math.nan

    low_frequency_inductance = float(reactance[0] / (2.0 * math.pi * frequency[0]))
    equivalent_capacitance = (
        1.0 / ((2.0 * math.pi * resonance) ** 2 * low_frequency_inductance)
        if resonance > 0.0 and low_frequency_inductance > 0.0
        else math.nan
    )
    peak_index = int(np.argmax(resistance))
    phase = np.degrees(np.arctan2(reactance, resistance))
    dataset_error = max(
        _finite(summary.get("dataset_frequency_relative_error"), "dataset_frequency_relative_error"),
        _finite(summary.get("dataset_resistance_relative_error"), "dataset_resistance_relative_error"),
        _finite(summary.get("dataset_reactance_relative_error"), "dataset_reactance_relative_error"),
    )
    replay_error = max(
        _finite(summary.get("replay_frequency_relative_error"), "replay_frequency_relative_error"),
        _finite(summary.get("replay_resistance_relative_error"), "replay_resistance_relative_error"),
        _finite(summary.get("replay_reactance_relative_error"), "replay_reactance_relative_error"),
    )

    checks = {
        "uniform_frequency_grid": bool(
            np.allclose(np.diff(frequency), np.diff(frequency)[0], rtol=1.0e-12, atol=0.0)
        ),
        "positive_series_resistance_is_passive": bool(np.all(resistance > 0.0)),
        "single_inductive_to_capacitive_transition": len(crossing_indices) == 1
        and reactance[0] > 0.0
        and reactance[-1] < 0.0,
        "interpolated_resonance_is_inside_bracket": len(crossing_indices) == 1
        and lower_frequency < resonance < upper_frequency,
        "resistance_peak_is_adjacent_to_resonance": math.isfinite(resonance)
        and abs(float(frequency[peak_index]) - resonance) <= float(np.diff(frequency)[0]),
        "positive_low_frequency_l_and_equivalent_c": low_frequency_inductance > 0.0
        and equivalent_capacitance > 0.0,
        "phase_changes_sign_with_reactance": phase[0] > 0.0 and phase[-1] < 0.0,
        "dataset_aliases_are_equivalent": dataset_error <= 1.0e-14,
        "fresh_sweep_replay_is_stable": replay_error <= 1.0e-6,
    }
    checks = {name: bool(value) for name, value in checks.items()}
    return {
        "policy": "coil_self_resonance_sweep_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "point_count": int(frequency.size),
            "frequency_start_hz": float(frequency[0]),
            "frequency_stop_hz": float(frequency[-1]),
            "interpolated_self_resonance_hz": resonance,
            "low_frequency_series_inductance_h": low_frequency_inductance,
            "equivalent_self_capacitance_f": equivalent_capacitance,
            "resistance_peak_frequency_hz": float(frequency[peak_index]),
            "maximum_dataset_relative_error": dataset_error,
            "maximum_replay_relative_error": replay_error,
        },
        "lesson": (
            "Validate a coil sweep as complex impedance, not as a constant inductance. "
            "Positive R is the passivity gate; the sign of X separates inductive and "
            "capacitive regimes. Interpolate the X=0 crossing for self-resonance and use "
            "a declared low-frequency L only to derive an equivalent self-capacitance."
        ),
    }
