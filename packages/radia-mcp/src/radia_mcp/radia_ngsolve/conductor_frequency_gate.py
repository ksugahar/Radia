"""Solver-neutral frequency-response gates for symmetric conductors."""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence


def twin_conductor_skin_effect_frequency_gate(
    frequencies_hz: Iterable[float],
    resistance_ohm: Iterable[Sequence[float]],
    inductance_h: Iterable[Sequence[float]],
    *,
    symmetry_rtol: float = 5.0e-4,
) -> dict[str, object]:
    """Gate passive, symmetric R/L frequency trends and derived impedance."""

    frequencies = [float(value) for value in frequencies_hz]
    resistance = [[float(value) for value in row] for row in resistance_ohm]
    inductance = [[float(value) for value in row] for row in inductance_h]
    if len(frequencies) < 5 or len(resistance) != len(frequencies) or len(inductance) != len(frequencies):
        raise ValueError("frequency, resistance, and inductance rows must have the same length >= 5")
    if any(len(row) != 2 for row in resistance + inductance):
        raise ValueError("each resistance and inductance row must contain two conductors")
    if not math.isfinite(float(symmetry_rtol)) or symmetry_rtol < 0.0:
        raise ValueError("symmetry_rtol must be finite and non-negative")
    scalars = frequencies + [value for row in resistance + inductance for value in row]
    if not all(math.isfinite(value) for value in scalars):
        raise ValueError("frequency-response values must be finite")

    symmetry_errors = [
        abs(a - b) / max(abs(a), abs(b), 1.0e-300)
        for row in resistance + inductance
        for a, b in [row]
    ]
    impedance = [
        [math.hypot(r, 2.0 * math.pi * frequency * l) for r, l in zip(r_row, l_row)]
        for frequency, r_row, l_row in zip(frequencies, resistance, inductance)
    ]
    checks = {
        "frequency_strictly_increasing_positive": all(value > 0.0 for value in frequencies)
        and all(a < b for a, b in zip(frequencies, frequencies[1:])),
        "resistance_positive": all(value > 0.0 for row in resistance for value in row),
        "inductance_positive": all(value > 0.0 for row in inductance for value in row),
        "resistance_non_decreasing": all(
            all(a <= b for a, b in zip(series, series[1:]))
            for series in zip(*resistance)
        ),
        "inductance_non_increasing": all(
            all(a >= b for a, b in zip(series, series[1:]))
            for series in zip(*inductance)
        ),
        "twin_conductor_symmetry": max(symmetry_errors) <= float(symmetry_rtol),
        "skin_effect_resistance_growth": all(resistance[-1][index] >= 2.0 * resistance[0][index] for index in range(2)),
        "impedance_magnitude_increasing": all(
            all(a < b for a, b in zip(series, series[1:]))
            for series in zip(*impedance)
        ),
    }
    return {
        "policy": "twin_conductor_skin_effect_frequency_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "frequency_count": len(frequencies),
            "maximum_twin_relative_mismatch": max(symmetry_errors),
            "resistance_growth_ratio": [resistance[-1][index] / resistance[0][index] for index in range(2)],
            "inductance_retention_ratio": [inductance[-1][index] / inductance[0][index] for index in range(2)],
            "impedance_magnitude_ohm": impedance,
        },
        "lesson": (
            "A passive conductor frequency sweep should keep R and L positive, show non-decreasing R "
            "and non-increasing L as skin/proximity effects develop, preserve geometric twin symmetry, "
            "and yield increasing |R+j omega L|."
        ),
    }
