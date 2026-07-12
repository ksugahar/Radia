"""Solver-neutral reciprocal two-port S-parameter and power-balance gate."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from typing import Any


def _complex(value: Any, name: str) -> complex:
    try:
        if isinstance(value, Mapping):
            result = complex(float(value["real"]), float(value["imag"]))
        else:
            result = complex(float(value[0]), float(value[1]))
    except (KeyError, IndexError, TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be [real, imag] or a real/imag mapping") from exc
    if not math.isfinite(result.real) or not math.isfinite(result.imag):
        raise ValueError(f"{name} must be finite")
    return result


def reciprocal_two_port_power_sweep_gate(
    rows: Iterable[Mapping[str, Any]],
    adaptive_power_rows: Iterable[Mapping[str, Any]],
    *,
    reference_impedance_ohm: float,
    reciprocity_atol: float = 5.0e-6,
    reflection_symmetry_atol: float = 2.0e-4,
    passivity_atol: float = 1.0e-10,
    power_balance_atol: float = 5.0e-8,
) -> dict[str, object]:
    """Gate a reciprocal symmetric two-port sweep and independent power rows."""

    sweep = [dict(row) for row in rows]
    power = [dict(row) for row in adaptive_power_rows]
    if len(sweep) < 5:
        raise ValueError("rows must contain at least five frequency samples")
    if len(power) < 3:
        raise ValueError("adaptive_power_rows must contain at least three samples")
    z0 = float(reference_impedance_ohm)
    if not math.isfinite(z0) or z0 <= 0.0:
        raise ValueError("reference_impedance_ohm must be finite and positive")

    frequencies: list[float] = []
    reciprocity: list[float] = []
    symmetry: list[float] = []
    passivity_1: list[float] = []
    passivity_2: list[float] = []
    transmission: list[float] = []
    for row in sweep:
        try:
            frequency = float(row["frequency_hz"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("each sweep row must record finite frequency_hz") from exc
        if not math.isfinite(frequency):
            raise ValueError("frequency_hz must be finite")
        s11 = _complex(row.get("s11"), "s11")
        s12 = _complex(row.get("s12"), "s12")
        s21 = _complex(row.get("s21"), "s21")
        s22 = _complex(row.get("s22"), "s22")
        frequencies.append(frequency)
        reciprocity.append(abs(s12 - s21))
        symmetry.append(abs(s11 - s22))
        passivity_1.append(1.0 - abs(s11) ** 2 - abs(s21) ** 2)
        passivity_2.append(1.0 - abs(s22) ** 2 - abs(s12) ** 2)
        transmission.append(abs(s21))

    closure: list[float] = []
    stimulated: dict[int, list[float]] = {1: [], 2: []}
    power_frequencies: dict[int, list[float]] = {1: [], 2: []}
    for row in power:
        try:
            port = int(row["excitation_port"])
            frequency = float(row["frequency_hz"])
            balance = float(row["balance"])
            accepted = float(row["accepted_power_w"])
            source = float(row["stimulated_power_w"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("invalid adaptive power row") from exc
        if port not in (1, 2) or not all(math.isfinite(value) for value in (frequency, balance, accepted, source)):
            raise ValueError("adaptive power rows require ports 1/2 and finite values")
        power_frequencies[port].append(frequency)
        stimulated[port].append(source)
        closure.append(abs((1.0 - balance) - accepted))

    peak = max(range(len(transmission)), key=transmission.__getitem__)
    checks = {
        "frequency_axis_strictly_increasing": all(b > a for a, b in zip(frequencies, frequencies[1:])),
        "reference_impedance_positive": z0 > 0.0,
        "reciprocity_closes": max(reciprocity) <= float(reciprocity_atol),
        "symmetric_port_reflections_close": max(symmetry) <= float(reflection_symmetry_atol),
        "passivity_margins_nonnegative": min(passivity_1 + passivity_2) >= -float(passivity_atol),
        "both_excitation_power_axes_recorded": all(len(power_frequencies[port]) >= 3 for port in (1, 2)),
        "adaptive_power_axes_increase": all(
            all(b > a for a, b in zip(power_frequencies[port], power_frequencies[port][1:]))
            for port in (1, 2)
        ),
        "accepted_power_closes_balance": max(closure) <= float(power_balance_atol),
        "stimulated_power_constant_per_port": all(
            max(stimulated[port]) - min(stimulated[port]) <= 1.0e-12 for port in (1, 2)
        ),
    }
    return {
        "policy": "reciprocal_two_port_power_sweep_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "frequency_count": len(frequencies),
            "frequency_range_hz": [frequencies[0], frequencies[-1]],
            "maximum_reciprocity_absolute_error": max(reciprocity),
            "maximum_reflection_symmetry_absolute_error": max(symmetry),
            "passivity_margin_range": [min(passivity_1 + passivity_2), max(passivity_1 + passivity_2)],
            "maximum_power_balance_absolute_error_w": max(closure),
            "peak_transmission_frequency_hz": frequencies[peak],
            "peak_s21_magnitude": transmission[peak],
        },
        "lesson": (
            "A reciprocal two-port sweep should be checked as a complex matrix, not only as an S21 plot. "
            "Bind reciprocity, port symmetry, passivity, and independent accepted-power closure to aligned frequency axes."
        ),
    }
