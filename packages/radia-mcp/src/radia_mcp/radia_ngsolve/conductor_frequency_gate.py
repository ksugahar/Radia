"""Solver-neutral frequency-response gates for symmetric conductors."""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from typing import Any, Mapping


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


def homogenized_bundle_impedance_comparison_gate(
    rows: Iterable[Mapping[str, Any]],
    *,
    resistance_rtol: float = 0.03,
    inductance_rtol: float = 0.005,
    impedance_rtol: float = 0.01,
    observable_rtol: float = 1.0e-10,
    minimum_element_reduction: float = 5.0,
    minimum_speedup: float = 5.0,
) -> dict[str, object]:
    """Compare a homogenized stranded bundle with an explicit reference.

    The gate verifies passive complex impedance, recomputes ``Z=V/I`` and
    ``L=Im(Z)/omega``, then balances approximation error against mesh and solve
    cost. It is solver-neutral and suitable for round-wire homogenization,
    litz-wire surrogates, and explicit-strand reference models.
    """

    items = [dict(row) for row in rows]
    if len(items) != 2:
        raise ValueError("rows must contain one homogenized and one explicit_reference model")
    by_role = {str(row.get("model_role") or "").strip().lower(): row for row in items}
    if set(by_role) != {"homogenized", "explicit_reference"}:
        raise ValueError("model_role values must be homogenized and explicit_reference")

    normalized: dict[str, dict[str, Any]] = {}
    for role, row in by_role.items():
        try:
            frequency = float(row["frequency_hz"])
            current_pair = [float(value) for value in row["current_a_complex"]]
            voltage_pair = [float(value) for value in row["voltage_v_complex"]]
            resistance = float(row["resistance_ohm"])
            inductance = float(row["inductance_h"])
            elements = int(row["element_count"])
            solve_time = float(row["solve_time_s"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{role} row has missing or invalid observables") from exc
        if len(current_pair) != 2 or len(voltage_pair) != 2:
            raise ValueError("complex current and voltage must be [real, imag]")
        scalars = [frequency, *current_pair, *voltage_pair, resistance, inductance, solve_time]
        if not all(math.isfinite(value) for value in scalars):
            raise ValueError("all observables must be finite")
        current = complex(*current_pair)
        voltage = complex(*voltage_pair)
        if frequency <= 0.0 or abs(current) == 0.0 or elements <= 0 or solve_time <= 0.0:
            raise ValueError("frequency, current magnitude, element count, and solve time must be positive")
        impedance = voltage / current
        derived_l = impedance.imag / (2.0 * math.pi * frequency)
        normalized[role] = {
            "frequency_hz": frequency,
            "current": current,
            "impedance": impedance,
            "resistance_ohm": resistance,
            "inductance_h": inductance,
            "derived_resistance_ohm": impedance.real,
            "derived_inductance_h": derived_l,
            "element_count": elements,
            "solve_time_s": solve_time,
        }

    approximate = normalized["homogenized"]
    reference = normalized["explicit_reference"]
    resistance_error = abs(approximate["resistance_ohm"] - reference["resistance_ohm"]) / abs(reference["resistance_ohm"])
    inductance_error = abs(approximate["inductance_h"] - reference["inductance_h"]) / abs(reference["inductance_h"])
    impedance_error = abs(approximate["impedance"] - reference["impedance"]) / abs(reference["impedance"])
    element_reduction = reference["element_count"] / approximate["element_count"]
    speedup = reference["solve_time_s"] / approximate["solve_time_s"]

    observable_errors = {}
    for role, row in normalized.items():
        observable_errors[role] = {
            "resistance_relative": abs(row["derived_resistance_ohm"] - row["resistance_ohm"]) / max(abs(row["resistance_ohm"]), 1.0e-300),
            "inductance_relative": abs(row["derived_inductance_h"] - row["inductance_h"]) / max(abs(row["inductance_h"]), 1.0e-300),
        }

    checks = {
        "frequency_matches": approximate["frequency_hz"] == reference["frequency_hz"],
        "current_phasor_matches": approximate["current"] == reference["current"],
        "passive_positive_resistance": all(row["resistance_ohm"] > 0.0 for row in normalized.values()),
        "positive_series_inductance": all(row["inductance_h"] > 0.0 for row in normalized.values()),
        "reported_observables_match_voltage_current": all(
            error <= float(observable_rtol)
            for errors in observable_errors.values()
            for error in errors.values()
        ),
        "homogenized_resistance_accurate": resistance_error <= float(resistance_rtol),
        "homogenized_inductance_accurate": inductance_error <= float(inductance_rtol),
        "homogenized_complex_impedance_accurate": impedance_error <= float(impedance_rtol),
        "explicit_reference_has_more_elements": element_reduction >= float(minimum_element_reduction),
        "homogenized_model_is_faster": speedup >= float(minimum_speedup),
    }
    return {
        "policy": "homogenized_bundle_impedance_comparison_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "resistance_relative_error": resistance_error,
            "inductance_relative_error": inductance_error,
            "complex_impedance_relative_error": impedance_error,
            "element_count_reduction": element_reduction,
            "solve_time_speedup": speedup,
            "observable_reconstruction_errors": observable_errors,
        },
        "lesson": (
            "A strand homogenization is useful only when passive complex impedance agrees with an explicit "
            "reference and the saved element/time reduction is demonstrated. Compare R, L, and complex Z; "
            "a small |Z| error can otherwise hide a material resistance error."
        ),
    }
