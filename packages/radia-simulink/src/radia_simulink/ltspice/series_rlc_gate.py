"""Solver-neutral complex-impedance gate for a current-driven series RLC."""
from __future__ import annotations

import math
from collections.abc import Mapping


def _positive_number(value: object, name: str) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise TypeError(f"{name} must be numeric")
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return number


def _nonnegative_number(value: object, name: str) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise TypeError(f"{name} must be numeric")
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return number


def series_rlc_complex_impedance_gate(
    summary: Mapping[str, object],
    *,
    max_complex_relative_l2: float = 5.0e-7,
    max_pointwise_relative_error: float = 2.0e-5,
    max_minimum_impedance_relative_error: float = 2.0e-3,
) -> dict[str, object]:
    """Gate full-complex AC response, resonance, and conversion semantics.

    The current source is oriented from the driven node toward ground.  Its
    one-ampere AC convention therefore gives
    ``V(in)=-(R + j*w*L + 1/(j*w*C))``.  Magnitude-only evidence cannot pass
    this gate because it loses both source orientation and phase.
    """

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")

    resistance = _positive_number(summary.get("resistance_ohm"), "resistance_ohm")
    inductance = _positive_number(summary.get("inductance_h"), "inductance_h")
    capacitance = _positive_number(summary.get("capacitance_f"), "capacitance_f")
    points = int(_positive_number(summary.get("points"), "points"))
    reported_f0 = _positive_number(
        summary.get("analytic_resonance_frequency_hz"),
        "analytic_resonance_frequency_hz",
    )
    minimum_impedance = _positive_number(
        summary.get("minimum_impedance_ohm"), "minimum_impedance_ohm"
    )
    pointwise_error = _nonnegative_number(
        summary.get("maximum_input_pointwise_relative_error"),
        "maximum_input_pointwise_relative_error",
    )

    errors = summary.get("complex_relative_l2")
    if not isinstance(errors, Mapping) or len(errors) < 3:
        raise ValueError("complex_relative_l2 must contain at least three traces")
    complex_errors = {
        str(name): _nonnegative_number(value, f"complex_relative_l2.{name}")
        for name, value in errors.items()
    }
    bracket = summary.get("resonance_bracket_hz")
    if not isinstance(bracket, (list, tuple)) or len(bracket) != 2:
        raise ValueError("resonance_bracket_hz must contain two frequencies")
    low = _positive_number(bracket[0], "resonance_bracket_hz[0]")
    high = _positive_number(bracket[1], "resonance_bracket_hz[1]")
    if low > high:
        raise ValueError("resonance_bracket_hz must be ordered")

    computed_f0 = 1.0 / (2.0 * math.pi * math.sqrt(inductance * capacitance))
    f0_relative_error = abs(reported_f0 - computed_f0) / computed_f0
    minimum_impedance_relative_error = abs(minimum_impedance - resistance) / resistance
    checks = {
        "ac_analysis_recorded": summary.get("analysis") == "ac",
        "source_orientation_preserved": summary.get("source_orientation")
        == "current_leaves_driven_node",
        "frequency_sweep_is_resolved": points >= 100,
        "three_full_complex_traces_match": max(complex_errors.values())
        <= max_complex_relative_l2,
        "pointwise_input_impedance_matches": pointwise_error
        <= max_pointwise_relative_error,
        "reported_resonance_matches_components": f0_relative_error <= 1.0e-12,
        "resonance_is_bracketed": low <= computed_f0 <= high,
        "minimum_impedance_matches_resistance": minimum_impedance_relative_error
        <= max_minimum_impedance_relative_error,
        "converted_and_reference_raw_are_equivalent": summary.get(
            "converted_and_reference_raw_are_equivalent"
        )
        is True,
        "converted_netlist_semantics_verified": summary.get(
            "converted_netlist_semantics_verified"
        )
        is True,
        "known_broken_schematic_is_rejected": summary.get(
            "known_broken_schematic_is_rejected"
        )
        is True,
    }
    return {
        "schema": "radia-spice-lab.series-rlc-complex-impedance.v1",
        "policy": "series_rlc_complex_impedance_and_conversion_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "computed_resonance_frequency_hz": computed_f0,
            "resonance_frequency_relative_error": f0_relative_error,
            "minimum_impedance_relative_error": minimum_impedance_relative_error,
            "maximum_complex_relative_l2": max(complex_errors.values()),
            "maximum_input_pointwise_relative_error": pointwise_error,
            "points": points,
        },
        "lesson": (
            "Accept schematic conversion only after the generated connectivity "
            "and full-complex AC response agree; component inventory alone is insufficient."
        ),
    }
