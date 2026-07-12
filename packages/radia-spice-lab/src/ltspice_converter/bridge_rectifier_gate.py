"""Solver-neutral single-phase full-wave bridge rectifier gate."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


def _number(value: object, name: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    return parsed


def bridge_rectifier_gate(summary: Mapping[str, object]) -> dict:
    """Gate frequency doubling, diagonal conduction, and current balance."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    tolerances = summary.get("tolerances")
    units = summary.get("units")
    window = summary.get("analysis_window_s")
    if not isinstance(tolerances, Mapping) or not isinstance(units, Mapping):
        raise ValueError("tolerances and units must be objects")
    if not isinstance(window, Sequence) or isinstance(window, (str, bytes)) or len(window) != 2:
        raise ValueError("analysis_window_s must contain [start, stop]")

    input_frequency = _number(summary.get("input_frequency_hz"), "input_frequency_hz")
    ripple_frequency = _number(summary.get("ripple_frequency_hz"), "ripple_frequency_hz")
    output_average = _number(summary.get("vout_average_v"), "vout_average_v")
    output_minimum = _number(summary.get("vout_min_v"), "vout_min_v")
    load_average = _number(summary.get("load_average_a"), "load_average_a")
    capacitor_average = _number(summary.get("capacitor_average_a"), "capacitor_average_a")
    diode_sum = _number(summary.get("diode_average_sum_a"), "diode_average_sum_a")
    pair_a_error = _number(
        summary.get("diagonal_pair_a_waveform_relative_error"),
        "diagonal_pair_a_waveform_relative_error",
    )
    pair_b_error = _number(
        summary.get("diagonal_pair_b_waveform_relative_error"),
        "diagonal_pair_b_waveform_relative_error",
    )
    alternate_overlap = _number(
        summary.get("alternate_pair_overlap_fraction"),
        "alternate_pair_overlap_fraction",
    )
    kcl_error = _number(summary.get("kcl_max_relative_error"), "kcl_max_relative_error")
    start = _number(window[0], "analysis_window_s[0]")
    stop = _number(window[1], "analysis_window_s[1]")
    if min(input_frequency, ripple_frequency, output_average, load_average, diode_sum) <= 0.0:
        raise ValueError("frequencies, output average, and average currents must be positive")
    if stop <= start:
        raise ValueError("analysis window stop must exceed start")

    frequency_tolerance = _number(
        tolerances.get("ripple_frequency_ratio_relative_error"),
        "tolerances.ripple_frequency_ratio_relative_error",
    )
    current_tolerance = _number(
        tolerances.get("diode_sum_to_twice_load_relative_error"),
        "tolerances.diode_sum_to_twice_load_relative_error",
    )
    pair_tolerance = _number(
        tolerances.get("diagonal_pair_waveform_relative_error"),
        "tolerances.diagonal_pair_waveform_relative_error",
    )
    overlap_tolerance = _number(
        tolerances.get("alternate_pair_overlap_fraction"),
        "tolerances.alternate_pair_overlap_fraction",
    )
    kcl_tolerance = _number(
        tolerances.get("kcl_max_relative_error"),
        "tolerances.kcl_max_relative_error",
    )
    capacitor_tolerance = _number(
        tolerances.get("capacitor_average_to_load_relative_error"),
        "tolerances.capacitor_average_to_load_relative_error",
    )
    if min(
        frequency_tolerance,
        current_tolerance,
        pair_tolerance,
        overlap_tolerance,
        kcl_tolerance,
        capacitor_tolerance,
    ) < 0.0:
        raise ValueError("tolerances must be nonnegative")

    ripple_ratio = ripple_frequency / input_frequency
    ripple_ratio_error = abs(ripple_ratio - 2.0) / 2.0
    current_sum_error = abs(diode_sum - 2.0 * load_average) / max(
        diode_sum, 2.0 * load_average
    )
    capacitor_average_error = abs(capacitor_average) / load_average
    checks = {
        "units_explicit": units.get("voltage") == "V"
        and units.get("current") == "A"
        and units.get("frequency") == "Hz"
        and units.get("time") == "s",
        "single_phase_full_wave_bridge_declared": summary.get("topology")
        == "single_phase_full_wave_bridge",
        "analysis_covers_two_input_cycles": (stop - start) * input_frequency >= 2.0 - 1.0e-9,
        "output_is_nonnegative": output_minimum >= 0.0,
        "ripple_frequency_is_twice_input": ripple_ratio_error <= frequency_tolerance,
        "four_diode_average_sum_is_twice_load": current_sum_error <= current_tolerance,
        "diagonal_pair_waveforms_match": max(pair_a_error, pair_b_error) <= pair_tolerance,
        "opposite_pairs_do_not_overlap": alternate_overlap <= overlap_tolerance,
        "output_node_kcl_closes": kcl_error <= kcl_tolerance,
        "capacitor_average_current_closes": capacitor_average_error <= capacitor_tolerance,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "schema": "radia-spice-bridge-rectifier/v1",
        "policy": "full_wave_bridge_frequency_pair_and_kcl_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "ripple_frequency_ratio": ripple_ratio,
            "ripple_frequency_ratio_relative_error": ripple_ratio_error,
            "diode_sum_to_twice_load_relative_error": current_sum_error,
            "diagonal_pair_max_waveform_relative_error": max(pair_a_error, pair_b_error),
            "alternate_pair_overlap_fraction": alternate_overlap,
            "kcl_max_relative_error": kcl_error,
            "capacitor_average_to_load_relative_error": capacitor_average_error,
        },
        "notes": [
            "A single-phase bridge has two output peaks per input cycle.",
            "Each load-current path contains two diodes, so the sum of four diode average currents is twice the load average current.",
            "Do not apply the reservoir-capacitor I/(2*f*C) ripple estimate unless the waveform remains in its small-ripple regime.",
        ],
    }
