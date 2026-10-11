"""Solver-neutral gate for a balanced three-phase delta-connected RL load."""

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


def _triple(
    value: object,
    name: str,
    *,
    positive: bool = False,
) -> list[float]:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != 3
    ):
        raise ValueError(f"{name} must contain exactly three values")
    return [
        _finite(item, f"{name}[{index}]", positive=positive)
        for index, item in enumerate(value)
    ]


def _complex_pair(value: object, name: str) -> complex:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != 2
    ):
        raise ValueError(f"{name} must contain [real, imag]")
    return complex(_finite(value[0], name), _finite(value[1], name))


def _relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(abs(expected), 1.0e-300)


def balanced_three_phase_delta_rl_gate(summary: Mapping[str, object]) -> dict[str, Any]:
    """Gate sequence, delta identities, complex power, and deterministic replay."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    model = summary.get("model_contract")
    metrics = summary.get("metrics")
    timing = summary.get("timing_breakdown_s")
    if not isinstance(model, Mapping) or not isinstance(metrics, Mapping):
        raise ValueError("model_contract and metrics must be objects")
    positive = metrics.get("positive")
    if not isinstance(positive, Mapping):
        raise ValueError("metrics.positive must be an object")

    frequency = _finite(model.get("frequency_hz"), "frequency_hz", positive=True)
    phase_voltage = _finite(
        model.get("phase_voltage_rms_v"), "phase_voltage_rms_v", positive=True
    )
    resistance = _finite(
        model.get("branch_resistance_ohm"), "branch_resistance_ohm", positive=True
    )
    inductance = _finite(
        model.get("branch_inductance_h"), "branch_inductance_h", positive=True
    )
    resistances = _triple(
        model.get("delta_resistances_ohm"),
        "delta_resistances_ohm",
        positive=True,
    )
    inductances = _triple(
        model.get("delta_inductances_h"),
        "delta_inductances_h",
        positive=True,
    )
    source_rows = model.get("source_rows")
    if (
        not isinstance(source_rows, Sequence)
        or isinstance(source_rows, (str, bytes))
        or len(source_rows) != 3
        or any(not isinstance(row, Mapping) for row in source_rows)
    ):
        raise ValueError("source_rows must contain exactly three objects")

    source_peaks = [
        _finite(row.get("peak_voltage_v"), "source peak", positive=True)
        for row in source_rows
    ]
    source_frequencies = [
        _finite(row.get("frequency_hz"), "source frequency", positive=True)
        for row in source_rows
    ]
    source_phases = [
        _finite(row.get("phase_deg"), "source phase") for row in source_rows
    ]

    reactance = 2.0 * math.pi * frequency * inductance
    impedance = complex(resistance, reactance)
    expected_line_voltage = math.sqrt(3.0) * phase_voltage
    expected_branch_current = expected_line_voltage / abs(impedance)
    expected_line_current = math.sqrt(3.0) * expected_branch_current
    expected_power = 3.0 * expected_line_voltage**2 * resistance / abs(impedance) ** 2
    expected_reactive = (
        3.0 * expected_line_voltage**2 * reactance / abs(impedance) ** 2
    )
    expected_pf = resistance / abs(impedance)

    model_impedance = _complex_pair(
        model.get("expected_branch_impedance_ohm"),
        "expected_branch_impedance_ohm",
    )
    expected_fields = {
        "line_voltage": _finite(
            model.get("expected_line_voltage_rms_v"),
            "expected_line_voltage_rms_v",
            positive=True,
        ),
        "branch_current": _finite(
            model.get("expected_branch_current_rms_a"),
            "expected_branch_current_rms_a",
            positive=True,
        ),
        "line_current": _finite(
            model.get("expected_line_current_rms_a"),
            "expected_line_current_rms_a",
            positive=True,
        ),
        "power": _finite(
            model.get("expected_active_power_w"),
            "expected_active_power_w",
            positive=True,
        ),
        "reactive": _finite(
            model.get("expected_reactive_power_var"),
            "expected_reactive_power_var",
            positive=True,
        ),
        "power_factor": _finite(
            model.get("expected_power_factor"),
            "expected_power_factor",
            positive=True,
        ),
    }

    phase_rms = _triple(
        positive.get("phase_voltage_rms_v"), "phase_voltage_rms_v", positive=True
    )
    line_rms = _triple(
        positive.get("line_voltage_rms_v"), "line_voltage_rms_v", positive=True
    )
    branch_rms = _triple(
        positive.get("branch_current_rms_a"), "branch_current_rms_a", positive=True
    )
    line_current_rms = _triple(
        positive.get("line_current_rms_a"), "line_current_rms_a", positive=True
    )
    branch_impedances_raw = positive.get("branch_impedance_ohm")
    if (
        not isinstance(branch_impedances_raw, Sequence)
        or isinstance(branch_impedances_raw, (str, bytes))
        or len(branch_impedances_raw) != 3
    ):
        raise ValueError("branch_impedance_ohm must contain three complex pairs")
    branch_impedances = [
        _complex_pair(value, "branch_impedance_ohm")
        for value in branch_impedances_raw
    ]
    source_power = _complex_pair(
        positive.get("source_complex_power_va"), "source_complex_power_va"
    )
    load_power = _complex_pair(
        positive.get("load_complex_power_va"), "load_complex_power_va"
    )

    point_count = int(
        _finite(positive.get("point_count"), "point_count", positive=True)
    )
    fit_start = _finite(positive.get("fit_window_start_s"), "fit_window_start_s")
    fit_stop = _finite(
        positive.get("fit_window_stop_s"), "fit_window_stop_s", positive=True
    )
    replay_error = _finite(
        metrics.get("maximum_phasor_replay_relative_error"),
        "maximum_phasor_replay_relative_error",
    )

    metric = {
        name: _finite(positive.get(name), name)
        for name in (
            "maximum_phase_voltage_relative_error",
            "maximum_line_voltage_relative_error",
            "maximum_branch_current_relative_error",
            "maximum_line_current_relative_error",
            "maximum_branch_impedance_relative_error",
            "line_current_kcl_relative_error",
            "phase_voltage_positive_to_negative_sequence_ratio",
            "phase_voltage_zero_sequence_ratio",
            "branch_current_magnitude_spread_relative",
            "line_current_magnitude_spread_relative",
            "source_load_complex_power_relative_error",
            "active_power_relative_error",
            "reactive_power_relative_error",
            "power_factor_absolute_error",
            "instantaneous_power_mean_relative_error",
            "instantaneous_power_ripple_relative",
            "maximum_phasor_fit_relative_residual",
        )
    }
    measured_pf = _finite(positive.get("power_factor"), "power_factor")

    expected_model_errors = [
        abs(model_impedance - impedance) / abs(impedance),
        _relative_error(expected_fields["line_voltage"], expected_line_voltage),
        _relative_error(expected_fields["branch_current"], expected_branch_current),
        _relative_error(expected_fields["line_current"], expected_line_current),
        _relative_error(expected_fields["power"], expected_power),
        _relative_error(expected_fields["reactive"], expected_reactive),
        _relative_error(expected_fields["power_factor"], expected_pf),
    ]

    timing_ok = False
    if isinstance(timing, Mapping) and len(timing) == 4:
        try:
            timing_ok = all(_finite(value, "timing") >= 0.0 for value in timing.values())
        except ValueError:
            timing_ok = False

    checks = {
        "balanced_y_source_delta_rl_model_contract": model.get("topology")
        == "balanced_y_source_delta_rl_load"
        and model.get("phase_sequence") == "abc",
        "source_amplitudes_frequency_and_abc_phases_close": max(
            _relative_error(value, phase_voltage * math.sqrt(2.0))
            for value in source_peaks
        )
        <= 1.0e-12
        and max(_relative_error(value, frequency) for value in source_frequencies)
        <= 1.0e-12
        and all(
            abs(((actual - expected + 180.0) % 360.0) - 180.0) <= 1.0e-12
            for actual, expected in zip(source_phases, (0.0, -120.0, 120.0))
        ),
        "three_equal_delta_branches_and_derived_expectations_close": max(
            _relative_error(value, resistance) for value in resistances
        )
        <= 1.0e-12
        and max(_relative_error(value, inductance) for value in inductances)
        <= 1.0e-12
        and max(expected_model_errors) <= 1.0e-12,
        "phase_sequence_and_sqrt3_line_voltage_identity_close": max(
            metric["maximum_phase_voltage_relative_error"],
            metric["maximum_line_voltage_relative_error"],
            max(_relative_error(value, phase_voltage) for value in phase_rms),
            max(_relative_error(value, expected_line_voltage) for value in line_rms),
        )
        <= 2.0e-5
        and max(
            metric["phase_voltage_positive_to_negative_sequence_ratio"],
            metric["phase_voltage_zero_sequence_ratio"],
        )
        <= 1.0e-8,
        "delta_impedance_branch_and_line_current_identities_close": max(
            metric["maximum_branch_impedance_relative_error"],
            metric["maximum_branch_current_relative_error"],
            metric["maximum_line_current_relative_error"],
            metric["line_current_kcl_relative_error"],
            max(abs(value - impedance) / abs(impedance) for value in branch_impedances),
            max(_relative_error(value, expected_branch_current) for value in branch_rms),
            max(_relative_error(value, expected_line_current) for value in line_current_rms),
        )
        <= 2.0e-4,
        "balanced_branch_and_line_current_magnitudes_close": max(
            metric["branch_current_magnitude_spread_relative"],
            metric["line_current_magnitude_spread_relative"],
        )
        <= 1.0e-8,
        "complex_power_power_factor_and_conservation_close": max(
            metric["source_load_complex_power_relative_error"],
            metric["active_power_relative_error"],
            metric["reactive_power_relative_error"],
            metric["power_factor_absolute_error"],
            abs(source_power - load_power) / abs(load_power),
            _relative_error(load_power.real, expected_power),
            _relative_error(load_power.imag, expected_reactive),
            abs(measured_pf - expected_pf),
        )
        <= 2.0e-4,
        "balanced_instantaneous_three_phase_power_is_constant": max(
            metric["instantaneous_power_mean_relative_error"],
            metric["instantaneous_power_ripple_relative"],
        )
        <= 2.0e-4,
        "steady_state_phasor_fit_has_two_periods_and_small_residual": point_count
        >= 500
        and fit_stop - fit_start >= 1.9 / frequency
        and metric["maximum_phasor_fit_relative_residual"] <= 2.0e-4,
        "positive_phasor_replay_is_deterministic": replay_error <= 1.0e-12,
        "exactly_four_timing_stages": timing_ok,
    }
    return {
        "policy": "balanced_three_phase_delta_rl_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "derived_line_voltage_rms_v": expected_line_voltage,
            "derived_branch_current_rms_a": expected_branch_current,
            "derived_line_current_rms_a": expected_line_current,
            "derived_active_power_w": expected_power,
            "derived_reactive_power_var": expected_reactive,
            "derived_power_factor": expected_pf,
            "maximum_model_expectation_relative_error": max(expected_model_errors),
            "maximum_phasor_replay_relative_error": replay_error,
        },
        "lesson": (
            "For a balanced Y source feeding an equal delta RL load, verify the "
            "ABC sequence, sqrt(3) voltage and current identities, branch impedance, "
            "complex-power conservation, and the cancellation of twice-line-frequency "
            "instantaneous power. A single unequal branch must break the balance gate."
        ),
    }
