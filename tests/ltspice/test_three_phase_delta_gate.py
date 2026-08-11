from __future__ import annotations

import copy
import math

from radia.ltspice.mcp_server import balanced_three_phase_delta_load_gate
from radia.ltspice.three_phase_delta_gate import (
    balanced_three_phase_delta_rl_gate,
)


def _summary() -> dict:
    phase_voltage = 100.0
    frequency = 50.0
    resistance = 24.0
    inductance = 0.05729
    reactance = 2.0 * math.pi * frequency * inductance
    impedance = math.hypot(resistance, reactance)
    line_voltage = math.sqrt(3.0) * phase_voltage
    branch_current = line_voltage / impedance
    line_current = math.sqrt(3.0) * branch_current
    power = 3.0 * line_voltage**2 * resistance / impedance**2
    reactive = 3.0 * line_voltage**2 * reactance / impedance**2
    power_factor = resistance / impedance
    positive = {
        "point_count": 1042,
        "fit_window_start_s": 0.06,
        "fit_window_stop_s": 0.1,
        "phase_voltage_rms_v": [phase_voltage] * 3,
        "line_voltage_rms_v": [line_voltage] * 3,
        "branch_current_rms_a": [branch_current] * 3,
        "line_current_rms_a": [line_current] * 3,
        "branch_impedance_ohm": [[resistance, reactance]] * 3,
        "maximum_phase_voltage_relative_error": 1.0e-9,
        "maximum_line_voltage_relative_error": 1.0e-9,
        "maximum_branch_current_relative_error": 6.0e-5,
        "maximum_line_current_relative_error": 6.0e-5,
        "maximum_branch_impedance_relative_error": 6.0e-5,
        "line_current_kcl_relative_error": 1.0e-9,
        "phase_voltage_positive_to_negative_sequence_ratio": 1.0e-9,
        "phase_voltage_zero_sequence_ratio": 1.0e-9,
        "branch_current_magnitude_spread_relative": 1.0e-9,
        "line_current_magnitude_spread_relative": 1.0e-9,
        "source_complex_power_va": [power, reactive],
        "load_complex_power_va": [power, reactive],
        "source_load_complex_power_relative_error": 1.0e-9,
        "active_power_relative_error": 7.0e-5,
        "reactive_power_relative_error": 4.0e-5,
        "power_factor": power_factor,
        "power_factor_absolute_error": 1.0e-5,
        "instantaneous_power_mean_relative_error": 7.0e-5,
        "instantaneous_power_ripple_relative": 1.0e-6,
        "maximum_phasor_fit_relative_residual": 1.0e-7,
    }
    return {
        "model_contract": {
            "topology": "balanced_y_source_delta_rl_load",
            "source_rows": [
                {"peak_voltage_v": phase_voltage * math.sqrt(2), "frequency_hz": frequency, "phase_deg": 0.0},
                {"peak_voltage_v": phase_voltage * math.sqrt(2), "frequency_hz": frequency, "phase_deg": -120.0},
                {"peak_voltage_v": phase_voltage * math.sqrt(2), "frequency_hz": frequency, "phase_deg": 120.0},
            ],
            "phase_voltage_rms_v": phase_voltage,
            "frequency_hz": frequency,
            "phase_sequence": "abc",
            "delta_resistances_ohm": [resistance] * 3,
            "delta_inductances_h": [inductance] * 3,
            "branch_resistance_ohm": resistance,
            "branch_inductance_h": inductance,
            "expected_branch_impedance_ohm": [resistance, reactance],
            "expected_line_voltage_rms_v": line_voltage,
            "expected_branch_current_rms_a": branch_current,
            "expected_line_current_rms_a": line_current,
            "expected_active_power_w": power,
            "expected_reactive_power_var": reactive,
            "expected_power_factor": power_factor,
        },
        "metrics": {
            "maximum_phasor_replay_relative_error": 0.0,
            "positive": positive,
        },
        "timing_breakdown_s": {
            "preflight": 0.01,
            "solve": 0.1,
            "analyze": 0.01,
            "serialize": 0.01,
        },
    }


def test_accepts_balanced_delta_sequence_currents_and_power() -> None:
    result = balanced_three_phase_delta_rl_gate(_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())
    assert balanced_three_phase_delta_load_gate(_summary())["status"] == "ok"


def test_rejects_unbalanced_branch_current_and_power_ripple() -> None:
    bad = copy.deepcopy(_summary())
    metrics = bad["metrics"]["positive"]
    metrics["branch_current_rms_a"][1] *= 0.86
    metrics["branch_impedance_ohm"][1][0] = 30.0
    metrics["maximum_branch_current_relative_error"] = 0.14
    metrics["maximum_line_current_relative_error"] = 0.10
    metrics["maximum_branch_impedance_relative_error"] = 0.20
    metrics["branch_current_magnitude_spread_relative"] = 0.15
    metrics["line_current_magnitude_spread_relative"] = 0.10
    metrics["instantaneous_power_ripple_relative"] = 0.15
    result = balanced_three_phase_delta_rl_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["delta_impedance_branch_and_line_current_identities_close"] is False
    assert result["checks"]["balanced_branch_and_line_current_magnitudes_close"] is False
    assert result["checks"]["balanced_instantaneous_three_phase_power_is_constant"] is False


def test_rejects_wrong_phase_sequence_and_short_fit_window() -> None:
    bad = copy.deepcopy(_summary())
    bad["model_contract"]["source_rows"][1]["phase_deg"] = 120.0
    bad["model_contract"]["source_rows"][2]["phase_deg"] = -120.0
    bad["metrics"]["positive"]["fit_window_start_s"] = 0.08
    result = balanced_three_phase_delta_rl_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["source_amplitudes_frequency_and_abc_phases_close"] is False
    assert result["checks"]["steady_state_phasor_fit_has_two_periods_and_small_residual"] is False
