from __future__ import annotations

import copy
import json
import math

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from radia_mcp.radia_ngsolve.server import pwm_controlled_motor_loss_gate as mcp_gate


def _payload() -> dict:
    count = 101
    time_s = [index * 1.0e-4 for index in range(count)]
    phase_currents = []
    for index in range(count):
        angle = 2.0 * math.pi * index / (count - 1)
        phase_currents.append(
            [
                10.0 * math.sin(angle),
                10.0 * math.sin(angle - 2.0 * math.pi / 3.0),
                10.0 * math.sin(angle + 2.0 * math.pi / 3.0),
            ]
        )
    power_components = [[10.0 + index, 2.0 - 0.5 * index] for index in range(count)]
    reported_power = [sum(row) for row in power_components]
    eddy_bins = [[3.0, 1.0], [2.0, 0.5], [1.0, 0.5], [0.0, 0.0], [0.0, 0.0]]
    hysteresis = [[2.0, 0.5]] + [[0.0, 0.0] for _ in range(4)]
    iron = [[5.0, 1.5]] + [[0.0, 0.0] for _ in range(4)]
    return {
        "time_series": {
            "time_s": time_s,
            "angle_deg": [3600.0 * value for value in time_s],
            "speed_rpm": [600.0] * count,
            "torque_nm": [20.0 + math.sin(index / 5.0) for index in range(count)],
            "phase_currents_a": phase_currents,
            "control_command_a": [[-10.0, 20.0] for _ in range(count)],
            "control_feedback_a": [[-9.8, 19.6] for _ in range(count)],
            "power_components_w": power_components,
            "reported_total_power_w": reported_power,
        },
        "loss_spectrum": {
            "frequency_hz": [0.0, 100.0, 200.0, 300.0, 400.0],
            "eddy_components_w": eddy_bins,
            "hysteresis_components_w": hysteresis,
            "iron_components_w": iron,
            "reported_eddy_total_w": [sum(row) for row in eddy_bins],
            "reported_hysteresis_total_w": [sum(row) for row in hysteresis],
            "reported_iron_total_w": [sum(row) for row in iron],
        },
        "ratio_diagnostic_policy": "exclude_ratio_when_denominator_current_is_below_floor",
    }


def test_pwm_motor_loss_gate_accepts_balanced_control_and_summary_row_semantics() -> None:
    result = pwm_controlled_motor_loss_gate(_payload())
    assert result["status"] == "ok"
    assert result["metrics"]["three_phase_kcl_global_relative_error"] < 1.0e-12
    assert result["checks"]["eddy_aggregate_matches_harmonic_bin_sum"] is True


def test_pwm_motor_loss_gate_rejects_control_tracking_failure() -> None:
    payload = copy.deepcopy(_payload())
    for row in payload["time_series"]["control_feedback_a"]:
        row[1] = 10.0
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["tail_current_control_tracks_commands"] is False


def test_pwm_motor_loss_gate_rejects_wrong_iron_summary() -> None:
    payload = copy.deepcopy(_payload())
    payload["loss_spectrum"]["iron_components_w"][0][0] += 1.0
    payload["loss_spectrum"]["reported_iron_total_w"][0] += 1.0
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["aggregate_iron_equals_eddy_plus_hysteresis"] is False


def test_pwm_motor_loss_gate_rejects_unguarded_ratio_outputs() -> None:
    payload = _payload()
    payload["ratio_diagnostic_policy"] = "accept_all_resistance_and_inductance_ratios"
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["near_zero_current_ratio_outputs_are_diagnostic_only"] is False


def test_pwm_motor_loss_mcp_tool_dispatches() -> None:
    result = json.loads(mcp_gate(_payload()))
    assert result["status"] == "ok"
