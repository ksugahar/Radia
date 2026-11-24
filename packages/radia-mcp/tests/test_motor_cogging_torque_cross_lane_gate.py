"""Surrogates for cross-lane cogging torque comparison."""

import json

import numpy as np
import pytest

from radia_mcp.motor.cogging_torque_cross_lane_gate import cogging_torque_cross_lane_gate


_LIMITS = {"peak_to_peak_nm": 1e-9, "waveform_rms_nm": 1e-9, "harmonic_nm": 1e-9}


def _lane(*, symmetry=True, closed=False, samples=None):
    count = samples or (32 if symmetry else 128)
    period = 30 if symmetry else 360
    angle = 7 + np.arange(count + int(closed)) * period / count
    phase = np.deg2rad(12 * (angle - 7))
    return {
        "slots": 6,
        "poles": 4,
        "slot_pole_symmetry": symmetry,
        "mechanical_angles_deg": angle.tolist(),
        "torque_nm": (0.3 * np.cos(phase) + 0.1 * np.sin(2 * phase)).tolist(),
        "extraction_method": "maxwell_stress",
        "skew": "none",
    }


def test_symmetry_and_mechanical_harmonic_orders():
    """Guards: correct period, harmonic scaling and advisory method differences.
    Assumes: a two-tone waveform sampled over its primitive period.
    """
    first, second = _lane(), _lane()
    second["extraction_method"] = "arkkio"
    result = cogging_torque_cross_lane_gate(first, second, _LIMITS)
    assert result["status"] == "pass"
    assert result["findings"] == ["extraction_methods_differ"]
    assert result["comparison_period_deg_mechanical"] == 30
    assert all(h["order_per_mechanical_revolution"] % 12 == 0 for h in result["harmonics"])
    assert result["harmonics"][0]["order_per_mechanical_revolution"] == 12
    assert result["harmonics"][0]["amplitudes_nm"] == pytest.approx([0.3, 0.3])
    assert result["harmonics"][1]["amplitudes_nm"] == pytest.approx([0.1, 0.1])
    assert result["peak_to_peak_nm"] == pytest.approx([np.ptp(first["torque_nm"])] * 2)


def test_full_turn_when_either_lane_does_not_declare_symmetry():
    """Guards: no sector reduction and full-waveform comparison of off-order content.
    Assumes: a once-per-revolution term in the asymmetric lane.
    """
    first, second = _lane(symmetry=False), _lane(symmetry=False)
    first["slot_pole_symmetry"] = True

    second["torque_nm"] = (np.asarray(second["torque_nm"]) + 0.2 * np.cos(
        np.deg2rad(second["mechanical_angles_deg"])
    )).tolist()
    result = cogging_torque_cross_lane_gate(first, second, _LIMITS)
    assert result["comparison_period_deg_mechanical"] == 360
    assert result["checks"]["slot_pole_harmonics"]
    assert not result["checks"]["signed_waveform_rms"]
    assert result["status"] == "mismatch"
    assert cogging_torque_cross_lane_gate(_lane(), second, _LIMITS)["status"] == "invalid_input"


def test_signed_torque_and_skew_mismatch():
    """Guards: equal peak-to-peak cannot hide inverted torque or skew differences.
    Assumes: identical angular datum and opposite torque sign.
    """
    first, second = _lane(), _lane()
    second["torque_nm"] = [-value for value in second["torque_nm"]]
    second["skew"] = "linear, 1 mechanical degree, axial average"
    result = cogging_torque_cross_lane_gate(first, second, _LIMITS)
    assert result["status"] == "mismatch"
    assert result["findings"] == ["skew_declarations_differ"]
    assert result["checks"]["peak_to_peak"]
    assert not result["checks"]["signed_waveform_rms"]
    assert result["harmonics"][0]["complex_difference_nm"] == pytest.approx(0.6)


def test_endpoint_and_nyquist():
    """Guards: explicit endpoint removal and rejection of unresolved slot/pole order.
    Assumes: uniformly spaced synthetic samples.
    """
    first, second = _lane(closed=True), _lane()
    assert "endpoint" in cogging_torque_cross_lane_gate(first, second, _LIMITS)["error"]
    first["drop_duplicate_endpoint"] = True
    assert cogging_torque_cross_lane_gate(first, second, _LIMITS)["status"] == "pass"
    coarse = _lane(symmetry=False, samples=24)
    assert "Nyquist" in cogging_torque_cross_lane_gate(coarse, coarse, _LIMITS)["error"]


def test_tool_json_and_grid_alignment():
    """Guards: structured JSON results and grid alignment in the public entry point.
    Assumes: callers align the mechanical angle datum before invoking the gate.
    """
    from radia_mcp.motor.server import motor_cogging_torque_cross_lane_gate

    first = json.dumps(_lane())
    limits = json.dumps(_LIMITS)
    assert json.loads(motor_cogging_torque_cross_lane_gate(first, first, limits))["status"] == "pass"
    assert json.loads(motor_cogging_torque_cross_lane_gate("[", first, limits))["checks"] == {"json_inputs": False}
    shifted = _lane()
    shifted["mechanical_angles_deg"] = [value + 0.2 for value in shifted["mechanical_angles_deg"]]
    result = json.loads(motor_cogging_torque_cross_lane_gate(first, json.dumps(shifted), limits))
    assert "common mechanical angle grid" in result["error"]
