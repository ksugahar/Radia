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
        "pole_count": 4,
        "torque_basis": "total",
        "stack_length_m": 0.05,
        "model_dimension": "2d",
        "model_fraction": 1,
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


def test_equal_dc_bias_does_not_pass_cogging():
    """Guards: detect nonzero mean even when the two traces match.
    Assumes: equal constant torque offsets on uniform periodic traces.
    """
    first = _lane()
    first["torque_nm"] = [value + 0.5 for value in first["torque_nm"]]
    result = cogging_torque_cross_lane_gate(first, first, _LIMITS)
    assert result["status"] == "mismatch"
    assert result["mean_torque_nm"] == pytest.approx([0.5, 0.5])
    assert not result["checks"]["zero_mean_torque"]
    assert "nonzero_mean_torque" in result["findings"]
    first["excitation"] = {"phase_current_A": 1}
    assert cogging_torque_cross_lane_gate(first, first, _LIMITS)["status"] == "invalid_input"


@pytest.mark.parametrize("fraction", [1, 0.25])
def test_convert_per_metre_and_model_fraction_once(fraction):
    """Guards: normalize length and sector factors without folding either twice.
    Assumes: equivalent full-machine torque with declared 2D and 3D bases.
    """
    first, second = _lane(), _lane()
    first.update(torque_basis="per_metre", model_fraction=fraction)
    first["torque_nm"] = [value * fraction / 0.05 for value in first["torque_nm"]]
    second["model_dimension"] = "3d"
    result = cogging_torque_cross_lane_gate(first, second, _LIMITS)
    assert result["status"] == "pass"
    assert result["findings"] == ["model_difference_declared_required"]
    assert result["torque_to_full_machine_factors"] == pytest.approx([0.05 / fraction, 1])


def test_stack_length_scaling_diagnosis_in_both_lane_orders():
    """Guards: report a likely misdeclared stack conversion without correcting it.
    Assumes: same signed waveform scaled by the declared stack length.
    """
    first, second = _lane(), _lane()
    second["torque_nm"] = [value * 0.05 for value in second["torque_nm"]]
    for left, right in ((first, second), (second, first)):
        result = cogging_torque_cross_lane_gate(left, right, _LIMITS)
        assert result["status"] == "mismatch"
        assert result["findings"] == ["stack_scaling_suspected"]


def test_pole_count_factor_two_is_diagnosed_before_grid_validation():
    """Guards: reject pole-count versus pole-pair ambiguity.
    Assumes: otherwise matching lane metadata and declared grids.
    """
    first, second = _lane(), _lane()
    second["pole_count"] = 2
    result = cogging_torque_cross_lane_gate(first, second, _LIMITS)
    assert result["status"] == "invalid_input"
    assert result["findings"] == ["pole_count_vs_pole_pairs_ambiguity"]


def test_none_skew_declarations_are_case_insensitive():
    """Guards: equivalent no-skew declarations compare equally.
    Assumes: case and surrounding whitespace differ.
    """
    first, second = _lane(), _lane()
    second["skew"] = "  None  "
    assert cogging_torque_cross_lane_gate(first, second, _LIMITS)["status"] == "pass"
