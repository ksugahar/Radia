from __future__ import annotations

import copy

from ltspice_converter.complex_zero_gate import second_order_complex_zero_gate
from ltspice_converter.mcp_server import second_order_complex_zero_transfer_gate


def _summary() -> dict:
    run = {
        "point_count": 203,
        "frequency_min_hz": 100.0,
        "frequency_max_hz": 10000.0,
        "minimum_magnitude": 0.23854646,
        "minimum_magnitude_frequency_hz": 3584.7217,
        "analytic_minimum_magnitude": 0.23854647,
        "analytic_minimum_frequency_hz": 3584.7217,
    }
    return {
        "model_contract": {
            "topology": "second_order_complex_zero",
            "pole_natural_frequency_hz": 1000.0,
            "pole_quality_factor": 2.5,
            "zero_natural_frequency_hz": 2000.0,
            "zero_quality_factor": 1.0,
            "dc_gain": 1.0,
            "high_frequency_gain": 0.25,
        },
        "runs": [copy.deepcopy(run), copy.deepcopy(run)],
        "pole_zero": {
            "poles": [
                [-1256.637061435917, 6156.239184776947],
                [-1256.637061435917, -6156.239184776947],
            ],
            "zeros": [
                [-6283.185307179586, 10882.796185405306],
                [-6283.185307179586, -10882.796185405306],
            ],
        },
        "metrics": {
            "maximum_analytic_complex_relative_l2": 6.1e-7,
            "maximum_analytic_complex_point_relative_error": 1.1e-6,
            "maximum_complex_replay_relative_error": 0.0,
        },
        "timing_breakdown_s": {
            "preflight": 0.1,
            "solve": 1.0,
            "parse": 0.1,
            "serialize": 0.1,
        },
    }


def test_accepts_complex_zero_transfer_roots_asymptotes_and_replay() -> None:
    result = second_order_complex_zero_gate(_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())
    assert second_order_complex_zero_transfer_gate(_summary())["status"] == "ok"


def test_rejects_real_axis_zero_claim_and_wrong_high_frequency_gain() -> None:
    bad = copy.deepcopy(_summary())
    bad["model_contract"]["high_frequency_gain"] = 1.0
    bad["runs"][1]["minimum_magnitude"] = 1.0e-9
    result = second_order_complex_zero_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["dc_and_high_frequency_gain_follow_frequency_ratio"] is False
    assert result["checks"]["complex_zero_has_finite_analytic_real_axis_dip"] is False


def test_rejects_unstable_zero_jacobian_drift_and_replay_drift() -> None:
    bad = copy.deepcopy(_summary())
    bad["pole_zero"]["zeros"][0][0] *= -1.0
    bad["metrics"]["maximum_analytic_complex_relative_l2"] = 1.0e-2
    bad["metrics"]["maximum_complex_replay_relative_error"] = 1.0e-3
    result = second_order_complex_zero_gate(bad)
    assert result["status"] == "needs_attention"
    assert (
        result["checks"]["minimum_phase_conjugate_zeros_recover_fn_and_qn"]
        is False
    )
    assert (
        result["checks"]["complex_transfer_matches_two_pole_two_zero_identity"]
        is False
    )
    assert result["checks"]["complex_observable_replay_is_deterministic"] is False
