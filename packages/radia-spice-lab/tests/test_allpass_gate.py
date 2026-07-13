from __future__ import annotations

import copy

from ltspice_converter.allpass_gate import second_order_allpass_gate
from ltspice_converter.mcp_server import second_order_allpass_phase_group_delay_gate as mcp_gate


def _summary() -> dict:
    run = {
        "point_count": 203,
        "frequency_min_hz": 100.0,
        "frequency_max_hz": 10000.0,
        "minimum_magnitude": 99.99979,
        "maximum_magnitude": 99.99999,
        "center_frequency_sample_hz": 1000.0,
        "center_phase_error_deg": 4.0e-11,
        "low_high_phase_sum_error_deg": 8.4e-12,
        "group_delay_at_center_s": 0.000445598,
        "phase_monotonic_violation_rad": 0.0,
    }
    return {
        "model_contract": {"topology": "second_order_allpass", "center_frequency_hz": 1000.0, "quality_factor": 0.7, "nominal_gain": 100.0},
        "runs": [copy.deepcopy(run), copy.deepcopy(run)],
        "pole_zero": {"poles": [[-4487.99, 4487.99], [-4487.99, -4487.99]], "zeros": [[4487.99, 4487.99], [4487.99, -4487.99]], "mirror_relative_error": 0.0},
        "metrics": {"maximum_analytic_complex_relative_l2": 1.41e-6, "maximum_complex_replay_relative_error": 0.0},
        "timing_breakdown_s": {"preflight": 0.1, "solve": 1.0, "parse": 0.1, "serialize": 0.1},
    }


def test_accepts_complex_allpass_phase_and_group_delay() -> None:
    result = second_order_allpass_gate(_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())
    assert mcp_gate(_summary())["status"] == "ok"


def test_rejects_flat_magnitude_only_false_allpass() -> None:
    bad = copy.deepcopy(_summary())
    bad["runs"][1]["center_phase_error_deg"] = 90.0
    bad["runs"][1]["group_delay_at_center_s"] = 1.0e-5
    bad["pole_zero"]["zeros"][0][0] *= -1.0
    result = second_order_allpass_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["center_sample_and_minus_180_phase_are_resolved"] is False
    assert result["checks"]["group_delay_matches_four_q_over_omega0"] is False
    assert result["checks"]["stable_poles_and_right_half_plane_zeros_are_mirrored"] is False


def test_rejects_inadequate_frequency_grid_replay_and_timing() -> None:
    bad = _summary()
    bad["runs"][0]["point_count"] = 20
    bad["metrics"]["maximum_complex_replay_relative_error"] = 0.1
    bad["timing_breakdown_s"].pop("serialize")
    result = second_order_allpass_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["two_dense_ac_replays_cover_two_decades"] is False
    assert result["checks"]["complex_observable_replay_is_deterministic"] is False
    assert result["checks"]["exactly_four_timing_stages"] is False
