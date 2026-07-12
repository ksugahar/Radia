import copy
import json

from radia_mcp.litz_transmission.proximity_pair_gate import (
    litz_proximity_approximation_pair_gate,
)


def _payload():
    return {
        "schema": "litz.proximity-approximation-pair.v1",
        "frequency_hz": 50000.0,
        "models": {
            "approximate": {
                "current_a": 1.0,
                "total_current_a_turn": 114.0000000000002,
                "impedance_ohm": {"real": 8.36452174059618, "imag": 46.64697852963744},
                "total_loss_w": 4.181591548575002,
                "half_real_vi_star_w": 4.18226087029809,
                "field_samples": [
                    {"id": "inner_axis", "b_abs_t": 0.006522048901208048},
                    {"id": "upper_bundle", "b_abs_t": 0.0024080865517217805},
                    {"id": "outer_axis", "b_abs_t": 0.0009567790961699914},
                ],
                "element_count": 12953,
                "solve_time_s": 0.49786599999060854,
            },
            "exact": {
                "current_a": 1.0,
                "total_current_a_turn": 113.9999999957913,
                "impedance_ohm": {"real": 8.579232608733186, "imag": 46.79852460556663},
                "total_loss_w": 4.289616303413933,
                "half_real_vi_star_w": 4.289616304366593,
                "field_samples": [
                    {"id": "inner_axis", "b_abs_t": 0.006606819966940693},
                    {"id": "upper_bundle", "b_abs_t": 0.002377697467416484},
                    {"id": "outer_axis", "b_abs_t": 0.000991030485664384},
                ],
                "element_count": 90604,
                "solve_time_s": 5.666773400007514,
            },
        },
    }


def test_litz_proximity_pair_accepts_loss_impedance_field_and_cost_evidence():
    result = litz_proximity_approximation_pair_gate(json.dumps(_payload()))
    assert result["status"] == "ok"
    assert result["metrics"]["total_loss_relative_error"] < 0.026
    assert result["metrics"]["max_field_sample_relative_error"] < 0.035
    assert result["checks"]["complex_power_closes_total_loss"] is True


def test_litz_proximity_pair_rejects_loss_only_fit_with_wrong_field_and_power():
    payload = copy.deepcopy(_payload())
    payload["models"]["approximate"]["half_real_vi_star_w"] *= 0.8
    payload["models"]["approximate"]["field_samples"][2]["b_abs_t"] *= 0.7
    result = litz_proximity_approximation_pair_gate(json.dumps(payload))
    assert result["status"] == "needs_attention"
    assert result["checks"]["complex_power_closes_total_loss"] is False
    assert result["checks"]["sampled_field_within_4pct"] is False
