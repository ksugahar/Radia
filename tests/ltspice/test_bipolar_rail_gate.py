import copy

import pytest

from radia.ltspice.bipolar_rail_gate import bipolar_rail_power_quality_gate
from radia.ltspice.mcp_server import bipolar_rail_power_quality_gate as mcp_gate


def good_summary():
    return {
        "units": {"voltage": "V", "power": "W", "efficiency": "percent", "time": "s"},
        "measure_window_s": [0.025, 0.030],
        "target_rail_voltage_v": 12.0,
        "positive_output_voltage_v": 12.204829375,
        "negative_output_voltage_v": -12.2048290959,
        "positive_output_ripple_pp_v": 0.014965057373,
        "negative_output_ripple_pp_v": 0.0149660110474,
        "delivered_input_power_w": 13.5698828388,
        "positive_output_power_w": 6.10241468752,
        "negative_output_power_w": 6.10241454795,
        "reported_efficiency_percent": 89.9405645609,
    }


def test_accepts_balanced_bipolar_rail_power_quality():
    result = bipolar_rail_power_quality_gate(good_summary())
    assert result["status"] == "ok"
    assert result["checks"]["bipolar_polarity_is_signed"] is True
    assert result["metrics"]["rail_imbalance_relative"] < 1.0e-7
    assert mcp_gate(good_summary())["status"] == "ok"


def test_rejects_solver_complete_but_unregulated_variant():
    row = good_summary()
    row.update(
        positive_output_voltage_v=0.0425099700206,
        negative_output_voltage_v=-0.0425099700206,
        positive_output_ripple_pp_v=0.284551682766,
        negative_output_ripple_pp_v=0.284551682766,
        delivered_input_power_w=0.0994088947388,
        positive_output_power_w=0.00870372775311,
        negative_output_power_w=0.00870372775311,
        reported_efficiency_percent=17.5109637341,
    )
    result = bipolar_rail_power_quality_gate(row)
    assert result["status"] == "needs_attention"
    assert result["checks"]["both_rails_meet_target_regulation"] is False
    assert result["checks"]["both_rail_ripples_are_small"] is False


def test_rejects_wrong_polarity_and_nonpassive_power():
    row = copy.deepcopy(good_summary())
    row["negative_output_voltage_v"] *= -1.0
    row["negative_output_power_w"] = 8.0
    result = bipolar_rail_power_quality_gate(row)
    assert result["checks"]["bipolar_polarity_is_signed"] is False
    assert result["checks"]["power_balance_is_passive"] is False


def test_rejects_relaxed_limits_and_negative_ripple():
    relaxed = good_summary()
    relaxed["max_ripple_fraction"] = 1.0
    with pytest.raises(ValueError, match="policy maxima"):
        bipolar_rail_power_quality_gate(relaxed)

    negative = good_summary()
    negative["positive_output_ripple_pp_v"] = -1.0e-3
    with pytest.raises(ValueError, match="nonnegative"):
        bipolar_rail_power_quality_gate(negative)
