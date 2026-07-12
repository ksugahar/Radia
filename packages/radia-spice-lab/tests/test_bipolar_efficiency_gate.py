import copy

from ltspice_converter.bipolar_efficiency_gate import bipolar_converter_efficiency_gate
from ltspice_converter.mcp_server import bipolar_converter_efficiency_gate as mcp_gate


def good_summary():
    return {
        "units": {"power": "W", "efficiency": "percent"},
        "measure_window_s": [0.003, 0.0035],
        "input_source_power_w": -26.623942597,
        "output_1_power_w": 11.996061959,
        "output_2_power_w": 11.991883345,
        "reported_efficiency_percent": 90.0991474747,
        "max_efficiency_closure_percent": 1.0e-6,
        "max_output_imbalance_relative": 0.02,
    }


def test_accepts_signed_balanced_efficiency_measurements():
    result = bipolar_converter_efficiency_gate(good_summary())
    assert result["status"] == "ok"
    assert result["checks"]["reported_efficiency_recomputes"] is True
    assert result["metrics"]["power_loss_w"] > 0.0
    assert mcp_gate(good_summary())["status"] == "ok"


def test_rejects_source_power_with_load_sign():
    payload = copy.deepcopy(good_summary())
    payload["input_source_power_w"] *= -1.0
    result = bipolar_converter_efficiency_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["source_branch_sign_means_delivered_power"] is False
    assert result["checks"]["power_balance_is_passive"] is False


def test_rejects_efficiency_drift_and_unbalanced_outputs():
    payload = copy.deepcopy(good_summary())
    payload["output_2_power_w"] = 8.0
    payload["reported_efficiency_percent"] = 95.0
    result = bipolar_converter_efficiency_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["reported_efficiency_recomputes"] is False
    assert result["checks"]["dual_outputs_are_balanced"] is False
