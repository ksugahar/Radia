import copy

from radia.ltspice.bipolar_startup_gate import bipolar_supply_startup_gate
from radia.ltspice.mcp_server import bipolar_supply_startup_gate as mcp_gate


def good():
    time_s = [index * 1.0e-5 for index in range(201)]
    positive = [12.0 * min(value / 8.0e-4, 1.0) for value in time_s]
    negative = [-12.0 * min(value / 1.2e-3, 1.0) for value in time_s]
    pg_positive = [0.0 if value < 8.0e-4 else 5.0 for value in time_s]
    pg_negative = [0.0 if value < 1.2e-3 else 5.0 for value in time_s]
    return {
        "time_unit": "s",
        "voltage_unit": "V",
        "target_magnitude_v": 12.0,
        "time_s": time_s,
        "positive_v": positive,
        "negative_v": negative,
        "power_good_positive_v": pg_positive,
        "power_good_negative_v": pg_negative,
    }


def test_accepts_ordered_bipolar_startup():
    result = bipolar_supply_startup_gate(good())
    assert result["status"] == "ok"
    assert result["checks"]["power_good_asserts_after_rail_t90"] is True
    assert mcp_gate(good())["status"] == "ok"


def test_rejects_wrong_negative_polarity_and_rail_imbalance():
    payload = copy.deepcopy(good())
    payload["negative_v"] = [-value * 0.7 for value in payload["negative_v"]]
    result = bipolar_supply_startup_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["final_negative_regulated"] is False
    assert result["checks"]["rail_magnitudes_balanced"] is False


def test_rejects_power_good_before_regulation():
    payload = copy.deepcopy(good())
    payload["power_good_negative_v"] = [0.0 if value < 2.0e-4 else 5.0 for value in payload["time_s"]]
    result = bipolar_supply_startup_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["power_good_asserts_after_rail_t90"] is False
