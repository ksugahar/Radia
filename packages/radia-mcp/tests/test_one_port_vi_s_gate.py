import copy
import json
import math

from radia_mcp.radia_ngsolve.server import one_port_vi_s_impedance_gate


def c(value):
    return {"real": value.real, "imag": value.imag}


def good():
    rows = []
    power_rows = []
    for index in range(7):
        frequency = 1.0e6 + index * 1.0e6
        zref = 50.0 + 0.0j
        impedance = 2.0 + 1j * (10.0 + index)
        s11 = (impedance - zref) / (impedance + zref)
        scale = math.sqrt(zref.real)
        voltage = scale * (1.0 + s11)
        current = (1.0 - s11) / scale
        rows.append({"frequency_hz": frequency, "s11": c(s11), "zref_ohm": c(zref), "voltage_v": c(voltage), "current_a": c(current)})
        power_rows.append({"frequency_hz": frequency, "accepted_power_w": 0.5 * (1.0 - abs(s11) ** 2), "stimulated_power_w": 0.5, "balance_magnitude": abs(s11)})
    return {"stimulated_power_w": 0.5, "rows": rows, "power_rows": power_rows}


def test_accepts_impedance_and_power_representations():
    result = json.loads(one_port_vi_s_impedance_gate(json.dumps(good())))
    assert result["status"] == "ok"
    assert result["metrics"]["maximum_impedance_relative_error"] < 1.0e-14


def test_rejects_stale_voltage_and_power_rows():
    payload = copy.deepcopy(good())
    payload["rows"][3]["voltage_v"]["real"] += 1.0
    payload["power_rows"][4]["accepted_power_w"] *= 0.5
    result = json.loads(one_port_vi_s_impedance_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["vi_matches_s_impedance_transform"] is False
    assert result["checks"]["accepted_power_matches_s11"] is False


def test_rejects_wrong_reference_impedance_and_balance():
    payload = good()
    payload["rows"][0]["zref_ohm"] = {"real": -50.0, "imag": 0.0}
    payload["power_rows"][0]["balance_magnitude"] += 0.1
    result = json.loads(one_port_vi_s_impedance_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["reference_impedance_positive_real"] is False
    assert result["checks"]["balance_matches_s11_magnitude"] is False


def test_accepts_roundoff_in_sparse_frequency_keys():
    payload = good()
    payload["power_rows"][3]["frequency_hz"] += 1.0e-8
    result = json.loads(one_port_vi_s_impedance_gate(json.dumps(payload)))
    assert result["status"] == "ok"
    assert result["metrics"]["maximum_frequency_match_relative_error"] > 0.0
