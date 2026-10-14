import copy
import json

from radia_mcp.radia_ngsolve.server import reciprocal_two_port_power_sweep_gate as mcp_gate
from radia_mcp.radia_ngsolve.two_port_power_gate import reciprocal_two_port_power_sweep_gate


def _rows():
    return [
        {
            "frequency_hz": frequency,
            "s11": [reflection, 0.0],
            "s12": [transmission, 0.0],
            "s21": [transmission, 0.0],
            "s22": [reflection, 0.0],
        }
        for frequency, reflection, transmission in (
            (10.0, 0.99, 0.10),
            (20.0, 0.95, 0.20),
            (30.0, 0.90, 0.30),
            (40.0, 0.85, 0.40),
            (50.0, 0.80, 0.50),
        )
    ]


def _power():
    return [
        {
            "excitation_port": port,
            "frequency_hz": frequency,
            "balance": 0.9999,
            "accepted_power_w": 0.0001,
            "stimulated_power_w": 0.5,
        }
        for port in (1, 2)
        for frequency in (10.0, 30.0, 50.0)
    ]


def test_two_port_power_gate_accepts_reciprocal_passive_sweep_and_dispatches():
    result = reciprocal_two_port_power_sweep_gate(_rows(), _power(), reference_impedance_ohm=50.0)
    assert result["status"] == "ok"
    assert json.loads(mcp_gate(_rows(), _power(), 50.0))["status"] == "ok"


def test_two_port_power_gate_rejects_nonreciprocal_and_bad_balance():
    rows = copy.deepcopy(_rows())
    rows[2]["s12"] = [0.1, 0.0]
    power = copy.deepcopy(_power())
    power[0]["accepted_power_w"] = 0.01
    result = reciprocal_two_port_power_sweep_gate(rows, power, reference_impedance_ohm=50.0)
    assert result["status"] == "needs_attention"
    assert result["checks"]["reciprocity_closes"] is False
    assert result["checks"]["accepted_power_closes_balance"] is False
