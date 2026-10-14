import json

from radia_mcp.radia_ngsolve.one_port_power_gate import one_port_power_balance_gate
from radia_mcp.radia_ngsolve.server import one_port_power_balance_sweep_gate


def _summary():
    rows = []
    for frequency, s11 in [(50.0, 0.3 + 0.4j), (100.0, 0.1 - 0.2j), (150.0, -0.2 + 0.1j)]:
        stimulated = 0.5
        rows.append(
            {
                "frequency": frequency,
                "s11_real": s11.real,
                "s11_imag": s11.imag,
                "stimulated_power_w": stimulated,
                "accepted_power_w": stimulated * (1.0 - abs(s11) ** 2),
                "balance_magnitude": abs(s11),
                "zref_real_ohm": 50.0,
                "zref_imag_ohm": 0.0,
            }
        )
    return {
        "frequency_unit": "MHz",
        "power_unit": "W",
        "reference_impedance_unit": "ohm",
        "sparameter_basis": "power_wave",
        "rows": rows,
    }


def test_one_port_power_gate_accepts_passive_power_closure():
    result = one_port_power_balance_gate(json.dumps(_summary()))
    assert result["status"] == "ok"
    assert result["metrics"]["maximum_power_relative_residual"] == 0.0
    assert json.loads(one_port_power_balance_sweep_gate(json.dumps(_summary())))["status"] == "ok"


def test_one_port_power_gate_rejects_power_and_balance_drift():
    summary = _summary()
    summary["rows"][1]["accepted_power_w"] += 0.02
    summary["rows"][2]["balance_magnitude"] += 0.1
    result = one_port_power_balance_gate(json.dumps(summary))
    assert result["status"] == "needs_attention"
    assert result["checks"]["accepted_power_closes_from_s11"] is False
    assert result["checks"]["reported_balance_matches_s11"] is False


def test_one_port_power_gate_rejects_active_s11_and_zref_drift():
    summary = _summary()
    summary["rows"][0]["s11_real"] = 1.01
    summary["rows"][0]["s11_imag"] = 0.0
    summary["rows"][0]["accepted_power_w"] = 0.0
    summary["rows"][1]["zref_real_ohm"] = 75.0
    result = one_port_power_balance_gate(json.dumps(summary))
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_rows_passive"] is False
    assert result["checks"]["reference_impedance_is_stable"] is False
