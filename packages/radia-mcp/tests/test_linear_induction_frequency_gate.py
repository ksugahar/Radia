import copy
import json

from radia_mcp.radia_ngsolve.linear_induction_gate import linear_induction_frequency_sweep_gate
from radia_mcp.radia_ngsolve.server import linear_induction_frequency_sweep_gate as mcp_gate


ROWS = [
    {"frequency_hz": 0.25, "lorentz_thrust_n": 3975.980335, "weighted_stress_thrust_n": 3976.591450, "plate_resistive_loss_w": 848.241451, "phase_current_sum_a": {"real": 0.0, "imag": 0.0}, "three_phase_complex_power_va": {"real": 10563.477970, "imag": 4350.368889}, "node_count": 7559, "element_count": 14661},
    {"frequency_hz": 0.5, "lorentz_thrust_n": 5230.224167, "weighted_stress_thrust_n": 5230.760694, "plate_resistive_loss_w": 2239.992087, "phase_current_sum_a": {"real": 0.0, "imag": 0.0}, "three_phase_complex_power_va": {"real": 13346.979235, "imag": 6177.620654}, "node_count": 7559, "element_count": 14661},
    {"frequency_hz": 1.0, "lorentz_thrust_n": 4359.945615, "weighted_stress_thrust_n": 4360.518003, "plate_resistive_loss_w": 3792.030169, "phase_current_sum_a": {"real": 0.0, "imag": 0.0}, "three_phase_complex_power_va": {"real": 16451.055398, "imag": 6735.569602}, "node_count": 7559, "element_count": 14661},
    {"frequency_hz": 2.0, "lorentz_thrust_n": 2510.192332, "weighted_stress_thrust_n": 2510.836152, "plate_resistive_loss_w": 4642.188580, "phase_current_sum_a": {"real": 0.0, "imag": 0.0}, "three_phase_complex_power_va": {"real": 18151.372250, "imag": 7749.177044}, "node_count": 7559, "element_count": 14661},
    {"frequency_hz": 5.0, "lorentz_thrust_n": 767.371869, "weighted_stress_thrust_n": 768.070439, "plate_resistive_loss_w": 5383.066494, "phase_current_sum_a": {"real": 0.0, "imag": 0.0}, "three_phase_complex_power_va": {"real": 19633.128029, "imag": 13943.647505}, "node_count": 7559, "element_count": 14661},
    {"frequency_hz": 10.0, "lorentz_thrust_n": 90.475018, "weighted_stress_thrust_n": 91.169312, "plate_resistive_loss_w": 6548.548234, "phase_current_sum_a": {"real": 0.0, "imag": 0.0}, "three_phase_complex_power_va": {"real": 21964.091677, "imag": 24780.759856}, "node_count": 7559, "element_count": 14661},
    {"frequency_hz": 20.0, "lorentz_thrust_n": -59.945865, "weighted_stress_thrust_n": -59.312586, "plate_resistive_loss_w": 8389.168251, "phase_current_sum_a": {"real": 0.0, "imag": 0.0}, "three_phase_complex_power_va": {"real": 25645.331465, "imag": 44164.613223}, "node_count": 7559, "element_count": 14661},
    {"frequency_hz": 50.0, "lorentz_thrust_n": -3.207994, "weighted_stress_thrust_n": -2.625919, "plate_resistive_loss_w": 10608.163400, "phase_current_sum_a": {"real": 0.0, "imag": 0.0}, "three_phase_complex_power_va": {"real": 30083.321899, "imag": 99376.906896}, "node_count": 7559, "element_count": 14661},
]


def test_linear_induction_live_shape_passes_and_dispatches():
    result = linear_induction_frequency_sweep_gate(ROWS)
    assert result["status"] == "ok"
    assert result["metrics"]["peak_frequency_hz"] == 0.5
    assert json.loads(mcp_gate(ROWS))["status"] == "ok"


def test_linear_induction_gate_rejects_phase_mesh_and_force_method_drift():
    bad = copy.deepcopy(ROWS)
    bad[3]["phase_current_sum_a"]["real"] = 1.0
    bad[4]["element_count"] += 1
    bad[5]["weighted_stress_thrust_n"] += 10.0
    result = linear_induction_frequency_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["three_phase_currents_balanced"] is False
    assert result["checks"]["mesh_inventory_invariant"] is False
    assert result["checks"]["lorentz_and_weighted_thrust_agree"] is False
