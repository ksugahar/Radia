import copy
import json
import math

from radia_mcp.radia_ngsolve.conductor_frequency_gate import opposed_busbar_skin_force_gate
from radia_mcp.radia_ngsolve.server import opposed_busbar_skin_force_gate as mcp_gate


def rows():
    result = []
    frequencies = (1.0e3, 1.0e4, 3.0e4, 1.0e5)
    resistance = (1.0e-5, 1.1e-5, 1.5e-5, 2.2e-5)
    for replay in (1, 2):
        for index, (frequency, r_ac) in enumerate(zip(frequencies, resistance)):
            skin = math.sqrt(2.0 / (2.0 * math.pi * frequency * 4.0e-7 * math.pi * 58.0e6)) * 1000.0
            current = 500.0
            loss = r_ac * current**2
            flux = 1.0e-7
            inner, center, outer = ((100.0, 100.0, 100.0) if index == 0 else (150.0, 100.0, 90.0) if index < 3 else (300.0, 90.0, 55.0))
            result.append({"replay": replay, "frequency_hz": frequency, "skin_depth_mm": skin, "ac_resistance_ohm": r_ac, "circuit_current_re_a": current, "circuit_current_im_a": 0.0, "voltage_re_v": 2.0 * loss / current, "voltage_im_v": 2.0 * math.pi * frequency * flux, "flux_re_wb_turn": flux, "total_loss_w": loss, "top_loss_w": 0.5 * loss, "bottom_loss_w": 0.5 * loss, "top_current_re_a": current, "top_current_im_a": 0.0, "bottom_current_re_a": -current, "bottom_current_im_a": 0.0, "top_lorentz_y_re_n": 0.01, "bottom_lorentz_y_re_n": -0.01, "energy_j": 1.0e-5, "coenergy_j": 1.0e-5, "inner_j_magnitude_ma_m2": inner, "center_j_magnitude_ma_m2": center, "outer_j_magnitude_ma_m2": outer, "node_count": 1000, "element_count": 1900})
    return result


def test_gate_accepts_phasor_skin_and_force_closure():
    result = opposed_busbar_skin_force_gate(rows(), conductor_thickness_mm=0.5, conductivity_s_per_m=58.0e6, commanded_current_a=500.0)
    assert result["status"] == "ok"
    assert result["checks"]["high_frequency_inner_face_crowding"] is True
    wrapped = json.loads(mcp_gate(rows(), 0.5, 58.0e6, 500.0))
    assert wrapped["status"] == "ok"


def test_gate_rejects_missing_current_and_force_closure():
    bad = copy.deepcopy(rows())
    bad[0]["bottom_current_re_a"] = -450.0
    bad[0]["bottom_lorentz_y_re_n"] = -0.005
    result = opposed_busbar_skin_force_gate(bad, conductor_thickness_mm=0.5, conductivity_s_per_m=58.0e6, commanded_current_a=500.0)
    assert result["status"] == "needs_attention"
    assert result["checks"]["opposed_conductor_current_closure"] is False
    assert result["checks"]["lorentz_action_reaction_closure"] is False
