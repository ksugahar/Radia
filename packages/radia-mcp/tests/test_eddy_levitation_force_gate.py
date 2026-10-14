import copy
import json

from radia_mcp.radia_ngsolve.server import linear_eddy_levitation_force_gate


def good() -> dict:
    rows = []
    for current in (10.0, 20.0, 30.0):
        rows.append(
            {
                "current_a": current,
                "lorentz_dc_z_n": 0.00255 * current**2,
                "weighted_stress_dc_z_n": 0.00249 * current**2,
                "resistive_loss_w": 0.03378 * current**2,
                "node_count": 16563,
                "element_count": 32610,
            }
        )
    return {
        "contract": {
            "frequency_hz": 50.0,
            "linear_materials": True,
            "force_component": "time_average_dc",
            "target_kind": "conducting_body",
            "weighted_stress_mask": "target_surrounded_by_air",
        },
        "rows": rows,
    }


def test_accepts_two_force_routes_and_i2_scaling():
    result = json.loads(linear_eddy_levitation_force_gate(json.dumps(good())))
    assert result["status"] == "ok"
    assert result["checks"]["lorentz_force_obeys_i2"] is True
    assert result["checks"]["force_methods_agree"] is True


def test_rejects_2x_component_and_invalid_stress_mask():
    payload = good()
    payload["contract"]["force_component"] = "two_x_phasor"
    payload["contract"]["weighted_stress_mask"] = "target_abuts_non_air"
    result = json.loads(linear_eddy_levitation_force_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["dc_time_average_force_component"] is False
    assert result["checks"]["weighted_stress_mask_valid"] is False


def test_rejects_broken_i2_law_force_disagreement_and_remesh():
    payload = copy.deepcopy(good())
    payload["rows"][2]["lorentz_dc_z_n"] *= 1.2
    payload["rows"][2]["weighted_stress_dc_z_n"] *= -1.0
    payload["rows"][2]["node_count"] += 1
    result = json.loads(linear_eddy_levitation_force_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["lorentz_force_obeys_i2"] is False
    assert result["checks"]["force_methods_have_same_sign"] is False
    assert result["checks"]["mesh_inventory_positive_and_stable"] is False
