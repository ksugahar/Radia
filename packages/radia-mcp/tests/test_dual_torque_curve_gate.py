import copy
import json
import math

from radia_mcp.motor.server import motor_dual_torque_method_curve_gate


def good():
    angles = [index * 5.0 for index in range(19)]
    primary = [5.0 * math.sin(math.radians(2.0 * angle)) for angle in angles]
    secondary = [value * 1.004 for value in primary]
    return {
        "angles_deg": angles,
        "primary_torque_nm": primary,
        "secondary_torque_nm": secondary,
        "angle_unit": "deg",
        "torque_unit": "N*m",
    }


def test_accepts_two_consistent_static_torque_curves():
    result = json.loads(motor_dual_torque_method_curve_gate(json.dumps(good())))
    assert result["status"] == "ok"
    assert result["checks"]["peak_angles_agree"] is True


def test_rejects_stale_curve_and_large_endpoint_residual():
    payload = copy.deepcopy(good())
    payload["secondary_torque_nm"][9] += 1.0
    payload["secondary_torque_nm"][-1] = 0.5
    result = json.loads(motor_dual_torque_method_curve_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["pointwise_method_agreement"] is False
    assert result["checks"]["low_torque_endpoints"] is False


def test_rejects_shifted_peak_and_wrong_units():
    payload = good()
    payload["secondary_torque_nm"] = payload["secondary_torque_nm"][3:] + [0.0, 0.0, 0.0]
    payload["torque_unit"] = "N"
    result = json.loads(motor_dual_torque_method_curve_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["peak_angles_agree"] is False
    assert result["checks"]["torque_unit_newton_metre"] is False
