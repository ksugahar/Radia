from __future__ import annotations

import copy

from radia_mcp.motor.hdiv_hex_torque_gate import hdiv_hex_motor_torque_gate


SPACE = {
    "family": "HDiv",
    "order": 1,
    "cell_family": "HEX",
    "project_lane": "BDM1",
    "strict_name": "tensor_product_hdiv_order1",
    "simplex_analogue": "BDM1",
}


def _artifact():
    angles = [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0]
    torque = [0.0, -1.18, -2.04, -2.36, -2.04, -1.18, 0.0]
    levels = []
    for elements, ndof, scale in ((40, 1136, 0.995), (320, 8384, 0.999), (1080, 27504, 1.0)):
        moment = [scale * value for value in torque]
        levels.append(
            {
                "angles_deg": angles,
                "hex_element_count": elements,
                "ndof": ndof,
                "torque_moment_Nm": moment,
                "torque_virtual_work_Nm": [value * (1.0 + 1.0e-7) for value in moment],
                "discrete_space": copy.deepcopy(SPACE),
                "operator_build_count": 1,
            }
        )
    identity = "a" * 64
    return {
        "schema": "radia.hdiv-hex-motor-torque-evidence.v1",
        "lane": "hdiv_mmm",
        "physics_identity_sha256": identity,
        "levels": levels,
        "reference": {
            "executed": True,
            "physics_identity_sha256": identity,
            "torque_Nm": [value * 1.002 for value in torque],
        },
    }


def test_hdiv_hex_motor_torque_gate_accepts_converged_same_identity_evidence():
    result = hdiv_hex_motor_torque_gate(_artifact())

    assert result["status"] == "ok"
    assert all(result["checks"].values())
    assert result["metrics"]["final_hex_element_count"] == 1080


def test_hdiv_hex_motor_torque_gate_rejects_non_hex_or_order_drift():
    artifact = _artifact()
    artifact["levels"][1]["discrete_space"]["cell_family"] = "TRI"
    artifact["levels"][2]["discrete_space"]["order"] = 2

    result = hdiv_hex_motor_torque_gate(artifact)

    assert result["status"] == "needs_attention"
    assert result["checks"]["all_levels_use_bdm1_project_lane_on_hex"] is False


def test_hdiv_hex_motor_torque_gate_rejects_identity_and_torque_route_breaks():
    artifact = _artifact()
    artifact["reference"]["physics_identity_sha256"] = "b" * 64
    artifact["levels"][-1]["torque_virtual_work_Nm"][3] = -1.0

    result = hdiv_hex_motor_torque_gate(artifact)

    assert result["status"] == "needs_attention"
    assert result["checks"]["reference_is_executed_same_identity"] is False
    assert result["checks"]["moment_and_virtual_work_torque_agree"] is False
