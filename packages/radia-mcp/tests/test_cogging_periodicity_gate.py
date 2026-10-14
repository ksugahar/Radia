import copy
import json
import math

from radia_mcp.radia_ngsolve.server import cogging_torque_periodicity_gate


def good() -> dict:
    angles = [0.5 * index for index in range(21)]
    torques = [2.0 * math.sin(2.0 * math.pi * angle / 10.0) for angle in angles]
    return {
        "machine": {"slots": 36, "poles": 4, "phase_currents_a": [0.0, 0.0, 0.0]},
        "torque_observable": {
            "family": "weighted_stress_body_torque",
            "selected_body": "complete_rotor",
        },
        "expected_dominant_period_harmonic": 1,
        "rows": [
            {"angle_mech_deg": angle, "torque_nm": torque}
            for angle, torque in zip(angles, torques)
        ],
    }


def test_accepts_complete_rotor_zero_current_lcm_period():
    result = json.loads(cogging_torque_periodicity_gate(json.dumps(good())))
    assert result["status"] == "ok"
    assert result["metrics"]["expected_period_mech_deg"] == 10.0
    assert result["metrics"]["dominant_period_harmonic"] == 1


def test_rejects_partial_body_and_nonzero_current():
    payload = good()
    payload["torque_observable"]["selected_body"] = "magnets_only"
    payload["machine"]["phase_currents_a"][1] = 1.0
    result = json.loads(cogging_torque_periodicity_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["weighted_stress_selects_complete_rotor"] is False
    assert result["checks"]["zero_current_cogging_condition"] is False


def test_rejects_wrong_period_and_endpoint_drift():
    payload = copy.deepcopy(good())
    payload["rows"][-1]["angle_mech_deg"] = 9.5
    payload["rows"][-1]["torque_nm"] = 1.0
    result = json.loads(cogging_torque_periodicity_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["angle_step_uniform"] is False
    assert result["checks"]["one_lcm_period_covered"] is False
    assert result["checks"]["periodic_endpoint_closure"] is False
