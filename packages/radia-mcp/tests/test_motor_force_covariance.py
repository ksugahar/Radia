import json

from radia_mcp.motor.force_covariance import evaluate_force_rotation_covariance
from radia_mcp.motor.server import motor_force_rotation_covariance_gate


TOP = {"Fx": 0.1830219608280039, "Fy": 399.8711360045095}
RIGHT = {"Fx": 399.866493590553, "Fy": -0.2079453356343497}


def test_force_rotation_covariance_accepts_live_bearing_pair():
    result = evaluate_force_rotation_covariance(TOP, RIGHT, -90.0, 1.0e-3)
    assert result["status"] == "ok"
    assert result["vector_relative_error"] < 7.0e-5
    assert result["magnitude_relative_error"] < 2.0e-5

    wrapped = json.loads(
        motor_force_rotation_covariance_gate(
            json.dumps(TOP), json.dumps(RIGHT), -90.0, 1.0e-3
        )
    )
    assert wrapped["status"] == "ok"


def test_force_rotation_covariance_rejects_wrong_sign_and_zero_force():
    wrong = evaluate_force_rotation_covariance(TOP, RIGHT, 90.0, 1.0e-3)
    assert wrong["status"] == "needs_attention"
    assert wrong["checks"]["vector_rotates_covariantly"] is False
    zero = evaluate_force_rotation_covariance({"Fx": 0, "Fy": 0}, {"Fx": 0, "Fy": 0}, -90)
    assert zero["status"] == "needs_attention"
    assert zero["checks"]["reference_force_nonzero"] is False
