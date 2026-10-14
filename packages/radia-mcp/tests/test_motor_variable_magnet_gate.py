import json

from radia_mcp.motor.server import motor_variable_magnet_material_gate
from radia_mcp.motor.variable_magnet_gate import variable_magnet_material_parameter_gate


PARAMETERS = {"iHc": -100000.0, "Br": 1.5, "mur": 1.2, "mug": 100.0, "Br0": 0.5}


def test_live_shape_variable_magnet_parameters_pass_and_dispatch():
    result = variable_magnet_material_parameter_gate(
        PARAMETERS,
        parameter_authority="saved_solver_report_message",
        study_label_is_parameter_authority=False,
    )
    assert result["status"] == "ok"
    assert result["metrics"]["initial_to_full_remanence_ratio"] == 1.0 / 3.0
    dispatched = json.loads(motor_variable_magnet_material_gate(PARAMETERS, "saved_solver_report_message"))
    assert dispatched["status"] == "ok"


def test_variable_magnet_gate_rejects_study_label_inference_and_bad_sign():
    bad = dict(PARAMETERS)
    bad["iHc"] = 100000.0
    result = variable_magnet_material_parameter_gate(
        bad,
        parameter_authority="study_label",
        study_label_is_parameter_authority=True,
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["coercive_field_uses_negative_internal_sign"] is False
    assert result["checks"]["study_label_not_used_as_parameter_value"] is False
