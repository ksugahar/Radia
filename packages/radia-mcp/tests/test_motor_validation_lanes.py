# -*- coding: utf-8 -*-
"""Fast tests for the radia-motor dual validation-lane contract."""

import json
import os
import sys


_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.motor.validation_lanes_knowledge import (  # noqa: E402
    format_artifact_gate_result,
    format_motor_validation_lanes,
    lane_template,
    validate_motor_validation_artifact,
)


def _base_artifact(lane: str, observable: str) -> dict:
    data = {
        "schema_version": "radia-motor-validation-artifact/v1",
        "timestamp_utc": "2026-07-03T00:00:00Z",
        "radia_version": "test",
        "motor_validation_lane": lane,
        "reference_source_class": "analytic_reference",
        "observable_family": observable,
        "case_count": 3,
        "status": "pass",
        "tolerances": {"relative_error": 1.0e-2},
        "metrics": {},
        "timing_breakdown_s": {"setup": 0.001, "solve": 0.002},
        "artifact_feedback": {
            "status": "candidate",
            "public_lesson": "dual-lane artifact contract has complete metadata",
        },
    }
    if lane == "hdiv_vim_reduced_fem":
        data["metrics"] = {"max_abs_relative_error": 1.0e-3}
        data["reduced_fem_contract"] = {"basis": "P1", "quantity": observable}
        data["vim_operator_contract"] = {"space": "HDiv", "quantity": observable}
    else:
        data["metrics"] = {"torque_relative_error": 1.0e-3}
        data["age_gate_ids"] = ["age_rotation_torque"]
        data["pytest_targets"] = ["tests/test_airgap_machine_rotation.py"]
    return data


def test_dual_lane_report_names_independent_radia_paths():
    report = format_motor_validation_lanes("lane_matrix")
    assert "HDiv-VIM + reduced FEM" in report
    assert "NGSolve+AGE" in report
    assert "pickup_flux" in report
    assert "cogging_torque" in report


def test_lane_templates_expose_required_artifact_contracts():
    hdiv = lane_template("hdiv_vim_reduced_fem")
    age = lane_template("ngsolve_age")
    assert "vim_operator_contract" in hdiv["required_fields"]
    assert "reduced_fem_contract" in hdiv["required_fields"]
    assert "age_gate_ids" in age["required_fields"]
    assert "pytest_targets" in age["required_fields"]


def test_hdiv_vim_artifact_gate_accepts_pickup_flux_contract():
    artifact = _base_artifact("hdiv_vim_reduced_fem", "pickup_flux")
    result = validate_motor_validation_artifact(artifact, "hdiv_vim_reduced_fem")
    assert result["status"] == "pass"
    assert result["accepted_for_mcp_learning"] is True
    text = format_artifact_gate_result(result)
    assert "accepted for MCP learning: `True`" in text


def test_ngsolve_age_artifact_gate_accepts_torque_contract_json():
    artifact = _base_artifact("ngsolve_age", "torque")
    result = validate_motor_validation_artifact(json.dumps(artifact), "ngsolve_age")
    assert result["status"] == "pass"
    assert result["accepted_for_mcp_learning"] is True


def test_artifact_gate_rejects_mixed_lane_observable():
    artifact = _base_artifact("hdiv_vim_reduced_fem", "torque")
    result = validate_motor_validation_artifact(artifact, "hdiv_vim_reduced_fem")
    assert result["status"] == "fail"
    assert any("observable_family" in err for err in result["errors"])


def test_validation_lane_docs_do_not_embed_private_absolute_paths():
    docs = format_motor_validation_lanes("all")
    for drive in ("S", "W", "C"):
        assert f"{drive}:" + "\\" not in docs
    for private_source in ("COMSOL", "FEMM", "JMAG", "CST", "ELF/MAGIC"):
        assert private_source not in docs
    assert "product_local_reference" in docs
