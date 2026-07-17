# -*- coding: utf-8 -*-
"""Fast tests for the radia-motor validation-lane contract."""

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
        "shared_mesh_material_identity": {
            "geometry_sha256": "1" * 64,
            "material_sha256": "2" * 64,
            "excitation_sha256": "3" * 64,
        },
        "solver_ready_artifact": {
            "artifact_id": f"{lane}_fixture_v1",
            "verification": ["pytest solver-ready fixture"],
        },
    }
    if lane == "hdiv_mmm_hcurl_eddy_bubble":
        data["metrics"] = {"max_abs_relative_error": 1.0e-3}
        data["coupling_design_status"] = "solver_validated"
        data["hdiv_mmm_operator_contract"] = {
            "space": "HDiv",
            "quantity": observable,
        }
        data["hcurl_eddy_bubble_contract"] = {
            "space": "HCurl",
            "basis": "eddy_bubble",
        }
        data["coupling_operator_contract"] = {
            "blocks": ["hdiv_mmm", "hcurl_eddy_bubble"],
        }
        data["solver_ready_artifact"] = {
            "artifact_id": "hdiv_mmm_hcurl_eddy_bubble_fixture_v1",
            "verification": ["pytest tests/test_vim_eddy_hybrid.py -q"],
        }
    else:
        data["metrics"] = {"torque_relative_error": 1.0e-3}
        data["age_gate_ids"] = ["age_rotation_torque"]
        data["pytest_targets"] = ["tests/test_airgap_machine_rotation.py"]
    return data


def test_validation_lane_report_names_independent_radia_paths():
    report = format_motor_validation_lanes("lane_matrix")
    assert "HDiv-MMM + HCurl eddy-bubble" in report
    assert "NGSolve+AGE" in report
    assert "pickup_flux" in report
    assert "linear_pm_flux" in report
    assert "linear_thrust" in report
    assert "motor_family_sweep" in report
    assert "rotary_flux_linkage" in report
    assert "cogging_torque" in report


def test_lane_templates_expose_required_artifact_contracts():
    hdiv = lane_template("hdiv_mmm_hcurl_eddy_bubble")
    age = lane_template("ngsolve_age")
    assert "hdiv_mmm_operator_contract" in hdiv["required_fields"]
    assert "hcurl_eddy_bubble_contract" in hdiv["required_fields"]
    assert "coupling_operator_contract" in hdiv["required_fields"]
    assert "shared_mesh_material_identity" in hdiv["required_fields"]
    assert hdiv["support_status"] == "required_validation_path"
    assert "age_gate_ids" in age["required_fields"]
    assert "pytest_targets" in age["required_fields"]
    assert "solver_ready_artifact" in age["required_fields"]
    assert "shared_mesh_material_identity" in age["required_fields"]
    assert age["support_status"] == "supported_validation_path"
    assert (
        "tests/test_vim_eddy_hybrid.py::test_eddy_bubble_hcurl_basis_is_vim_and_hdiv_mmm_ready"
        in hdiv["public_evidence"]
    )
    assert "validation_test/cln/hcurl_vim_hdiv_mmm_end_to_end.py" in hdiv["public_evidence"]


def test_hdiv_mmm_eddy_artifact_gate_accepts_solver_ready_contract():
    artifact = _base_artifact("hdiv_mmm_hcurl_eddy_bubble", "pickup_flux")
    result = validate_motor_validation_artifact(
        artifact, "hdiv_mmm_hcurl_eddy_bubble"
    )
    assert result["status"] == "pass"
    assert result["accepted_for_mcp_learning"] is True
    assert result["accepted_for_mcp_rfc_learning"] is False
    assert result["support_status"] == "required_validation_path"
    assert result["validated_solver_path"] is False
    assert result["validated_experimental_solver_path"] is False
    assert result["validated_required_solver_path"] is True
    text = format_artifact_gate_result(result)
    assert "accepted for MCP learning: `True`" in text
    assert "validated required solver path: `True`" in text


def test_legacy_hdiv_vim_alias_is_parseable_but_cannot_train():
    artifact = _base_artifact("hdiv_mmm_hcurl_eddy_bubble", "pickup_flux")
    artifact["motor_validation_lane"] = "hdiv_vim_reduced_fem"
    result = validate_motor_validation_artifact(artifact, "hdiv_vim_reduced_fem")
    assert result["status"] == "pass"
    assert result["lane"] == "hdiv_mmm_hcurl_eddy_bubble"
    assert result["accepted_for_mcp_learning"] is False
    assert result["warnings"]


def test_linear_motor_observables_are_dual_lane_training_targets():
    age = _base_artifact("ngsolve_age", "linear_thrust")
    age["metrics"] = {"field_relative_error": 1.0e-3}
    age_result = validate_motor_validation_artifact(age, "ngsolve_age")

    hdiv = _base_artifact("hdiv_mmm_hcurl_eddy_bubble", "linear_pm_flux")
    hdiv_result = validate_motor_validation_artifact(
        hdiv, "hdiv_mmm_hcurl_eddy_bubble"
    )

    assert age_result["status"] == "pass"
    assert age_result["accepted_for_mcp_learning"] is True
    assert hdiv_result["status"] == "pass"
    assert hdiv_result["accepted_for_mcp_learning"] is True


def test_rotary_motor_family_sweep_is_dual_lane_training_target():
    age = _base_artifact("ngsolve_age", "motor_family_sweep")
    age["metrics"] = {"quantity_specific_residual": 2.0e-3}
    age_result = validate_motor_validation_artifact(age, "ngsolve_age")

    hdiv = _base_artifact("hdiv_mmm_hcurl_eddy_bubble", "rotary_flux_linkage")
    hdiv_result = validate_motor_validation_artifact(
        hdiv, "hdiv_mmm_hcurl_eddy_bubble"
    )

    assert age_result["status"] == "pass"
    assert age_result["accepted_for_mcp_learning"] is True
    assert hdiv_result["status"] == "pass"
    assert hdiv_result["accepted_for_mcp_learning"] is True


def test_hdiv_mmm_eddy_saliency_contract_can_train_mixed_lane():
    artifact = _base_artifact(
        "hdiv_mmm_hcurl_eddy_bubble", "force_or_torque_trend"
    )
    artifact["metrics"] = {
        "signed_agreement_count": 2,
        "max_abs_relative_error": 2.0e-2,
    }
    result = validate_motor_validation_artifact(
        artifact, "hdiv_mmm_hcurl_eddy_bubble"
    )

    assert result["status"] == "pass"
    assert result["validated_required_solver_path"] is True
    assert result["accepted_for_mcp_learning"] is True


def test_hdiv_mmm_eddy_artifact_requires_verification_list():
    artifact = _base_artifact("hdiv_mmm_hcurl_eddy_bubble", "pickup_flux")
    artifact["solver_ready_artifact"] = {
        "artifact_id": "hdiv_mmm_eddy_pickup_flux_solver_ready_v1",
        "verification": "todo",
    }

    result = validate_motor_validation_artifact(
        artifact, "hdiv_mmm_hcurl_eddy_bubble"
    )

    assert result["status"] == "fail"
    assert result["accepted_for_mcp_learning"] is False
    assert any("verification list" in item for item in result["errors"])


def test_ngsolve_age_artifact_gate_accepts_torque_contract_json():
    artifact = _base_artifact("ngsolve_age", "torque")
    result = validate_motor_validation_artifact(json.dumps(artifact), "ngsolve_age")
    assert result["status"] == "pass"
    assert result["accepted_for_mcp_learning"] is True
    assert result["validated_solver_path"] is True
    assert result["validated_supported_path"] is True


def test_ngsolve_age_artifact_gate_requires_solver_ready_execution():
    artifact = _base_artifact("ngsolve_age", "torque")
    artifact.pop("solver_ready_artifact")

    result = validate_motor_validation_artifact(artifact, "ngsolve_age")

    assert result["status"] == "fail"
    assert result["validated_solver_path"] is False
    assert result["accepted_for_mcp_learning"] is False
    assert any("solver_ready_artifact" in error for error in result["errors"])


def test_artifact_gate_rejects_mixed_lane_observable():
    artifact = _base_artifact("hdiv_mmm_hcurl_eddy_bubble", "torque")
    result = validate_motor_validation_artifact(
        artifact, "hdiv_mmm_hcurl_eddy_bubble"
    )
    assert result["status"] == "fail"
    assert any("observable_family" in err for err in result["errors"])


def test_validation_lane_docs_do_not_embed_private_absolute_paths():
    docs = format_motor_validation_lanes("all")
    for drive in ("S", "W", "C"):
        assert f"{drive}:" + "\\" not in docs
    for private_source in ("COMSOL", "FEMM", "JMAG", "CST", "ELF/MAGIC"):
        assert private_source not in docs
    assert "product_local_reference" in docs
    assert "hcurl_vim_hdiv_mmm_end_to_end.py" in docs
