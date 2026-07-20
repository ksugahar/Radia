# -*- coding: utf-8 -*-
"""Fast tests for ELF-seeded radia-motor triple-check planning."""

import json
import os
import sys


_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.motor.triple_check_knowledge import (  # noqa: E402
    format_motor_triple_check_plan,
    format_triple_check_gate_result,
    route_motor_triple_check,
    validate_motor_source_deck_review_packet,
    validate_motor_triple_check_artifact,
)


def _lane_artifact(lane: str, observable: str) -> dict:
    artifact = {
        "schema_version": "radia-motor-validation-artifact/v1",
        "timestamp_utc": "2026-07-04T00:00:00Z",
        "radia_version": "test",
        "motor_validation_lane": lane,
        "reference_source_class": "stored_regression",
        "observable_family": observable,
        "case_count": 2,
        "status": "pass",
        "tolerances": {"relative_error": 1.0e-2},
        "metrics": {},
        "timing_breakdown_s": {"route": 0.001, "gate": 0.001},
        "artifact_feedback": {
            "status": "candidate",
            "public_lesson": "triple-check metadata is complete",
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
        artifact["metrics"] = {"max_abs_relative_error": 1.0e-3}
        artifact["coupling_design_status"] = "solver_validated"
        artifact["hdiv_mmm_operator_contract"] = {
            "space": "HDiv",
            "observable": observable,
        }
        artifact["hcurl_eddy_bubble_contract"] = {
            "space": "HCurl",
            "basis": "eddy_bubble",
        }
        artifact["coupling_operator_contract"] = {
            "blocks": ["hdiv_mmm", "hcurl_eddy_bubble"],
        }
        artifact["solver_ready_artifact"] = {
            "artifact_id": "hdiv_mmm_hcurl_eddy_bubble_fixture_v1",
            "verification": ["pytest tests/test_vim_eddy_hybrid.py -q"],
        }
    else:
        artifact["metrics"] = {"torque_relative_error": 1.0e-3}
        artifact["age_gate_ids"] = ["age_rotation_torque"]
        artifact["pytest_targets"] = ["validation_test/radia_mcp/test_airgap_machine_rotation.py"]
    return artifact


def test_triple_check_plan_requires_hdiv_mmm_hcurl_eddy_bubble():
    plan = route_motor_triple_check("IPM hairpin motor flux linkage and MTPA")
    assert plan["inferred_family"] == "ipm"
    assert plan["standard_comparison"]["primary_required_lanes"] == [
        "ngsolve_age",
        "hdiv_mmm_hcurl_eddy_bubble",
    ]
    assert plan["standard_comparison"]["optional_auxiliary_lanes"] == []
    assert "elf_motor_hybrid_router" in plan["source_mcp_seed"]["calls"][0]
    assert "application/motor/emdlab_ipm_hairpin_10/eip001/eip001.mai" in (
        plan["source_mcp_seed"]["representative_public_decks"]
    )
    assert (
        plan["radia_lanes"]["ngsolve_age"]["support_status"]
        == "supported_validation_path"
    )
    assert (
        plan["radia_lanes"]["hdiv_mmm_hcurl_eddy_bubble"]["support_status"]
        == "required_validation_path"
    )
    text = format_motor_triple_check_plan(plan)
    assert "HDiv-MMM + HCurl Eddy-Bubble" in text
    assert "primary required lanes" in text
    assert "optional auxiliary lanes: `none`" in text


def test_source_deck_review_packet_gate_requires_both_motor_lanes():
    packet = {
        "schema_version": "motor-source-deck-review-packet/v1",
        "observable_id": "torque",
        "selected_decks": [
            {"family": "motor/ipm", "case": "case001", "mai_path": "motor/ipm/case001.mai"}
        ],
        "required_lanes": ["ngsolve_age", "hdiv_mmm_hcurl_eddy_bubble"],
        "required_result_fields": [
            "observable_id", "observable_unit", "coordinate_frame", "sign_convention",
            "solver_version", "run_date_utc", "timing_breakdown_s",
            "shared_mesh_material_identity", "solver_ready_artifact",
        ],
        "publication_boundary": "metadata only",
    }
    ok = validate_motor_source_deck_review_packet(packet)
    assert ok["status"] == "ok"

    age_only = {**packet, "required_lanes": ["ngsolve_age"]}
    bad = validate_motor_source_deck_review_packet(age_only)
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["both_primary_lanes_required"] is False


def test_triple_check_gate_accepts_both_solver_ready_lanes_with_shared_identity():
    artifact = {
        "schema_version": "radia-motor-triple-check-artifact/v1",
        "goal": "IPM hairpin motor flux linkage and MTPA",
        "source_mcp_seed": {
            "source_mcp_calls": [
                "elf_motor_hybrid_router('IPM hairpin motor flux linkage and MTPA')"
            ],
            "representative_public_decks": [
                "application/motor/emdlab_ipm_hairpin_10/eip001/eip001.mai"
            ],
        },
        "lane_artifacts": {
            "hdiv_mmm_hcurl_eddy_bubble": _lane_artifact(
                "hdiv_mmm_hcurl_eddy_bubble", "pickup_flux"
            ),
            "ngsolve_age": _lane_artifact("ngsolve_age", "torque"),
        },
        "mcp_feedback": {
            "public_status": "verified",
            "public_summary": (
                "AGE and the HDiv-MMM/HCurl mixed system share one model identity."
            ),
            "learning_targets": ["radia_mcp.motor.triple_check_knowledge"],
            "verification": ["pytest tests/test_motor_triple_check.py"],
        },
    }
    result = validate_motor_triple_check_artifact(json.dumps(artifact))
    assert result["status"] == "pass"
    assert result["accepted_for_supported_mcp_learning"] is True
    assert result["accepted_for_mcp_rfc_learning"] is False
    assert result["accepted_for_mcp_learning"] is True
    assert result["research_triple_check_ready"] is True
    assert result["validated_supported_solver_check"] is True
    assert result["validated_dual_solver_check"] is True
    assert result["shared_model_identity_matches"] is True
    assert result["accepted_for_primary_dual_learning"] is True
    text = format_triple_check_gate_result(result)
    assert "validated supported solver check: `True`" in text
    assert "validated dual solver check: `True`" in text
    assert "accepted for MCP learning: `True`" in text


def test_triple_check_gate_rejects_mismatched_model_identity():
    hdiv = _lane_artifact(
        "hdiv_mmm_hcurl_eddy_bubble", "force_or_torque_trend"
    )
    age = _lane_artifact("ngsolve_age", "torque")
    age["shared_mesh_material_identity"]["geometry_sha256"] = "9" * 64
    artifact = {
        "schema_version": "radia-motor-triple-check-artifact/v1",
        "goal": "SPM shared-identity rejection",
        "source_mcp_seed": {
            "source_mcp_calls": [
                "elf_motor_hybrid_router('SPM shared-identity rejection')"
            ],
            "representative_public_decks": [
                "application/motor/spm_surface_pm_10/spm001/spm001.mai"
            ],
        },
        "lane_artifacts": {
            "hdiv_mmm_hcurl_eddy_bubble": hdiv,
            "ngsolve_age": age,
        },
        "mcp_feedback": {
            "public_status": "verified",
            "public_summary": (
                "Both lanes are solver-ready but intentionally use different geometry."
            ),
            "learning_targets": ["radia_mcp.motor.triple_check_knowledge"],
            "verification": ["pytest tests/test_motor_triple_check.py"],
        },
    }
    result = validate_motor_triple_check_artifact(json.dumps(artifact))
    assert result["status"] == "fail"
    assert result["shared_model_identity_matches"] is False
    assert result["accepted_for_mcp_learning"] is False
    assert any("identities do not match" in error for error in result["errors"])


def test_triple_check_gate_rejects_non_hex_age_identity():
    age = _lane_artifact("ngsolve_age", "torque")
    age["shared_mesh_material_identity"]["geometry_sha256"] = "z" * 64
    artifact = {
        "schema_version": "radia-motor-triple-check-artifact/v1",
        "goal": "reject malformed identity",
        "source_mcp_seed": {
            "source_mcp_calls": ["elf_motor_hybrid_router('identity')"],
            "representative_public_decks": ["application/motor/case.mai"],
        },
        "lane_artifacts": {
            "ngsolve_age": age,
            "hdiv_mmm_hcurl_eddy_bubble": _lane_artifact(
                "hdiv_mmm_hcurl_eddy_bubble", "pickup_flux"
            ),
        },
        "mcp_feedback": {
            "public_status": "verified",
            "public_summary": "malformed digest must not train the MCP",
            "learning_targets": ["radia_mcp.motor.triple_check_knowledge"],
            "verification": ["pytest tests/test_motor_triple_check.py"],
        },
    }

    result = validate_motor_triple_check_artifact(artifact)

    assert result["status"] == "fail"
    assert result["shared_model_identity_matches"] is False
    assert result["accepted_for_mcp_learning"] is False


def test_primary_dual_learning_requires_only_age_and_hdiv_lanes():
    hdiv = _lane_artifact(
        "hdiv_mmm_hcurl_eddy_bubble", "force_or_torque_trend"
    )
    artifact = {
        "schema_version": "radia-motor-triple-check-artifact/v1",
        "goal": "SPM primary AGE and HDiv-MMM/HCurl comparison",
        "source_mcp_seed": {
            "source_mcp_calls": [
                "elf_motor_hybrid_router('SPM primary mixed-system comparison')"
            ],
            "representative_public_decks": [
                "application/motor/spm_surface_pm_10/spm001/spm001.mai"
            ],
        },
        "lane_artifacts": {
            "hdiv_mmm_hcurl_eddy_bubble": hdiv,
            "ngsolve_age": _lane_artifact("ngsolve_age", "torque"),
        },
        "mcp_feedback": {
            "public_status": "verified",
            "public_summary": "AGE and HDiv-MMM/HCurl are the mandatory pair.",
            "learning_targets": ["radia_mcp.motor.triple_check_knowledge"],
            "verification": ["pytest tests/test_motor_triple_check.py"],
        },
    }

    result = validate_motor_triple_check_artifact(json.dumps(artifact))

    assert result["status"] == "pass"
    assert result["validated_dual_solver_check"] is True
    assert result["accepted_for_mcp_learning"] is True


def test_primary_dual_learning_rejects_age_only_artifact():
    artifact = {
        "schema_version": "radia-motor-triple-check-artifact/v1",
        "goal": "AGE only must not train radia-motor",
        "source_mcp_seed": {
            "source_mcp_calls": ["elf_motor_hybrid_router('AGE only')"],
            "representative_public_decks": [
                "application/motor/spm_surface_pm_10/spm001/spm001.mai"
            ],
        },
        "lane_artifacts": {
            "ngsolve_age": _lane_artifact("ngsolve_age", "torque"),
        },
        "mcp_feedback": {
            "public_status": "verified",
            "public_summary": "single-lane AGE artifacts are not enough",
            "learning_targets": ["radia_mcp.motor.triple_check_knowledge"],
            "verification": ["pytest tests/test_motor_triple_check.py"],
        },
    }

    result = validate_motor_triple_check_artifact(json.dumps(artifact))

    assert result["status"] == "fail"
    assert result["accepted_for_mcp_learning"] is False
    assert "missing lane artifact: hdiv_mmm_hcurl_eddy_bubble" in result["errors"]
