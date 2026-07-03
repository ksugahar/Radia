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
    }
    if lane == "hdiv_vim_reduced_fem":
        artifact["metrics"] = {"max_abs_relative_error": 1.0e-3}
        artifact["coupling_design_status"] = "experimental_rfc"
        artifact["interface_operator_contract"] = {
            "rotor_side": "HDiv-VIM source field",
            "stator_side": "fixed-stator reduced FEM",
            "status": "experimental RFC",
        }
        artifact["reduced_fem_contract"] = {"basis": "P1", "observable": observable}
        artifact["vim_operator_contract"] = {"space": "HDiv", "observable": observable}
    elif lane == "mmmm2d_coarse":
        artifact["metrics"] = {"torque_relative_error": 1.0e-3}
        artifact["mmmm2d_contract"] = {
            "solver": "radia.mmmm2d",
            "material_input": "per-region dict",
            "sweep": "factor-once linear torque_angle_sweep",
        }
        artifact["region_material_contract"] = {
            "regions": ["inner", "outer"],
            "missing_region_policy": "raise",
        }
        artifact["pytest_targets"] = ["validation_test/feec/test_moment2d_perregion.py"]
    else:
        artifact["metrics"] = {"torque_relative_error": 1.0e-3}
        artifact["age_gate_ids"] = ["age_rotation_torque"]
        artifact["pytest_targets"] = ["tests/test_airgap_machine_rotation.py"]
    return artifact


def test_triple_check_plan_marks_hdiv_reduced_fem_as_experimental():
    plan = route_motor_triple_check("IPM hairpin motor flux linkage and MTPA")
    assert plan["inferred_family"] == "ipm"
    assert "elf_motor_hybrid_router" in plan["source_mcp_seed"]["calls"][0]
    assert "application/motor/emdlab_ipm_hairpin_10/eip001/eip001.mai" in (
        plan["source_mcp_seed"]["representative_public_decks"]
    )
    assert (
        plan["radia_lanes"]["ngsolve_age"]["support_status"]
        == "supported_validation_path"
    )
    assert (
        plan["radia_lanes"]["mmmm2d_coarse"]["support_status"]
        == "supported_coarse_path"
    )
    assert (
        plan["radia_lanes"]["hdiv_vim_reduced_fem"]["support_status"]
        == "experimental_rfc"
    )
    text = format_motor_triple_check_plan(plan)
    assert "experimental_rfc" in text
    assert "mmmm2d_coarse" in text
    assert "supported_coarse_path" in text
    assert "validated solver path" in text


def test_triple_check_gate_accepts_learning_with_mmmm_but_not_hdiv_solver_validation():
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
            "hdiv_vim_reduced_fem": _lane_artifact(
                "hdiv_vim_reduced_fem", "pickup_flux"
            ),
            "mmmm2d_coarse": _lane_artifact("mmmm2d_coarse", "torque"),
            "ngsolve_age": _lane_artifact("ngsolve_age", "torque"),
        },
        "mcp_feedback": {
            "public_status": "verified",
            "public_summary": (
                "AGE is the supported full path, MMMM is the supported coarse "
                "path, and HDiv-VIM plus reduced FEM is an RFC."
            ),
            "learning_targets": ["radia_mcp.motor.triple_check_knowledge"],
            "verification": ["pytest tests/test_motor_triple_check.py"],
        },
    }
    result = validate_motor_triple_check_artifact(json.dumps(artifact))
    assert result["status"] == "pass"
    assert result["accepted_for_mcp_learning"] is True
    assert result["research_triple_check_ready"] is True
    assert result["validated_supported_solver_check"] is True
    assert result["validated_dual_solver_check"] is False
    text = format_triple_check_gate_result(result)
    assert "validated supported solver check: `True`" in text
    assert "validated dual solver check: `False`" in text
