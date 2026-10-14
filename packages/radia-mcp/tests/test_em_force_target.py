from radia_mcp.radia_ngsolve.em_force_target import (
    build_em_force_target_artifact,
    em_force_public_row_for_slot,
    select_em_force_slots,
)
from radia_mcp.radia_ngsolve.force import electromagnetic_force_method_selection_gate


def _slot(tool, slot_id, lesson_axis="force calculation examples"):
    return {
        "index": slot_id,
        "slot_id": slot_id,
        "lap": 1,
        "slot_index_in_lap": slot_id,
        "tool": tool,
        "lesson_family": "force_torque_motor",
        "source_type": "upstream_example",
        "source_present": True,
        "required_fields_present": True,
        "lesson_axis": lesson_axis,
        "intended_validation": "compare EM force identity with radia-ngsolve",
        "learning_lanes": {"public": "verified", "source_tool": "candidate"},
    }


def test_em_force_target_extracts_force_slots_without_claiming_live_solvers():
    source_artifact = {
        "artifact_role": "autonomous_basic_learning",
        "result_artifact_id": "autonomous_basic_learning_fixture",
        "slots": [
            _slot("COMSOL", 1, "air-gap actuator force"),
            _slot("FEMM", 2, "parallel wire force"),
            _slot("JMAG", 3, "IPM motor torque"),
            _slot("ELF(MAGIC product)", 4, "PM gap force sweep"),
            {
                "slot_id": 5,
                "tool": "build123d",
                "lesson_family": "geometry_mesh",
                "learning_lanes": {"public": "verified", "source_tool": "none"},
            },
        ],
    }

    artifact = build_em_force_target_artifact(
        source_artifact,
        artifact_id="em_force_target_test",
        source_artifact_id="autonomous_basic_learning_fixture",
        run_date_utc="2026-07-03T00:10:00Z",
        radia_mcp_version="test",
        command="python validation/force/electromagnetic_force_target.py --source-json fixture.json --out-dir out",
    )

    assert artifact["pass"] is True
    assert artifact["summary"]["force_slot_count"] == 4
    assert artifact["summary"]["row_count"] == 4
    assert artifact["summary"]["source_tool_candidate_count"] == 4
    assert artifact["summary"]["target_counts"] == {
        "ipm_dq_torque_components": 1,
        "magnetic_air_gap_pressure": 1,
        "parallel_wire_lorentz": 1,
        "pm_force_gap_sweep": 1,
    }
    assert artifact["row_gate"]["status"] == "ok"
    assert artifact["mcp_feedback"]["feedback_gate"]["status"] == "ok"
    assert artifact["learning_lanes"] == {"public": "verified", "source_tool": "candidate"}
    assert artifact["checks"]["source_tool_learning_not_overclaimed"] is True
    assert all(row["pass"] for row in artifact["rows"])


def test_select_em_force_slots_can_classify_queue_like_records():
    source_artifact = {
        "slots": [
            {"tool": "JMAG", "lesson_axis": "motor torque table"},
            {"tool": "MATLAB", "lesson_axis": "optimization objective"},
        ],
    }

    slots = select_em_force_slots(source_artifact)

    assert len(slots) == 1
    assert slots[0]["tool"] == "JMAG"


def test_em_force_public_row_for_slot_uses_tool_specific_gate():
    assert em_force_public_row_for_slot(_slot("FEMM", 10))["target_kind"] == "parallel_wire_lorentz"
    assert em_force_public_row_for_slot(_slot("JMAG", 11))["target_kind"] == "ipm_dq_torque_components"
    assert em_force_public_row_for_slot(_slot("ELF(MAGIC product)", 12))["target_kind"] == "pm_force_gap_sweep"
    assert em_force_public_row_for_slot(_slot("COMSOL", 13))["target_kind"] == "magnetic_air_gap_pressure"


def test_electromagnetic_force_method_selection_prefers_robust_primary_evidence():
    conductor = electromagnetic_force_method_selection_gate(
        "current_conductor",
        "lorentz_body_force",
        relative_permeability=1.0,
    )
    assert conductor["status"] == "ok"
    assert conductor["recommended_method"] == "lorentz_body_force"

    magnetic = electromagnetic_force_method_selection_gate(
        "magnetic_body",
        "weighted_stress_volume",
        relative_permeability=1000.0,
        weighted_stress_available=True,
        virtual_work_samples_available=True,
    )
    assert magnetic["status"] == "ok"
    assert magnetic["recommended_method"] == "weighted_stress_volume"

    contour_only = electromagnetic_force_method_selection_gate(
        "magnetic_body",
        "contour_maxwell_stress",
        relative_permeability=1000.0,
        weighted_stress_available=True,
        contour_clearance_mesh_layers=2,
    )
    assert contour_only["status"] == "needs_attention"
    assert contour_only["checks"]["contour_not_used_as_primary"] is False
