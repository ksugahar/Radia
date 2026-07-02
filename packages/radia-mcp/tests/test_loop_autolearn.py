from radia_mcp.radia_ngsolve.loop_autolearn import (
    build_autonomous_basic_learning_artifact,
    classify_slot_family,
)


def test_autonomous_basic_learning_artifact_processes_all_slots_without_overclaiming():
    queue = {
        "created_at": "2026-07-03T00:00:00Z",
        "rounds": 1,
        "total_slots": 4,
        "rotation": ["COMSOL", "Coreform(Cubit)", "FEMM", "MATLAB"],
        "slots": [
            {
                "tool": "COMSOL",
                "source_native_example": "public-example:acdc/coaxial-cable",
                "source_type": "public_doc",
                "lesson_axis": "model-page replay and LiveLink attach-only discipline",
                "intended_validation": "extract metadata and map physics to public analogues",
                "status": "queued_source_native_preflight",
                "lap": 1,
                "slot_id": 1,
            },
            {
                "tool": "Coreform(Cubit)",
                "source_native_example": "upstream-example:cubit/sphere",
                "source_type": "upstream_example",
                "lesson_axis": "sphere geometry and mesh reference",
                "intended_validation": "analytic sphere volume/area and export contract",
                "status": "queued_source_native_preflight",
                "lap": 1,
                "slot_id": 2,
            },
            {
                "tool": "FEMM",
                "source_native_example": "upstream-example:femm/coilgun",
                "source_type": "upstream_example",
                "lesson_axis": "force calculation examples",
                "intended_validation": "compare force identity with radia-ngsolve",
                "status": "queued_source_native_preflight",
                "lap": 1,
                "slot_id": 3,
            },
            {
                "tool": "MATLAB",
                "source_native_example": "public-example:mathworks/least-squares",
                "source_type": "public_doc",
                "lesson_axis": "optimization teaching lane",
                "intended_validation": "objective/constraint/sensitivity examples with reference gates",
                "status": "queued_source_native_preflight",
                "lap": 1,
                "slot_id": 4,
            },
        ],
    }

    artifact = build_autonomous_basic_learning_artifact(
        queue,
        artifact_id="autonomous_basic_learning_test",
        queue_id="fixture_queue",
        run_date_utc="2026-07-03T00:01:00Z",
        radia_mcp_version="test",
        command="python validation/loop_learning/autonomous_basic_learning.py --queue-json fixture.json --out-dir out",
        check_local_sources=False,
    )

    assert artifact["pass"] is True
    assert artifact["summary"]["slot_count"] == 4
    assert artifact["summary"]["row_count"] == 4
    assert artifact["queue_gate"]["status"] == "ok"
    assert artifact["row_gate"]["status"] == "ok"
    assert artifact["mcp_feedback"]["feedback_gate"]["status"] == "ok"
    assert artifact["learning_lanes"]["public"] == "verified"
    assert artifact["learning_lanes"]["source_tool"] == "candidate"
    assert artifact["summary"]["source_tool_candidate_count"] == 3
    assert artifact["checks"]["source_tool_learning_not_overclaimed"] is True
    assert all(row["pass"] for row in artifact["rows"])


def test_autonomous_basic_learning_family_classification():
    assert classify_slot_family({"tool": "MATLAB", "lesson_axis": "optimization teaching lane"}) == "matlab_optimization"
    assert classify_slot_family({"tool": "CST", "lesson_axis": "Touchstone S-parameter export"}) == "rf_acoustic"
    assert classify_slot_family({"tool": "build123d", "lesson_axis": "CAD volume mesh route"}) == "geometry_mesh"
    assert classify_slot_family({"tool": "JMAG", "lesson_axis": "force calculation examples"}) == "force_torque_motor"
    assert classify_slot_family({"tool": "MATLAB", "lesson_axis": "Gypsilab FEM/BEM coupling"}) == "fem_bem"
