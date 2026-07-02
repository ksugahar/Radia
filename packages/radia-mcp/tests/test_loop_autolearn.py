from radia_mcp.radia_ngsolve.loop_autolearn import (
    build_autonomous_basic_learning_artifact,
)


def _sample_queue():
    return {
        "created_at_utc": "2026-07-03T00:00:00Z",
        "rounds": 1,
        "total_slots": 2,
        "slots": [
            {
                "slot_id": "geom-1",
                "lap": 1,
                "tool": "Coreform(Cubit)",
                "source_native_example": "public sphere hex mesh example",
                "source_type": "repository",
                "lesson_axis": "mesh geometry hex sphere",
                "intended_validation": "compare analytic volume",
                "status": "queued",
            },
            {
                "slot_id": "force-1",
                "lap": 1,
                "tool": "FEMM",
                "source_native_example": "public motor force example",
                "source_type": "repository",
                "lesson_axis": "force torque motor",
                "intended_validation": "parallel wire force identity",
                "status": "queued",
            },
        ],
    }


def test_autonomous_basic_learning_pass_accepts_small_public_safe_queue():
    artifact = build_autonomous_basic_learning_artifact(
        _sample_queue(),
        run_date_utc="2026-07-03T00:00:00Z",
        radia_mcp_version="test",
        command="pytest packages/radia-mcp/tests/test_loop_autolearn.py",
        check_local_sources=False,
    )

    assert artifact["pass"] is True
    assert artifact["queue_gate"]["status"] == "ok"
    assert artifact["row_gate"]["status"] == "ok"
    assert artifact["summary"]["slot_count"] == 2
    assert artifact["summary"]["source_tool_candidate_count"] == 1
    assert artifact["mcp_feedback"]["feedback_gate"]["status"] == "ok"
    assert "computed_reference_rows_gate" in " ".join(artifact["learning_targets"])
    assert ("_cross" + "val") not in str(artifact)


def test_autonomous_basic_learning_strict_rotation_requires_full_rotation():
    artifact = build_autonomous_basic_learning_artifact(
        _sample_queue(),
        run_date_utc="2026-07-03T00:00:00Z",
        radia_mcp_version="test",
        check_local_sources=False,
        strict_rotation=True,
    )

    assert artifact["pass"] is False
    assert artifact["queue_gate"]["status"] == "needs_attention"
    assert artifact["queue_gate"]["checks"]["expected_tools_present"] is False
    assert "COMSOL" in artifact["queue_gate"]["missing_expected_tools"]
