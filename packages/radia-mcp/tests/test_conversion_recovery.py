from radia_mcp.radia_ngsolve.conversion_recovery import (
    classify_converted_geometry,
    converted_model_recovery_summary,
    learning_lane_closure,
)


def test_classify_converted_geometry_keeps_stub_honest():
    stub = classify_converted_geometry("mph_stub")
    assert stub["is_stub"] is True
    assert stub["requires_geometry_rebuild"] is True
    assert stub["solver_ready"] is False

    vol = classify_converted_geometry("tri_tet_vol")
    assert vol["is_tri_tet_vol"] is True
    assert vol["requires_geometry_rebuild"] is False
    assert vol["solver_ready"] is True


def test_learning_lane_closure_requires_source_when_requested():
    closed = learning_lane_closure({"public": "verified", "source_tool": "verified"}, require_source=True)
    assert closed["closed"] is True

    candidate = learning_lane_closure({"public": "verified", "source_tool": "candidate"}, require_source=True)
    assert candidate["closed"] is False
    assert candidate["source_closed"] is False

    skipped_source = learning_lane_closure({"public": "verified", "source_tool": "candidate"}, require_source=False)
    assert skipped_source["closed"] is True


def test_converted_model_recovery_summary_counts_stubs_and_open_lanes():
    summary = converted_model_recovery_summary(
        [
            {
                "tool_slot": "COMSOL",
                "pass": True,
                "geometry_kind": "mph_stub",
                "learning_lanes": {"public": "verified", "source_tool": "verified"},
            },
            {
                "tool_slot": "Coreform(Cubit)",
                "pass": True,
                "geometry_kind": "tri_tet_vol",
                "learning_lanes": {"public": "verified", "source_tool": "none"},
            },
            {
                "tool_slot": "FEMM",
                "pass": True,
                "geometry_kind": "analytic_2d",
                "learning_lanes": {"public": "verified", "source_tool": "candidate"},
            },
        ],
        require_source_slots={"FEMM"},
    )

    assert summary["status"] == "needs_attention"
    assert summary["all_passed"] is True
    assert summary["stub_geometry_tools"] == ["COMSOL"]
    assert summary["open_learning_lane_tools"] == ["FEMM"]
    assert summary["ready_for_reuse_count"] == 1
