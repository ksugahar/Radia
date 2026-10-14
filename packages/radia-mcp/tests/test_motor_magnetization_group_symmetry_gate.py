import copy
import json
import math

from radia_mcp.motor.magnetization_group_symmetry_gate import (
    mirror_symmetric_three_magnet_handoff_gate,
)
from radia_mcp.motor.server import motor_mirror_symmetric_three_magnet_handoff_gate


def _run() -> dict:
    vectors = []
    for index in range(4):
        angle = 0.1 + index * 0.08
        left = [math.sin(angle), math.cos(angle), 0.0]
        center = [0.0, math.cos(angle), math.sin(angle)]
        right = [-left[0], left[1], left[2]]
        for group, direction, offset in ((1, left, 0), (2, center, 100), (3, right, 200)):
            vectors.append(
                {
                    "element_id": offset + index + 1,
                    "group": group,
                    "direction": direction,
                    "magnitude": 0.4 - 0.02 * index,
                }
            )
    return {
        "solver_returncode": 0,
        "end_of_moment_solution": True,
        "error_count": 0,
        "magnetization_vectors": vectors,
        "bmax_t_by_material": {"1": 0.18, "2": 0.16, "3": 0.18},
    }


def summary() -> dict:
    first = _run()
    return {
        "runs": [first, copy.deepcopy(first)],
        "result_contract": {
            "group_ids": [1, 2, 3],
            "vectors_per_group": 4,
            "mirror_axis": "x",
        },
    }


def test_accepts_unit_mirror_vectors_and_repeated_outer_response():
    result = mirror_symmetric_three_magnet_handoff_gate(summary())
    assert result["status"] == "ok"
    assert result["handoff_ready"] is True
    assert json.loads(
        motor_mirror_symmetric_three_magnet_handoff_gate(json.dumps(summary()))
    )["status"] == "ok"


def test_rejects_nonunit_direction_and_outer_response_drift():
    bad = summary()
    bad["runs"][1]["magnetization_vectors"][0]["direction"][0] *= 2.0
    bad["runs"][1]["bmax_t_by_material"]["3"] = 0.14
    result = mirror_symmetric_three_magnet_handoff_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_magnetization_directions_are_unit_vectors"] is False
    assert result["checks"]["outer_group_responses_are_symmetric"] is False
    assert result["checks"]["fresh_runs_repeat_vector_package"] is False


def test_rejects_missing_group_vector_and_incomplete_solver_run():
    bad = summary()
    bad["runs"][0]["magnetization_vectors"].pop()
    bad["runs"][0]["end_of_moment_solution"] = False
    result = mirror_symmetric_three_magnet_handoff_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["three_groups_have_expected_vector_count"] is False
    assert result["checks"]["two_or_more_complete_solver_runs"] is False
