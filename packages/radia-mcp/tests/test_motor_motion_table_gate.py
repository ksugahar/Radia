import json

from radia_mcp.motor.motion_table_gate import motion_table_coordinate_gate
from radia_mcp.motor.server import motor_motion_table_coordinate_gate


def test_motion_table_accepts_independent_translation_and_rotation_axes():
    result = motion_table_coordinate_gate(
        [0.0, 3.0], [[0, 0, 0], [0.3, 0.3, 0.3]],
        [0.0, 4.0], [[0, 0, 0], [0, 90, 0]],
        coordinate_frame_id="body_frame",
    )
    assert result["status"] == "ok"
    assert result["independent_time_axes"] is True
    assert result["combined_motion_end_time_s"] == 4.0


def test_motion_table_rejects_nonmonotone_time_and_ambiguous_units():
    result = motion_table_coordinate_gate(
        [0.0, 3.0, 2.0], [[0, 0, 0], [1, 0, 0], [2, 0, 0]],
        [0.0, 4.0], [[0, 0, 0], [0, 90, 0]],
        coordinate_frame_id="",
        rotation_unit="rpm",
        motion_semantics="incremental",
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["translation_time_strictly_increases"] is False
    assert result["checks"]["coordinate_frame_recorded"] is False


def test_motion_table_mcp_dispatches_json():
    result = json.loads(motor_motion_table_coordinate_gate(
        [0.0, 3.0], [[0, 0, 0], [0.3, 0.3, 0.3]],
        [0.0, 4.0], [[0, 0, 0], [0, 90, 0]],
        "body_frame",
    ))
    assert result["status"] == "ok"
    assert result["policy"] == "motion_table_coordinate_gate_v1"
