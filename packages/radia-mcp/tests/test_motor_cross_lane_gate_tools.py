"""MCP adapter contracts for the motor cross-lane gates.

The adapters live in radia_mcp.motor.server, whose import graph reaches
optional dependencies; the pure gates are covered by their own files, which
stay collectable in the minimal matrix.
"""

import json

from test_motor_back_emf_cross_lane_gate import LIMITS as _BACK_EMF_LIMITS
from test_motor_back_emf_cross_lane_gate import _lane as _back_emf_lane
from test_motor_cogging_torque_cross_lane_gate import _LIMITS as _COGGING_LIMITS
from test_motor_cogging_torque_cross_lane_gate import _lane as _cogging_lane


def test_cogging_tool_json_and_grid_alignment():
    """Guards: structured JSON results and grid alignment in the public entry point.
    Assumes: callers align the mechanical angle datum before invoking the gate.
    """
    from radia_mcp.motor.server import motor_cogging_torque_cross_lane_gate

    first = json.dumps(_cogging_lane())
    limits = json.dumps(_COGGING_LIMITS)
    assert json.loads(motor_cogging_torque_cross_lane_gate(first, first, limits))["status"] == "pass"
    assert json.loads(motor_cogging_torque_cross_lane_gate("[", first, limits))["checks"] == {"json_inputs": False}
    shifted = _cogging_lane()
    shifted["mechanical_angles_deg"] = [value + 0.2 for value in shifted["mechanical_angles_deg"]]
    result = json.loads(motor_cogging_torque_cross_lane_gate(first, json.dumps(shifted), limits))
    assert "common mechanical angle grid" in result["error"]


def test_back_emf_tool_uses_named_validation_group_and_json_strings():
    """Guards: the MCP adapter uses JSON strings and the validation registry.

    Assumes: the server has installed its validation tools at import time.
    """
    from radia_mcp.motor import server

    lane = json.dumps(_back_emf_lane())
    command = 'motor_back_emf_cross_lane_gate'
    handler = getattr(server, command)
    assert server._validation._entries[command].function is handler
    assert json.loads(handler(lane, lane, json.dumps(_BACK_EMF_LIMITS)))['status'] == 'pass'
    assert json.loads(handler('{', lane, json.dumps(_BACK_EMF_LIMITS)))['status'] == 'invalid_input'
