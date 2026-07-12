import copy
import json

from radia_mcp.maglev.periodic_settling_gate import rotating_conductor_periodic_settling_gate
from radia_mcp.maglev.server import rotating_conductor_periodic_settling_gate as mcp_gate


def _response():
    final = [10.0, 11.0, 10.5, 9.5]
    return [
        value
        for scale in (0.70, 0.93, 0.985, 0.997, 1.0)
        for value in (scale * item for item in final)
    ]


def test_periodic_settling_accepts_contracted_replay_and_dispatches():
    response = _response()
    result = rotating_conductor_periodic_settling_gate(
        response,
        steps_per_period=4,
        angle_step_deg=90.0,
        reference_response=copy.deepcopy(response),
        maximum_final_period_relative_l2=5.0e-3,
    )
    assert result["status"] == "ok"
    assert result["metrics"]["period_count"] == 5
    assert json.loads(
        mcp_gate(response, 4, 90.0, response, 5.0e-3)
    )["status"] == "ok"


def test_periodic_settling_rejects_unsettled_last_turn():
    response = _response()
    response[-4:] = [8.0, 12.0, 8.0, 12.0]
    result = rotating_conductor_periodic_settling_gate(
        response,
        steps_per_period=4,
        angle_step_deg=90.0,
        maximum_final_period_relative_l2=5.0e-3,
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["final_period_is_settled"] is False
