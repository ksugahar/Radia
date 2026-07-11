import json

import pytest

from radia_mcp.radia_ngsolve.rotational_time_axis_gate import rotational_kinematics_time_axis_gate
from radia_mcp.radia_ngsolve.server import rotational_kinematics_time_axis_gate as mcp_gate


def _case():
    time = [0.1 * index for index in range(10)]
    angle = [0.6 * index for index in range(10)]
    speed = [0.0] + [1.0] * 9
    return time, angle, speed


def test_rotational_time_axis_accepts_si_values_with_display_unit_metadata():
    time, angle, speed = _case()
    result = rotational_kinematics_time_axis_gate(
        time, angle, speed, reported_time_unit="ms", time_value_basis="si_seconds"
    )
    assert result["status"] == "ok"
    assert result["max_central_relative_error"] < 1.0e-12
    assert result["display_unit_label_is_metadata_only"] is True
    assert result["display_unit_interpretation_relative_error"] == pytest.approx(0.999)


def test_rotational_time_axis_rejects_scaling_si_values_by_display_label():
    time, angle, speed = _case()
    result = rotational_kinematics_time_axis_gate(
        time, angle, speed, reported_time_unit="ms", time_value_basis="display_unit"
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["rotational_kinematics_match"] is False


def test_rotational_time_axis_mcp_dispatches_and_reports_invalid_shape():
    time, angle, speed = _case()
    result = json.loads(mcp_gate(time, angle, speed, "ms", "si_seconds"))
    assert result["status"] == "ok"
    bad = json.loads(mcp_gate(time, angle[:-1], speed, "ms", "si_seconds"))
    assert bad["status"] == "invalid_input"
