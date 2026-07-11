import json
import pytest

from radia_mcp.radia_ngsolve.force_position_profile_gate import force_position_profile_gate
from radia_mcp.radia_ngsolve.server import force_position_profile_gate as mcp_force_position_profile_gate


POSITIONS = [0.1 * index for index in range(21)]
FORCES = [
    27.2971, 24.5771, 26.7999, 30.6236, 35.4418, 42.1264, 47.9895,
    42.1516, 37.8650, 37.5447, 38.5627, 40.0615, 40.6014, 36.4427,
    29.8449, 24.3283, 20.3980, 17.5330, 15.4716, 13.8515, 12.5382,
]
NODES = [26225, 26112, 26215, 26230, 26193, 26154, 26237, 26112, 26067, 26152,
         26081, 26134, 26097, 26121, 26232, 26025, 26101, 26099, 26059, 26133, 26054]
ELEMENTS = [51854, 51629, 51840, 51881, 51803, 51727, 51890, 51643, 51550, 51717,
            51578, 51684, 51607, 51654, 51872, 51465, 51624, 51618, 51536, 51679, 51522]


def test_force_position_profile_gate_accepts_nonmonotone_interior_peak_and_mesh_drift():
    result = force_position_profile_gate(
        POSITIONS,
        FORCES,
        node_counts=NODES,
        element_counts=ELEMENTS,
        require_interior_peak=True,
        require_nonnegative=True,
    )
    assert result["status"] == "ok"
    assert result["peak_index"] == 6
    assert result["peak_position"] == pytest.approx(0.6)
    assert result["node_count_relative_span"] < 0.02
    assert result["element_count_relative_span"] < 0.02


def test_force_position_profile_gate_rejects_flat_or_large_mesh_drift():
    result = force_position_profile_gate(
        POSITIONS,
        [1.0] * len(POSITIONS),
        node_counts=[100] * 20 + [150],
        element_counts=[200] * len(POSITIONS),
        require_interior_peak=True,
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["force_profile_nontrivial"] is False
    assert result["checks"]["node_count_drift_ok"] is False
    assert result["checks"]["interior_peak_present_when_required"] is False


def test_force_position_profile_mcp_tool_dispatches_json():
    result = json.loads(mcp_force_position_profile_gate(
        POSITIONS,
        FORCES,
        NODES,
        ELEMENTS,
        require_interior_peak=True,
        require_nonnegative=True,
    ))
    assert result["status"] == "ok"
    assert result["interior_peak"] is True
