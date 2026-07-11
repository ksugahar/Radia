import json

from radia_mcp.radia_ngsolve.radar_range_angle_gate import radar_range_angle_localization_gate
from radia_mcp.radia_ngsolve.server import radar_range_angle_localization_gate as mcp_gate


def _frequency():
    return [76.0e9 + i * 5.0e6 for i in range(401)]


def _targets():
    return [
        {"target_id": "a", "expected_range_m": 5.10, "estimated_range_m": 5.145, "expected_angle_deg": -11.38, "estimated_angle_deg": -12.41},
        {"target_id": "b", "expected_range_m": 4.24, "estimated_range_m": 4.289, "expected_angle_deg": 44.94, "estimated_angle_deg": 45.75},
    ]


def test_range_angle_gate_accepts_two_localized_targets():
    result = radar_range_angle_localization_gate(_frequency(), _targets())
    assert result["status"] == "ok"
    assert json.loads(mcp_gate(_frequency(), json.dumps(_targets())))["status"] == "ok"


def test_range_angle_gate_rejects_range_and_angle_drift():
    rows = _targets()
    rows[1]["estimated_range_m"] = 4.7
    rows[1]["estimated_angle_deg"] = 50.0
    result = radar_range_angle_localization_gate(_frequency(), rows)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_targets_localized"] is False


def test_range_angle_gate_rejects_duplicate_target_assignment():
    rows = _targets()
    rows[1]["target_id"] = "a"
    result = radar_range_angle_localization_gate(_frequency(), rows)
    assert result["checks"]["target_ids_unique"] is False
