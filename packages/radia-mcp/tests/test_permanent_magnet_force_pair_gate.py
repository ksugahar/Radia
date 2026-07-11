import json

from radia_mcp.motor.permanent_magnet_force_pair_gate import permanent_magnet_force_pair_gate
from radia_mcp.motor.server import motor_permanent_magnet_force_pair_gate


def _summary():
    return {
        "interaction_axis": "x",
        "positive_axis_interaction": "repulsion",
        "force_unit": "N",
        "torque_unit": "N*m",
        "component_frame": "global_cartesian",
        "reference_length_m": 0.1,
        "cases": [
            {"pole_relation": "like", "force_N": [2.0, 1.0e-6, 0.0], "torque_Nm": [0.0, 0.0, 2.0e-7]},
            {"pole_relation": "opposite", "force_N": [-1.99, -1.0e-6, 0.0], "torque_Nm": [0.0, 0.0, -1.0e-7]},
        ],
    }


def test_force_pair_gate_accepts_reversal_and_nearly_equal_magnitude():
    result = permanent_magnet_force_pair_gate(json.dumps(_summary()))
    assert result["status"] == "ok"
    assert all(result["checks"].values())
    assert json.loads(motor_permanent_magnet_force_pair_gate(json.dumps(_summary())))["status"] == "ok"


def test_force_pair_gate_rejects_lost_reversal():
    summary = _summary()
    summary["cases"][1]["force_N"][0] = 1.99
    result = permanent_magnet_force_pair_gate(json.dumps(summary))
    assert result["status"] == "needs_attention"
    assert result["checks"]["opposite_poles_attract"] is False


def test_force_pair_gate_rejects_selection_or_frame_drift():
    summary = _summary()
    summary["cases"][1]["force_N"] = [-1.2, 0.2, 0.0]
    result = permanent_magnet_force_pair_gate(json.dumps(summary))
    assert result["status"] == "needs_attention"
    assert result["checks"]["axial_magnitudes_match"] is False
    assert result["checks"]["off_axis_force_small"] is False
