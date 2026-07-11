import json

from radia_mcp.motor.server import motor_thermal_handoff_gate
from radia_mcp.motor.thermal_handoff import evaluate_motor_thermal_handoff


LOSSES = {"winding": 120.0, "stator_iron": 45.0, "rotor": 8.0}
NETWORK = {
    "ambient_node": "ambient",
    "nodes": [
        {"id": "winding", "source_regions": ["winding"], "capacitance_J_per_K": 620.0},
        {"id": "stator", "source_regions": ["stator_iron"], "capacitance_J_per_K": 1800.0},
        {"id": "rotor", "source_regions": ["rotor"], "capacitance_J_per_K": 950.0},
        {"id": "housing", "source_regions": [], "capacitance_J_per_K": 2400.0},
        {"id": "ambient", "source_regions": []},
    ],
    "branches": [
        {"from": "winding", "to": "stator", "resistance_K_per_W": 0.18},
        {"from": "rotor", "to": "stator", "resistance_K_per_W": 0.32},
        {"from": "stator", "to": "housing", "resistance_K_per_W": 0.11},
        {"from": "housing", "to": "ambient", "resistance_K_per_W": 0.24},
    ],
}
MESH = [
    {"region": "winding", "cell_type": "hex8", "cell_count": 3200, "loss_W": 120.0},
    {"region": "stator_iron", "cell_type": "hex20", "cell_count": 4800, "loss_W": 45.0},
    {"region": "rotor", "cell_type": "hex8", "cell_count": 2100, "loss_W": 8.0},
]


def test_motor_thermal_handoff_accepts_shared_loss_contract():
    result = evaluate_motor_thermal_handoff(LOSSES, NETWORK, MESH)
    assert result["status"] == "ok"
    assert result["source_total_loss_W"] == 173.0
    assert all(result["checks"].values())

    wrapped = json.loads(
        motor_thermal_handoff_gate(
            json.dumps(LOSSES), json.dumps(NETWORK), json.dumps(MESH)
        )
    )
    assert wrapped["status"] == "ok"
    assert wrapped["policy"] == "loss_to_lptn_and_hex_thermal_mesh"


def test_motor_thermal_handoff_rejects_non_hex_and_loss_mismatch():
    mesh = [dict(row) for row in MESH]
    mesh[1]["cell_type"] = "tet4"
    mesh[2]["loss_W"] = 9.0
    result = evaluate_motor_thermal_handoff(LOSSES, NETWORK, mesh)
    assert result["status"] == "needs_attention"
    assert result["checks"]["mesh_all_hexahedral"] is False
    assert result["checks"]["regional_loss_values_match"] is False
    assert result["checks"]["total_heat_conserved"] is False


def test_motor_thermal_handoff_rejects_missing_and_disconnected_network_regions():
    network = {
        **NETWORK,
        "nodes": [dict(row) for row in NETWORK["nodes"]],
        "branches": [dict(row) for row in NETWORK["branches"][:-1]],
    }
    network["nodes"][2]["source_regions"] = []
    result = evaluate_motor_thermal_handoff(LOSSES, network, MESH)
    assert result["status"] == "needs_attention"
    assert result["checks"]["network_regions_one_to_one"] is False
    assert result["checks"]["network_connected_to_ambient"] is False


def test_motor_thermal_handoff_rejects_invalid_json_at_mcp_boundary():
    try:
        motor_thermal_handoff_gate("not-json", json.dumps(NETWORK), json.dumps(MESH))
    except ValueError as exc:
        assert "loss_buckets_json must be valid JSON" in str(exc)
    else:
        raise AssertionError("invalid JSON was accepted")
