import copy
import json

from radia_mcp.radia_ngsolve.server import voice_coil_force_flux_sweep_gate as mcp_gate
from radia_mcp.radia_ngsolve.voice_coil_gate import voice_coil_force_flux_sweep_gate


ROWS = [
    {"current_a": -2.0, "circuit_current_a": -2.0, "axial_force_n": -9.1289652746, "flux_linkage_wb_turn": [-0.0176567683, 0.0], "node_count": 10321, "element_count": 20280},
    {"current_a": -1.0, "circuit_current_a": -1.0, "axial_force_n": -4.4613303962, "flux_linkage_wb_turn": [-0.0168761056, 0.0], "node_count": 10321, "element_count": 20280},
    {"current_a": 0.0, "circuit_current_a": 0.0, "axial_force_n": 0.0221434354, "flux_linkage_wb_turn": [-0.0160798515, 0.0], "node_count": 10321, "element_count": 20280},
    {"current_a": 1.0, "circuit_current_a": 1.0, "axial_force_n": 4.3093265809, "flux_linkage_wb_turn": [-0.0152673166, 0.0], "node_count": 10321, "element_count": 20280},
    {"current_a": 2.0, "circuit_current_a": 2.0, "axial_force_n": 8.3878807256, "flux_linkage_wb_turn": [-0.0144384824, 0.0], "node_count": 10321, "element_count": 20280},
]


def test_live_shape_voice_coil_sweep_passes_and_dispatches():
    result = voice_coil_force_flux_sweep_gate(ROWS)
    assert result["status"] == "ok"
    assert result["metrics"]["maximum_odd_residual_relative"] < 0.10
    assert json.loads(mcp_gate(ROWS))["status"] == "ok"


def test_voice_coil_gate_rejects_stale_mesh_and_negative_incremental_flux():
    bad = copy.deepcopy(ROWS)
    bad[-1]["node_count"] += 1
    bad[-1]["flux_linkage_wb_turn"][0] = -0.020
    result = voice_coil_force_flux_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["mesh_inventory_invariant"] is False
    assert result["checks"]["incremental_flux_linkage_positive"] is False


def test_voice_coil_gate_rejects_idealized_flat_force():
    bad = copy.deepcopy(ROWS)
    for row in bad:
        row["axial_force_n"] = 0.0
    result = voice_coil_force_flux_sweep_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["axial_force_strictly_increases"] is False
