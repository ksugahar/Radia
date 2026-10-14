from __future__ import annotations

import copy
import json

from radia_mcp.radia_ngsolve.linear_axisymmetric_circuit_energy_gate import (
    linear_axisymmetric_circuit_energy_gate as gate,
)
from radia_mcp.radia_ngsolve.server import linear_axisymmetric_circuit_energy_gate


def _summary() -> dict:
    inductance = 0.02292739419524418
    field_per_current = 0.012962730356142147
    energy_per_current_squared = 0.01146081375626481
    return {
        "model_contract": {
            "physics": "magnetostatic",
            "coordinate_system": "axisymmetric",
            "all_materials_linear": True,
            "same_mesh": True,
            "same_boundary_conditions": True,
            "only_circuit_current_changed": True,
        },
        "units": {
            "current": "A",
            "flux_linkage": "Wb-turn",
            "magnetic_flux_density": "T",
            "energy": "J",
            "coenergy": "J",
        },
        "rows": [
            {
                "current_a": current,
                "circuit_current_a": current,
                "flux_linkage_wb_turn": inductance * current,
                "magnetic_flux_density_t": field_per_current * current,
                "energy_j": energy_per_current_squared * current**2,
                "coenergy_j": energy_per_current_squared * current**2,
                "node_count": 9121,
                "element_count": 17760,
            }
            for current in (0.5, 1.0, 2.0)
        ],
    }


def test_gate_accepts_linear_axisymmetric_circuit_identities() -> None:
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["flux_per_current_relative_span"] == 0.0
    assert result["metrics"]["energy_identity_relative_error_max"] < 5.0e-4


def test_gate_rejects_field_drift_and_mesh_change() -> None:
    summary = copy.deepcopy(_summary())
    summary["rows"][2]["magnetic_flux_density_t"] *= 0.8
    summary["rows"][2]["element_count"] += 1
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "field_scales_with_current" in result["issues"]
    assert "fixed_positive_mesh_inventory" in result["issues"]


def test_gate_rejects_nonlinear_material_claim() -> None:
    summary = copy.deepcopy(_summary())
    summary["model_contract"]["all_materials_linear"] = False
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "linear_axisymmetric_winding_contract" in result["issues"]


def test_mcp_wrapper_returns_structured_invalid_input() -> None:
    result = json.loads(linear_axisymmetric_circuit_energy_gate('{"rows": []}'))
    assert result["status"] == "invalid_input"
