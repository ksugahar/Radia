import copy
import json

import pytest

from radia_mcp.radia_ngsolve.server import thermal_robin_boundary_balance_gate
from radia_mcp.radia_ngsolve.thermal_robin_balance_gate import (
    thermal_robin_boundary_balance_gate as evaluate_gate,
)


def _summary():
    return {
        "boundary_groups": [
            {"name": "heated", "heat_rate_W": -100.0},
            {"name": "cooled", "heat_rate_W": 100.005},
        ],
        "internal_cut_heat_rate_W": 100.001,
        "mesh_ladder": [
            {
                "role": "plateau_a",
                "node_count": 100,
                "element_count": 180,
                "average_temperature_K": 500.0,
                "robin_throughput_W": 100.0025,
                "balance_relative": 2.0e-4,
            },
            {
                "role": "plateau_b",
                "node_count": 100,
                "element_count": 180,
                "average_temperature_K": 500.0,
                "robin_throughput_W": 100.0025,
                "balance_relative": 2.0e-4,
            },
            {
                "role": "fine",
                "node_count": 160,
                "element_count": 300,
                "average_temperature_K": 500.01,
                "robin_throughput_W": 100.003,
                "balance_relative": 5.0e-5,
            },
            {
                "role": "fine_repeat",
                "node_count": 160,
                "element_count": 300,
                "average_temperature_K": 500.01,
                "robin_throughput_W": 100.003,
                "balance_relative": 5.0e-5,
            },
        ],
        "exact_boundary_flux_status": "rejected_nonfinite_boundary_flux",
        "symmetry_relative_error": 8.0e-6,
        "constitutive_relative_error": 3.0e-16,
        "temperature_reflection": {
            "temperature_relative_errors": {"field": 4.0e-9},
            "flux_relative_errors": {"cut": 1.6e-8},
        },
    }


def test_accepts_balance_plateau_refinement_replay_and_reflection():
    result = evaluate_gate(_summary())
    assert result["status"] == "ok"
    assert result["checks"]["adaptive_mesh_plateau_is_explicit"] is True
    assert result["metrics"]["balance_relative"] < 1.0e-4


def test_rejects_a_same_sign_boundary_flux_control():
    bad = copy.deepcopy(_summary())
    bad["boundary_groups"][0]["heat_rate_W"] = 100.0
    result = evaluate_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["opposed_boundary_heat_rate_signs_exist"] is False


def test_rejects_a_false_exact_boundary_flux_success():
    bad = _summary()
    bad["exact_boundary_flux_status"] = "accepted"
    result = evaluate_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["exact_boundary_nonfinite_flux_is_rejected"] is False


def test_stdio_wrapper_returns_json_and_invalid_input():
    good = json.loads(thermal_robin_boundary_balance_gate(json.dumps(_summary())))
    assert good["status"] == "ok"
    invalid = json.loads(thermal_robin_boundary_balance_gate("[]"))
    assert invalid["status"] == "invalid_input"
    assert "summary must be an object" in invalid["error"]


def test_requires_all_four_mesh_roles():
    bad = _summary()
    bad["mesh_ladder"].pop()
    with pytest.raises(ValueError, match="required role"):
        evaluate_gate(bad)
