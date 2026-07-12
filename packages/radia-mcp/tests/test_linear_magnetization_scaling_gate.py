from __future__ import annotations

import copy
import json

from radia_mcp.radia_ngsolve.linear_magnetization_scaling_gate import (
    linear_magnetization_scaling_gate as gate,
)
from radia_mcp.radia_ngsolve.server import linear_magnetization_scaling_gate


def _summary() -> dict:
    return {
        "model_contract": {
            "physics": "magnetostatic",
            "source": "prescribed_magnetization",
            "all_materials_linear": True,
            "same_mesh": True,
            "same_boundary_conditions": True,
            "only_source_scale_changed": True,
        },
        "units": {
            "magnetic_vector_potential": "Wb/m",
            "magnetic_flux_density": "T",
            "energy": "J",
            "coenergy": "J",
        },
        "full_scale": {
            "max_abs_a": 4.0,
            "max_b": 2.0,
            "energy": 8.0,
            "coenergy": 8.0,
        },
        "scaled": {
            "source_scale": 0.5,
            "max_abs_a": 2.0,
            "max_b": 1.0,
            "energy": 2.0,
            "coenergy": 2.0,
        },
        "fieldwise_errors": {"a_relative": 0.0, "b_relative": 0.0},
        "independent_reference": {
            "solver_family": "independent_fem",
            "space": "H1 P1",
            "separate_implementation": True,
            "energy_relative_error": 0.002,
            "max_abs_a_relative_error": 0.004,
            "maximum_refinement_relative_change": 0.003,
            "maximum_linear_residual": 7.0e-15,
        },
    }


def test_linear_magnetization_scaling_gate_accepts_exact_scaling():
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["energy_ratio"] == 0.25


def test_linear_magnetization_scaling_gate_rejects_mesh_change():
    summary = _summary()
    summary["model_contract"]["same_mesh"] = False
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "same_mesh_and_boundary_contract" in result["issues"]


def test_linear_magnetization_scaling_gate_rejects_scalar_only_false_positive():
    summary = copy.deepcopy(_summary())
    summary["fieldwise_errors"]["b_relative"] = 0.1
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "b_field_scales_linearly" in result["issues"]


def test_linear_magnetization_scaling_gate_rejects_wrong_energy_ratio():
    summary = copy.deepcopy(_summary())
    summary["scaled"]["energy"] = 4.0
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "energy_scales_quadratically" in result["issues"]


def test_linear_magnetization_scaling_gate_rejects_unconverged_independent_reference():
    summary = copy.deepcopy(_summary())
    summary["independent_reference"]["maximum_refinement_relative_change"] = 0.2
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "independent_refinement_is_stable" in result["issues"]


def test_linear_magnetization_scaling_gate_is_exposed_over_mcp_wrapper():
    result = json.loads(linear_magnetization_scaling_gate(json.dumps(_summary())))
    assert result["status"] == "ok"
