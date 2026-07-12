from __future__ import annotations

import copy
import json

import numpy as np

from radia_mcp.radia_ngsolve.cq_urn import bdf_delta
from radia_mcp.radia_ngsolve.coupled_cq_refinement_gate import (
    coupled_cq_refinement_gate as gate,
)
from radia_mcp.radia_ngsolve.server import coupled_cq_refinement_gate


def _summary() -> dict:
    runs = []
    for name, method, count, step in (
        ("bdf1_coarse", "BDF1", 8, 0.6),
        ("bdf2_coarse", "BDF2", 8, 0.6),
        ("bdf2_medium", "BDF2", 15, 0.3),
        ("bdf2_fine", "BDF2", 29, 0.15),
        ("bdf2_medium_replay", "BDF2", 15, 0.3),
    ):
        radius = 1.0e-4 ** (1.0 / count)
        zeta = radius * np.exp(-2j * np.pi * np.arange(count) / count)
        laplace = bdf_delta(zeta, method.lower()) / step
        runs.append({
            "name": name,
            "method": method,
            "num_time": count,
            "time_step_s": step,
            "cq_radius": radius,
            "cq_laplace_parameter_real": laplace.real.tolist(),
            "cq_laplace_parameter_imag": laplace.imag.tolist(),
            "max_relative_residual": 1.0e-16,
            "max_abs_interior_pressure": 1.0,
            "max_abs_exterior_pressure": 0.1,
            "max_imag_interior_before_real": 1.0e-12,
            "max_imag_exterior_before_real": 1.0e-13,
            "num_volume_nodes": 5,
            "num_boundary_nodes": 4,
            "num_interior_only_nodes": 1,
        })
    return {
        "model_contract": {
            "volume_element": "P1 tetrahedron",
            "boundary_element": "P1 triangle",
            "coupling": "JohnsonNedelecCalderon",
            "time_method_family": "Lubich CQ",
            "absorbing_boundary": "none",
        },
        "runs": runs,
        "comparison": {
            "bdf1_coarse_to_bdf2_fine_relative_error": 1.05,
            "bdf2_coarse_to_fine_relative_error": 0.87,
            "bdf2_medium_to_fine_relative_error": 0.67,
            "bdf2_refinement_error_ratio": 0.77,
            "bdf2_medium_replay_relative_error": 0.0,
            "initial_exterior_relative_amplitude": 3.7e-4,
            "trusted_window_end_s": 3.15,
            "contour_alias_target": 1.0e-4,
            "trusted_window_fraction": 0.75,
        },
    }


def test_gate_accepts_coupled_cq_refinement_contract() -> None:
    result = gate(_summary())
    assert result["status"] == "ok"
    assert result["metrics"]["maximum_bdf_symbol_relative_error"] < 1.0e-14


def test_gate_accepts_descriptive_solver_method_names() -> None:
    summary = _summary()
    for run in summary["runs"]:
        order = run["method"].lower()
        run["method"] = (
            f"lubich_{order}_cq_volume_p1_fem_"
            "johnson_nedelec_calderon_bem_coupling"
        )
    assert gate(summary)["status"] == "ok"


def test_gate_rejects_symbol_and_trusted_window_drift() -> None:
    summary = copy.deepcopy(_summary())
    summary["runs"][3]["cq_laplace_parameter_imag"][2] += 0.1
    summary["comparison"]["trusted_window_fraction"] = 1.0
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert "bdf_symbols_match_saved_laplace_nodes" in result["issues"]
    assert "trusted_window_is_declared" in result["issues"]


def test_mcp_wrapper_returns_structured_invalid_input() -> None:
    result = json.loads(coupled_cq_refinement_gate("{}"))
    assert result["status"] == "invalid_input"
