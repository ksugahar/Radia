from __future__ import annotations

import copy
import math
from pathlib import Path

import numpy as np
import pytest

from radia_mcp.radia_ngsolve.vol2d_circuit import write_structured_rect_vol
from radia_mcp.radia_ngsolve.vol2d_dynamics import (
    assemble_vol2d_dynamics,
    compile_vol2d_state_space,
    solve_vol2d_nonlinear_static,
    solve_vol2d_transient,
)


MU0 = 4.0e-7 * math.pi


def _request(tmp_path: Path, *, quads: bool = False, family: str = "P1", formulation: str = "planar") -> dict:
    path = tmp_path / ("dynamic_quad.vol" if quads else "dynamic_tri.vol")
    write_structured_rect_vol(
        path,
        x0=0.1 if formulation == "axisymmetric_henrotte" else 0.0,
        x1=1.0,
        nx=5,
        ny=5,
        quads=quads,
        material="domain",
    )
    return {
        "vol_text": path.read_text(encoding="utf-8"),
        "source_name": path.name,
        "element_family": family,
        "formulation": formulation,
        "dirichlet_boundaries": ["bottom", "right", "top", "left"],
        "branches": [{"name": "coil", "material": "domain", "turns": 12.0}],
        "materials": {
            "domain": {
                "permeability_h_per_m": 200.0 * MU0,
                "conductivity_s_per_m": 3.0,
            }
        },
    }


@pytest.mark.parametrize(
    "quads,family,formulation",
    [
        (False, "P1", "planar"),
        (True, "Q1", "planar"),
        (False, "P2", "axisymmetric_henrotte"),
        (True, "Q2", "axisymmetric_henrotte"),
    ],
)
def test_vol2d_k_m_source_assembly_is_symmetric_psd(
    tmp_path: Path, quads: bool, family: str, formulation: str
) -> None:
    result = assemble_vol2d_dynamics(
        _request(tmp_path, quads=quads, family=family, formulation=formulation)
    )
    stiffness = np.asarray(result["assembly"]["field_matrix"])
    mass = np.asarray(result["conductivity_mass_matrix"])

    assert result["status"] == "assembled"
    assert np.max(np.abs(stiffness - stiffness.T)) < 1.0e-8 * np.max(np.abs(stiffness))
    assert np.max(np.abs(mass - mass.T)) < 1.0e-12 * max(1.0, np.max(np.abs(mass)))
    assert result["mass_matrix_minimum_eigenvalue"] > -1.0e-10
    assert len(result["operator_sha256"]) == 64


def test_nonlinear_picard_is_sublinear_and_reports_history(tmp_path: Path) -> None:
    request = _request(tmp_path)
    request["materials"]["domain"] = {
        "bh_curve": [
            {"b_t": 0.0, "h_a_per_m": 0.0},
            {"b_t": 0.2, "h_a_per_m": 20.0},
            {"b_t": 0.8, "h_a_per_m": 200.0},
            {"b_t": 1.2, "h_a_per_m": 2000.0},
            {"b_t": 1.6, "h_a_per_m": 20000.0},
        ],
        "conductivity_s_per_m": 3.0,
    }
    request.update({"branch_current_a": [2.0], "relaxation": 0.4, "relative_tolerance": 1.0e-7})
    low = solve_vol2d_nonlinear_static(request)
    request["branch_current_a"] = [10.0]
    high = solve_vol2d_nonlinear_static(request)
    response_ratio = high["field_state_l2"] / low["field_state_l2"]

    assert low["converged"] and high["converged"]
    assert low["iterations"] >= 3
    assert high["residual"]["field_inf"] < 1.0e-5
    assert 1.0 < response_ratio < 5.0


def test_backward_euler_passive_decay_and_static_limit(tmp_path: Path) -> None:
    request = _request(tmp_path)
    assembly = assemble_vol2d_dynamics(request)
    stiffness = np.asarray(assembly["assembly"]["field_matrix"])
    sources = np.asarray(assembly["assembly"]["source_matrix"])
    initial = np.linalg.solve(stiffness, sources[:, 0]).tolist()
    request.update(
        {
            "time_s": [0.0, 0.02, 0.04, 0.06],
            "branch_current_history_a": [[0.0], [0.0], [0.0], [0.0]],
            "initial_state": initial,
            "theta": 1.0,
        }
    )
    decay = solve_vol2d_transient(request)
    assert decay["passive_energy_decay"] is True
    assert decay["magnetic_energy_history_j"][-1] < decay["magnetic_energy_history_j"][0]
    assert decay["maximum_step_residual_inf"] < 1.0e-10

    static_request = copy.deepcopy(request)
    static_request["materials"]["domain"]["conductivity_s_per_m"] = 0.0
    static_request["initial_state"] = "zero"
    static_request["branch_current_history_a"] = [[0.0], [1.0], [1.0], [1.0]]
    static = solve_vol2d_transient(static_request)
    expected = np.linalg.solve(
        np.asarray(assemble_vol2d_dynamics(static_request)["assembly"]["field_matrix"]),
        np.asarray(assemble_vol2d_dynamics(static_request)["assembly"]["source_matrix"])[:, 0],
    )
    assert np.allclose(static["field_state_history"][-1], expected, rtol=1.0e-11, atol=1.0e-12)


def test_vol_identity_reaches_native_mex_state_space_contract(tmp_path: Path) -> None:
    request = _request(tmp_path)
    request["state_space"] = {
        "branch_resistance_ohm": [0.4],
        "sample_time_s": 1.0e-3,
        "voltage_input_mode": "per_branch",
    }
    result = compile_vol2d_state_space(request)

    assert result["backend"] == "native-mex-sfunction"
    assert result["stable"] is True
    assert result["python_per_step"] is False
    assert set(result["artifact_identity"]) == {
        "mesh_contract_sha256",
        "material_contract_sha256",
        "operator_sha256",
    }


def test_dynamic_request_negative_controls_fail_closed(tmp_path: Path) -> None:
    request = _request(tmp_path)
    request.update(
        {
            "time_s": [0.0, 0.1],
            "branch_current_history_a": [[0.0], [1.0]],
            "theta": 0.49,
        }
    )
    with pytest.raises(ValueError, match="theta"):
        solve_vol2d_transient(request)

    request["theta"] = 1.0
    request["expected_operator_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="expected_operator_sha256"):
        solve_vol2d_transient(request)

    del request["expected_operator_sha256"]
    request["branch_current_history_a"] = [[0.0, 1.0], [1.0, 2.0]]
    with pytest.raises(ValueError, match="shape"):
        solve_vol2d_transient(request)
