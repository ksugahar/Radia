from __future__ import annotations

import json
import math

import numpy as np
import pytest

from radia_mcp.motor.server import motor_circuit_field_analysis
from radia_mcp.radia_ngsolve.circuit_system import analyze_circuit_field


def _parallel_payload(total_current=6.0, frequency_hz=0.0):
    return {
        "operation": "parallel",
        "field_matrix": [[4.0]],
        "source_matrix": [[0.0, 0.0]],
        "field_rhs": [0.0],
        "branch_impedance_ohm": [1.0, 2.0],
        "total_current_a": total_current,
        "frequency_hz": frequency_hz,
    }


def test_parallel_resistive_split_is_solved_not_assumed_equal():
    result = analyze_circuit_field(_parallel_payload())
    currents = np.array([complex(*value) for value in result["branch_current_a"]])

    assert result["status"] == "solved"
    assert result["equal_current_split_assumed"] is False
    assert currents == pytest.approx([4.0, 2.0])
    assert sum(currents) == pytest.approx(6.0)
    assert result["residual"]["maximum"] < 1.0e-12


def test_parallel_augmented_field_and_circuit_residuals_close():
    payload = _parallel_payload(total_current={"re": 2.0, "im": -0.5}, frequency_hz=50.0)
    payload.update(
        {
            "field_matrix": [[5.0, 1.0], [1.0, 3.0]],
            "source_matrix": [[1.0, -0.4], [0.25, 0.8]],
            "field_rhs": [0.2, -0.1],
            "branch_impedance_ohm": [{"re": 0.8, "im": 0.1}, 1.7],
        }
    )
    result = analyze_circuit_field(payload)

    assert result["augmented_shape"] == [5, 5]
    assert result["residual"]["field_inf"] < 1.0e-11
    assert result["residual"]["branch_voltage_inf"] < 1.0e-11
    assert result["residual"]["total_current_abs"] < 1.0e-11


def test_zero_total_current_retains_induced_circulating_branches():
    payload = _parallel_payload(total_current=0.0, frequency_hz=25.0)
    payload.update(
        {
            "field_matrix": [[4.0]],
            "source_matrix": [[1.0, 2.0]],
            "field_rhs": [1.0],
            "branch_impedance_ohm": [1.0, 2.0],
        }
    )
    result = analyze_circuit_field(payload)
    currents = np.array([complex(*value) for value in result["branch_current_a"]])
    flux = np.array([complex(*value) for value in result["flux_linkage_wb_turn"]])

    assert sum(currents) == pytest.approx(0.0, abs=1.0e-12)
    assert abs(currents[0]) > 1.0e-6
    assert currents[0] == pytest.approx(-currents[1])
    assert np.max(np.abs(flux)) > 0.0


def test_series_solution_preserves_signed_source_flux_and_voltage():
    result = analyze_circuit_field(
        {
            "operation": "series",
            "field_matrix": [[3.0]],
            "source_matrix": [[2.0, -1.0]],
            "field_rhs": [0.0],
            "branch_impedance_ohm": [0.5, 0.25],
            "circuit_current_a": 4.0,
            "frequency_hz": 10.0,
        }
    )
    field = complex(*result["field_state"][0])
    flux = [complex(*value) for value in result["flux_linkage_wb_turn"]]

    assert field == pytest.approx(4.0 / 3.0)
    assert flux == pytest.approx([8.0 / 3.0, -4.0 / 3.0])
    assert result["residual"]["field_inf"] < 1.0e-12
    assert abs(complex(*result["circuit_terminal_voltage_v"])) > 0.0


def test_annular_and_planar_age_sweeps_close_without_refactorization():
    annular = analyze_circuit_field(
        {
            "operation": "age_sweep",
            "kind": "annular_age",
            "positions": [0.0, math.pi / 3.0, 2.0 * math.pi],
            "harmonics": [1, 3],
        }
    )
    planar = analyze_circuit_field(
        {
            "operation": "age_sweep",
            "kind": "planar_age",
            "positions": [0.0, 0.01],
            "wavenumbers_per_m": [2.0 * math.pi / 0.01],
        }
    )

    for result in (annular, planar):
        assert result["factorization_count"] == 1
        assert result["operator_rebuild_count"] == 0
        assert result["mesh_rebuild_count"] == 0
        assert max(result["endpoint_closure_error"].values()) < 1.0e-12


def test_state_space_reuses_native_simulink_mex_contract():
    result = analyze_circuit_field(
        {
            "operation": "state_space",
            "field_matrix": [[4.0, 1.0], [1.0, 3.0]],
            "source_matrix": [[1.0, 0.0], [0.0, 1.0]],
            "field_rhs": [0.0, 0.0],
            "branch_resistance_ohm": [1.0, 2.0],
            "sample_time_s": 1.0e-4,
            "voltage_input_mode": "common",
        }
    )

    assert result["backend"] == "native-mex-sfunction"
    assert result["mex_s_function"] == "radia_state_space_mex_sfunction"
    assert result["state_order"] == 2
    assert result["input_count"] == 1
    assert result["output_count"] == 4
    assert np.asarray(result["Ad"]).shape == (2, 2)
    assert np.asarray(result["Bd"]).shape == (2, 1)
    assert result["stable"] is True
    assert result["python_per_step"] is False


def test_invalid_or_singular_systems_fail_closed_and_mcp_calls_real_kernel():
    singular = _parallel_payload()
    singular["field_matrix"] = [[0.0]]
    with pytest.raises(ValueError, match="singular"):
        analyze_circuit_field(singular)

    response = json.loads(motor_circuit_field_analysis(json.dumps(_parallel_payload())))
    assert response["status"] == "solved"
    invalid = json.loads(motor_circuit_field_analysis('{"operation":"parallel"}'))
    assert invalid["status"] == "invalid_input"

    with pytest.raises(ValueError, match="positive integers"):
        analyze_circuit_field(
            {
                "operation": "age_sweep",
                "kind": "annular_age",
                "positions": [0.0],
                "harmonics": [0],
            }
        )
    with pytest.raises(ValueError, match="positive and finite"):
        analyze_circuit_field(
            {
                "operation": "age_sweep",
                "kind": "planar_age",
                "positions": [0.0],
                "wavenumbers_per_m": [-1.0],
            }
        )
