"""Focused contract tests for the improved-DUCAS panel integration."""

from types import SimpleNamespace

import numpy as np
import pytest

from radia.panels.calc_streamfunction import (
    _aca_factorization_metadata,
    _aca_tolerance,
    _abe_result_metadata,
    _parse_abe_allowed_modes,
    _solve_and_metrics,
)
from radia.stream_function import aca_tsvd


def test_surface_aca_factorization_provenance_is_explicit():
    factor = SimpleNamespace(k_aca=7, modes=5)
    args = SimpleNamespace(aca_eps=2.5e-9)

    assert _aca_factorization_metadata({"base": factor}, args) == {
        "factorization": "aca_plus_qr_tsvd",
        "aca_eps": 2.5e-9,
        "aca_rank": 7,
        "tsvd_modes": 5,
    }
    assert _aca_tolerance(SimpleNamespace()) == 1.0e-10
    with pytest.raises(ValueError, match="finite and positive"):
        _aca_tolerance(SimpleNamespace(aca_eps=0.0))

    weighted_factor = SimpleNamespace(k_aca=4, modes=3)
    problem = {
        "base": factor,
        "_last_abe_solution": SimpleNamespace(factor=weighted_factor),
    }
    assert _aca_factorization_metadata(problem, args)["aca_rank"] == 4


def _args(**updates):
    values = {
        "inverse_method": "abe",
        "alpha": 0.0,
        "abe_node_weights": "uniform",
        "abe_initial_potential": "",
        "abe_allowed_modes": "",
        "abe_min_mode_strength": 0.0,
        "abe_min_target_correlation": 0.0,
        "abe_residual_peak_to_peak": None,
        "abe_residual_rms": None,
        "abe_distance_floor": None,
    }
    values.update(updates)
    return SimpleNamespace(**values)


def test_abe_cli_mode_numbers_are_one_based_and_non_contiguous():
    assert _parse_abe_allowed_modes("1,7,17,37").tolist() == [0, 6, 16, 36]
    assert _parse_abe_allowed_modes("1-3,7").tolist() == [0, 1, 2, 6]


def test_panel_abe_path_reports_selected_physical_modes():
    A = np.diag([9.0, 5.0, 2.0])
    target = np.array([1.0, 1.0e-3, 1.0])
    base = aca_tsvd(3, 3, lambda i, j: float(A[i, j]), method="dense")
    problem = {
        "Af": A, "Bc": target, "Am": A, "Bm": target,
        "reg": None, "md": 1.0, "n_free": 3, "base": base,
        "R": np.eye(3), "node_points": np.zeros((3, 3)),
        "cpts": np.ones((2, 3)),
    }
    potential, fit_rms, homogeneity = _solve_and_metrics(
        problem, 0.0, _args(abe_min_mode_strength=0.1))
    result = problem["_last_abe_solution"]
    assert result.selected_modes.tolist() == [0, 2]
    assert np.allclose(potential, [1.0 / 9.0, 0.0, 0.5])
    assert fit_rms == homogeneity
    assert fit_rms > 0.0
    metadata = _abe_result_metadata(problem)
    assert metadata["abe_solve_converged"]
    assert not metadata["abe_residual_target_specified"]
    assert metadata["abe_residual_target_met"] is None


def test_shield_weight_does_not_change_abe_residual_units():
    physical = np.array([[1.0, 0.2], [0.4, 1.0]])
    weights = np.array([1.0, 7.0])
    target = np.array([1.0, 0.0])
    weighted = weights[:, None] * physical
    base = aca_tsvd(2, 2, lambda i, j: float(weighted[i, j]),
                    method="dense")
    problem = {
        "Af": weighted, "Bc": target, "Am": weighted, "Bm": target,
        "abe_response": physical, "abe_target": target,
        "abe_field_weights": weights,
        "reg": None, "md": 1.0, "n_free": 2, "base": base,
        "R": np.eye(2), "node_points": np.zeros((2, 3)),
        "cpts": np.ones((2, 3)),
    }
    potential, _, _ = _solve_and_metrics(
        problem, 0.0, _args(abe_allowed_modes="1"))
    result = problem["_last_abe_solution"]
    physical_residual = target - physical @ potential
    np.testing.assert_allclose(result.residual_field, physical_residual)
    assert result.residual_rms == np.sqrt(np.mean(physical_residual**2))
