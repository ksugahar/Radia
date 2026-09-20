"""Protect the live reduced-operator probe without rerunning the FE solve."""
from types import SimpleNamespace

import numpy as np
import pytest

from validation_test.maglev.team28_hcurl_vim_force import _measured_interaction_diagnostics


def interaction(matrix, diagnostics=None):
    if diagnostics is None:
        diagnostics = {"minimum_eigenvalue_H": None}
    return SimpleNamespace(
        diagnostics=lambda: diagnostics,
        matrix=SimpleNamespace(shape=matrix.shape, to_dense=lambda: matrix),
    )


def test_live_probe_measures_eigenvalues():
    result = _measured_interaction_diagnostics(interaction(np.diag([2.0, 5.0])))
    assert result["minimum_eigenvalue_H"] == 2.0
    assert result["maximum_eigenvalue_H"] == 5.0
    assert result["minimum_to_maximum_eigenvalue_ratio"] == pytest.approx(0.4)


def test_live_probe_does_not_mutate_source_diagnostics():
    diagnostics = {"minimum_eigenvalue_H": None}
    _measured_interaction_diagnostics(interaction(np.eye(2), diagnostics))
    assert diagnostics == {"minimum_eigenvalue_H": None}


def test_live_probe_accepts_roundoff_at_structural_zero():
    matrix = np.array([[2.0, 1.0e-30], [0.0, 5.0]])
    result = _measured_interaction_diagnostics(interaction(matrix))
    assert result["hermitian_relative_error"] < 1.0e-10


@pytest.mark.parametrize("matrix", [np.array([[np.nan]]), np.array([[1., 1.], [0., 1.]])])
def test_live_probe_rejects_invalid_operator(matrix):
    with pytest.raises(ValueError, match="finite and Hermitian"):
        _measured_interaction_diagnostics(interaction(matrix))


def test_live_probe_is_bounded():
    with pytest.raises(ValueError, match="128 reduced modes"):
        _measured_interaction_diagnostics(interaction(np.eye(129)))
