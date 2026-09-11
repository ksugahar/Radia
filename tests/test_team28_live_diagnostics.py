"""Protect the live reduced-operator probe without rerunning the FE solve."""
from types import SimpleNamespace

import numpy as np
import pytest

from validation_test.maglev.team28_hcurl_vim_force import _measured_interaction_diagnostics


def interaction(matrix):
    return SimpleNamespace(
        diagnostics=lambda: {"minimum_eigenvalue_H": None},
        matrix=SimpleNamespace(shape=matrix.shape, to_dense=lambda: matrix),
    )


def test_live_probe_measures_eigenvalues():
    result = _measured_interaction_diagnostics(interaction(np.diag([2.0, 5.0])))
    assert result["minimum_eigenvalue_H"] == 2.0
    assert result["maximum_eigenvalue_H"] == 5.0


@pytest.mark.parametrize("matrix", [np.array([[np.nan]]), np.array([[1., 1.], [0., 1.]])])
def test_live_probe_rejects_invalid_operator(matrix):
    with pytest.raises(ValueError, match="finite and Hermitian"):
        _measured_interaction_diagnostics(interaction(matrix))


def test_live_probe_is_bounded():
    with pytest.raises(ValueError, match="128 reduced modes"):
        _measured_interaction_diagnostics(interaction(np.eye(129)))
