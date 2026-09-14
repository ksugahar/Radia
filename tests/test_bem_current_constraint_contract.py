"""Small algebraic tests of the physical-current validation boundary."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'torus_current', ROOT / 'validation_test/bem/validate_torus_current_constraint.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_conserved_unit_current_energy_and_basis_covariance():
    matrix = np.diag([2., 2., 3.])
    divergence = np.array([[1., -1., 0.]])
    current = np.array([0., 0., 1.])
    solution, energy, nullity = MODULE.constrained_energy(matrix, divergence, current)
    np.testing.assert_allclose(solution, [0, 0, 1], atol=1e-14)
    assert energy == pytest.approx(3)
    assert nullity == 2
    transform = np.array([[0., 0., -1.], [1., 0., 0.], [0., -1., 0.]])
    changed, changed_energy, _ = MODULE.constrained_energy(
        transform.T@matrix@transform, divergence@transform, transform.T@current)
    np.testing.assert_allclose(transform@changed, solution, atol=1e-14)
    assert changed_energy == pytest.approx(energy)


def test_unresolved_loop_current_fails_instead_of_normalizing_roundoff():
    with pytest.raises(ValueError, match='no resolved'):
        MODULE.constrained_energy(np.eye(3), np.array([[1., -1., 0.]]),
                                  np.array([1., -1., 0.]))
