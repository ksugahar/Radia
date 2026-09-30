"""BEM SIBC coefficient support and augmented-equation acceptance contracts."""
import numpy as np
import pytest

from radia import bem_sibc_solver as module


def _two_mode_surface():
    solver = module.ScalarBIESIBCSolver.__new__(module.ScalarBIESIBCSolver)
    solver._intree_lagrange_p2 = True
    solver.ndof = 2
    solver.M = np.eye(2)
    solver.M_inv = np.eye(2)
    solver.K = np.array([[1., -1.], [-1., 1.]])
    solver.SL = np.eye(2)
    solver.DL = np.zeros((2, 2))
    solver._c_gauge = np.ones(2)
    return solver


@pytest.mark.parametrize("impedance", [0j, 0.01+0.02j])
def test_uniform_surface_matches_independent_eigenmode_solution(impedance):
    solver = _two_mode_surface()
    omega = 1000.
    source = np.array([1., -1.])
    expected = source / (0.5 + 2*impedance/(1j*omega*module.MU_0))
    result = solver.solve(source, impedance, omega)
    repeated = solver.solve(source, np.full(2, impedance), omega)
    np.testing.assert_allclose(result['phi_vec'], expected, rtol=1e-13)
    np.testing.assert_array_equal(result['phi_vec'], repeated['phi_vec'])
    assert result['linear_residual_rel'] <= module.RELATIVE_LIMIT
    assert result['linear_residual_limit'] == module.RELATIVE_LIMIT


def test_variable_surface_is_rejected_before_wrong_observation_scaling():
    solver = _two_mode_surface()
    with pytest.raises(NotImplementedError, match='source-side coefficient-weighted'):
        solver.solve(np.array([1., -1.]), np.array([0.01, 0.02]), 1000.)


@pytest.mark.parametrize("info", [1, -1])
def test_gmres_failure_is_never_accepted_even_with_zero_residual(info):
    with pytest.raises(RuntimeError, match='GMRES failed'):
        module._check_bie_solution(np.eye(2), np.ones(2), np.ones(2), info=info)


@pytest.mark.parametrize("scale", [1., 1e-20])
def test_true_residual_uses_fixed_load_scale(scale):
    rhs = scale*np.array([1., 0.])
    with pytest.raises(RuntimeError, match='true relative residual'):
        module._check_bie_solution(np.eye(2), rhs*1.000002, rhs)
    assert module._check_bie_solution(np.eye(2), rhs, rhs) == 0.


def test_dense_bad_solution_is_rejected(monkeypatch):
    solver = _two_mode_surface()
    monkeypatch.setattr(module, 'scipy_solve', lambda matrix, rhs: np.zeros_like(rhs))
    with pytest.raises(RuntimeError, match='true relative residual'):
        solver._solve_with_gauge(np.eye(2), np.array([1., -1.]))
