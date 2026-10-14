"""Fast algebraic contracts for complete (not loss-inferred) reciprocity."""
import numpy as np
import pytest
from radia.workpiece_surface import _complete_sibc_reaction
from radia.workpiece_surface import _check_sibc_reaction_power


@pytest.mark.parametrize('surface,reaction', [(6.,1.), (1.,-1.), (float('nan'),1.), (1.,float('inf'))])
def test_inconsistent_power_is_not_publishable(surface,reaction):
    with pytest.raises(RuntimeError, match='SIBC reciprocal power'):
        _check_sibc_reaction_power(surface,reaction)


def test_consistent_and_zero_power():
    assert _check_sibc_reaction_power(10,10.01) == pytest.approx(.001)
    assert _check_sibc_reaction_power(0,0) == 0


def test_electric_term_is_bilinear_and_current_normalized():
    k = np.array([[1., -1.], [-1., 1.]])
    incident = np.array([1+2j, 0j])
    total = np.array([3-4j, 0j])
    z, omega, current = 2+3j, 10., 2.
    expected = .7j + z/(1j*omega*current**2)*(1+2j)*(3-4j)
    assert _complete_sibc_reaction(.7j, incident, total, k, z, omega, current) == pytest.approx(expected)
    assert _complete_sibc_reaction(.7j, incident+9, total-7j, k, z, omega, current) == pytest.approx(expected)


def test_current_rescaling_preserves_inductance():
    k = np.array([[1., -1.], [-1., 1.]])
    p, q = np.array([1., 0]), np.array([1-2j, 0])
    a = _complete_sibc_reaction(1j, p, q, k, 1+1j, 4., 1.)
    b = _complete_sibc_reaction(1j, 3*p, 3*q, k, 1+1j, 4., 3.)
    assert a == pytest.approx(b)


@pytest.mark.parametrize('z,w,i', [(np.ones(2), 1, 1), (-1+1j,1,1), (1,0,1), (1,1,0), (complex('nan'),1,1)])
def test_invalid_reciprocity_rejected(z,w,i):
    with pytest.raises(ValueError):
        _complete_sibc_reaction(0, np.ones(2), np.ones(2), np.eye(2), z,w,i)


def test_mismatched_basis_rejected():
    with pytest.raises(ValueError, match='matching'):
        _complete_sibc_reaction(0, np.ones(2), np.ones(3), np.eye(2), 1,1,1)
