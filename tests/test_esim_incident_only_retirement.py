import numpy as np
import pytest
from radia.esim_coupled_solver import InductionHeatingCoil, ESIMCoupledSolver


def test_incident_only_iteration_is_not_a_coupled_solver():
    solver = ESIMCoupledSolver(InductionHeatingCoil(coil_type='loop'), None, 1000)
    with pytest.raises(NotImplementedError, match='reaction-field'):
        solver.solve()
    assert not solver.converged


@pytest.mark.parametrize('phase', [1j, (1+1j)/np.sqrt(2), -1])
def test_coil_preserves_peak_phasor_phase(phase):
    coil = InductionHeatingCoil(coil_type='loop', radius=.05)
    reference = coil.compute_field_at_point([0,0,.02])
    coil.set_current(phase)
    np.testing.assert_allclose(coil.compute_field_at_point([0,0,.02]),phase*reference,rtol=1e-14)


@pytest.mark.parametrize('current',[2,2j, np.sqrt(2)*(1+1j)])
def test_reflected_impedance_uses_average_power_and_peak_current(current):
    coil=InductionHeatingCoil(coil_type='loop');coil.set_current(current)
    solver=ESIMCoupledSolver(coil,None,1000)
    impedance=3+4j
    power=.5*abs(current)**2*impedance
    assert solver.compute_reflected_impedance(power.real,power.imag)==pytest.approx(impedance)
    assert solver.M_mutual is None and solver.coupling_factor is None


def test_zero_current_cannot_identify_impedance():
    coil=InductionHeatingCoil(coil_type='loop');coil.set_current(0)
    solver=ESIMCoupledSolver(coil,None,1000)
    with pytest.raises(ValueError,match='nonzero'):
        solver.compute_reflected_impedance(0,0)
