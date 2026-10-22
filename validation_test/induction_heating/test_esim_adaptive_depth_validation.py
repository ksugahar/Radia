"""ESIM physical/analytic and mesh convergence validation; no private inputs.

Analytical references certify the declared linear cell formulation. Nonlinear
seed agreement is discretization consistency, not independent FEM validation.
"""
import numpy as np
import pytest
from scipy.special import ive

from radia.esim_cell_problem import ESIMCellProblemSolver, ESIMFiniteSlabSolver

MU0 = 4e-7*np.pi
LINEAR = [[h, 100*MU0*h] for h in (0., 10., 100., 1e3, 1e4, 1e5, 5e5)]


@pytest.mark.parametrize('geometry', ['slab', 'cylinder'])
def test_thin_skin_finite_cells_match_analytic(geometry):
    solver = ESIMFiniteSlabSolver(.005, sigma=5e6, frequency=222e3,
                                 mu_r=100, geometry=geometry, n_nodes=24)
    seed = solver.mesh_points.copy()
    result = solver.solve(1000.)
    gamma = (1+1j)/solver.delta
    expected = solver.rho*gamma
    expected *= (ive(1, gamma*.005)/ive(0, gamma*.005) if geometry == 'cylinder'
                 else np.tanh(gamma*.005))
    assert abs(result['Z']/expected-1) < .003
    assert result['mesh_converged']
    assert result['linear_residual'] < 1e-6
    assert result['P_prime'] == pytest.approx(.5*result['Z'].real*1e6)
    assert result['mesh_history'][-1]['nodes'] > len(seed)



@pytest.mark.parametrize('nodes', [100, 800, 1000])
def test_halfspace_refinement_does_not_exponentially_collapse(nodes):
    solver = ESIMCellProblemSolver(bh_curve=LINEAR, sigma=5e6,
                                  frequency=222e3, n_nodes=nodes)
    result = solver.solve(1000.)
    expected = (1+1j)*np.sqrt(2*np.pi*222e3*100*MU0/(2*5e6))
    assert abs(result['Z']/expected-1) < .003
    assert result['depth_converged']
    assert result['Z'].real > 0



def test_nonlinear_synthetic_cell_seed_convergence():
    fields = np.r_[0., np.geomspace(1., 1e6, 35)]
    bh = [[h, MU0*h + 1.2*np.tanh(h/1500.)] for h in fields]
    answers = []
    for count in (24, 60):
        solver = ESIMFiniteSlabSolver(.005, bh_curve=bh, sigma=5e6,
                    frequency=222e3, geometry='cylinder', n_nodes=count)
        answers.append(solver.solve(4000.))
    assert abs(answers[0]['Z']/answers[1]['Z']-1) < .003
    assert abs(answers[0]['P_prime']/answers[1]['P_prime']-1) < .003


def test_lossy_material_sign_and_halfspace_depth_are_checked():
    cell = ESIMCellProblemSolver(complex_mu=(100., 20.), sigma=5e6,
                                frequency=222e3, n_nodes=40, num_skin_depths=1)
    result = cell.solve(1000.)
    expected = np.sqrt(1j*cell.omega*MU0*(100-20j)*cell.rho)
    assert abs(result['Z']/expected-1) < .003
    assert len(result['depth_history']) > 1
    assert result['P_magnetic'] > 0
