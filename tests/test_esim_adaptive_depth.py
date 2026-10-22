"""ESIM source contracts: state preservation, errors, cache and CLI wiring."""
import numpy as np
import pytest

from radia.esim_cell_problem import ESIMFiniteSlabSolver
from radia.esim_panel_evaluator import PanelESIMEvaluator

MU0 = 4e-7*np.pi


@pytest.mark.parametrize('geometry', ['slab', 'cylinder'])
def test_adaptive_solve_preserves_mesh_seed(geometry):
    solver = ESIMFiniteSlabSolver(.005, sigma=5e6, frequency=222e3,
                                 mu_r=100, geometry=geometry, n_nodes=24)
    seed = solver.mesh_points.copy()
    solver.solve(1000.)
    np.testing.assert_array_equal(seed, solver.mesh_points)


def test_nonlinear_iteration_limit_fails_loudly():
    fields = np.r_[0., np.geomspace(1., 1e6, 35)]
    bh = [[h, MU0*h + 1.2*np.tanh(h/1500.)] for h in fields]
    solver = ESIMFiniteSlabSolver(.005, bh_curve=bh, sigma=5e6,
                frequency=222e3, geometry='cylinder', n_nodes=60)
    with pytest.raises(RuntimeError, match='iteration limit'):
        solver.solve(4000., max_iter=1)


def test_unresolved_mesh_budget_fails_and_cache_controls_are_hashed():
    solver = ESIMFiniteSlabSolver(.005, mu_r=100, sigma=5e6,
                                 frequency=222e3, n_nodes=24, max_nodes=24)
    with pytest.raises(RuntimeError, match='node budget'):
        solver.solve(1000.)
    evaluator = PanelESIMEvaluator(solver)
    solver.mesh_tol *= .5
    with pytest.raises(RuntimeError, match='configuration changed'):
        evaluator._fields(np.array([1000.]))



def test_panel_evaluator_identity_stays_stable_after_adaptive_solve():
    cell = ESIMFiniteSlabSolver(.005, mu_r=100, sigma=5e6,
                               frequency=222e3, n_nodes=24)
    evaluator = PanelESIMEvaluator(cell)
    first = evaluator._direct(1000., 'iteration_cell_calls')
    evaluator._fields(np.array([1000.]))
    second = evaluator._direct(1000., 'certification_cell_calls')
    assert first == second



@pytest.mark.parametrize('kelvin', [False, True])
def test_designspec_cell_controls_reach_the_cli_parser(kelvin):
    import argparse
    from radia.ih_design import IHDesignSpec
    from radia.esim_depth import add_depth_arguments, depth_controls
    spec = IHDesignSpec(impedance_model='Nonlinear ESIM', bh_file='synthetic.txt',
        esim_cell_max_iter=157, esim_cell_mesh_tol=.0004, esim_cell_max_nodes=5000)
    command = []
    spec._append_esim_args(command, include_anderson=False, kelvin=kelvin)
    parser = argparse.ArgumentParser()
    add_depth_arguments(parser)
    args, _ = parser.parse_known_args(command)
    assert depth_controls(args) == dict(cell_max_iter=157, mesh_tol=.0004, max_nodes=5000)
