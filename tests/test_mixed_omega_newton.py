"""Nonlinear residual, field and failure contracts of mixed Omega Newton."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

ng = pytest.importorskip('ngsolve')
from radia.mixed_omega_newton import (
    solve_magnetostatic_mixed_total_reduced_omega_newton_kelvin as solve,
    MixedOmegaNewtonNotConverged)


@pytest.fixture(scope='module')
def case():
    spec = importlib.util.spec_from_file_location('mixed_helpers', Path(__file__).with_name('test_kelvin_mixed_omega.py'))
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    return helper._picard_case(maxh=0.65)


def run(case, **options):
    mesh, source, potential, table = case
    args = dict(bh_table=table, nonlinear_materials=('total',),
                reduced_materials=('reduced',), total_materials=('total',),
                interface_boundary='source_total_interface', dirichlet_bbbnd='outer',
                kelvin_mats=(), tolerance=1e-6, residual_tolerance=1e-9)
    args.update(options)
    with ng.TaskManager():
        return solve(mesh, source, potential, 1., (3., 0., 0.), **args)


@pytest.mark.parametrize('order', [1, 2])
def test_newton_solves_the_quadrature_law_and_reduces_the_residual(case, order):
    result = run(case, order=order)
    stats = result['nonlinear_stats']
    assert stats['converged']
    assert stats['residual_relative'] <= 1e-9
    assert result['constitutive_field_audit']['relative_B_constitutive_L2'] < 1e-8
    assert all(row['line_search_accepted'] for row in stats['history'])
    assert all(row['residual_relative'] <= row['residual_relative_before']
               for row in stats['history'])


def test_newton_linear_law_agrees_with_linear_mixed_solver(case):
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin
    mesh, source, potential, _ = case
    mu0 = 4e-7 * np.pi
    table = [[0., 0.], [1000., mu0 * 500. * 1000.]]
    result = run(case, bh_table=table)
    with ng.TaskManager():
        reference = solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, source, potential, 1., (3., 0., 0.), mu_r_by_material={'total':500.},
            reduced_materials=('reduced',), total_materials=('total',),
            interface_boundary='source_total_interface', dirichlet_bbbnd='outer', kelvin_mats=())
    for p in [(0.5, 0.1, 0.2),(-0.5, 0.1, 0.2)]:
        np.testing.assert_allclose(result['B_cf'](mesh(*p)), reference['B_cf'](mesh(*p)), rtol=1e-7, atol=1e-13)


def test_iteration_limit_never_returns_an_accepted_field(case):
    with pytest.raises(MixedOmegaNewtonNotConverged) as exc:
        run(case, max_iterations=1)
    assert not exc.value.state['nonlinear_stats']['converged']


def test_invalid_law_is_rejected_before_assembly(case):
    with pytest.raises(ValueError, match=r'B\(H\)'):
        run(case, bh_table=[[0.,0.],[1.,-1.]])


@pytest.mark.parametrize('order', [1, 2])
def test_matching_trace_newton_preserves_harmonic_source_and_field(case, order):
    from radia.kelvin_solver import project_source_interface_potential
    mesh, source, _, table = case
    with ng.TaskManager():
        trace = project_source_interface_potential(mesh, source, 'source_total_interface',
                                                  order=order)['potential']
    consistent = mesh, source, trace, table
    options = dict(order=order, total_source_h=ng.CF((0.01, 0., 0.02)),
                   total_source_materials=('total',))
    reference = run(consistent, **options)
    condensed = run(consistent, condense_matching_trace=True, linear_solver='cg', **options)
    assert condensed['nonlinear_stats']['converged']
    assert condensed['fes'].ndof < reference['fes'].ndof
    for point in [(0.5,0.1,0.2),(-0.5,0.1,0.2)]:
        np.testing.assert_allclose(condensed['B_cf'](mesh(*point)),
                                   reference['B_cf'](mesh(*point)), rtol=2e-6, atol=1e-12)
