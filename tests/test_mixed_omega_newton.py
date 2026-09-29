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
    from netgen.meshing import Element0D
    from radia.kelvin_solver import project_source_interface_potential
    mesh, source, _, table = helper._picard_case(maxh=0.65)
    material = mesh.GetMaterials().index('total') + 1
    vertex = next(v for e in mesh.ngmesh.Elements3D() if e.index == material
                  for v in e.vertices if mesh.ngmesh.Points()[v].p[0] > 0.)
    gauge_index = len(mesh.GetBBBoundaries()) + 1
    mesh.ngmesh.Add(Element0D(vertex, index=gauge_index))
    mesh.ngmesh.SetCD3Name(gauge_index, 'GND')
    with ng.TaskManager():
        mesh = ng.Mesh(mesh.ngmesh)
        potential = project_source_interface_potential(
            mesh, source, 'source_total_interface', order=2)['potential']
    return mesh, source, potential, table


def run(case, **options):
    mesh, source, potential, table = case
    args = dict(bh_table=table, nonlinear_materials=('total',),
                reduced_materials=('reduced',), total_materials=('total',),
                interface_boundary='source_total_interface', dirichlet_bbbnd='GND',
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
            interface_boundary='source_total_interface', dirichlet_bbbnd='GND', kelvin_mats=())
    for p in [(0.5, 0.1, 0.2),(-0.5, 0.1, 0.2)]:
        np.testing.assert_allclose(result['B_cf'](mesh(*p)), reference['B_cf'](mesh(*p)), rtol=1e-7, atol=1e-13)


def test_material_quadrature_is_set_apart_from_the_assembly_bonus(case):
    mesh = case[0]
    default = run(case, order=2, bonus_intorder=6)
    same = run(case, order=2, bonus_intorder=6, material_bonus_intorder=6)
    lower = run(case, order=2, bonus_intorder=6, material_bonus_intorder=2)
    assert default['nonlinear_stats']['material_bonus_intorder'] == 6
    assert lower['nonlinear_stats']['material_bonus_intorder'] == 2
    assert lower['nonlinear_stats']['bonus_intorder'] == 6
    assert lower['nonlinear_stats']['converged']
    for point in [(0.5, 0.1, 0.2), (-0.5, 0.1, 0.2)]:
        # None keeps the assembly bonus (threaded assembly rounds run to run).
        np.testing.assert_allclose(same['B_cf'](mesh(*point)), default['B_cf'](mesh(*point)),
                                   rtol=1e-12, atol=1e-18)
        # A lower co-energy rule changes only the quadrature of a smooth law.
        np.testing.assert_allclose(lower['B_cf'](mesh(*point)), default['B_cf'](mesh(*point)),
                                   rtol=1e-3, atol=1e-10)
    for bad in (-1, 2.0):
        with pytest.raises(ValueError, match='material_bonus_intorder'):
            run(case, material_bonus_intorder=bad)


def test_spline_law_is_the_pchip_law_and_its_coenergy():
    from radia.esrf_examples import get_esrf_bh_table
    from radia.scalar_potential_solver import (
        _build_bh_coefficient_function, _build_bh_coenergy_coefficient_function,
        _build_bh_spline_law)
    table = np.asarray(get_esrf_bh_table(6), dtype=float)
    b_of, coenergy_of, limit = _build_bh_spline_law(table)
    mesh = ng.Mesh(ng.unit_cube.GenerateMesh(maxh=0.5))
    point = mesh(0.5, 0.5, 0.5)
    h = ng.Parameter(0.0)
    pairs = ((b_of(h), _build_bh_coefficient_function(h, table)),
             (coenergy_of(h), _build_bh_coenergy_coefficient_function(h, table)))
    # Every knot, the tail piece, H_limit itself (where an NGSolve BSpline is
    # zero) and the analytic continuation beyond it.
    for value in np.concatenate([table[:, 0], np.geomspace(1e-6, 5 * limit, 400), [limit]]):
        h.Set(float(value))
        for spline, reference in pairs:
            expected = reference(point)
            assert spline(point) == pytest.approx(expected, rel=1e-13, abs=1e-300)
    with pytest.raises(ValueError, match=r'\[0, 0\]'):
        _build_bh_spline_law([[1., 0.], [2., 1.]])


def test_iron_integrals_never_evaluate_the_air_source(case):
    """A compiled MaterialCF evaluates every entry; the coil lives in the air one."""
    from radia.esrf_examples import get_esrf_bh_table
    from radia.scalar_potential_solver import _build_bh_coefficient_function
    mesh, source, potential, table = case
    # A costly air-only addition (221-interval nested IfPos) of negligible size.
    costly = 1e-30 * _build_bh_coefficient_function(
        ng.sqrt(ng.x * ng.x + 1.0), get_esrf_bh_table(6))
    cheap = run(case, order=2)
    heavy = run((mesh, source + ng.CF((costly, costly, costly)), potential, table), order=2)
    cheap_s = sum(row['field_change_s'] for row in cheap['nonlinear_stats']['history'])
    heavy_s = sum(row['field_change_s'] for row in heavy['nonlinear_stats']['history'])
    # Evaluating the air entry at iron points costs about 80 times more here.
    assert heavy_s < 5.0 * cheap_s + 0.05
    point = mesh(0.5, 0.1, 0.2)
    np.testing.assert_allclose(heavy['B_cf'](point), cheap['B_cf'](point), rtol=1e-9, atol=1e-16)


def test_public_workflow_rejects_material_bonus_without_bh_iron():
    from radia.static_electromagnet import solve_static_electromagnet_mixed_total_reduced_omega
    with pytest.raises(ValueError, match='nonlinear \\(bh_table\\) solve only'):
        solve_static_electromagnet_mixed_total_reduced_omega(
            None, None, None, 1., (0., 0., 0.), order=1,
            linear_mu_r_by_material={'iron': 100.}, nonlinear_material_bonus_intorder=2)


@pytest.mark.parametrize('order', [1, 2])
def test_hybrid_linear_total_material_keeps_its_permeability(order):
    """A linear total material needs its own mu_r; a nonlinear one cannot be overridden."""
    from netgen.occ import Box, Pnt, Glue, OCCGeometry, X
    from radia.kelvin_material import MU_0
    air = Box(Pnt(-1, -1, -1), Pnt(0, 1, 1))
    iron = Box(Pnt(0, -1, -1), Pnt(1, 1, 1))
    pole = Box(Pnt(1, -1, -1), Pnt(2, 1, 1))
    air.mat('reduced'); iron.mat('total'); pole.mat('pole')
    air.faces.name = iron.faces.name = pole.faces.name = 'outer'
    air.faces.Max(X).name = iron.faces.Min(X).name = 'source_total_interface'
    mesh = ng.Mesh(OCCGeometry(Glue([air, iron, pole])).GenerateMesh(maxh=1))
    options = dict(bh_table=[(0, 0), (1, 100*MU_0), (2, 150*MU_0)], nonlinear_materials=('total',),
                   reduced_materials=('reduced',), total_materials=('total', 'pole'),
                   interface_boundary='source_total_interface', dirichlet_bbbnd='outer',
                   kelvin_mats=(), order=order, tolerance=1e-6)
    zero = ng.CF((0., 0., 0.))
    with ng.TaskManager():
        # Before the check existed, the pole was silently solved as air.
        with pytest.raises(ValueError, match='linear total materials'):
            solve(mesh, zero, ng.CF(0.), 1., (3., 0., 0.), **options)
        result = solve(mesh, zero, ng.CF(0.), 1., (3., 0., 0.),
                       mu_r_by_material={'pole': 5000.}, **options)
        with pytest.raises(ValueError, match='override nonlinear'):
            solve(mesh, zero, ng.CF(0.), 1., (3., 0., 0.),
                  mu_r_by_material={'pole': 5000., 'total': 1.}, **options)
    assert result['mu_cf'](mesh(1.5, .13, .17)) / MU_0 == pytest.approx(5000.)


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


def test_matching_trace_rejects_a_missing_point_gauge(case):
    with pytest.raises(ValueError, match="requires a gauge"):
        run(case, order=2, condense_matching_trace=True, linear_solver='cg',
            dirichlet_bbbnd='missing')
