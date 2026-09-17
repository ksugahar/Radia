"""Physical boundary identity under an exactly representable source lift."""
import numpy as np
import pytest


def _case():
    import ngsolve as ng
    from netgen.occ import Box, Pnt, Glue, OCCGeometry
    core = Box(Pnt(-.3, -.3, -.3), Pnt(.3, .3, .3))
    core.mat('iron')
    core.faces.name = 'interface'
    outer = Box(Pnt(-1, -1, -1), Pnt(1, 1, 1))
    outer.faces.name = 'outer'
    air = outer - core
    air.mat('air')
    mesh = ng.Mesh(OCCGeometry(Glue([core, air])).GenerateMesh(maxh=.8))
    options = dict(reduced_materials=('air',), total_materials=('iron',),
                   interface_boundary='interface', kelvin_mats=(),
                   surface_dirichlet={'reduced': {'outer': 0.}})
    return mesh, ng.CoefficientFunction((0., 0., 100.)), options


def test_surface_dirichlet_fixes_physical_field_without_point_gauge():
    import ngsolve as ng
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin as solve
    mesh, source, options = _case()
    with ng.TaskManager():
        result = solve(mesh, source, -100*ng.z, 1., (0, 0, 0),
                       mu_r_by_material={'air': 1., 'iron': 1.}, **options)
    for point in ((.1, .1, .1), (.6, .1, .1)):
        assert np.asarray(result['H_cf'](mesh(*point))) == pytest.approx([0, 0, 100], abs=1e-8)


def test_vacuum_hodge_split_requires_response_representable_lift():
    import ngsolve as ng
    from radia.kelvin_solver import (
        solve_magnetostatic_mixed_total_reduced_omega_kelvin as solve,
        project_source_total_hodge,
    )
    mesh, _, options = _case()
    source = ng.CF((2*ng.x, 0., -2*ng.z))
    selector = mesh.Materials('iron')
    errors = {}
    with ng.TaskManager():
        for response_order, projection_order in ((1,1), (1,2), (2,2)):
            split = project_source_total_hodge(mesh, source, ('iron',), order=projection_order)
            result = solve(mesh, source, split['potential'], 1., (0,0,0),
                mu_r_by_material={'air':1., 'iron':1.}, order=response_order,
                total_source_h=split['harmonic_field'], total_source_materials=('iron',),
                **options)
            delta = result['H_cf']-source
            errors[response_order, projection_order] = float(ng.Integrate(
                ng.InnerProduct(delta,delta), mesh, definedon=selector, order=6))**.5
    assert errors[1,1] < 1e-10
    assert errors[1,2] > 1e-3
    assert errors[2,2] < 1e-10


def test_constitutive_audit_detects_centroid_secant_defect():
    import ngsolve as ng
    from radia.kelvin_solver import audit_mixed_omega_constitutive_field
    from radia.scalar_potential_solver import _build_bh_coefficient_function
    mesh, _, _ = _case()
    h = np.linspace(0., 10000., 101)
    table = np.column_stack((h, 4e-7*np.pi*h + 1.5*np.tanh(h/1000)))
    magnitude = 2000 + 1800*ng.x/.3
    field = ng.CF((magnitude, 0, 0))
    target = ng.CF((_build_bh_coefficient_function(magnitude, table), 0, 0))
    centroid_mu = (4e-7*np.pi*2000+1.5*np.tanh(2))/2000
    with ng.TaskManager():
        exact = audit_mixed_omega_constitutive_field(mesh, field, target, table, ('iron',))
        frozen = audit_mixed_omega_constitutive_field(mesh, field, centroid_mu*field,
                                                    table, ('iron',))
    assert exact['relative_B_constitutive_L2'] < 1e-12
    assert frozen['relative_B_constitutive_L2'] > .2
    assert frozen['accuracy_accepted'] is False


def test_pointwise_picard_solves_linear_material_control():
    import ngsolve as ng
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin
    mesh, source, options = _case()
    table = np.array([[0.,0.], [100., 100.*4e-7*np.pi*20],
                      [1000., 1000.*4e-7*np.pi*20]])
    with ng.TaskManager():
        result = solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin(
            mesh, source, -100*ng.z, 1., (0,0,0), bh_table=table,
            nonlinear_materials=('iron',), mu_r_initial=10.,
            material_sampling='integration_point', relaxation=1.,
            max_iterations=8, tolerance=1e-9, **options)
    assert result['nonlinear_stats']['converged']
    assert result['constitutive_field_audit']['relative_B_constitutive_L2'] < 1e-9


def test_pointwise_picard_nonlinear_material_has_quadrature_consistency():
    import ngsolve as ng
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin
    mesh, source, options = _case()
    h = np.linspace(0., 1000., 21)
    table = np.column_stack((h, 4e-7*np.pi*h + .01*np.tanh(h/100)))
    with ng.TaskManager():
        result = solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin(
            mesh, source, -100*ng.z, 1., (0,0,0), bh_table=table,
            nonlinear_materials=('iron',), mu_r_initial=10.,
            material_sampling='integration_point', relaxation=1.,
            max_iterations=40, tolerance=1e-5, **options)
    assert result['nonlinear_stats']['converged']
    assert result['constitutive_field_audit']['relative_B_constitutive_L2'] < 1e-5


def test_linear_spline_table_matches_piecewise_linear_curve():
    import ngsolve as ng
    from radia.scalar_potential_solver import (
        _build_bh_linear_spline_coefficient_function,
        _build_bh_linear_spline_coenergy_coefficient_function)
    mesh, _, _ = _case()
    table = np.array([[0.,0.], [100.,.2], [200.,.35]])
    field = _build_bh_linear_spline_coefficient_function(1000*ng.x, table)
    energy = _build_bh_linear_spline_coenergy_coefficient_function(1000*ng.x, table)
    for x in (.001, .02, .07, .12, .18, .25):
        actual = float(field(mesh(x, .0, .0)))
        expected = (np.interp(1000*x, table[:,0], table[:,1])
                    if x <= .2 else table[-1,1] + 4e-7*np.pi*(1000*x-table[-1,0]))
        assert actual == pytest.approx(expected, abs=1e-12)
        step = 1e-6
        derivative = (float(energy(mesh(x+step,0,0))) -
                      float(energy(mesh(x-step,0,0))))/(2000*step)
        assert derivative == pytest.approx(actual, rel=1e-8)


def test_pointwise_picard_linear_spline_nonlinear_control():
    import ngsolve as ng
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin
    mesh, source, options = _case()
    h = np.linspace(0., 1000., 21)
    table = np.column_stack((h, 4e-7*np.pi*h + .01*np.tanh(h/100)))
    with ng.TaskManager():
        result = solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin(
            mesh, source, -100*ng.z, 1., (0,0,0), bh_table=table,
            nonlinear_materials=('iron',), mu_r_initial=10.,
            material_sampling='integration_point', bh_interpolation='linear_spline',
            relaxation=1., max_iterations=40, tolerance=1e-5, **options)
    assert result['nonlinear_stats']['converged']
    assert result['constitutive_field_audit']['interpolation'] == 'linear_spline'
    assert result['constitutive_field_audit']['relative_B_constitutive_L2'] < 1e-5
    energy = result['energy_observables']
    assert energy['energy_J'] > 0
    assert energy['coenergy_J'] > 0
    assert energy['energy_J'] + energy['coenergy_J'] == pytest.approx(
        energy['h_dot_b_integral_J'], rel=1e-12)


def test_pointwise_picard_damped_material_update_converges():
    import ngsolve as ng
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin
    mesh, source, options = _case()
    h = np.linspace(0., 1000., 21)
    table = np.column_stack((h, 4e-7*np.pi*h + .01*np.tanh(h/100)))
    with ng.TaskManager():
        result = solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin(
            mesh, source, -100*ng.z, 1., (0,0,0), bh_table=table,
            nonlinear_materials=('iron',), mu_r_initial=10.,
            material_sampling='integration_point', bh_interpolation='linear_spline',
            relaxation=.5, max_iterations=60, tolerance=1e-5, **options)
    assert result['nonlinear_stats']['converged']
    assert result['constitutive_field_audit']['relative_B_constitutive_L2'] < 1e-5


def test_total_surface_requires_source_lift_not_reduced_zero():
    import ngsolve as ng
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin as solve
    mesh, source, options = _case()
    options.update(reduced_materials=('iron',), total_materials=('air',),
                   mu_r_by_material={'air': 1., 'iron': 1.})
    with ng.TaskManager():
        options['surface_dirichlet'] = {'total': {'outer': -100*ng.z}}
        lifted = solve(mesh, source, -100*ng.z, 1., (0, 0, 0), **options)
        options['surface_dirichlet'] = {'total': {'outer': 0.}}
        zero = solve(mesh, source, -100*ng.z, 1., (0, 0, 0), **options)
    for point in ((.1, .1, .1), (.6, .1, .1)):
        assert np.asarray(lifted['H_cf'](mesh(*point))) == pytest.approx([0, 0, 100], abs=1e-8)
        assert np.linalg.norm(zero['H_cf'](mesh(*point))) < 1e-8


@pytest.mark.parametrize('bad,match', [
    ({'surface_dirichlet': []}, 'mapping'),
    ({'surface_dirichlet': {'reduced': {'missing': 0}}}, 'exact'),
    ({'surface_dirichlet': {'total': {'outer': 0}}}, 'exterior'),
    ({'surface_dirichlet': {'reduced': {'interface': 0}}}, 'exterior'),
    ({'reduced_zero_normal_boundary': 'outer'}, 'overlap'),
    ({'total_dirichlet_cf': 0.}, 'gauge'),
])
def test_surface_dirichlet_rejects_ambiguous_boundary(bad, match):
    import ngsolve as ng
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin as solve
    mesh, source, options = _case()
    options.update(bad)
    with pytest.raises(ValueError, match=match):
        solve(mesh, source, -100*ng.z, 1., (0, 0, 0), **options)


def test_dirichlet_junction_is_rejected_before_singular_multiplier_solve():
    import ngsolve as ng
    from netgen.occ import Box, Pnt, Glue, OCCGeometry, X
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin as solve
    a = Box(Pnt(-1,-1,-1), Pnt(0,1,1))
    b = Box(Pnt(0,-1,-1), Pnt(1,1,1))
    a.mat('air'); b.mat('iron')
    a.faces.name = 'left'; b.faces.name = 'right'
    a.faces.Max(X).name = 'interface'; b.faces.Min(X).name = 'interface'
    mesh = ng.Mesh(OCCGeometry(Glue([a,b])).GenerateMesh(maxh=1.))
    with pytest.raises(ValueError, match='junction'):
        solve(mesh, ng.CoefficientFunction((0,0,1)), -ng.z, 1., (0,0,0),
              reduced_materials=('air',), total_materials=('iron',),
              interface_boundary='interface', kelvin_mats=(),
              surface_dirichlet={'reduced': {'left': 0}, 'total': {'right': -ng.z}})


@pytest.mark.parametrize('order', [1, 2])
def test_nonlinear_total_and_reduced_representations_match(order):
    import ngsolve as ng
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin as solve
    mesh, source, options = _case()
    options.update(bh_table=[(0., 0.), (50., .001), (100., .0015), (200., .002)],
                   nonlinear_materials=('iron',), mu_r_initial=10.,
                   tolerance=1e-7, max_iterations=150, relaxation=.5,
                   cache_fixed_rhs=False, order=order,
                   material_update_order=order-1)
    with ng.TaskManager():
        total = solve(mesh, source, -100*ng.z, 1., (0, 0, 0), **options)
        reduced = solve(mesh, source, ng.CoefficientFunction(0.), 1., (0, 0, 0),
                        total_source_h=source, total_source_materials=('iron',), **options)
    for key in ('H_cf', 'B_cf'):
        for point in ((.1, .1, .1), (.6, .1, .1)):
            a = np.asarray(total[key](mesh(*point)))
            b = np.asarray(reduced[key](mesh(*point)))
            assert np.linalg.norm(a-b) < 1e-8*max(np.linalg.norm(a), 1e-6)
