"""Exact vertex gauges must preserve the projected physical source field."""
import numpy as np
import pytest

import ngsolve as ng
from netgen.occ import Box, Glue, OCCGeometry, Pnt
from radia.kelvin_solver import project_source_total_hodge


@pytest.fixture
def disconnected_mesh():
    solids = []
    for shift in (0., .02):
        solid = Box(Pnt(shift, 0., 0.), Pnt(shift+.006, .006, .006))
        solid.mat('iron')
        solids.append(solid)
    return ng.Mesh(OCCGeometry(Glue(solids)).GenerateMesh(maxh=.004))


@pytest.mark.parametrize('order', [1, 2])
@pytest.mark.parametrize('load', ['volume', 'surface_flux'])
def test_vertex_gauge_is_exact_and_field_independent(disconnected_mesh, order, load):
    mesh = disconnected_mesh
    source = ng.CF((1., 2., -.3))
    vertices = [[v.nr for v in mesh.vertices if (v.point[0] < .01) == first]
                for first in (True, False)]
    choices = [tuple(min(v) for v in vertices), tuple(max(v) for v in vertices)]
    results = []
    with ng.TaskManager():
        for chosen in choices:
            result = project_source_total_hodge(
                mesh, source, ('iron',), order=order, source_load=load,
                tangential_tolerance=1e-8 if load == 'surface_flux' else None,
                gauge_vertices=chosen)
            assert result['gauge']['vertices'] == chosen
            for dof in result['gauge']['dofs']:
                assert result['potential'].vec[dof] == 0.
                assert not result['free_dofs'][dof]
            assert result['linear_true_relative_residual'] < 1e-10
            assert result['gauge']['full_relative_residual'] < 1e-9
            results.append(result)
        difference = ng.grad(results[0]['potential'])-ng.grad(results[1]['potential'])
        relative = np.sqrt(ng.Integrate(ng.InnerProduct(difference, difference), mesh)
                           / ng.Integrate(ng.InnerProduct(source, source), mesh))
        assert relative < 1e-10
        for result in results:
            error = source + ng.grad(result['potential'])
            assert np.sqrt(ng.Integrate(ng.InnerProduct(error, error), mesh)
                           / ng.Integrate(ng.InnerProduct(source, source), mesh)) < 1e-10


def test_vertex_gauge_requires_one_owned_vertex_per_component(disconnected_mesh):
    mesh = disconnected_mesh
    source = ng.CF((1., 0., 0.))
    with ng.TaskManager():
        default = project_source_total_hodge(mesh, source, ('iron',))
        assert len(default['gauge']['vertices']) == 2
        for chosen in ([], [0], [0, 0], [-1, 0], [True, 1]):
            with pytest.raises(ValueError, match='gauge_vertices'):
                project_source_total_hodge(mesh, source, ('iron',), gauge_vertices=chosen)


@pytest.mark.parametrize('order', [2, 3])
def test_nonpolynomial_flux_is_balanced_per_component(disconnected_mesh, order):
    """NumberSpace's mean does not balance higher-order H1 quadrature."""
    mesh = disconnected_mesh
    # Analytic harmonic potential: div(H)=curl(H)=0. Different disconnected
    # components have the same material name, so a global correction is insufficient.
    k = 200.0
    source = ng.CF((ng.exp(k*ng.x)*ng.cos(k*ng.y),
                    -ng.exp(k*ng.x)*ng.sin(k*ng.y), 0.0))
    components = [[v.nr for v in mesh.vertices if (v.point[0] < .01) == first]
                  for first in (True, False)]
    with ng.TaskManager():
        # This deliberately coarse test isolates compatibility, not source accuracy.
        results = [project_source_total_hodge(
            mesh, source, ('iron',), order=order, bonus_intorder=0,
            source_load='surface_flux', tangential_tolerance=1.0,
            gauge_vertices=tuple(choose(v) for v in components))
            for choose in (min, max)]
        for result in results:
            balance = result['flux_balance']['component_compatibility']
            assert len(balance) == 2
            assert max(abs(c['assembled_net_load_before']) / c['raw_load_norm']
                       for c in balance) > 1e-12
            assert max(abs(c['assembled_net_load_after']) / c['raw_load_norm']
                       for c in balance) < 1e-12
            assert result['gauge']['full_relative_residual'] < 1e-10
        delta = ng.grad(results[0]['potential'])-ng.grad(results[1]['potential'])
        relative = np.sqrt(ng.Integrate(ng.InnerProduct(delta, delta), mesh)
                           / ng.Integrate(ng.InnerProduct(source, source), mesh))
        assert relative < 1e-10


@pytest.mark.parametrize('second_net', [0.0, -1.0])
def test_incompatible_components_are_rejected(disconnected_mesh, second_net):
    from radia.kelvin_solver import (
        _balance_hodge_surface_load, _hodge_vertex_gauges,
        material_set_boundary_signs)
    mesh = disconnected_mesh
    fes = ng.H1(mesh, order=2)
    with ng.TaskManager():
        _, _, pins, components = _hodge_vertex_gauges(mesh, fes, {'iron'}, None)
        load = ng.LinearForm(fes)
        load += 0.0 * fes.TestFunction() * ng.dx
        load.Assemble()
        load.vec[pins[0]] = 1.0
        load.vec[pins[1]] = second_net
        with pytest.raises(RuntimeError, match='component compatibility failed'):
            _balance_hodge_surface_load(
                mesh, fes, load, fes.TestFunction(),
                material_set_boundary_signs(mesh, ('iron',)), components, 0,
                ng.CF((1.0, 0.0, 0.0)))


def test_non_solenoidal_source_is_not_hidden_by_mean_removal(disconnected_mesh):
    with ng.TaskManager():
        with pytest.raises(RuntimeError, match='component compatibility failed'):
            project_source_total_hodge(
                disconnected_mesh, ng.CF((ng.x, ng.y, ng.z)), ('iron',),
                order=2, source_load='surface_flux', bonus_intorder=16,
                tangential_tolerance=1.0)
