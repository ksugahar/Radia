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
