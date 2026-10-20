"""Preserve high-order boundary fields when copying a matching H1 trace."""
import numpy as np
import pytest
import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh
from radia.kelvin_solver import _copy_matching_h1_trace_to_global


@pytest.mark.parametrize('order', [2, 3])
def test_boundary_field_survives_matching_trace_copy(order):
    mesh = MakeStructured3DMesh(nx=2, ny=2, nz=2)
    space = ng.H1(mesh, order=order)
    source, target = ng.GridFunction(space), ng.GridFunction(space)
    boundary = mesh.GetBoundaries()[0]
    with ng.TaskManager():
        source.Set(ng.x**order + 2*ng.y**order + 3*ng.z**order + ng.x*ng.y)
        _copy_matching_h1_trace_to_global(mesh, source, target, boundary, order=order)
        error = ng.Integrate((source-target)**2 * ng.ds(definedon=mesh.Boundaries(boundary)), mesh)
    assert np.sqrt(abs(error)) < 1e-12


def test_unsupported_trace_order_raises():
    mesh = MakeStructured3DMesh(nx=1, ny=1, nz=1)
    field = ng.GridFunction(ng.H1(mesh, order=4))
    with pytest.raises(ValueError, match='order 1, 2 or 3'):
        _copy_matching_h1_trace_to_global(mesh, field, field, mesh.GetBoundaries()[0], order=4)
