"""Connected scalar Kelvin gauges; isolated CAD points are not FE constraints."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest
from ngsolve import Mesh, H1, BilinearForm, LinearForm, GridFunction, grad, dx, TaskManager
from netgen.occ import Sphere, Pnt, OCCGeometry

# Test the geometry module independently of the optional Radia native backend.
_module = Path(__file__).resolve().parents[2] / 'src/radia/kelvin_geometry.py'
spec = importlib.util.spec_from_file_location('kelvin_ground_geometry', _module)
geometry = importlib.util.module_from_spec(spec)
spec.loader.exec_module(geometry)


def _mesh():
    shape = Sphere(Pnt(0, 0, 0), 1)
    shape.name = 'kelvin'
    return Mesh(OCCGeometry(shape).GenerateMesh(maxh=0.7))


def test_missing_ground_rejected_and_explicit_ground_constrains_volume():
    mesh = _mesh()
    with pytest.raises(ValueError, match='connected vertex'):
        geometry.require_kelvin_scalar_ground(mesh, 'kelvin')
    mesh = geometry.ground_kelvin_scalar_mesh(mesh, (0, 0, 0))
    ground = geometry.require_kelvin_scalar_ground(mesh, 'kelvin')
    assert len(ground) == 1
    assert geometry.ground_kelvin_scalar_mesh(mesh, (0, 0, 0)) is mesh
    free = H1(mesh).FreeDofs()
    constrained = H1(mesh, dirichlet_bbbnd='GND').FreeDofs()
    assert sum(free) - sum(constrained) == 1


def test_disconnected_point_label_rejected():
    from netgen.meshing import MeshPoint, Element0D, Pnt as MeshPnt
    mesh = _mesh()
    point = mesh.ngmesh.Add(MeshPoint(MeshPnt(0, 0, 0)))
    mesh.ngmesh.Add(Element0D(point, index=1))
    mesh.ngmesh.SetCD3Name(1, 'GND')
    mesh = Mesh(mesh.ngmesh)
    with pytest.raises(ValueError, match='isolated OCC'):
        geometry.ground_kelvin_scalar_mesh(mesh, (0, 0, 0))


def test_scalar_gradient_independent_of_gauge_vertex():
    # Independent manufactured field u=x. The weak load has zero net source.
    # Anchoring two different existing vertices must give the same gradient.
    from ngsolve import CoefficientFunction, Integrate
    mesh = _mesh()
    mesh = geometry.ground_kelvin_scalar_mesh(mesh, (0, 0, 0), name='center_gauge')
    mesh = geometry.ground_kelvin_scalar_mesh(mesh, (1, 0, 0), name='edge_gauge')
    fields = []
    with TaskManager():
        for name in ['center_gauge', 'edge_gauge']:
            fes = H1(mesh, order=1, dirichlet_bbbnd=name)
            u, v = fes.TnT()
            a = BilinearForm(fes, symmetric=True)
            a += grad(u)*grad(v)*dx
            f = LinearForm(fes)
            f += CoefficientFunction((1, 0, 0))*grad(v)*dx
            a.Assemble(); f.Assemble()
            g = GridFunction(fes)
            g.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse='sparsecholesky') * f.vec
            fields.append(g)
        error = Integrate((grad(fields[0])-grad(fields[1]))**2, mesh)
        exact = Integrate((grad(fields[0])-CoefficientFunction((1, 0, 0)))**2, mesh)
    assert error < 1e-20
    assert exact < 1e-20


def test_geometry_builder_does_not_promise_an_isolated_ground():
    inner = Sphere(Pnt(0, 0, 0), 1)
    inner.name = 'air'
    inner.faces.name = 'kelvin_int'
    shape, info = geometry.add_kelvin_exterior_domain(inner, (3, 0, 0), 1)
    assert info['gnd_vertex'] is None
    mesh = Mesh(OCCGeometry(shape).GenerateMesh(maxh=0.8))
    with pytest.raises(ValueError, match='connected vertex'):
        geometry.require_kelvin_scalar_ground(mesh, 'kelvin')
    mesh = geometry.ground_kelvin_scalar_mesh(mesh, info['offset'])
    from ngsolve import Periodic
    ungrounded = Periodic(H1(mesh, order=1))
    grounded = Periodic(H1(mesh, order=1, dirichlet_bbbnd='GND'))
    assert sum(ungrounded.FreeDofs()) - sum(grounded.FreeDofs()) == 1
