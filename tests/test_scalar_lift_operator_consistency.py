"""Discrete lifting loads use the same curved-element operator as the solve.

These are assembly contracts with prescribed jumps, not Maxwell accuracy tests.
"""
import numpy as np
import pytest
import ngsolve as ng
from netgen.occ import Sphere, Pnt, Glue, OCCGeometry
from netgen.meshing import IdentificationType
from radia.kelvin_solver import (
    solve_magnetostatic_matching_trace_total_reduced_omega as matching,
    solve_magnetostatic_mixed_total_reduced_omega_kelvin as mixed,
)


def curved_mesh(kelvin=False):
    iron = Sphere(Pnt(0, 0, 0), .2)
    iron.mat('iron'); iron.faces.name = 'interface'; iron.maxh = .09
    outer = Sphere(Pnt(0, 0, 0), .5); outer.faces.name = 'kelvin_int' if kelvin else 'outer'
    air = outer - iron; air.mat('air')
    if not kelvin:
        mesh = ng.Mesh(OCCGeometry(Glue([iron, air])).GenerateMesh(maxh=.18, grading=.6))
    else:
        exterior = Sphere(Pnt(1.5, 0, 0), .5)
        exterior.mat('kelvin'); exterior.faces.name = 'kelvin_ext'
        for vertex in exterior.vertices:
            if vertex.p[2] > 0: vertex.name = 'GND'
        geo = Glue([iron, air, exterior])
        faces = {face.name: face for solid in geo.solids for face in solid.faces if face.name}
        faces['kelvin_int'].Identify(faces['kelvin_ext'], 'periodic', IdentificationType.PERIODIC)
        mesh = ng.Mesh(OCCGeometry(geo).GenerateMesh(maxh=.18, grading=.6))
    mesh.Curve(3)
    return mesh


def operator_action(result, lift, material):
    space = lift.space
    u, v = space.TnT()
    operator = ng.BilinearForm(space, symmetric=True)
    operator += result['mu_cf'] * ng.grad(u) * ng.grad(v) * ng.dx(
        definedon=space.mesh.Materials(material), bonus_intorder=10)
    operator.Assemble()
    action = lift.vec.CreateVector(); action.data = operator.mat * lift.vec
    return np.array(action.FV().NumPy(), copy=True)


@pytest.mark.parametrize('order', [2, 3])
def test_matching_reduced_lift_equals_original_operator_action(order):
    with ng.TaskManager():
        mesh = curved_mesh()
        source = ng.GridFunction(ng.H1(mesh, order=order))
        source.Set(-ng.z + .04*(ng.x**2-ng.y**2))
        result = matching(mesh, ng.CF((0, 0, 0)), source, order=order,
            mu_r_by_material={'air': 1., 'iron': 100.},
            reduced_materials=('air',), total_materials=('iron',),
            interface_boundary='interface', dirichlet_boundary='outer',
            bonus_intorder=10, return_system=True)
        expected = operator_action(result, result['source_lift'], 'air')
    actual = result['system']['linear_form'].vec.FV().NumPy()
    free = np.array(list(result['fes'].FreeDofs()), dtype=bool)
    assert np.linalg.norm(expected[free]) > 1e-12
    assert np.linalg.norm((actual-expected)[free])/np.linalg.norm(expected[free]) < 1e-12


@pytest.mark.parametrize('order', [2, 3])
def test_kelvin_lift_equals_original_operator_action(order):
    with ng.TaskManager():
        mesh = curved_mesh(kelvin=True)
        result = mixed(mesh, ng.CF((0, 0, 0)), ng.CF(0), .5, (1.5, 0, 0),
            mu_r_by_material={'air': 1., 'iron': 100.},
            reduced_materials=('air',), total_materials=('iron', 'kelvin'),
            interface_boundary='interface', kelvin_interface_boundary='kelvin_int',
            kelvin_source_potential=-ng.z + .04*(ng.x**2-ng.y**2),
            order=order, bonus_intorder=10, return_system=True)
        expected = -operator_action(result, result['kelvin_source_lift'], 'kelvin')
    actual = result['system']['linear_form'].vec.FV().NumPy()[:len(expected)]
    free = np.array(list(result['fes_reduced'].FreeDofs()), dtype=bool)
    assert np.linalg.norm(expected[free]) > 1e-12
    assert np.linalg.norm((actual-expected)[free])/np.linalg.norm(expected[free]) < 1e-12
