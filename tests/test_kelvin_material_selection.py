"""Geometry and metric must resolve the same Kelvin material set."""

import pytest
import ngsolve as ng
from netgen.occ import Box, Pnt, OCCGeometry
from radia.kelvin_material import make_kelvin_mu_cf, make_kelvin_nu_cf, MU_0


def _mesh(name):
    box = Box(Pnt(1, 1, 1), Pnt(2, 2, 2))
    box.mat(name)
    return ng.Mesh(OCCGeometry(box).GenerateMesh(maxh=1))


@pytest.mark.parametrize('key', ['kelvin', 'Kelvin', 'KELVIN'])
def test_case_insensitive_metric_and_reciprocal(key):
    mesh = _mesh('Kelvin_Exterior')
    mu = make_kelvin_mu_cf(mesh, 2, (0, 0, 0), kelvin_mats=(key,))
    nu = make_kelvin_nu_cf(mesh, 2, (0, 0, 0), kelvin_mats=(key,))
    point = mesh(1.5, 1.5, 1.5)
    assert mu(point) / MU_0 == pytest.approx(4 / 6.75)
    assert mu(point) * nu(point) == pytest.approx(1)


def test_exact_match_does_not_transform_similarly_named_physical_material():
    mesh = _mesh('Kelvin_air')
    options = dict(kelvin_mats=('Kelvin',), kelvin_match_exact=True,
                   mu_r_by_material={'Kelvin_air': 5})
    point = mesh(1.5, 1.5, 1.5)
    mu = make_kelvin_mu_cf(mesh, 2, (0, 0, 0), **options)
    nu = make_kelvin_nu_cf(mesh, 2, (0, 0, 0), **options)
    assert mu(point) / MU_0 == pytest.approx(5)
    assert mu(point) * nu(point) == pytest.approx(1)


def test_mixed_solver_rejects_metric_on_reduced_domain():
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_kelvin
    from netgen.occ import Glue
    air = Box(Pnt(-1, -1, -1), Pnt(0, 1, 1))
    iron = Box(Pnt(0, -1, -1), Pnt(1, 1, 1))
    air.mat('reduced')
    iron.mat('total')
    air.faces.name = iron.faces.name = 'source_total_interface'
    mesh = ng.Mesh(OCCGeometry(Glue([air, iron])).GenerateMesh(maxh=1))
    with pytest.raises(ValueError, match='must not match reduced'):
        solve_magnetostatic_mixed_total_reduced_omega_kelvin(
            mesh, ng.CoefficientFunction((0, 0, 1)), -ng.z, 2, (0, 0, 0),
            reduced_materials=('reduced',), total_materials=('total',),
            interface_boundary='source_total_interface', kelvin_mats=('REDUCED',))


def test_empty_selector_is_rejected():
    with pytest.raises(ValueError, match='must not be empty'):
        make_kelvin_mu_cf(_mesh('air'), 2, (0, 0, 0), kelvin_mats=('',))
