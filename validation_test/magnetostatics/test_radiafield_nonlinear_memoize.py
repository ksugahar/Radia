"""Nonlinear cached/uncached source-field equivalence; expensive validation.

The uncached reference intentionally repeats native source evaluation.
Pointwise cache and concurrency contracts live in tests/test_radiafield_memoize.py.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest
import ngsolve as ng
import radia as rad


def _coil():
    rad.UtiDelAll()
    arc = rad.ObjArcCur([-2.5, 0.0, 0.0], [.3, .6], [0.0, 1.2], .5, 40, 'man', 'z', 2.0e5)
    bar = rad.ObjRecCur([-2.2, 0.8, 0.0], [0.3, 0.8, 0.4], [0, 2.0e5, 0])
    return rad.ObjCnt([arc, bar])


@pytest.fixture(scope="module")
def picard_case():
    return _picard_case_data()


def _picard_case_data():
    spec = importlib.util.spec_from_file_location(
        "source_load_tests", Path(__file__).resolve().parents[2] / "tests" / "test_mixed_omega_source_load.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    from radia.kelvin_solver import project_source_total_hodge, project_source_interface_potential

    mesh = module._mesh(0.5)
    from netgen.meshing import Element0D
    material = mesh.GetMaterials().index("total") + 1
    vertex = next(v for e in mesh.ngmesh.Elements3D() if e.index == material
                  for v in e.vertices if mesh.ngmesh.Points()[v].p[0] > 0.)
    gauge = len(mesh.GetBBBoundaries()) + 1
    mesh.ngmesh.Add(Element0D(vertex, index=gauge))
    mesh.ngmesh.SetCD3Name(gauge, "GND")
    mesh = ng.Mesh(mesh.ngmesh)
    coil = _coil()
    source = rad.RadiaField(coil, "h")
    mu0 = 4.0e-7 * np.pi
    table = ((0.0, 0.0), (500.0, 500.0 * mu0 * 2000.0), (2.0e3, 1.4), (8.0e3, 2.0), (4.0e4, 2.4))
    with ng.TaskManager():
        hodge = project_source_total_hodge(mesh, source, ("total",), order=2)
        potential = hodge["potential"]
    return mesh, source, potential, hodge["harmonic_field"], table


@pytest.mark.parametrize("lane", ["picard_projected_p2", "newton_p2"])
def test_nonlinear_memoisation_leaves_the_solution_unchanged(picard_case, lane):
    from radia.kelvin_solver import solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin as picard
    from radia.mixed_omega_newton import solve_magnetostatic_mixed_total_reduced_omega_newton_kelvin as newton

    mesh, source, potential, harmonic, table = picard_case
    common = dict(bh_table=table, nonlinear_materials=("total",), reduced_materials=("reduced",),
                  total_materials=("total",), interface_boundary="source_total_interface",
                  dirichlet_bbbnd="GND", kelvin_mats=(), order=2,
                  total_source_h=harmonic, total_source_materials=("total",))
    solve = newton if lane.startswith("newton") else picard
    options = (dict(tolerance=1e-8, residual_tolerance=1e-10) if lane.startswith("newton")
               else dict(tolerance=1e-8, max_iterations=200, material_update_order=1, relaxation=0.3))
    # Serial, so the assembly sums are reproducible and any difference would
    # come from the memoised values; concurrency is covered by the test above.
    memoised = solve(mesh, source, potential, 1.0, (3.0, 0.0, 0.0), **options, **common)
    assert not source.memoize and source.GetCacheStats()["size"] == 0
    plain = solve.__wrapped__(mesh, source, potential, 1.0, (3.0, 0.0, 0.0), **options, **common)
    assert memoised["nonlinear_stats"]["converged"]
    for point in ((-0.5, 0.15, -0.1), (0.5, 0.1, 0.2), (0.2, -0.3, -0.1)):
        np.testing.assert_array_equal(np.asarray(memoised["B_cf"](mesh(*point))),
                                      np.asarray(plain["B_cf"](mesh(*point))))
