"""Physical energy diagnostics use shared volume weights, not gap samples."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")
spec = importlib.util.spec_from_file_location(
    "energy_validation", Path(__file__).parents[1] / "validation_test" /
    "c_type_three_engine" / "magnetic_energy.py")
energy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(energy)


def test_linear_energy_and_coenergy_include_volume_and_direction():
    h = np.array([[1., 2., 3.], [2., 0., 0.]])
    mu = energy.MU0 * 7
    weights = np.array([2., 3.])
    r = energy.integrate_energy(mu*h, h, weights, ["iron"]*2, {"iron": 7.})
    expected = .5 * mu * (2*14 + 3*4)
    assert r["energy_J"] == pytest.approx(expected)
    assert r["coenergy_J"] == pytest.approx(expected)
    assert abs(r["fenchel_gap_J"]) < 1e-18
    inconsistent = energy.integrate_energy(-mu*h, h, weights, ["iron"]*2, {"iron": 7.})
    assert inconsistent["fenchel_gap_J"] == pytest.approx(4*expected)


def test_nonlinear_energy_matches_independent_integral_and_saturation_tail():
    from scipy.integrate import quad
    from radia.scalar_potential_solver import _build_bh_interpolator
    table = np.array([[0., 0.], [100., .5], [400., 1.1], [1000., 1.5]])
    law = _build_bh_interpolator(table)
    for h in [0., 55., 700., 1300.]:
        b = law(h)
        r = energy.integrate_energy([[b, 0, 0]], [[h, 0, 0]], [.03], ["iron"], {"iron": table})
        expected = .03 * quad(law, 0, h, points=[x for x in table[:, 0] if 0 < x < h], epsabs=1e-10)[0]
        assert r["coenergy_J"] == pytest.approx(expected, abs=1e-10)
        assert r["energy_J"] == pytest.approx(.03*h*b - expected, abs=1e-10)
        assert abs(r["fenchel_gap_J"]) < 1e-9
    assert r["energy_J"] != pytest.approx(r["coenergy_J"], rel=.01)


@pytest.mark.parametrize("weights, law", [([-1.], 1.), ([float('nan')], 1.), ([1.], [[1., .1], [2., .2]])])
def test_invalid_volume_or_remanent_law_is_rejected(weights, law):
    with pytest.raises(ValueError):
        energy.integrate_energy([[1., 0, 0]], [[1., 0, 0]], weights, ["iron"], {"iron": law})


def test_common_mesh_quadrature_agrees_with_ngsolve_integral():
    from netgen.occ import unit_cube
    netmesh = unit_cube.GenerateMesh(maxh=.5)
    netmesh.SetMaterial(1, "iron")
    mesh = ng.Mesh(netmesh)
    table = [[0., 0.], [10., 10.*energy.MU0*5]]
    observer = energy.PhysicalVolumeEnergy(mesh, table, 4)
    H = ng.CF((ng.x, 2*ng.y, 3*ng.z))
    B = energy.MU0*5*H
    r = observer.fem(mesh, B, H)
    expected = .5 * energy.MU0*5 * 14/3
    assert r["volume_m3"] == pytest.approx(1.)
    assert r["energy_J"] == pytest.approx(expected, rel=1e-10)
    assert r["coenergy_J"] == pytest.approx(expected, rel=1e-10)
    assert r["contract"]["all_space_energy"] is False


def test_scalar_ordering_requires_the_same_domain_contract():
    names = ["hdiv_mmm", "reduced_a", "mixed_total_reduced_omega"]
    records = {name: {"contract": {"quadrature_sha256": "shared"}, "energy_J": value, "coenergy_J": value}
               for name, value in zip(names, [2., 3., 1.])}
    assert energy.compare_energy(records)["energy_J"]["hdiv_between"]
    records[names[0]]["contract"] = {"quadrature_sha256": "different"}
    with pytest.raises(ValueError, match="contracts differ"):
        energy.compare_energy(records)


def test_kelvin_exterior_dipole_energy_and_field_pullback():
    from netgen.occ import OCCGeometry, Sphere, Pnt
    from radia.kelvin_source import kelvin_solution_to_computational
    shape = Sphere(Pnt(0, 0, 0), 1).mat("kelvin")
    mesh = ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=.3))
    mesh.Curve(2)
    observer = energy.PhysicalVolumeEnergy(
        mesh, [[0., 0.], [1., energy.MU0]], 6,
        kelvin_center=(0., 0., 0.), kelvin_radius=1.)
    radii = np.linalg.norm(observer.points, axis=1)
    n = observer.points / radii[:, None]
    H = (3*n*n[:, 2, None] - np.array([0., 0., 1.])) / (4*np.pi*radii[:, None]**3)
    B = energy.MU0*H
    direct = observer._record(B, H)
    # Integral of the dipole's vacuum energy outside radius R=1, m=1.
    assert direct["energy_J"] == pytest.approx(energy.MU0/(12*np.pi), rel=.005)
    def computational(values, form):
        value = kelvin_solution_to_computational(
            values, observer.computational_points, kelvin_center=(0., 0., 0.),
            physical_center=(0., 0., 0.), radius=1., form=form)
        return lambda mapped: value
    transformed = observer.omega(mesh, computational(B, "flux_density"), computational(H, "field_strength"))
    assert transformed["energy_J"] == pytest.approx(direct["energy_J"], rel=1e-12)
    assert transformed["coenergy_J"] == pytest.approx(direct["energy_J"], rel=1e-12)
    assert transformed["contract"]["scope"] == "kelvin_exterior_only"
    assert transformed["contract"]["all_space_energy"] is False


@pytest.mark.parametrize("method", ["linear", "picard"])
def test_reduced_a_explicit_sparsecholesky_matches_direct_field(method):
    from netgen.occ import unit_cube
    from radia.vector_potential_solver import VectorPotentialSolver
    netmesh = unit_cube.GenerateMesh(maxh=.5)
    netmesh.SetMaterial(1, "iron")
    mesh = ng.Mesh(netmesh)
    fields = []
    for backend in ("sparsecholesky", "direct"):
        solver = VectorPotentialSolver(mesh, iron_domains="iron", mu_r=10, order=1)
        solver.set_source_cf(ng.CF((0, 0, ng.x)))
        if method == "linear":
            solver.solve_linear(dirichlet=".*", solver=backend, eps=1e-3)
            stats = solver._last_linear_stats
            assert stats["relative_residual"] < 1e-6
        else:
            solver.solve_nonlinear(np.array([[0., 0.], [1e6, energy.MU0*10*1e6]]),
                                   dirichlet=".*", solver=backend, eps=1e-3,
                                   relax=1., maxiter=4, tol=1e-8, verbose=False)
            stats = solver._last_nonlinear_stats
            assert stats["converged"]
        if backend == "sparsecholesky":
            assert stats["direct_inverse"] == backend
        fields.append(np.array(solver.get_B()(mesh(.3, .4, .5))))
    assert np.linalg.norm(fields[0] - fields[1]) < 1e-8

