"""The star-shaped Kelvin tensors against the Jacobian of the map, on the assembly and the checking expressions.

Guards: `radia.kelvin_star`: the closed forms of J, nu', mu' and the source pullback, the coefficient function that the
solver assembles, and the identities of the boundary conditions.
Assumes: finite differences of the map as the reference for J; the limits leave room for their truncation. This is not
independent validation: `validation_test/magnetostatics/mixed_omega_star_kelvin_sphere_probe.py` carries the physics.
"""
import numpy as np
import pytest

from radia.kelvin_star import (
    RoundedDiskBoundary, SphereBoundary, star_kelvin_jacobian, star_kelvin_map, star_kelvin_mu, star_kelvin_nu, star_kelvin_physical_h,
)

CENTER = np.array([0.8, 0.1, -0.2])
BOUNDARIES = {"sphere": SphereBoundary(0.3), "disk": RoundedDiskBoundary(0.2, 0.1), "thin disk": RoundedDiskBoundary(0.25, 0.05)}


def _points(boundary, count=60, seed=3):
    """Interior points of the computational ball, kept off the rim of the disk where rho has a curvature jump."""
    rng = np.random.default_rng(seed)
    directions = rng.standard_normal((4*count, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    if isinstance(boundary, RoundedDiskBoundary):
        rim = boundary.flat_radius*np.abs(directions[:, 2]) - boundary.edge_radius*np.hypot(directions[:, 0], directions[:, 1])
        directions = directions[np.abs(rim) > 0.02][:count]
    rho, _ = boundary.radius_and_normal(directions)
    return CENTER + directions*(rho*rng.uniform(0.3, 0.95, len(directions)))[:, None]


@pytest.mark.parametrize("name", BOUNDARIES)
def test_the_map_is_an_involution_fixing_the_boundary(name):
    boundary = BOUNDARIES[name]
    points = _points(boundary)
    back = star_kelvin_map(star_kelvin_map(points, boundary, CENTER), boundary, CENTER)
    assert np.max(np.abs(back - points)) <= 1e-12, np.max(np.abs(back - points))
    directions = (points - CENTER)/np.linalg.norm(points - CENTER, axis=1, keepdims=True)
    on_boundary = CENTER + directions*boundary.radius_and_normal(directions)[0][:, None]
    assert np.max(np.abs(star_kelvin_map(on_boundary, boundary, CENTER) - on_boundary)) <= 1e-12


@pytest.mark.parametrize("name", BOUNDARIES)
def test_the_jacobian_matches_finite_differences_of_the_map(name):
    boundary = BOUNDARIES[name]
    points = _points(boundary)
    h = 1e-6
    numeric = np.stack([(star_kelvin_map(points + h*e, boundary, CENTER) - star_kelvin_map(points - h*e, boundary, CENTER))/(2*h) for e in np.eye(3)], axis=-1)
    closed = star_kelvin_jacobian(points, boundary, CENTER)
    relative = np.max(np.linalg.norm(numeric - closed, axis=(1, 2))/np.linalg.norm(closed, axis=(1, 2)))
    assert relative <= 1e-5, f"{name}: Jacobian relative error {relative:.3e}"
    assert np.all(np.linalg.det(closed) < 0.0), "the map must reverse orientation"


@pytest.mark.parametrize("name", BOUNDARIES)
def test_the_tensors_are_the_hodge_stars_of_the_jacobian_and_inverse_to_each_other(name):
    boundary = BOUNDARIES[name]
    points = _points(boundary)
    J = star_kelvin_jacobian(points, boundary, CENTER)
    det = np.abs(np.linalg.det(J))[:, None, None]
    nu, mu = star_kelvin_nu(points, boundary, CENTER, nu_0=1.0), star_kelvin_mu(points, boundary, CENTER, mu_0=1.0)
    expected_nu = np.swapaxes(J, 1, 2)@J/det
    Ji = np.linalg.inv(J)
    expected_mu = Ji@np.swapaxes(Ji, 1, 2)*det
    for actual, expected, label in ((nu, expected_nu, "nu"), (mu, expected_mu, "mu"), (mu@nu, np.broadcast_to(np.eye(3), mu.shape), "mu nu")):
        worst = np.max(np.abs(actual - expected)/np.max(np.abs(expected), axis=(1, 2), keepdims=True))
        assert worst <= 1e-12, f"{name} {label}: relative {worst:.3e}"


@pytest.mark.parametrize("name", BOUNDARIES)
def test_the_boundary_conditions_hold_for_the_pulled_back_field(name):
    boundary = BOUNDARIES[name]
    rng = np.random.default_rng(5)
    directions = (_points(boundary) - CENTER)
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    rho, normal = boundary.radius_and_normal(directions)
    points = CENTER + directions*rho[:, None]
    H = rng.standard_normal(points.shape)
    H_prime = -np.einsum("nji,nj->ni", star_kelvin_jacobian(points, boundary, CENTER), H)
    tangent = np.cross(normal, rng.standard_normal(points.shape))
    tangential = np.max(np.abs(np.sum(tangent*(H_prime + H), axis=1)))
    flux = np.einsum("ni,nij,nj->n", normal, star_kelvin_mu(points, boundary, CENTER, mu_0=1.0), H_prime) - np.sum(normal*H, axis=1)
    assert tangential <= 1e-10, f"{name}: tangential jump {tangential:.3e}"
    assert np.max(np.abs(flux)) <= 1e-10, f"{name}: normal flux jump {np.max(np.abs(flux)):.3e}"
    round_trip = star_kelvin_physical_h(H_prime, points, boundary, CENTER)
    assert np.max(np.abs(round_trip - H)) <= 1e-12, np.max(np.abs(round_trip - H))


def test_the_sphere_is_the_isotropic_factor():
    boundary = BOUNDARIES["sphere"]
    points = _points(boundary)
    r = np.linalg.norm(points - CENTER, axis=1)[:, None, None]
    assert np.max(np.abs(star_kelvin_mu(points, boundary, CENTER, mu_0=1.0) - (0.3/r)**2*np.eye(3))) <= 1e-12


def test_the_coefficient_function_is_the_checked_tensor(tmp_path):
    ng = pytest.importorskip("ngsolve")
    from netgen.occ import Box, OCCGeometry, Pnt
    from radia.kelvin_star import star_kelvin_mu_cf, star_kelvin_source_h_cf
    block = Box(Pnt(0.65, -0.05, -0.3), Pnt(0.95, 0.25, -0.1))
    block.mat("kelvin")
    mesh = ng.Mesh(OCCGeometry(block).GenerateMesh(maxh=0.2))
    boundary = BOUNDARIES["disk"]
    assembled = star_kelvin_mu_cf(mesh, boundary, CENTER)
    probes = np.array([[0.70, 0.02, -0.15], [0.90, 0.20, -0.25], [0.81, 0.12, -0.28], [0.72, 0.21, -0.18]])
    uniform = ng.CF((0.0, 0.0, 1.0))
    source = star_kelvin_source_h_cf(lambda X, Y, Z: ng.CF((X*0.0 + 0.0, Y*0.0 + 0.0, Z*0.0 + 1.0)), boundary, CENTER)
    assert uniform.dim == source.dim == 3
    for point in probes:
        value = np.array(assembled(mesh(*point))).reshape(3, 3)
        expected = star_kelvin_mu(point, boundary, CENTER)
        assert np.max(np.abs(value - expected)) <= 1e-12*np.max(np.abs(expected)), (point, value, expected)
        pulled = np.array(source(mesh(*point)))
        expected_h = -star_kelvin_jacobian(point, boundary, CENTER).T@np.array([0.0, 0.0, 1.0])
        assert np.max(np.abs(pulled - expected_h)) <= 1e-12*np.max(np.abs(expected_h)), (point, pulled, expected_h)


def test_a_relative_permeability_on_the_exterior_refuses(tmp_path):
    ng = pytest.importorskip("ngsolve")
    from netgen.occ import Box, OCCGeometry, Pnt
    from radia.kelvin_star import star_kelvin_mu_cf
    block = Box(Pnt(0.65, -0.05, -0.3), Pnt(0.95, 0.25, -0.1))
    block.mat("kelvin")
    mesh = ng.Mesh(OCCGeometry(block).GenerateMesh(maxh=0.2))
    with pytest.raises(ValueError, match="override Kelvin"):
        star_kelvin_mu_cf(mesh, BOUNDARIES["disk"], CENTER, mu_r_by_material={"kelvin": 2.0})
    with pytest.raises(ValueError, match="not mesh materials"):
        star_kelvin_mu_cf(mesh, BOUNDARIES["disk"], CENTER, mu_r_by_material={"iron": 2.0})
    with pytest.raises(ValueError, match="positive"):
        RoundedDiskBoundary(-1.0, 0.1)
