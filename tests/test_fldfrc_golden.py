"""Golden tests for the C++ Maxwell-stress force API (rad.FldFrc).

``rad.FldFrc(obj, shape)`` integrates the magnetic Maxwell stress over the
rectangle created by ``rad.FldFrcShpRtg(center, [Lx, Ly])``.  The rectangle
lies in the plane z = center[2] and carries the normal +z, so the returned
three-vector is

    F = oint_S (1/mu0) [ (B.n) B - 0.5 |B|^2 n ] dS ,    n = +z ,

which is the force on the material on the -z side of the rectangle (close the
surface at infinity below the plane).  The API returns force ONLY -- three
components, no torque.

These tests lock the SI normalisation.  Before 2026-09-10 the integrand used
1/mu0 as its prefactor while the field handed to it is H in A/m, so the whole
API returned (1/mu0^2) ~ 6.33e11 times the correct force with no error raised
and no test covering it.
"""

import math

import numpy as np
import pytest

import radia as rad
from radia.force import integrate_maxwell_surface_force

MU0 = 4.0e-7 * math.pi


@pytest.fixture(autouse=True)
def _clean_radia():
    rad.UtiDelAll()
    yield
    rad.UtiDelAll()


def test_uniform_background_matches_analytic_maxwell_pressure():
    """A uniform B0 z-hat over a flat rectangle gives exactly B0^2 A / (2 mu0).

    The field is constant, so the surface quadrature is exact and this pins the
    SI prefactor to seven digits rather than to an integration tolerance.
    """
    b0 = 1.0
    side = 0.1
    obj = rad.ObjBckg(lambda p: [0.0, 0.0, b0])
    shape = rad.FldFrcShpRtg([0.0, 0.0, 0.0], [side, side])

    force = np.asarray(rad.FldFrc(obj, shape))

    expected_fz = b0 * b0 * side * side / (2.0 * MU0)  # 3978.8736 N
    assert force.shape == (3,)
    assert force[2] == pytest.approx(expected_fz, rel=1e-6)
    assert abs(force[0]) < 1e-6 * expected_fz
    assert abs(force[1]) < 1e-6 * expected_fz


def test_fldfrc_returns_force_only_three_components():
    """The API returns [Fx, Fy, Fz] -- Torque_ is never set by ComputeFieldForce."""
    obj = rad.ObjBckg(lambda p: [0.0, 0.0, 0.5])
    shape = rad.FldFrcShpRtg([0.0, 0.0, 0.0], [0.05, 0.05])
    force = np.asarray(rad.FldFrc(obj, shape))
    assert force.shape == (3,), "FldFrc returns force only, not force+torque"


def _box_surface_quadrature(lower, upper, order):
    """Gauss-Legendre tensor quadrature on the six faces of an axis-aligned box."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    points, normals, areas = [], [], []
    for axis in range(3):
        first, second = [k for k in range(3) if k != axis]
        c1 = 0.5 * (upper[first] + lower[first])
        h1 = 0.5 * (upper[first] - lower[first])
        c2 = 0.5 * (upper[second] + lower[second])
        h2 = 0.5 * (upper[second] - lower[second])
        for side, plane in ((+1, upper[axis]), (-1, lower[axis])):
            for g1, w1 in zip(nodes, weights):
                for g2, w2 in zip(nodes, weights):
                    point = np.zeros(3)
                    point[axis] = plane
                    point[first] = c1 + h1 * g1
                    point[second] = c2 + h2 * g2
                    normal = np.zeros(3)
                    normal[axis] = side
                    points.append(point)
                    normals.append(normal)
                    areas.append(h1 * h2 * w1 * w2)
    return np.array(points), np.array(normals), np.array(areas)


def test_magnet_pair_matches_independent_box_surface_integral():
    """Cross-engine check: C++ FldFrc vs the numpy radia.force kernel.

    Two coaxial cuboid magnets magnetised +z attract.  A closed box drawn in the
    air around the UPPER magnet gives the force on it; the mid-plane rectangle
    handed to FldFrc carries n = +z and therefore reports the force on the LOWER
    magnet, which is the same magnitude with the opposite sign.
    """
    magnetization = [0.0, 0.0, 954929.6585]  # A/m  (Br = 1.2 T)
    sides = [0.02, 0.02, 0.01]
    gap = 0.01
    offset = 0.5 * (sides[2] + gap)

    lower = rad.ObjRecMag([0.0, 0.0, -offset], sides, magnetization)
    upper = rad.ObjRecMag([0.0, 0.0, +offset], sides, magnetization)
    system = rad.ObjCnt([lower, upper])

    # Independent reference: closed box in the air around the upper magnet.
    points, normals, areas = _box_surface_quadrature(
        (-0.03, -0.03, 0.0), (0.03, 0.03, 0.03), 48
    )
    field = np.asarray(rad.Fld(system, "b", points.tolist())).reshape(-1, 3)
    reference = integrate_maxwell_surface_force(field, normals, areas)

    assert reference[2] < 0.0, "coaxial like-magnetised magnets must attract"
    assert reference[2] == pytest.approx(-21.9079, rel=2e-4)
    assert np.hypot(reference[0], reference[1]) < 1e-9 * abs(reference[2])

    shape = rad.FldFrcShpRtg([0.0, 0.0, 0.0], [0.2, 0.2])
    force = np.asarray(rad.FldFrc(system, shape))

    # Force on the -z side of the rectangle = force on the lower magnet.
    assert force[2] == pytest.approx(-reference[2], rel=5e-3)
