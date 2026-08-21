"""Lock the native 3D orbit tracker against a scipy DOP853 reference.

The composite field here is a pair of permanent-magnet blocks (pure C++
sources through the RadFldBatchSerial term), so the native tracker and the
scipy reference integrate the same full-3D Lorentz force with no Python
callbacks inside the released-GIL region.
"""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

import radia as rad
from radia import _radia_pybind as _native

RIGIDITY = 0.5
ENTRANCE = np.array([-0.06, 0.0, 0.0])
DIRECTION = np.array([1.0, 0.0, 0.0])
EXIT_X = 0.06
MAGNETIZATION = 9.5e5


def _mirrored_pair():
    upper = rad.ObjRecMag([0.0, 0.0, 0.03], [0.08, 0.08, 0.02],
                          [0.0, 0.0, MAGNETIZATION])
    lower = rad.ObjRecMag([0.0, 0.0, -0.03], [0.08, 0.08, 0.02],
                          [0.0, 0.0, MAGNETIZATION])
    return rad.ObjCnt([upper, lower])


def _scipy_reference(container, rigidity):
    def rhs(_s, state):
        field = np.asarray(rad.Fld(container, "b", state[:3].tolist()))
        tangent = state[3:]
        return np.r_[tangent, np.cross(tangent, field) / rigidity]

    def exit_plane(_s, state):
        return state[0] - EXIT_X

    exit_plane.terminal = True
    exit_plane.direction = 1.0
    solution = solve_ivp(
        rhs, (0.0, 0.5), np.r_[ENTRANCE, DIRECTION], method="DOP853",
        rtol=1.0e-10, atol=1.0e-12, max_step=1.0e-3, events=exit_plane)
    assert solution.success and len(solution.t_events[0]) == 1
    final = solution.y_events[0][0]
    return final[:3], final[3:] / np.linalg.norm(final[3:]), float(
        solution.t_events[0][0])


def test_native_tracker_matches_scipy_reference():
    rad.UtiDelAll()
    container = _mirrored_pair()
    positions, tangents, stations, curvature, length, oop_m, oop_t = (
        _native.track_reference_orbit_native(
            None, 0.0, int(container), False, RIGIDITY, ENTRANCE, DIRECTION,
            EXIT_X, 5.0e-4, 0.5, 1.0e-6, 33))
    ref_position, ref_tangent, ref_length = _scipy_reference(
        container, RIGIDITY)
    assert abs(length - ref_length) < 5.0e-7
    assert float(np.max(np.abs(positions[-1] - ref_position))) < 5.0e-7
    assert float(np.max(np.abs(tangents[-1] - ref_tangent))) < 5.0e-6
    # The symmetric pair keeps the measured planarity at integrator noise.
    assert oop_m < 1.0e-9
    assert oop_t < 1.0e-9
    # Midpoint curvature must collocate -B_z/(B rho) of the same field.
    middle = len(stations) // 2
    midpoint = 0.5 * (positions[middle] + positions[middle + 1])
    field = np.asarray(rad.Fld(container, "b", midpoint.tolist()))
    assert curvature[middle] == pytest.approx(-field[2] / RIGIDITY,
                                              rel=1.0e-3)
    rad.UtiDelAll()


def test_native_tracker_mirror_matches_explicit_pair():
    rad.UtiDelAll()
    pair = _mirrored_pair()
    upper_only = rad.ObjRecMag([0.0, 0.0, 0.03], [0.08, 0.08, 0.02],
                               [0.0, 0.0, MAGNETIZATION])
    explicit = _native.track_reference_orbit_native(
        None, 0.0, int(pair), False, RIGIDITY, ENTRANCE, DIRECTION,
        EXIT_X, 5.0e-4, 0.5, 1.0e-6, 33)
    # 0.5 * (B_up(r) + M B_up(Mr)) is exactly half the explicit pair field,
    # so halving the rigidity reproduces the same trajectory.
    symmetrized = _native.track_reference_orbit_native(
        None, 0.0, int(upper_only), True, 0.5 * RIGIDITY, ENTRANCE,
        DIRECTION, EXIT_X, 5.0e-4, 0.5, 1.0e-6, 33)
    assert float(np.max(np.abs(explicit[0] - symmetrized[0]))) < 1.0e-12
    assert float(np.max(np.abs(explicit[3] - symmetrized[3]))) < 1.0e-9
    rad.UtiDelAll()


def test_native_tracker_planarity_gate_trips_on_asymmetric_field():
    rad.UtiDelAll()
    upper_only = rad.ObjRecMag([0.0, 0.0, 0.03], [0.08, 0.08, 0.02],
                               [0.0, 0.0, MAGNETIZATION])
    with pytest.raises(RuntimeError, match="left the bend plane"):
        _native.track_reference_orbit_native(
            None, 0.0, int(upper_only), False, RIGIDITY, ENTRANCE,
            DIRECTION, EXIT_X, 5.0e-4, 0.5, 1.0e-6, 33)
    rad.UtiDelAll()


def test_native_general_plane_tracker_recovers_uniform_field_sector():
    rigidity = 1.7
    radius = 3.2
    angle = 0.25
    bending_field = rigidity / radius
    radial_0 = np.array([0.0, -1.0, 0.0])
    tangent_0 = np.array([1.0, 0.0, 0.0])
    rotation = np.array([
        [np.cos(angle), -np.sin(angle), 0.0],
        [np.sin(angle), np.cos(angle), 0.0],
        [0.0, 0.0, 1.0],
    ])
    radial_1 = rotation @ radial_0
    tangent_1 = rotation @ tangent_0
    positions, tangents, stations, curvature, length, oop_m, oop_t = (
        _native.track_reference_orbit_to_plane_native(
            None, 0.0, -1, False,
            np.array([0.0, 0.0, bending_field]), -rigidity,
            radius * radial_0, tangent_0, tangent_1, 0.0,
            5.0e-4, 2.0, 1.0e-8, 33))

    np.testing.assert_allclose(
        positions[-1], radius * radial_1, atol=2.0e-9)
    np.testing.assert_allclose(tangents[-1], tangent_1, atol=2.0e-9)
    np.testing.assert_allclose(length, radius * angle, atol=2.0e-9)
    np.testing.assert_allclose(stations, np.linspace(0.0, length, 33))
    np.testing.assert_allclose(curvature, 1.0 / radius, atol=2.0e-12)
    assert oop_m == 0.0 and oop_t == 0.0


def test_hdiv_field_evaluator_point_gradient_is_analytic():
    evaluator = _native._HDivFieldEvaluator.from_cloud(
        np.array([[0.0, 0.0, 0.0]]), np.array([2.0]))
    point = np.array([[1.0, 2.0, 3.0]])
    displacement = point[0]
    radius2 = float(displacement@displacement)
    expected = 2.0*(
        np.eye(3)/radius2**1.5
        - 3.0*np.outer(displacement, displacement)/radius2**2.5)
    np.testing.assert_allclose(
        evaluator.field_gradient(point)[0], expected,
        rtol=2.0e-15, atol=2.0e-17)
