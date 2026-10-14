"""Planar optics routes share one right-handed curvilinear convention.

The local frame is ``x = bend_axis x tangent``, ``y = bend_axis``,
``s = tangent``.  Positive signed curvature bends the design orbit toward
``-x`` (metric ``1 + h*x``), and ``curvature_sign`` is the charge sign in
``h = q*B_y/p`` and ``k1 = q*(dB_y/dx)/p``.  Every map below is compared with
Cartesian Lorentz tracking of the same analytic field, projected on
``PlanarDesignOrbit.frame_at``, so an internally consistent but mirrored
transverse axis cannot pass.
"""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from radia.accelerator_lie_topopt import (
    third_order_lie_map_from_multipoles,
    track_canonical_hamiltonian_s,
)
from radia.accelerator_magnet_topopt import (
    PlanarDesignOrbit,
    PlanarTransferMatrixObjective,
    planar_orbit_field_observations,
)
from radia.accelerator_taylor_topopt import second_order_taylor_map_from_multipoles

RADIUS = 3.2
RIGIDITY = 1.7
ANGLE = 0.6
GRADIENT = 0.4
STATIONS = 25


def _orbit(axis_sign):
    """Counter-clockwise circle about +y; ``bend_axis`` selects the sign of h."""
    angles = np.linspace(0.0, ANGLE, STATIONS)
    return PlanarDesignOrbit(
        positions=np.column_stack(
            (RADIUS * np.sin(angles), 0.0 * angles, RADIUS * np.cos(angles))
        ),
        tangents=np.column_stack(
            (np.cos(angles), 0.0 * angles, -np.sin(angles))
        ),
        magnetic_rigidity=RIGIDITY,
        bend_axis=np.array([0.0, axis_sign, 0.0]),
        path_length_stations=RADIUS * angles,
    )


class _SectorField:
    """``B_y = B0 + g*x`` and ``B_x = g*y`` in the orbit's curvilinear frame."""

    def __init__(self, orbit, charge, gradient):
        self.axis = orbit.bend_axis
        self.h = float(orbit.signed_curvature[0])
        self.b0 = self.h * RIGIDITY / charge
        self.gradient = gradient

    def __call__(self, point):
        radial = point.copy()
        radial[1] = 0.0
        distance = np.linalg.norm(radial)
        outward = radial / distance
        # h > 0 bends toward -x, so +x points away from the centre.
        e_x = np.sign(self.h) * outward
        x = np.sign(self.h) * (distance - RADIUS)
        y = float(point @ self.axis)
        return (self.b0 + self.gradient * x) * self.axis + self.gradient * y * e_x


def _exact_exit(orbit, field, charge, state):
    """Track ``(x, px, y, py, delta)`` to the design exit plane."""
    x, px, y, py, delta = state
    h0, v0, t0 = (value[0] for value in orbit.frame_at(np.array([0.0])))
    s_end = orbit.arc_length_stations[-1]
    origin = orbit.position_at(s_end)
    h1, v1, t1 = orbit.frame_at(s_end)
    momentum = 1.0 + delta
    direction = (
        px * h0 + py * v0 + np.sqrt(momentum**2 - px**2 - py**2) * t0
    ) / momentum
    scale = charge / (RIGIDITY * momentum)

    def rhs(_path, value):
        return np.concatenate(
            (value[3:], scale * np.cross(value[3:], field(value[:3])))
        )

    def exit_plane(_path, value):
        return float((value[:3] - origin) @ t1)

    exit_plane.terminal = True
    exit_plane.direction = 1.0
    start = orbit.position_at(0.0) + x * h0 + y * v0
    solution = solve_ivp(
        rhs,
        (0.0, 2.0 * s_end),
        np.concatenate((start, direction)),
        events=exit_plane,
        method="DOP853",
        rtol=1.0e-13,
        atol=1.0e-16,
    )
    final = solution.y_events[0][0]
    offset = final[:3] - origin
    return np.array(
        [
            offset @ h1,
            momentum * (final[3:] @ h1),
            offset @ v1,
            momentum * (final[3:] @ v1),
            delta,
        ]
    ), final[:3]


def _exact_linear_map(orbit, field, charge):
    """Central-difference ``R`` in ``(x, px, y, py, delta)``."""
    step = 1.0e-6
    matrix = np.zeros((4, 5))
    for column in range(5):
        seed = np.zeros(5)
        seed[column] = step
        plus, _ = _exact_exit(orbit, field, charge, seed)
        minus, _ = _exact_exit(orbit, field, charge, -seed)
        matrix[:, column] = (plus[:4] - minus[:4]) / (2.0 * step)
    return matrix


def _multipole_response(orbit, field, gradient):
    coefficients = np.zeros((7, len(orbit.segment_lengths)))
    coefficients[0] = field.b0
    coefficients[1] = gradient
    return coefficients.reshape(-1)


CASES = [
    pytest.param(+1.0, +1.0, id="h-positive-proton"),
    pytest.param(-1.0, +1.0, id="h-negative-proton"),
    pytest.param(+1.0, -1.0, id="h-positive-negative-charge"),
    pytest.param(-1.0, -1.0, id="h-negative-negative-charge"),
]
PHASE_SPACE = [0, 1, 2, 3, 5]


def test_design_orbit_frame_is_right_handed_and_bends_toward_negative_x():
    for axis_sign in (+1.0, -1.0):
        orbit = _orbit(axis_sign)
        horizontal, vertical, tangent = orbit.frame_at(orbit.arc_length_stations)
        np.testing.assert_allclose(
            np.einsum("ij,ij->i", np.cross(horizontal, vertical), tangent),
            1.0,
            atol=2.0e-15,
        )
        # The orbit turns counter-clockwise about +y; with bend_axis=+y the
        # centre therefore lies on +x and h is negative.
        np.testing.assert_allclose(
            orbit.signed_curvature, -axis_sign / RADIUS, rtol=2.0e-12
        )
        centre_direction = -orbit.positions[0] / RADIUS
        assert np.sign(centre_direction @ horizontal[0]) == -np.sign(
            orbit.signed_curvature[0]
        )


@pytest.mark.parametrize(("axis_sign", "charge"), CASES)
def test_linear_combined_function_map_matches_cartesian_tracking(axis_sign, charge):
    orbit = _orbit(axis_sign)
    field = _SectorField(orbit, charge, GRADIENT)
    points, weights = planar_orbit_field_observations(orbit, gradient_offset=1.0e-4)
    rows = np.einsum(
        "rpk,pk->r", weights, np.asarray([field(point) for point in points])
    )
    objective = PlanarTransferMatrixObjective(
        orbit, np.eye(6), 1.0, 1.0, curvature_sign=charge
    )
    # Hermite midpoints sit ~1e-9 relative off the circle, inside the
    # gradient; the sign is what this contract protects.
    np.testing.assert_allclose(
        objective.required_bend_field, rows[: len(orbit.segment_lengths)],
        rtol=1.0e-8,
    )
    model = objective.evaluate_transfer_map(rows).matrix
    exact = _exact_linear_map(orbit, field, charge)
    np.testing.assert_allclose(
        model[np.ix_([0, 1, 2, 3], PHASE_SPACE)], exact, rtol=0.0, atol=2.0e-6
    )


@pytest.mark.parametrize(("axis_sign", "charge"), CASES)
def test_taylor_and_lie_multipole_maps_match_cartesian_tracking(axis_sign, charge):
    orbit = _orbit(axis_sign)
    field = _SectorField(orbit, charge, GRADIENT)
    response = _multipole_response(orbit, field, GRADIENT)
    exact = _exact_linear_map(orbit, field, charge)
    taylor = second_order_taylor_map_from_multipoles(
        response[: 5 * len(orbit.segment_lengths)],
        orbit.segment_lengths,
        RIGIDITY,
        curvature_sign=charge,
        maximum_step_m=0.004,
    )
    lie = third_order_lie_map_from_multipoles(
        response,
        orbit.segment_lengths,
        RIGIDITY,
        curvature_sign=charge,
        reference_curvature_per_m=orbit.signed_curvature,
        maximum_step_m=0.02,
    )
    for matrix in (taylor.R, lie.R):
        np.testing.assert_allclose(
            matrix[np.ix_([0, 1, 2, 3], PHASE_SPACE)], exact, rtol=0.0, atol=2.0e-6
        )


@pytest.mark.parametrize(("axis_sign", "charge"), CASES)
def test_canonical_s_tracking_reconstructs_the_cartesian_trajectory(axis_sign, charge):
    orbit = _orbit(axis_sign)
    field = _SectorField(orbit, charge, 0.0)
    initial = np.array([0.05, 0.004, 0.0, 0.0, 0.0, 0.002])
    tracked = track_canonical_hamiltonian_s(
        orbit,
        _multipole_response(orbit, field, 0.0),
        initial,
        curvature_sign=charge,
        integrator="DOP853",
        maximum_step_m=0.004,
    )
    exact, exit_point = _exact_exit(orbit, field, charge, initial[PHASE_SPACE])
    np.testing.assert_allclose(
        tracked.final_state[PHASE_SPACE], exact, rtol=0.0, atol=2.0e-10
    )
    np.testing.assert_allclose(
        tracked.global_positions_m[-1], exit_point, rtol=0.0, atol=2.0e-10
    )
