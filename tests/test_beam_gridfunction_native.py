from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

ngsolve = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

from netgen.occ import Box, Pnt
from ngsolve import CF, GridFunction, Mesh, TaskManager, VectorH1, x, y, z
from scipy.linalg import expm

from radia.beam import (
    CartesianState,
    ClassicalRK4,
    GridFunctionField,
    LorentzEquation,
    ParticleSpecies,
    Tracker,
    TrackPlan,
    propagate_grid_function_linear_map,
    propagate_grid_function_multipole_map,
)
from radia.accelerator_lie_topopt import fourth_order_lie_map_from_tracked_orbit
from radia.accelerator_magnet_topopt import PlanarDesignOrbit


@pytest.fixture(scope="module")
def affine_combined_function_field():
    geometry = Box(Pnt(-0.1, -0.1, 0.0), Pnt(0.1, 0.1, 1.0))
    mesh = geometry.GenerateMesh(maxh=0.12)
    space = VectorH1(mesh, order=1)
    field = GridFunction(space)
    bend_t = 1.2
    normal_gradient_t_per_m = 2.4
    skew_gradient_t_per_m = -0.6
    coefficient = CF(
        (
            normal_gradient_t_per_m * y + skew_gradient_t_per_m * x,
            bend_t + normal_gradient_t_per_m * x
            - skew_gradient_t_per_m * y,
            0.0,
        )
    )
    with TaskManager():
        field.Set(coefficient)
    return field, bend_t, normal_gradient_t_per_m, skew_gradient_t_per_m


@pytest.fixture(scope="module")
def cubic_multipole_field():
    geometry = Box(Pnt(-0.05, -0.05, 0.0), Pnt(0.05, 0.05, 1.0))
    mesh = geometry.GenerateMesh(maxh=0.1)
    field = GridFunction(VectorH1(mesh, order=3))
    normal = np.array([0.0, 2.4, 15.0, -80.0])
    skew = np.array([0.0, -0.6, -4.0, 20.0])
    by = (
        normal[1] * x
        - skew[1] * y
        + normal[2] * (x * x - y * y)
        - 2 * skew[2] * x * y
        + normal[3] * (x**3 - 3 * x * y * y)
        - skew[3] * (3 * x * x * y - y**3)
    )
    bx = (
        skew[1] * x
        + normal[1] * y
        + skew[2] * (x * x - y * y)
        + 2 * normal[2] * x * y
        + skew[3] * (x**3 - 3 * x * y * y)
        + normal[3] * (3 * x * x * y - y**3)
    )
    with TaskManager():
        field.Set(CF((bx, by, 0.0)))
    return field, normal, skew


def test_grid_function_field_is_linearized_and_accumulated_in_cpp(
    affine_combined_function_field,
):
    field, bend_t, normal_gradient, skew_gradient = (
        affine_combined_function_field
    )
    lengths = np.array([0.2, 0.3, 0.25])
    positions = np.column_stack(
        (np.zeros(3), np.zeros(3), np.array([0.15, 0.45, 0.8]))
    )
    tangents = np.tile([0.0, 0.0, 1.0], (3, 1))
    rigidity = 3.0

    result = propagate_grid_function_linear_map(
        field,
        lengths,
        positions,
        tangents,
        rigidity,
        sample_radius_m=0.01,
        names=["entrance", "body", "exit"],
        maximum_step_m=0.01,
    )

    curvature = bend_t / rigidity
    k1 = normal_gradient / rigidity
    k1s = skew_gradient / rigidity
    expected_a = np.zeros((6, 6))
    expected_a[0, 1] = 1.0
    expected_a[1, 0] = -(curvature**2 + k1)
    expected_a[1, 2] = k1s
    expected_a[1, 5] = curvature
    expected_a[2, 3] = 1.0
    expected_a[3, 0] = k1s
    expected_a[3, 2] = k1
    expected_a[4, 0] = curvature

    assert result["schema"] == "radia.beam.grid-function-linear-map.result.v1"
    assert result["backend"] == "native-cpp-ngsolve-gridfunction"
    assert result["maximum_order"] == 1
    np.testing.assert_allclose(result["frame_horizontal"], [[1, 0, 0]] * 3)
    np.testing.assert_allclose(result["frame_vertical"], [[0, 1, 0]] * 3)
    np.testing.assert_allclose(result["curvature_per_m"], curvature, atol=1e-12)
    np.testing.assert_allclose(result["normal_gradient_per_m2"], k1, atol=1e-11)
    np.testing.assert_allclose(result["skew_gradient_per_m2"], k1s, atol=1e-11)
    np.testing.assert_allclose(
        result["local_A_per_m"], np.tile(expected_a, (3, 1, 1)), atol=1e-11
    )
    np.testing.assert_allclose(
        result["R"], expm(expected_a * lengths.sum()), atol=2e-11
    )
    np.testing.assert_allclose(result["T"], 0.0, atol=1e-14)
    np.testing.assert_allclose(result["U"], 0.0, atol=1e-14)
    np.testing.assert_allclose(
        result["transverse_divergence_t_per_m"], 0.0, atol=1e-11
    )
    np.testing.assert_allclose(
        result["transverse_curl_mismatch_t_per_m"], 0.0, atol=1e-11
    )
    np.testing.assert_array_equal(result["fit_rank"], 3)
    np.testing.assert_allclose(result["scaled_design_condition"], 1.5)
    np.testing.assert_allclose(result["center_fit_bias_t"], 0.0, atol=1e-12)
    assert np.max(result["maximum_fit_residual_t"]) < 1e-11


def test_grid_function_sampling_fails_loudly_outside_mesh(
    affine_combined_function_field,
):
    field = affine_combined_function_field[0]
    with pytest.raises(RuntimeError, match="outside the volume mesh"):
        propagate_grid_function_linear_map(
            field,
            [0.1],
            [[0.0, 0.0, 1.2]],
            [[0.0, 0.0, 1.0]],
            3.0,
            sample_radius_m=0.01,
        )


def test_direct_grid_function_field_rejects_nonfinite_position(
    affine_combined_function_field,
):
    field = GridFunctionField(affine_combined_function_field[0])
    with pytest.raises(ValueError, match="position_m must contain finite"):
        field.evaluate([np.nan, 0.0, 0.0])


def test_python_binding_matches_the_matlab_fixture_contract():
    fixture = (
        Path(__file__).parent / "fixtures" / "beam" / "affine_field_tetra.vol"
    )
    mesh = Mesh(str(fixture))
    field = GridFunction(VectorH1(mesh, order=1))
    with TaskManager():
        field.Set(CF((x, y, z)))

    lengths = np.array([0.2, 0.3])
    result = propagate_grid_function_linear_map(
        field,
        lengths,
        [[0.0, 0.2, 0.25], [0.0, 0.2, 0.75]],
        [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]],
        2.0,
        sample_radius_m=0.02,
        names=["entrance", "exit_quadrupole"],
    )

    generator = np.zeros((6, 6))
    generator[0, 1] = 1.0
    generator[1, 0] = -(0.1**2)
    generator[1, 5] = 0.1
    generator[2, 3] = 1.0
    generator[4, 0] = 0.1
    np.testing.assert_allclose(result["curvature_per_m"], 0.1, atol=1e-12)
    np.testing.assert_allclose(
        result["normal_gradient_per_m2"], 0.0, atol=1e-12
    )
    np.testing.assert_allclose(
        result["skew_gradient_per_m2"], 0.0, atol=1e-12
    )
    np.testing.assert_allclose(
        result["transverse_divergence_t_per_m"], 2.0, atol=1e-11
    )
    np.testing.assert_allclose(
        result["local_A_per_m"], np.tile(generator, (2, 1, 1)), atol=1e-12
    )
    np.testing.assert_allclose(
        result["R"], expm(generator * lengths.sum()), atol=2e-11
    )


def test_cubic_grid_function_builds_multipoles_and_nonlinear_map(
    cubic_multipole_field,
):
    field, normal, skew = cubic_multipole_field
    lengths = np.array([0.03, 0.04])
    result = propagate_grid_function_multipole_map(
        field,
        lengths,
        [[0.0, 0.0, 0.3], [0.0, 0.0, 0.7]],
        [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]],
        3.0,
        sample_radius_m=0.01,
        names=["sextupole_body", "octupole_body"],
        maximum_step_m=5e-4,
    )

    assert result["schema"] == (
        "radia.beam.grid-function-multipole-map.result.v1"
    )
    assert result["maximum_order"] == 3
    assert result["linearization_order"] == 3
    np.testing.assert_allclose(
        result["multipole_normal_t_per_m_power"],
        np.tile(normal, (2, 1)),
        rtol=2e-10,
        atol=2e-10,
    )
    np.testing.assert_allclose(
        result["multipole_skew_t_per_m_power"],
        np.tile(skew, (2, 1)),
        rtol=2e-10,
        atol=2e-10,
    )
    assert np.max(result["multipole_maximum_fit_residual_t"]) < 2e-11
    np.testing.assert_allclose(result["local_A_per_m"][:, 1, 0], -0.8)
    np.testing.assert_allclose(
        result["local_F2_per_m"][:, 1, 0, 0], -10.0
    )
    np.testing.assert_allclose(
        result["local_F3_per_m"][:, 1, 0, 0, 0], 160.0
    )
    assert np.max(np.abs(result["T"])) > 0.0
    assert np.max(np.abs(result["U"])) > 0.0
    assert result["diagnostics"]["T_reconstruction_error"] < 2e-11
    assert result["diagnostics"]["U_reconstruction_error"] < 2e-10


def test_rk_orbit_moving_frame_feeds_fourth_order_lie_map(cubic_multipole_field):
    field, normal, skew = cubic_multipole_field
    orbit = PlanarDesignOrbit(
        positions=np.array(
            [[0.0, 0.0, 0.3], [0.0, 0.0, 0.5], [0.0, 0.0, 0.7]]
        ),
        tangents=np.tile([0.0, 0.0, 1.0], (3, 1)),
        magnetic_rigidity=3.0,
        bend_axis=np.array([0.0, 1.0, 0.0]),
    )
    with TaskManager():
        result = fourth_order_lie_map_from_tracked_orbit(
            field,
            orbit,
            sample_radius_m=0.01,
            maximum_step_m=0.01,
        )

    expected = np.vstack(
        (
            np.full(2, normal[0]),
            np.full(2, normal[1]),
            np.full(2, skew[1]),
            np.full(2, normal[2]),
            np.full(2, skew[2]),
            np.full(2, normal[3]),
            np.full(2, skew[3]),
        )
    ).reshape(-1)
    np.testing.assert_allclose(result.multipole_response, expected, atol=2e-9)
    assert result.transfer.V.shape == (6, 6, 6, 6, 6)
    assert result.transfer.f5.shape == (6, 6, 6, 6, 6)
    assert np.max(np.abs(result.transfer.V)) > 0.0
    assert result.transfer.factorization.reconstructed_symplectic_residual.maximum < 1e-12


def test_multipole_schema_and_fit_model_do_not_depend_on_map_order(
    cubic_multipole_field,
):
    result = propagate_grid_function_multipole_map(
        cubic_multipole_field[0],
        [0.02],
        [[0.0, 0.0, 0.3]],
        [[0.0, 0.0, 1.0]],
        3.0,
        sample_radius_m=0.01,
        multipole_order=2,
        maximum_map_order=1,
    )
    assert result["schema"] == (
        "radia.beam.grid-function-multipole-map.result.v1"
    )
    assert result["maximum_order"] == 1
    assert result["linearization_order"] == 2
    assert result["fit_model"].endswith("through order 2")


def test_multipole_map_agrees_with_direct_grid_function_tracking(
    cubic_multipole_field,
):
    field = cubic_multipole_field[0]
    rigidity = 3.0
    length = 0.02
    map_result = propagate_grid_function_multipole_map(
        field,
        [length],
        [[0.0, 0.0, 0.3]],
        [[0.0, 0.0, 1.0]],
        rigidity,
        sample_radius_m=0.01,
        maximum_step_m=1e-4,
    )
    initial = np.array([2e-4, 3e-4, -1e-4, -2e-4, 0.0, 0.0])
    predicted = (
        map_result["R"] @ initial
        + 0.5 * np.einsum("ijk,j,k->i", map_result["T"], initial, initial)
        + (1.0 / 6.0)
        * np.einsum("ijkl,j,k,l->i", map_result["U"], initial, initial, initial)
    )

    species = ParticleSpecies.proton()
    momentum = species.charge_c * rigidity
    px = momentum * initial[1]
    py = momentum * initial[3]
    pz = np.sqrt(momentum * momentum - px * px - py * py)
    state = CartesianState([initial[0], initial[2], 0.3], [px, py, pz])
    equation = LorentzEquation(
        species, GridFunctionField(field), independent="path_length"
    )
    plan = TrackPlan()
    plan.start = 0.0
    plan.stop = length
    plan.maximum_step = 1e-4
    with TaskManager():
        trajectory = Tracker(equation, ClassicalRK4()).track(state, plan)
    final = trajectory.samples[-1]
    direct = np.array(
        [
            final.position_m[0],
            final.kinetic_momentum_kg_m_s[0] / momentum,
            final.position_m[1],
            final.kinetic_momentum_kg_m_s[1] / momentum,
        ]
    )
    np.testing.assert_allclose(predicted[:4], direct, atol=3e-8, rtol=2e-5)
    assert trajectory.summary.momentum_conservation_applicable
