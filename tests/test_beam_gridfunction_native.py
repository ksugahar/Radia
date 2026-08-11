from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

ngsolve = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

from netgen.occ import Box, Pnt
from ngsolve import CF, GridFunction, Mesh, TaskManager, VectorH1, x, y, z
from scipy.linalg import expm

from radia.beam import propagate_grid_function_linear_map


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
