from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

import radia.accelerator_lie_topopt as lie_module

ngsolve = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

from netgen.occ import Box, Pnt
from ngsolve import CF, GridFunction, HCurl, Mesh, TaskManager, VectorH1, x, y, z
from scipy.linalg import expm

from radia.accelerator_lie_topopt import (
    _fourth_order_lie_map_from_vector_potential_polynomials as _internal_lie_from_jet,
)
from radia.accelerator_lie_topopt import (
    certify_p5_lie_aperture_against_b_map,
    compare_hcurl_a_map_to_b_coefficient_rk,
    compare_hcurl_lie_map_to_direct_rk,
    compare_tracked_lie_map_to_rk,
    differentiate_hcurl_transverse_lie_map,
    fourth_order_lie_map_from_hcurl_transverse,
    fourth_order_lie_map_from_tracked_orbit,
    fourth_order_lie_map_p_convergence,
    project_earlytimes_grid_function_maps,
    track_hcurl_vector_potential_canonical_s,
)
from radia.accelerator_magnet_topopt import PlanarDesignOrbit
from radia.beam import (
    CartesianState,
    ClassicalRK4,
    GridFunctionField,
    LorentzEquation,
    ParticleSpecies,
    Tracker,
    TrackPlan,
    build_curvilinear_beam_mesh,
    fit_transverse_vector_potential_polynomials,
    project_design_orbit_gauge,
    propagate_grid_function_linear_map,
    propagate_grid_function_multipole_map,
    propagate_hcurl_grid_function_multipole_map,
    sample_transverse_vector_potential,
)


def _double_reflection_step(p0, p1, t0, t1, horizontal):
    chord = p1 - p0
    reflected_horizontal = horizontal - (
        2.0 * np.dot(chord, horizontal) / np.dot(chord, chord)
    ) * chord
    reflected_tangent = t0 - (
        2.0 * np.dot(chord, t0) / np.dot(chord, chord)
    ) * chord
    difference = t1 - reflected_tangent
    if np.dot(difference, difference) > 1.0e-28:
        reflected_horizontal -= (
            2.0
            * np.dot(difference, reflected_horizontal)
            / np.dot(difference, difference)
        ) * difference
    reflected_horizontal -= np.dot(reflected_horizontal, t1) * t1
    return reflected_horizontal / np.linalg.norm(reflected_horizontal)


def _double_reflection_rmf(positions, tangents, initial_horizontal):
    """Independent transcription of Wang et al. (2008), Algorithm 1."""
    positions = np.asarray(positions, dtype=float)
    tangents = np.asarray(tangents, dtype=float)
    tangents = tangents / np.linalg.norm(tangents, axis=1)[:, None]
    horizontal = np.empty_like(tangents)
    horizontal[0] = initial_horizontal - (
        np.dot(initial_horizontal, tangents[0]) * tangents[0]
    )
    horizontal[0] /= np.linalg.norm(horizontal[0])
    for index in range(len(positions) - 1):
        horizontal[index + 1] = _double_reflection_step(
            positions[index],
            positions[index + 1],
            tangents[index],
            tangents[index + 1],
            horizontal[index],
        )
    return tangents, horizontal, np.cross(tangents, horizontal)


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


@pytest.fixture(scope="module")
def hcurl_cubic_vector_potential():
    geometry = Box(Pnt(-0.05, -0.05, 0.0), Pnt(0.05, 0.05, 1.0))
    mesh = geometry.GenerateMesh(maxh=0.1)
    vector_potential = GridFunction(HCurl(mesh, order=4))
    normal = np.array([0.0, 2.4, 15.0, -80.0])
    skew = np.array([0.0, -0.6, -4.0, 20.0])
    az = (
        -0.5 * normal[1] * x**2
        + 0.5 * normal[1] * y**2
        + skew[1] * x * y
        - normal[2] * x**3 / 3.0
        + normal[2] * x * y**2
        + skew[2] * x**2 * y
        - skew[2] * y**3 / 3.0
        - normal[3] * x**4 / 4.0
        + 1.5 * normal[3] * x**2 * y**2
        - normal[3] * y**4 / 4.0
        + skew[3] * x**3 * y
        - skew[3] * x * y**3
    )
    with TaskManager():
        vector_potential.Set(CF((0.0, 0.0, az)))
    return vector_potential, normal, skew


@pytest.fixture(scope="module")
def hcurl_p5_xy_vector_potential():
    geometry = Box(Pnt(-0.05, -0.05, 0.0), Pnt(0.05, 0.05, 1.0))
    mesh = geometry.GenerateMesh(maxh=0.1)
    vector_potential = GridFunction(HCurl(mesh, order=5))
    Ay = (
        -0.06 * y
        + 0.7 * x * y
        - 1.0e4 * x**2 * y**3
    )
    As = (
        -0.4 * x**2
        + 0.3 * y**2
        + 2.0 * x**3
        - 1.5 * x * y**2
        + 2.0e4 * x**3 * y**2
    )
    with TaskManager():
        vector_potential.Set(CF((0.0, Ay, As)))
    return vector_potential


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
        field_representation="magnetic_flux_density",
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
    assert result["frame_convention"] == (
        "right-handed Bishop/RMF double reflection seeded by "
        "initial_horizontal"
    )
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


def test_grid_function_frame_matches_bishop_double_reflection_on_helix(
    affine_combined_function_field,
):
    field = affine_combined_function_field[0]
    angle = np.linspace(0.0, 1.2, 6)
    positions = np.column_stack(
        (0.02 * np.cos(angle), 0.02 * np.sin(angle), 0.15 + 0.5 * angle)
    )
    tangents = np.column_stack(
        (-0.02 * np.sin(angle), 0.02 * np.cos(angle), np.full(6, 0.5))
    )
    initial_horizontal = np.array([1.0, 0.0, 0.0])
    expected_t, expected_h, expected_v = _double_reflection_rmf(
        positions, tangents, initial_horizontal
    )

    with TaskManager():
        result = propagate_grid_function_linear_map(
            field,
            np.full(6, 0.02),
            positions,
            tangents,
            3.0,
            sample_radius_m=0.002,
            initial_horizontal=initial_horizontal,
            maximum_step_m=0.01,
            field_representation="magnetic_flux_density",
        )

    np.testing.assert_allclose(result["frame_tangent"], expected_t, atol=2e-14)
    np.testing.assert_allclose(
        result["frame_horizontal"], expected_h, atol=2e-14
    )
    np.testing.assert_allclose(result["frame_vertical"], expected_v, atol=2e-14)
    np.testing.assert_allclose(
        np.einsum(
            "ij,ij->i", result["frame_horizontal"], result["frame_tangent"]
        ),
        0.0,
        atol=2e-14,
    )


def test_bishop_frame_rejects_duplicate_reference_stations(
    affine_combined_function_field,
):
    field = affine_combined_function_field[0]
    with pytest.raises(ValueError, match="consecutive reference positions"):
        propagate_grid_function_linear_map(
            field,
            [0.1, 0.1],
            [[0.0, 0.0, 0.4], [0.0, 0.0, 0.4]],
            [[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]],
            3.0,
            sample_radius_m=0.002,
            field_representation="magnetic_flux_density",
        )


def test_periodic_minimal_twist_frame_distributes_closed_orbit_holonomy(
    affine_combined_function_field,
):
    field = affine_combined_function_field[0]
    angle = np.linspace(0.0, 2.0 * np.pi, 24, endpoint=False)
    positions = np.column_stack(
        (
            0.03 * np.cos(angle),
            0.025 * np.sin(angle),
            0.5 + 0.008 * np.sin(2.0 * angle) + 0.006 * np.cos(3.0 * angle),
        )
    )
    tangents = np.column_stack(
        (
            -0.03 * np.sin(angle),
            0.025 * np.cos(angle),
            0.016 * np.cos(2.0 * angle) - 0.018 * np.sin(3.0 * angle),
        )
    )
    tangents /= np.linalg.norm(tangents, axis=1)[:, None]
    common = {
        "sample_radius_m": 0.001,
        "initial_horizontal": [1.0, 0.0, 0.0],
        "maximum_step_m": 0.01,
        "field_representation": "magnetic_flux_density",
    }
    with TaskManager():
        open_frame = propagate_grid_function_linear_map(
            field,
            np.full(len(angle), 0.01),
            positions,
            tangents,
            3.0,
            **common,
        )
        periodic = propagate_grid_function_linear_map(
            field,
            np.full(len(angle), 0.01),
            positions,
            tangents,
            3.0,
            periodic_frame=True,
            **common,
        )

    raw_h = open_frame["frame_horizontal"]
    closure_h = _double_reflection_step(
        positions[-1], positions[0], tangents[-1], tangents[0], raw_h[-1]
    )
    correction = np.arctan2(
        np.dot(tangents[0], np.cross(closure_h, raw_h[0])),
        np.dot(closure_h, raw_h[0]),
    )
    chord_lengths = np.linalg.norm(
        np.roll(positions, -1, axis=0) - positions, axis=1
    )
    cumulative = np.r_[0.0, np.cumsum(chord_lengths[:-1])]
    roll = correction * cumulative / chord_lengths.sum()
    expected_h = (
        np.cos(roll)[:, None] * raw_h
        + np.sin(roll)[:, None] * np.cross(tangents, raw_h)
    )

    assert abs(correction) > 0.1
    assert periodic["periodic_frame"]
    assert periodic["frame_convention"].startswith(
        "right-handed periodic minimal-twist"
    )
    np.testing.assert_allclose(
        periodic["frame_holonomy_correction_rad"], correction, atol=2e-14
    )
    np.testing.assert_allclose(
        periodic["frame_horizontal"], expected_h, atol=3e-14
    )

    closure_from_periodic = _double_reflection_step(
        positions[-1],
        positions[0],
        tangents[-1],
        tangents[0],
        periodic["frame_horizontal"][-1],
    )
    closing_roll = correction * chord_lengths[-1] / chord_lengths.sum()
    closure_from_periodic = (
        np.cos(closing_roll) * closure_from_periodic
        + np.sin(closing_roll) * np.cross(tangents[0], closure_from_periodic)
    )
    np.testing.assert_allclose(
        closure_from_periodic, periodic["frame_horizontal"][0], atol=4e-14
    )


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
            field_representation="magnetic_flux_density",
        )


def test_direct_grid_function_field_rejects_nonfinite_position(
    affine_combined_function_field,
):
    field = GridFunctionField(
        affine_combined_function_field[0], "magnetic_flux_density"
    )
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
        field_representation="magnetic_flux_density",
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
        field_representation="magnetic_flux_density",
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


def test_hcurl_vector_potential_is_the_default_and_curled_before_map_fit(
    hcurl_cubic_vector_potential,
):
    vector_potential, normal, skew = hcurl_cubic_vector_potential
    result = propagate_grid_function_multipole_map(
        vector_potential,
        [0.03],
        [[0.0, 0.0, 0.5]],
        [[0.0, 0.0, 1.0]],
        3.0,
        sample_radius_m=0.01,
        multipole_order=3,
        maximum_step_m=5e-4,
    )

    assert result["field_representation"] == "hcurl_vector_potential"
    assert result["magnetic_evaluation"] == "ngsolve-native-curl(A)"
    assert result["grid_function_space_class"] == "HCurlHighOrderFESpace"
    assert result["grid_function_space_order"] == 4
    np.testing.assert_allclose(
        result["multipole_normal_t_per_m_power"][0],
        normal,
        rtol=2e-9,
        atol=2e-9,
    )
    np.testing.assert_allclose(
        result["multipole_skew_t_per_m_power"][0],
        skew,
        rtol=2e-9,
        atol=2e-9,
    )
    assert result["multipole_maximum_fit_residual_t"][0] < 2e-10


def test_hcurl_vector_potential_mode_rejects_non_hcurl_grid_function(
    cubic_multipole_field,
):
    with pytest.raises(ValueError, match="requires a GridFunction.*HCurl"):
        propagate_hcurl_grid_function_multipole_map(
            cubic_multipole_field[0],
            [0.03],
            [[0.0, 0.0, 0.5]],
            [[0.0, 0.0, 1.0]],
            3.0,
        )


def test_direct_tracker_field_reads_native_hcurl_curl(
    hcurl_cubic_vector_potential,
):
    vector_potential, normal, skew = hcurl_cubic_vector_potential
    point = np.array([0.008, -0.006, 0.5])
    x_value, y_value = point[:2]
    expected_by = (
        normal[1] * x_value
        - skew[1] * y_value
        + normal[2] * (x_value**2 - y_value**2)
        - 2.0 * skew[2] * x_value * y_value
        + normal[3] * (x_value**3 - 3.0 * x_value * y_value**2)
        - skew[3] * (3.0 * x_value**2 * y_value - y_value**3)
    )
    expected_bx = (
        skew[1] * x_value
        + normal[1] * y_value
        + skew[2] * (x_value**2 - y_value**2)
        + 2.0 * normal[2] * x_value * y_value
        + skew[3] * (x_value**3 - 3.0 * x_value * y_value**2)
        + normal[3] * (3.0 * x_value**2 * y_value - y_value**3)
    )
    field = GridFunctionField(vector_potential)
    with TaskManager():
        sample = field.evaluate(point)
    np.testing.assert_allclose(
        sample.magnetic_t,
        [expected_bx, expected_by, 0.0],
        rtol=2e-10,
        atol=2e-10,
    )


def test_rk_orbit_moving_frame_feeds_fourth_order_lie_map(
    hcurl_cubic_vector_potential,
):
    vector_potential, normal, skew = hcurl_cubic_vector_potential
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
            vector_potential,
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
            np.zeros(2),
            np.zeros(2),
        )
    ).reshape(-1)
    np.testing.assert_allclose(
        result.multipole_response[:-4], expected[:-4], atol=2e-9
    )
    # Curl evaluation is exact to pointwise roundoff, but the absent degree-4
    # coefficient divides the ring moment by radius**4 and therefore amplifies
    # that roundoff.  Keep the physical zero gate absolute and explicit.
    assert np.max(np.abs(result.multipole_response[-4:])) < 5e-7
    assert result.field_fit["field_representation"] == "hcurl_vector_potential"
    assert result.field_fit["grid_function_space_order"] == 4
    assert result.field_fit["reference_curvature_source"] == (
        "design-orbit-tangent-turning"
    )
    np.testing.assert_array_equal(
        result.field_fit["reference_curvature_per_m"], 0.0
    )
    assert result.transfer.V.shape == (6, 6, 6, 6, 6)
    assert result.transfer.f5.shape == (6, 6, 6, 6, 6)
    assert np.max(np.abs(result.transfer.V)) > 0.0
    assert (
        result.transfer.factorization.reconstructed_symplectic_residual.maximum
        < 1e-12
    )


def test_tracked_lie_map_self_contained_difference_against_field_rk(
    hcurl_cubic_vector_potential,
):
    vector_potential, _, _ = hcurl_cubic_vector_potential
    orbit = PlanarDesignOrbit(
        positions=np.array(
            [[0.0, 0.0, 0.3], [0.0, 0.0, 0.5], [0.0, 0.0, 0.7]]
        ),
        tangents=np.tile([0.0, 0.0, 1.0], (3, 1)),
        magnetic_rigidity=3.0,
        bend_axis=np.array([0.0, 1.0, 0.0]),
        path_length_stations=np.array([0.0, 0.2, 0.4]),
    )
    initial = np.array(
        [
            [2.0e-4, 3.0e-4, -1.0e-4, -2.0e-4, 0.0, 0.0],
            [-1.0e-4, -2.0e-4, 1.5e-4, 1.0e-4, 0.0, 2.0e-4],
        ]
    )
    with TaskManager():
        lie_map = fourth_order_lie_map_from_tracked_orbit(
            vector_potential,
            orbit,
            sample_radius_m=0.01,
            maximum_step_m=0.002,
        )
        comparison = compare_tracked_lie_map_to_rk(
            vector_potential,
            lie_map,
            initial,
            maximum_step_m=2.0e-4,
            exit_plane_tolerance_m=1.0e-12,
        )

    assert comparison.comparison_indices == (0, 1, 2, 3, 5)
    np.testing.assert_array_equal(comparison.reference_curvature_per_m, 0.0)
    np.testing.assert_allclose(
        comparison.reference_field_rk_state[[0, 1, 2, 3, 5]],
        0.0,
        atol=2.0e-12,
    )
    assert comparison.reference_exit_plane_residual_m < 1.0e-12
    assert comparison.maximum_lie_vs_field_rk_error < 2.0e-7
    for case in comparison.cases:
        assert case.exit_plane_residual_m < 1.0e-12
        assert case.field_rk_relative_momentum_error < 2.0e-13
        assert np.isnan(case.lie_minus_field_rk[4])
        np.testing.assert_allclose(
            case.lie_state,
            case.polynomial_state,
            atol=2.0e-13,
            rtol=0.0,
        )


def test_curvilinear_s_rk_uses_direct_hcurl_a_and_independent_b_source(
    hcurl_cubic_vector_potential,
):
    vector_potential, normal, skew = hcurl_cubic_vector_potential
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
    source_magnetic_flux_density = CF((bx, by, 0.0))
    orbit = PlanarDesignOrbit(
        positions=np.array(
            [[0.0, 0.0, 0.3], [0.0, 0.0, 0.5], [0.0, 0.0, 0.7]]
        ),
        tangents=np.tile([0.0, 0.0, 1.0], (3, 1)),
        magnetic_rigidity=3.0,
        bend_axis=np.array([0.0, 1.0, 0.0]),
        path_length_stations=np.array([0.0, 0.2, 0.4]),
    )
    initial = np.array([2.0e-4, 3.0e-4, -1.0e-4, -2.0e-4, 0.0, 2.0e-4])
    with TaskManager():
        comparison = compare_hcurl_a_map_to_b_coefficient_rk(
            vector_potential,
            source_magnetic_flux_density,
            orbit,
            initial,
            field_mesh=vector_potential.space.mesh,
            integrator="DOP853",
            maximum_step_m=2.0e-4,
        )

    assert comparison.a_map.field_representation == (
        "hcurl-vector-potential-direct-A"
    )
    assert comparison.b_map.field_representation == (
        "hdiv-mmm-B-coefficient-direct"
    )
    assert comparison.a_map.field_evaluations > 0
    assert comparison.b_map.field_evaluations > 0
    assert abs(comparison.b_map.exit_plane_residual_m) < 1.0e-11
    assert comparison.maximum_final_state_difference < 2.0e-9


def test_earlytimes_boundary_projects_hdiv_mmm_coefficients_to_gridfunctions(
    hcurl_cubic_vector_potential,
):
    mesh = hcurl_cubic_vector_potential[0].space.mesh
    vector_potential_coefficient = CF((0.0, 0.0, -1.2 * x**2))
    magnetic_flux_density_coefficient = CF((0.0, 2.4 * x, 0.0))
    with TaskManager():
        maps = project_earlytimes_grid_function_maps(
            vector_potential_coefficient,
            mesh,
            magnetic_flux_density_coefficient=magnetic_flux_density_coefficient,
            project_magnetic_flux_density=True,
        )
        point = mesh(0.01, 0.0, 0.5)
        sampled_a = np.asarray(maps.vector_potential(point), dtype=float)
        sampled_b = np.asarray(maps.magnetic_flux_density(point), dtype=float)

    assert "HCurl" in type(maps.vector_potential.space).__name__
    assert "HDiv" in type(maps.magnetic_flux_density.space).__name__
    assert maps.vector_potential_order == 5
    assert maps.magnetic_flux_density_order == 4
    np.testing.assert_allclose(sampled_a, [0.0, 0.0, -1.2e-4], atol=2.0e-12)
    np.testing.assert_allclose(sampled_b, [0.0, 0.024, 0.0], atol=2.0e-12)


def test_earlytimes_lie_and_a_rk_reject_unprojected_a_coefficient_function():
    orbit = PlanarDesignOrbit(
        positions=np.array([[0.0, 0.0, 0.3], [0.0, 0.0, 0.7]]),
        tangents=np.tile([0.0, 0.0, 1.0], (2, 1)),
        magnetic_rigidity=3.0,
        bend_axis=np.array([0.0, 1.0, 0.0]),
        path_length_stations=np.array([0.0, 0.4]),
    )
    with pytest.raises(TypeError, match="GridFunction on HCurl"):
        track_hcurl_vector_potential_canonical_s(
            CF((0.0, 0.0, -1.2 * x**2)),
            orbit,
            np.zeros(6),
        )


def test_hcurl_volume_recovers_full_p5_xy_jet_and_lie_gradient(
    hcurl_p5_xy_vector_potential,
):
    orbit = PlanarDesignOrbit(
        positions=np.array(
            [[0.0, 0.0, 0.3], [0.0, 0.0, 0.5], [0.0, 0.0, 0.7]]
        ),
        tangents=np.tile([0.0, 0.0, 1.0], (3, 1)),
        magnetic_rigidity=3.0,
        bend_axis=np.array([0.0, 1.0, 0.0]),
        path_length_stations=np.array([0.0, 0.2, 0.4]),
    )
    x_offsets = np.linspace(-0.02, 0.02, 13)
    y_offsets = np.linspace(-0.02, 0.02, 9)
    with TaskManager():
        result = fourth_order_lie_map_from_hcurl_transverse(
            hcurl_p5_xy_vector_potential,
            orbit,
            x_offsets,
            y_offsets,
            fit_tolerance_t_m=2.0e-9,
            left_right_tolerance_t_m=2.0e-9,
            reference_orbit_tolerance=2.0e-7,
            maximum_step_m=0.01,
        )

    fit = result.polynomial_fit
    x_radius = np.max(np.abs(x_offsets))
    y_radius = np.max(np.abs(y_offsets))
    np.testing.assert_allclose(
        fit.Ay_coefficients_t_m[:, 1, 1] * x_radius * y_radius,
        0.7 * x_radius * y_radius,
        atol=3.0e-10,
    )
    np.testing.assert_allclose(
        fit.Ay_coefficients_t_m[:, 2, 3] * x_radius**2 * y_radius**3,
        -1.0e4 * x_radius**2 * y_radius**3,
        atol=3.0e-10,
    )
    np.testing.assert_allclose(
        fit.As_coefficients_t_m[:, 3, 2] * x_radius**3 * y_radius**2,
        2.0e4 * x_radius**3 * y_radius**2,
        atol=3.0e-10,
    )
    assert fit.maximum_Ax_t_m < 1.0e-10
    assert fit.maximum_orbit_Ay_As_t_m < 1.0e-10
    assert result.field_certificate.symmetry_class == "normal"
    assert result.field_certificate.maximum_symmetry_defect_t_m < 2.0e-9
    assert fit.maximum_left_right_scaled_coefficient_discrepancy_t_m < 2.0e-9
    assert np.max(np.abs(result.lie_map.hamiltonian_linear)) < 2.0e-7
    assert np.max(np.abs(result.lie_map.transfer.V)) > 0.0

    transfer = result.lie_map.transfer
    vertical_parity = np.asarray([0, 0, 1, 1, 0, 0], dtype=np.int64)
    for tensor_name in ("R", "T", "U", "V"):
        tensor = getattr(transfer, tensor_name)
        for index in np.ndindex(tensor.shape):
            output_parity = vertical_parity[index[0]]
            input_parity = sum(vertical_parity[item] for item in index[1:]) % 2
            if output_parity != input_parity:
                assert abs(tensor[index]) < 2.0e-12
    assert abs(transfer.T[3, 0, 2]) > 0.1  # allowed second-order x-y term
    assert abs(transfer.T[0, 0, 5]) > 1.0e-3  # allowed x-delta term

    x_grid, y_grid = np.meshgrid(x_offsets, y_offsets, indexing="ij")
    response_shape = result.samples.Ay_t_m.shape + (1,)
    Ay_response = np.zeros(response_shape)
    As_response = np.zeros(response_shape)
    Ay_response[..., 0] = x_grid * y_grid
    differentiated = differentiate_hcurl_transverse_lie_map(
        result,
        Ay_response,
        As_response,
    )
    assert differentiated.derivative_backend == (
        "sample-fit-linear-chain-plus-forward-hamiltonian-lie-ad"
    )
    np.testing.assert_allclose(
        differentiated.hamiltonian_linear_response,
        0.0,
        atol=2.0e-12,
    )

    step = 1.0e-4
    plus_Ay = fit.Ay_coefficients_t_m.copy()
    minus_Ay = fit.Ay_coefficients_t_m.copy()
    plus_Ay[:, 1, 1] += step
    minus_Ay[:, 1, 1] -= step
    options = {
        "segment_lengths": orbit.segment_lengths,
        "magnetic_rigidity": orbit.magnetic_rigidity,
        "reference_curvature_per_m": orbit.signed_curvature,
        "maximum_step_m": 0.01,
        "reference_orbit_tolerance": 2.0e-7,
    }
    plus = _internal_lie_from_jet(
        plus_Ay,
        fit.As_coefficients_t_m,
        **options,
    )
    minus = _internal_lie_from_jet(
        minus_Ay,
        fit.As_coefficients_t_m,
        **options,
    )
    for name in ("R", "T", "U", "V", "f3", "f4", "f5"):
        finite_difference = (
            getattr(plus.transfer, name) - getattr(minus.transfer, name)
        ) / (2.0 * step)
        np.testing.assert_allclose(
            getattr(differentiated, name)[0],
            finite_difference,
            rtol=2.0e-4,
            atol=2.0e-8,
        )
    objective_gradient = differentiated.objective_gradient(
        V=result.lie_map.transfer.V
    )
    finite_difference_objective = (
        0.5 * np.sum(plus.transfer.V**2)
        - 0.5 * np.sum(minus.transfer.V**2)
    ) / (2.0 * step)
    np.testing.assert_allclose(
        objective_gradient[0],
        finite_difference_objective,
        rtol=2.0e-4,
        atol=2.0e-8,
    )
    orbit_breaking_response = np.zeros(response_shape)
    orbit_breaking_response[..., 0] = x_grid
    orbit_constraint = differentiate_hcurl_transverse_lie_map(
        result,
        np.zeros(response_shape),
        orbit_breaking_response,
    ).hamiltonian_linear_response
    assert np.max(np.abs(orbit_constraint)) > 0.1

    Bx = 0.6 * y - 3.0 * x * y + 4.0e4 * x**3 * y
    By = 0.8 * x - 6.0 * x**2 + 1.5 * y**2 - 6.0e4 * x**2 * y**2
    Bs = 0.7 * y - 2.0e4 * x * y**3
    B_coefficient = CF((Bx, By, Bs))
    initial = np.array(
        [1.0e-4, 2.0e-4, -8.0e-5, -1.0e-4, 0.0, 1.0e-4]
    )
    with TaskManager():
        comparison = compare_hcurl_lie_map_to_direct_rk(
            hcurl_p5_xy_vector_potential,
            B_coefficient,
            result,
            initial,
            field_mesh=hcurl_p5_xy_vector_potential.space.mesh,
            maximum_step_m=2.0e-4,
        )
    assert len(comparison.cases) == 1
    assert comparison.cases[0].a_field_evaluations > 0
    assert comparison.cases[0].b_field_evaluations > 0
    assert abs(comparison.cases[0].b_exit_plane_residual_m) < 1.0e-11
    assert comparison.maximum_lie_truncation_error < 2.0e-9
    assert comparison.maximum_a_b_field_route_error < 2.0e-9
    assert comparison.maximum_total_lie_b_error < 3.0e-9

    with TaskManager():
        certificate = certify_p5_lie_aperture_against_b_map(
            hcurl_p5_xy_vector_potential,
            B_coefficient,
            result,
            [1.0e-4],
            field_mesh=hcurl_p5_xy_vector_potential.space.mesh,
            spatial_angles=4,
            normalized_momentum_radii=(0.0,),
            momentum_angles=4,
            delta_values=(0.0,),
            position_tolerance_m=5.0e-7,
            normalized_momentum_tolerance=5.0e-7,
            delta_tolerance=1.0e-10,
            maximum_step_m=2.0e-3,
            exit_plane_tolerance_m=1.0e-11,
        )
    assert certificate.field_representation == (
        "hdiv-mmm-B-coefficient-direct"
    )
    assert certificate.certified_radius_m == pytest.approx(1.0e-4)
    assert certificate.checks[0].case_count == 4
    assert certificate.checks[0].passed


def test_hcurl_volume_reports_incompatible_left_and_right_jets(
    hcurl_p5_xy_vector_potential,
):
    orbit = PlanarDesignOrbit(
        positions=np.array([[0.0, 0.0, 0.3], [0.0, 0.0, 0.7]]),
        tangents=np.tile([0.0, 0.0, 1.0], (2, 1)),
        magnetic_rigidity=3.0,
        bend_axis=np.array([0.0, 1.0, 0.0]),
        path_length_stations=np.array([0.0, 0.4]),
    )
    x_offsets = np.linspace(-0.02, 0.02, 13)
    y_offsets = np.linspace(-0.02, 0.02, 9)
    with TaskManager():
        samples = sample_transverse_vector_potential(
            hcurl_p5_xy_vector_potential,
            orbit,
            x_offsets,
            y_offsets,
        )
    incompatible_As = samples.As_t_m.copy()
    incompatible_As[:, x_offsets > 0.0, :] += (
        0.02 * x_offsets[x_offsets > 0.0, None] ** 2
    )
    incompatible = replace(samples, As_t_m=incompatible_As)
    with pytest.raises(RuntimeError, match="left/right transverse A jets"):
        fit_transverse_vector_potential_polynomials(
            incompatible,
            degree=5,
            left_right_tolerance_t_m=1.0e-7,
        )

    parity_broken = replace(
        samples,
        Ay_t_m=samples.Ay_t_m + 0.01 * samples.x_m[None, :, None],
    )
    with pytest.raises(ValueError, match="median-plane symmetry"):
        fit_transverse_vector_potential_polynomials(
            parity_broken,
            degree=5,
            symmetry_tolerance_t_m=1.0e-8,
        )


def test_hcurl_lie_public_api_does_not_accept_arbitrary_polynomial_maps():
    assert "canonical_vector_potential_hamiltonian_jet" not in lie_module.__all__
    assert (
        "fourth_order_lie_map_from_vector_potential_polynomials"
        not in lie_module.__all__
    )


def test_lie_map_p_convergence_includes_odd_orders_and_selects_p5(
    hcurl_cubic_vector_potential,
):
    vector_potential, _, _ = hcurl_cubic_vector_potential
    orbit = PlanarDesignOrbit(
        positions=np.array(
            [[0.0, 0.0, 0.3], [0.0, 0.0, 0.5], [0.0, 0.0, 0.7]]
        ),
        tangents=np.tile([0.0, 0.0, 1.0], (3, 1)),
        magnetic_rigidity=3.0,
        bend_axis=np.array([0.0, 1.0, 0.0]),
    )
    with TaskManager():
        study = fourth_order_lie_map_p_convergence(
            vector_potential,
            vector_potential.space.mesh,
            orbit,
            orders=(3, 4, 5),
            minimum_order=5,
            sample_radius_m=0.01,
            maximum_step_m=0.01,
        )

    assert tuple(step.order for step in study.steps) == (3, 4, 5)
    assert study.converged
    assert study.selected_order == 5
    assert study.selected_step.order == 5
    assert study.steps[-1].maximum_normalized_change < 1.0
    assert study.steps[-1].symplectic_residual < 1.0e-12
    np.testing.assert_array_equal(
        study.selected_step.result.multipole_response[-4:], 0.0
    )
    assert np.any(
        study.selected_step.result.field_fit[
            "lie_normal_multipole_noise_clipped"
        ][:, 4]
    )


def test_design_orbit_gauge_preserves_the_earlytimes_lie_map():
    radius = 0.25
    angles = np.linspace(0.0, 0.4, 5)
    positions = np.column_stack(
        (
            radius * np.sin(angles),
            np.zeros_like(angles),
            radius * np.cos(angles),
        )
    )
    tangents = np.column_stack(
        (np.cos(angles), np.zeros_like(angles), -np.sin(angles))
    )
    orbit = PlanarDesignOrbit(
        positions=positions,
        tangents=tangents,
        magnetic_rigidity=1.5,
        bend_axis=np.array([0.0, 1.0, 0.0]),
    )
    tube = build_curvilinear_beam_mesh(
        orbit,
        half_width_m=0.008,
        half_height_m=0.004,
        maxh_m=0.008,
        curve_order=2,
    )
    coefficient = CF(
        (0.2 + 0.1 * x, 0.3 + 0.2 * z, -0.1 + 0.05 * y)
    )
    with TaskManager():
        gauged = project_design_orbit_gauge(
            coefficient,
            tube,
            order=5,
            gauge_tolerance=2.0e-6,
        )
        raw_map = fourth_order_lie_map_from_tracked_orbit(
            gauged.ungauged_vector_potential,
            orbit,
            sample_radius_m=0.002,
            maximum_step_m=0.005,
        )
        gauged_map = fourth_order_lie_map_from_tracked_orbit(
            gauged.vector_potential,
            orbit,
            sample_radius_m=0.002,
            maximum_step_m=0.005,
        )

    assert gauged.maximum_orbit_gauge_residual_t_m < 2.0e-6
    assert gauged.curl_change_l2_t_m32 < 1.0e-9
    np.testing.assert_allclose(
        gauged_map.field_fit["reference_curvature_per_m"],
        orbit.signed_curvature,
        rtol=0.0,
        atol=0.0,
    )
    assert gauged_map.field_fit["reference_curvature_source"] == (
        "design-orbit-tangent-turning"
    )
    for name in ("R", "T", "U", "V", "f3", "f4", "f5"):
        np.testing.assert_allclose(
            getattr(gauged_map.transfer, name),
            getattr(raw_map.transfer, name),
            atol=1.0e-12,
            rtol=0.0,
        )


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
        field_representation="magnetic_flux_density",
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
        field_representation="magnetic_flux_density",
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
        species,
        GridFunctionField(field, "magnetic_flux_density"),
        independent="path_length",
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
