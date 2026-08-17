from types import SimpleNamespace

import numpy as np

import radia.accelerator_magnet_topopt as magnet_topopt
from radia.accelerator_magnet_topopt import (
    MeasuredMedianPlaneFieldTarget,
    PlanarDesignOrbit,
    build_measured_median_plane_field_response_matrix,
    optimize_hdiv_mmm_magnet_to_measured_median_plane,
)


def _straight_orbit():
    return PlanarDesignOrbit(
        positions=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
        tangents=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]),
        magnetic_rigidity=2.0,
        bend_axis=np.array([0.0, 1.0, 0.0]),
        path_length_stations=np.array([0.0, 1.0]),
    )


def test_measured_median_plane_target_builds_physical_local_b_rows():
    measured = np.arange(12, dtype=float).reshape(2, 2, 3) * 1.0e-3
    target = MeasuredMedianPlaneFieldTarget(
        orbit=_straight_orbit(),
        s_m=np.array([0.25, 0.75]),
        x_m=np.array([-0.1, 0.1]),
        measured_B_local_t=measured,
        measurement_band_t=2.0e-5,
        components=("x", "y", "s"),
    )

    np.testing.assert_allclose(
        target.observation_points_m,
        [
            [-0.1, 0.0, 0.25],
            [0.1, 0.0, 0.25],
            [-0.1, 0.0, 0.75],
            [0.1, 0.0, 0.75],
        ],
    )
    expected_basis = np.eye(3)
    for point in range(4):
        np.testing.assert_allclose(
            target.observation_weights[3 * point : 3 * point + 3, point],
            expected_basis,
        )
    np.testing.assert_allclose(target.response_target, measured.reshape(-1))
    np.testing.assert_allclose(target.response_band, 2.0e-5)


def test_measured_response_matrix_uses_native_hdiv_functional_rows():
    target = MeasuredMedianPlaneFieldTarget(
        orbit=_straight_orbit(),
        s_m=np.array([0.25, 0.75]),
        x_m=np.array([-0.1, 0.1]),
        measured_B_local_t=np.zeros((2, 2, 1)),
        measurement_band_t=1.0e-4,
        components=("y",),
    )

    class ChargeGram:
        @staticmethod
        def configured_field_functional_rows(points, weights):
            assert points.shape == (4, 3)
            return weights[:, :, 1]

    rows = build_measured_median_plane_field_response_matrix(
        ChargeGram(), target, field_scale=2.0
    )
    np.testing.assert_allclose(rows, 2.0 * np.eye(4))


def test_measured_field_topology_objective_never_reconstructs_off_plane(
    monkeypatch,
):
    measured = np.zeros((2, 2, 3))
    measured[:, :, 1] = np.array([[1.0, 2.0], [3.0, 4.0]])
    target = MeasuredMedianPlaneFieldTarget(
        orbit=_straight_orbit(),
        s_m=np.array([0.25, 0.75]),
        x_m=np.array([-0.1, 0.1]),
        measured_B_local_t=measured,
        measurement_band_t=0.1,
        components=("y",),
    )
    response_matrix = np.array(
        [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [2.0, 1.0]]
    )
    generation = SimpleNamespace(
        state=np.array([1.0, 2.0]),
        active_elements=np.array([True, False]),
        converged=True,
    )
    captured = {}

    def fake_grow(**options):
        captured.update(options)
        return generation

    topology = object()
    monkeypatch.setattr(magnet_topopt, "grow_hdiv_mmm_by_superposition", fake_grow)
    monkeypatch.setattr(
        magnet_topopt, "ngsolve_growth_topology", lambda mesh, active: topology
    )
    result = optimize_hdiv_mmm_magnet_to_measured_median_plane(
        target,
        charge_gram=object(),
        fes=SimpleNamespace(ndof=2, mesh=object()),
        inv_chi=0.5,
        rhs=np.zeros(2),
        field_response_matrix=response_matrix,
        active_elements=np.array([True, False]),
        element_volumes=np.ones(2),
        volume_max=2.0,
    )

    np.testing.assert_allclose(captured["response_target"], [1.0, 2.0, 3.0, 4.0])
    np.testing.assert_allclose(captured["response_band"], 0.1)
    np.testing.assert_allclose(result.realized_field_response_t, [1.0, 2.0, 3.0, 4.0])
    assert result.maximum_measurement_band_ratio == 0.0
    assert result.topology is topology
    assert result.converged
