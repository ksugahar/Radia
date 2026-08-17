from types import SimpleNamespace

import numpy as np


def _one_segment_orbit(*, rigidity, radius=4.0, angle=0.04):
    from radia.accelerator_magnet_topopt import PlanarDesignOrbit

    positions = np.asarray([
        [radius, 0.0, 0.0],
        [radius * np.cos(angle), radius * np.sin(angle), 0.0],
    ])
    tangents = np.asarray([
        [0.0, 1.0, 0.0],
        [-np.sin(angle), np.cos(angle), 0.0],
    ])
    return PlanarDesignOrbit(
        positions, tangents, magnetic_rigidity=rigidity,
        bend_axis=np.asarray([0.0, 0.0, 1.0]))


def _single_taylor_objective(*, rigidity, raw):
    from radia.accelerator_taylor_topopt import (
        PlanarSecondOrderTaylorMapObjective,
        second_order_taylor_map_from_multipoles,
    )

    orbit = _one_segment_orbit(rigidity=rigidity)
    transfer = second_order_taylor_map_from_multipoles(
        raw, orbit.segment_lengths, rigidity, maximum_step_m=0.1)
    return PlanarSecondOrderTaylorMapObjective(
        orbit=orbit, target_R=transfer.R, target_T=transfer.T,
        R_band=1e-3, T_band=1e-2, normal_dipole_band=1e-3,
        R_entries=((0, 0), (0, 1), (1, 5)),
        T_entries=((0, 0, 5), (1, 0, 0)), maximum_step_m=0.1)


def test_multi_momentum_taylor_objective_has_block_ad_jacobian():
    from radia.ffag_taylor_topopt import (
        MultiMomentumSecondOrderTaylorMapObjective,
    )

    raw_a = np.asarray([0.21, 0.015, -0.004, 0.002, -0.001])
    raw_b = np.asarray([0.44, -0.012, 0.003, -0.001, 0.002])
    objective = MultiMomentumSecondOrderTaylorMapObjective((
        _single_taylor_objective(rigidity=0.8, raw=raw_a),
        _single_taylor_objective(rigidity=2.4, raw=raw_b),
    ))
    raw = np.r_[raw_a, raw_b]
    jacobian = objective.transform_jacobian(raw)
    first_rows = objective.objectives[0].response_target.size

    assert objective.raw_offsets.tolist() == [0, 5, 10]
    assert objective.source_calibration_rows.tolist() == [0, 5]
    assert objective.response_group_indices(("R",)).tolist() == [1, 2, 3, 7, 8, 9]
    assert objective.response_group_indices(("normal_dipole", "T")).tolist() == [
        0, 4, 5, 6, 10, 11]
    assert objective.group_max_band_ratio(objective.transform(raw), "R") == 0.0
    assert jacobian.shape == (2 * first_rows, 10)
    np.testing.assert_allclose(jacobian[:first_rows, 5:], 0.0)
    np.testing.assert_allclose(jacobian[first_rows:, :5], 0.0)
    direction = np.asarray([
        0.2, -0.1, 0.05, 0.08, -0.04,
        -0.1, 0.12, 0.03, -0.05, 0.07,
    ])
    step = 2e-7
    regression = (
        objective.transform(raw + step * direction)
        - objective.transform(raw - step * direction)) / (2.0 * step)
    np.testing.assert_allclose(
        jacobian @ direction, regression, rtol=2e-6, atol=2e-9)


def test_coilbuilder_incident_multipoles_use_same_observation_rows():
    from radia.accelerator_magnet_topopt import CoilBuilderHDivSource
    from radia.accelerator_taylor_topopt import (
        planar_orbit_multipole_observations,
    )
    from radia.ffag_taylor_topopt import (
        MultiMomentumSecondOrderTaylorMapObjective,
        incident_multi_orbit_multipole_response,
    )

    raw = np.asarray([0.2, 0.01, 0.0, 0.0, 0.0])
    objective = MultiMomentumSecondOrderTaylorMapObjective((
        _single_taylor_objective(rigidity=0.8, raw=raw),
        _single_taylor_objective(rigidity=2.4, raw=2.0 * raw),
    ))
    square = np.asarray([
        [[-1.0, -1.0, 0.5], [1.0, -1.0, 0.5]],
        [[1.0, -1.0, 0.5], [1.0, 1.0, 0.5]],
        [[1.0, 1.0, 0.5], [-1.0, 1.0, 0.5]],
        [[-1.0, 1.0, 0.5], [-1.0, -1.0, 0.5]],
    ])
    source = CoilBuilderHDivSource(((square, 1.2e5),))
    radius = 2e-3
    expected = []
    for item in objective.objectives:
        points, weights = planar_orbit_multipole_observations(
            item.orbit, sample_radius=radius)
        expected.append(np.einsum(
            "rpc,pc->r", weights, source.b_field(points)))

    actual = incident_multi_orbit_multipole_response(
        source, objective, sample_radius=radius)

    np.testing.assert_allclose(actual, np.concatenate(expected), rtol=1e-14)


def test_ffag_taylor_optimizer_eliminates_source_scale_on_every_solve(
        monkeypatch):
    import radia.accelerator_magnet_topopt as accelerator
    import radia.ffag_taylor_topopt as fusion
    import radia.topology_optimization as topopt

    raw_a = np.asarray([0.2, 0.01, 0.0, 0.0, 0.0])
    raw_b = np.asarray([0.4, -0.01, 0.0, 0.0, 0.0])
    objective = fusion.MultiMomentumSecondOrderTaylorMapObjective((
        _single_taylor_objective(rigidity=0.8, raw=raw_a),
        _single_taylor_objective(rigidity=2.4, raw=raw_b),
    ))
    square = np.asarray([
        [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0]],
        [[1.0, 0.0, 1.0], [0.0, 0.0, 1.0]],
    ])
    source = accelerator.CoilBuilderHDivSource(((square, 1.0),))
    monkeypatch.setattr(
        accelerator.CoilBuilderHDivSource, "assemble_hdiv_rhs",
        lambda self, fes: np.ones(fes.ndof))
    active = np.asarray([True, False, True])
    captured = {}

    def fake_grow(**kwargs):
        captured.update(kwargs)
        response = np.r_[raw_a, raw_b]
        return topopt.HDivMMMGenerationResult(
            active_elements=active.copy(), state=np.ones(4),
            response=response, history=tuple(), converged=False,
            source_scale=0.5,
            objective_response=kwargs["response_transform"](response),
            stop_reason="contract regression")

    monkeypatch.setattr(fusion, "grow_hdiv_mmm_by_superposition", fake_grow)
    monkeypatch.setattr(
        fusion, "ngsolve_growth_topology",
        lambda mesh, values: SimpleNamespace(valid=True))
    fes = SimpleNamespace(ndof=4, mesh=object())
    result = fusion.optimize_ffag_hdiv_mmm_from_second_order_taylor_maps(
        objective, source=source, charge_gram=object(), fes=fes,
        inv_chi=0.1, active_elements=active,
        element_volumes=np.ones(3), volume_max=3.0,
        sample_radius=1e-3, source_scale=2.0,
        multipole_response_matrix=np.zeros((10, 4)),
        incident_multipole_response=np.zeros(10),
        fixed_active_elements=np.asarray([True, False, False]))

    assert result.source_scale == 1.0
    np.testing.assert_array_equal(
        captured["source_calibration_rows"], [0, 5])
    np.testing.assert_allclose(
        captured["source_calibration_target"],
        objective.source_calibration_target)
    np.testing.assert_allclose(
        captured["source_calibration_band"],
        objective.source_calibration_band)
    assert captured["source_calibration_norm"] == "linf"
    assert captured["response_transform"].__self__ is objective
    assert captured["response_transform_jacobian"].__self__ is objective
    np.testing.assert_array_equal(result.active_elements, active)

    captured.clear()
    response = np.r_[raw_a, raw_b]
    dipole_limit = (
        objective.group_max_band_ratio(objective.transform(response),
                                       "normal_dipole") + 0.1)
    hierarchical = fusion.optimize_ffag_hdiv_mmm_from_second_order_taylor_maps(
        objective, source=source, charge_gram=object(), fes=fes,
        inv_chi=0.1, active_elements=active,
        element_volumes=np.ones(3), volume_max=3.0,
        sample_radius=1e-3, source_scale=2.0,
        multipole_response_matrix=np.zeros((10, 4)),
        incident_multipole_response=np.zeros(10),
        primary_response_groups=("R",),
        maximum_group_band_ratios={"normal_dipole": dipole_limit},
        fixed_active_elements=np.asarray([True, False, False]))
    R_rows = objective.response_group_indices(("R",))
    np.testing.assert_allclose(
        captured["response_target"], objective.response_target[R_rows])
    np.testing.assert_allclose(
        captured["response_band"], objective.response_band[R_rows])
    primary = captured["response_transform"](response)
    assert primary.shape == (len(R_rows),)
    assert captured["response_transform_jacobian"](response).shape == (
        len(R_rows), response.size)
    guard = captured["exact_response_validator"]
    assert guard(response, primary)
    bad = response.copy()
    bad[0] = (
        objective.objectives[0].required_normal_dipole[0]
        + (dipole_limit + 1.0)
        * objective.objectives[0].normal_dipole_band[0])
    assert not guard(bad, captured["response_transform"](bad))
    assert hierarchical.primary_response_groups == ("R",)
    assert hierarchical.primary_max_band_ratio == 0.0
    assert hierarchical.maximum_group_band_ratios == (
        ("normal_dipole", dipole_limit),)


def test_ffag_target_lift_preserves_requested_R_and_uses_loose_T_band():
    from radia.ffag_topopt import build_ffag_cell_target_family
    from radia.ffag_taylor_topopt import (
        build_ffag_second_order_taylor_objective,
    )

    family = build_ffag_cell_target_family(
        [31.0, 250.0], n_segments=16,
        transfer_matrix_band=0.2, bend_field_band=0.1)
    objective = build_ffag_second_order_taylor_objective(
        family, T_band=1e6, maximum_step_m=0.1)

    assert len(objective.objectives) == 2
    for index, item in enumerate(objective.objectives):
        np.testing.assert_allclose(
            item.target_R, family.objective.target_matrices[index])
        np.testing.assert_allclose(item.T_band, 1e6)
        assert item.raw_field_response_size == 5 * 16
