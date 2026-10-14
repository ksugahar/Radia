import numpy as np
import pytest

from radia.accelerator_magnet_topopt import PlanarDesignOrbit
from radia.accelerator_taylor_topopt import (
    SECOND_ORDER_MULTIPOLE_COMPONENTS,
    DecoupledFirstOrderTarget,
    PlanarSecondOrderTaylorMapObjective,
    Symplectic2x2KAN,
    certify_decoupled_first_order_reachability,
    optimize_hdiv_mmm_magnet_from_second_order_taylor_map,
    planar_orbit_multipole_observations,
    run_second_order_taylor_material_inverse_pipeline,
    second_order_taylor_map_from_multipoles,
    solve_decoupled_first_order_continuation,
)


def _straight_orbit(length=0.08, rigidity=3.0):
    return PlanarDesignOrbit(
        positions=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, length]]),
        tangents=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]),
        magnetic_rigidity=rigidity,
        bend_axis=np.array([0.0, 1.0, 0.0]),
    )


def _xy_objective(orbit, wanted):
    target = second_order_taylor_map_from_multipoles(
        wanted, orbit.segment_lengths, orbit.magnetic_rigidity, maximum_step_m=0.01
    )
    return PlanarSecondOrderTaylorMapObjective(
        orbit=orbit,
        target_R=target.R,
        target_T=target.T,
        R_band=1.0e-5,
        T_band=1.0e-5,
        normal_dipole_band=1.0e-6,
        R_entries=((1, 2), (3, 0)),
        T_entries=((1, 0, 2), (3, 0, 0), (3, 2, 2)),
        maximum_step_m=0.01,
    )


def test_nine_point_harmonic_rows_recover_normal_and_skew_multipoles():
    orbit = _straight_orbit()
    points, weights = planar_orbit_multipole_observations(orbit, sample_radius=0.01)
    coefficients = np.array([0.2, 1.1, -0.4, 2.3, 0.7])
    coordinate = points[:, 0] + 1.0j * points[:, 1]
    complex_field = (
        coefficients[0]
        + (coefficients[1] + 1.0j * coefficients[2]) * coordinate
        + (coefficients[3] + 1.0j * coefficients[4]) * coordinate**2
    )
    field = np.column_stack(
        (complex_field.imag, complex_field.real, np.zeros(points.shape[0]))
    )
    recovered = np.einsum("rpd,pd->r", weights, field)

    assert SECOND_ORDER_MULTIPOLE_COMPONENTS == (
        "normal_dipole",
        "normal_quadrupole",
        "skew_quadrupole",
        "normal_sextupole",
        "skew_sextupole",
    )
    np.testing.assert_allclose(recovered, coefficients, atol=2.0e-13)


def test_second_order_R_T_forward_ad_matches_centered_difference():
    raw = np.array(
        [
            0.20,
            0.22,
            1.00,
            0.80,
            0.30,
            -0.20,
            2.00,
            1.50,
            0.70,
            -0.50,
        ]
    )
    lengths = np.array([0.03, 0.04])
    differentiated = second_order_taylor_map_from_multipoles(
        raw, lengths, 3.0, maximum_step_m=0.01
    )

    assert differentiated.derivative_backend == ("forward-mode-rk4-taylor-ad")
    assert differentiated.value_backend == "native-cpp-variational-map"
    assert abs(differentiated.R[1, 2]) > 1.0e-5
    assert abs(differentiated.R[3, 0]) > 1.0e-5
    assert abs(differentiated.T[1, 0, 2]) > 1.0e-5
    np.testing.assert_allclose(
        differentiated.T, np.swapaxes(differentiated.T, 1, 2), atol=3.0e-15
    )

    for parameter in range(raw.size):
        step = 1.0e-6 * max(1.0, abs(raw[parameter]))
        direction = np.zeros_like(raw)
        direction[parameter] = step
        plus = second_order_taylor_map_from_multipoles(
            raw + direction, lengths, 3.0, maximum_step_m=0.01
        )
        minus = second_order_taylor_map_from_multipoles(
            raw - direction, lengths, 3.0, maximum_step_m=0.01
        )
        np.testing.assert_allclose(
            differentiated.R_jacobian[parameter],
            (plus.R - minus.R) / (2.0 * step),
            rtol=2.0e-7,
            atol=5.0e-10,
        )
        np.testing.assert_allclose(
            differentiated.T_jacobian[parameter],
            (plus.T - minus.T) / (2.0 * step),
            rtol=2.0e-7,
            atol=5.0e-10,
        )


def test_kan_parameterization_covers_random_symplectic_blocks_and_target_gate():
    rng = np.random.default_rng(20260816)
    for _ in range(32):
        parameters = rng.uniform(
            [-4.0 * np.pi, -1.5, -3.0],
            [4.0 * np.pi, 1.5, 3.0],
        )
        original = Symplectic2x2KAN(*parameters)
        recovered = Symplectic2x2KAN.from_matrix(original.matrix)
        assert np.linalg.det(original.matrix) == pytest.approx(1.0, abs=2.0e-14)
        np.testing.assert_allclose(recovered.matrix, original.matrix, atol=3.0e-14)

    target = DecoupledFirstOrderTarget(
        Symplectic2x2KAN(0.4, -0.2, 0.7),
        Symplectic2x2KAN(-0.3, 0.1, -0.5),
        0.08,
    )
    recovered = DecoupledFirstOrderTarget.from_matrix(target.matrix)
    np.testing.assert_allclose(recovered.matrix, target.matrix, atol=3.0e-15)

    nonsymplectic = target.matrix.copy()
    nonsymplectic[0, 0] *= 1.1
    with pytest.raises(ValueError, match="determinant one"):
        DecoupledFirstOrderTarget.from_matrix(nonsymplectic)
    coupled = target.matrix.copy()
    coupled[0, 2] = 0.01
    with pytest.raises(ValueError, match="must decouple"):
        DecoupledFirstOrderTarget.from_matrix(coupled)
    rf_block = target.matrix.copy()
    rf_block[5, 4] = 0.02
    with pytest.raises(ValueError, match="magnetostatic longitudinal"):
        DecoupledFirstOrderTarget.from_matrix(rf_block)


def test_decoupled_first_order_rank_certificate_rejects_drift_and_accepts_seed():
    count = 8
    lengths = np.full(count, 0.08)
    drift = np.zeros((5, count))
    drift_certificate = certify_decoupled_first_order_reachability(
        drift.reshape(-1),
        lengths,
        3.0,
        maximum_step_m=0.008,
    )
    assert drift_certificate.transverse_numerical_rank == 3
    assert not drift_certificate.full_transverse_rank

    seeded = drift.copy()
    seeded[1] = np.asarray(
        [70.5, -123.0, 0.0, -157.0, -8.33, 106.0, 131.0, 0.1]
    )
    seeded_certificate = certify_decoupled_first_order_reachability(
        seeded.reshape(-1),
        lengths,
        3.0,
        maximum_step_m=0.008,
    )
    assert seeded_certificate.transverse_numerical_rank == 6
    assert seeded_certificate.full_transverse_rank
    assert seeded_certificate.decoupling_residual == 0.0
    # The current piecewise RK4 value path is not a symplectic integrator;
    # this residual is its checked step-size error, not a target defect.
    assert seeded_certificate.symplectic_residual < 5.0e-9

    for direction in (1, 6):
        step = 1.0e-4
        plus = seeded.copy()
        minus = seeded.copy()
        plus[1, direction] += step
        minus[1, direction] -= step
        plus_parameters = certify_decoupled_first_order_reachability(
            plus.reshape(-1), lengths, 3.0, maximum_step_m=0.008
        ).parameters
        minus_parameters = certify_decoupled_first_order_reachability(
            minus.reshape(-1), lengths, 3.0, maximum_step_m=0.008
        ).parameters
        finite_difference = (plus_parameters - minus_parameters) / (2.0 * step)
        np.testing.assert_allclose(
            seeded_certificate.directional_parameter_jacobian[:, direction],
            finite_difference,
            rtol=8.0e-6,
            atol=3.0e-9,
        )

    drift_target = DecoupledFirstOrderTarget.from_matrix(
        drift_certificate.transfer_matrix
    )
    rejected = solve_decoupled_first_order_continuation(
        drift.reshape(-1),
        lengths,
        3.0,
        drift_target,
        maximum_step_m=0.008,
    )
    assert not rejected.converged
    assert rejected.final_transverse_rank == 3
    assert "rank 3" in rejected.status


def test_rank6_continuation_reaches_random_nearby_decoupled_targets():
    count = 8
    lengths = np.full(count, 0.08)
    seeded = np.zeros((5, count))
    seeded[1] = np.asarray(
        [70.5, -123.0, 0.0, -157.0, -8.33, 106.0, 131.0, 0.1]
    )
    certificate = certify_decoupled_first_order_reachability(
        seeded.reshape(-1),
        lengths,
        3.0,
        maximum_step_m=0.008,
    )
    assert certificate.full_transverse_rank

    rng = np.random.default_rng(20260817)
    for _ in range(2):
        wanted = certificate.parameters.copy()
        wanted[:6] += rng.uniform(-3.0e-4, 3.0e-4, 6)
        target = DecoupledFirstOrderTarget(
            Symplectic2x2KAN(*wanted[:3]),
            Symplectic2x2KAN(*wanted[3:6]),
            wanted[6],
        )
        result = solve_decoupled_first_order_continuation(
            seeded.reshape(-1),
            lengths,
            3.0,
            target,
            maximum_target_step=1.5e-4,
            minimum_target_step=1.0e-7,
            stage_tolerance=2.0e-8,
            final_tolerance=3.0e-8,
            maximum_stage_evaluations=100,
            maximum_step_m=0.008,
        )

        assert result.converged, result.status
        assert result.final_transverse_rank == 6
        assert result.stages
        assert all(stage.accepted for stage in result.stages)
        assert all(stage.transverse_numerical_rank == 6 for stage in result.stages)
        assert result.final_maximum_normalized_parameter_error < 3.0e-8
        assert result.final_matrix_error < 2.0e-7
        assert result.final_decoupling_residual < 2.0e-13
        assert result.final_symplectic_residual < 5.0e-9


def test_second_order_pipeline_selects_skew_xy_material_column():
    orbit = _straight_orbit()
    current = np.zeros(5)
    wanted = np.array([0.0, 0.0, 0.4, 0.0, 1.2])
    objective = _xy_objective(orbit, wanted)
    candidate_delta = np.column_stack(
        (
            wanted,
            np.array([0.0, 0.0, 0.0, 1.2, 0.0]),
        )
    )

    result = run_second_order_taylor_material_inverse_pipeline(
        objective,
        current,
        candidate_elements=np.array([10, 11]),
        candidate_multipole_response_delta=candidate_delta,
        candidate_volumes=np.ones(2),
        volume_budget=2.0,
        candidate_material_active=np.zeros(2, dtype=bool),
        maximum_changed_elements=2,
        field_inverse_relative_tolerance=1.0e-10,
        material_relative_tolerance=1.0e-10,
    )

    assert result.stage_order == (
        "normal-skew-multipole-distribution",
        "forward-ad-second-order-taylor-map",
        "target-R-T-difference",
        "tsvd-minimax-multipole-correction",
        "aca-thin-qr-tsvd-material-inverse",
        "native-exact-R-T-gate",
    )
    assert result.automatic_differentiation.backend == ("forward-mode-rk4-taylor-ad")
    assert result.multipole_correction.derivative_backend == (
        "forward-mode-rk4-taylor-ad"
    )
    np.testing.assert_array_equal(result.material_selection.selected_elements, [10])
    assert result.material_selection.aca_rank >= 1
    assert result.material_selection.numerical_rank >= 1
    assert result.proposed_exact_max_band_ratio < 1.0e-8
    np.testing.assert_allclose(result.proposed_R, objective.target_R)
    np.testing.assert_allclose(result.proposed_T, objective.target_T)


def test_whole_hex_topology_controls_R_and_xy_second_order_T():
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh

    from radia.topology_optimization import solve_hdiv_mmm_active_elements
    from radia.vim._vim import build_charge_gram

    mesh = MakeStructured3DMesh(hexes=True, nx=2, ny=1, nz=1)
    fes = ng.HDiv(mesh, order=1, discontinuous=True)
    with ng.TaskManager():
        _, gram, mass = build_charge_gram(
            fes, eps=1.0e-10, leafsize=256, eta=2.0, internal_interfaces=True
        )
    rng = np.random.default_rng(20260814)
    rhs = np.asarray(mass @ rng.normal(size=fes.ndof))
    initial = np.array([True, False])
    target_active = np.ones(2, dtype=bool)
    zero_response = np.zeros((1, fes.ndof))
    initial_state = solve_hdiv_mmm_active_elements(
        charge_gram=gram,
        fes=fes,
        inv_chi=0.2,
        rhs=rhs,
        response_matrix=zero_response,
        active_elements=initial,
        solve_tolerance=1.0e-11,
    )[0]
    target_state = solve_hdiv_mmm_active_elements(
        charge_gram=gram,
        fes=fes,
        inv_chi=0.2,
        rhs=rhs,
        response_matrix=zero_response,
        active_elements=target_active,
        solve_tolerance=1.0e-11,
    )[0]

    orbit = _straight_orbit()
    wanted = np.array([0.0, 0.0, 0.4, 0.0, 1.2])
    states = np.vstack((initial_state, target_state))
    response_matrix = np.vstack(
        [
            np.linalg.lstsq(states, np.array([0.0, value]), rcond=None)[0]
            for value in wanted
        ]
    )
    objective = _xy_objective(orbit, wanted)
    volumes = np.asarray(ng.Integrate(1.0, mesh, element_wise=True))

    result = optimize_hdiv_mmm_magnet_from_second_order_taylor_map(
        objective,
        charge_gram=gram,
        fes=fes,
        inv_chi=0.2,
        rhs=rhs,
        multipole_response_matrix=response_matrix,
        active_elements=initial,
        element_volumes=volumes,
        volume_max=float(np.sum(volumes)) + 1.0e-14,
        fixed_active_elements=initial,
        maximum_batch_elements=1,
        graph_front_proposal_limit=0,
        max_iterations=1,
        solve_tolerance=1.0e-11,
    )

    assert result.converged
    np.testing.assert_array_equal(result.active_elements, target_active)
    np.testing.assert_allclose(result.realized_multipole_response, wanted, atol=2.0e-13)
    np.testing.assert_allclose(result.realized_R, objective.target_R)
    np.testing.assert_allclose(result.realized_T, objective.target_T)
    assert result.normal_dipole_max_band_ratio < 1.0e-8
    assert result.taylor_map_max_band_ratio < 1.0e-7
    assert result.topology.valid
