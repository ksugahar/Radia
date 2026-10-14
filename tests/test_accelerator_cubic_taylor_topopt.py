from itertools import permutations

import numpy as np

from radia.accelerator_magnet_topopt import PlanarDesignOrbit
from radia.accelerator_taylor_topopt import (
    THIRD_ORDER_MULTIPOLE_COMPONENTS,
    PlanarThirdOrderTaylorMapObjective,
    certify_taylor_map_reachability,
    optimize_hdiv_mmm_magnet_from_third_order_taylor_map,
    planar_orbit_cubic_multipole_observations,
    run_third_order_taylor_material_inverse_pipeline,
    third_order_taylor_map_from_multipoles,
)


def _straight_orbit(length=0.08, rigidity=3.0):
    return PlanarDesignOrbit(
        positions=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, length]]),
        tangents=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]),
        magnetic_rigidity=rigidity,
        bend_axis=np.array([0.0, 1.0, 0.0]),
    )


def _skew_octupole_objective(orbit, wanted):
    target = third_order_taylor_map_from_multipoles(
        wanted,
        orbit.segment_lengths,
        orbit.magnetic_rigidity,
        maximum_step_m=0.01,
    )
    return PlanarThirdOrderTaylorMapObjective(
        orbit=orbit,
        target_R=target.R,
        target_T=target.T,
        target_U=target.U,
        R_band=1.0e-6,
        T_band=1.0e-6,
        U_band=1.0e-6,
        normal_dipole_band=1.0e-6,
        R_entries=((0, 0),),
        T_entries=((0, 1, 5),),
        U_entries=((1, 0, 0, 2),),
        maximum_step_m=0.01,
    )


def test_nine_point_harmonic_rows_recover_octupoles_too():
    orbit = _straight_orbit()
    points, weights = planar_orbit_cubic_multipole_observations(
        orbit, sample_radius=0.01
    )
    coefficients = np.array([0.2, 1.1, -0.4, 2.3, 0.7, -4.2, 1.6])
    coordinate = points[:, 0] + 1.0j * points[:, 1]
    complex_field = coefficients[0]
    for degree in range(1, 4):
        complex_field += (
            complex(coefficients[2 * degree - 1], coefficients[2 * degree])
            * coordinate**degree
        )
    field = np.column_stack(
        (complex_field.imag, complex_field.real, np.zeros(points.shape[0]))
    )
    recovered = np.einsum("rpd,pd->r", weights, field)

    assert THIRD_ORDER_MULTIPOLE_COMPONENTS[-2:] == (
        "normal_octupole",
        "skew_octupole",
    )
    np.testing.assert_allclose(recovered, coefficients, atol=2.0e-11)


def test_third_order_R_T_U_forward_ad_matches_centered_difference():
    raw = np.array([0.2, 1.0, 0.3, 2.0, 0.7, 4.0, -3.0])
    differentiated = third_order_taylor_map_from_multipoles(
        raw, [0.07], 3.0, maximum_step_m=0.01
    )

    assert differentiated.value_backend == "native-cpp-variational-map"
    assert abs(differentiated.U[1, 0, 0, 2]) > 1.0e-3
    assert abs(differentiated.U[0, 1, 5, 5]) > 1.0e-3
    np.testing.assert_allclose(
        differentiated.U,
        np.transpose(differentiated.U, (0, 2, 1, 3)),
        atol=3.0e-15,
    )
    np.testing.assert_allclose(
        differentiated.U,
        np.transpose(differentiated.U, (0, 1, 3, 2)),
        atol=3.0e-15,
    )

    for parameter in range(raw.size):
        step = 1.0e-6 * max(1.0, abs(raw[parameter]))
        direction = np.zeros_like(raw)
        direction[parameter] = step
        plus = third_order_taylor_map_from_multipoles(
            raw + direction, [0.07], 3.0, maximum_step_m=0.01
        )
        minus = third_order_taylor_map_from_multipoles(
            raw - direction, [0.07], 3.0, maximum_step_m=0.01
        )
        np.testing.assert_allclose(
            differentiated.R_jacobian[parameter],
            (plus.R - minus.R) / (2.0 * step),
            rtol=3.0e-7,
            atol=7.0e-10,
        )
        np.testing.assert_allclose(
            differentiated.T_jacobian[parameter],
            (plus.T - minus.T) / (2.0 * step),
            rtol=3.0e-7,
            atol=7.0e-10,
        )
        np.testing.assert_allclose(
            differentiated.U_jacobian[parameter],
            (plus.U - minus.U) / (2.0 * step),
            rtol=4.0e-7,
            atol=8.0e-10,
        )


def test_sextupole_cascade_generates_delta_free_cubic_map_term():
    sextupole_only = np.array([0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0])
    transfer = third_order_taylor_map_from_multipoles(
        sextupole_only, [0.08], 3.0, maximum_step_m=0.01
    )

    assert transfer.U[1, 0, 0, 0] > 4.0e-4


def test_reachability_certificate_identifies_uncontrollable_sigma_input_term():
    orbit = _straight_orbit()
    current = np.zeros(7)
    realized = third_order_taylor_map_from_multipoles(
        current, orbit.segment_lengths, orbit.magnetic_rigidity, maximum_step_m=0.01
    )
    target_U = realized.U.copy()
    for inputs in set(permutations((0, 0, 4))):
        target_U[(1,) + inputs] = 0.1
    objective = PlanarThirdOrderTaylorMapObjective(
        orbit=orbit,
        target_R=realized.R,
        target_T=realized.T,
        target_U=target_U,
        R_band=1.0e-3,
        T_band=1.0e-3,
        U_band=1.0e-3,
        normal_dipole_band=1.0e-3,
        R_entries=((0, 0),),
        T_entries=((0, 1, 5),),
        U_entries=((1, 0, 0, 4),),
        maximum_step_m=0.01,
    )

    certificate = certify_taylor_map_reachability(objective, current)

    assert not certificate.linearized_reachable
    assert certificate.max_unreachable_band_ratio == 100.0
    assert dict(certificate.component_max_unreachable_band_ratios)["U"] == 100.0
    np.testing.assert_allclose(certificate.parameter_step, 0.0, atol=1.0e-15)


def test_third_order_pipeline_selects_skew_octupole_material_column():
    orbit = _straight_orbit()
    current = np.zeros(7)
    wanted = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0])
    objective = _skew_octupole_objective(orbit, wanted)
    candidate_delta = np.column_stack(
        (wanted, np.array([0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0]))
    )

    result = run_third_order_taylor_material_inverse_pipeline(
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

    assert result.reachability.linearized_reachable
    np.testing.assert_array_equal(result.material_selection.selected_elements, [10])
    assert result.proposed_exact_max_band_ratio < 1.0e-8
    np.testing.assert_allclose(result.proposed_U, objective.target_U)


def test_whole_hex_topology_controls_skew_octupole_U_term():
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
    wanted = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0])
    states = np.vstack((initial_state, target_state))
    response_matrix = np.vstack(
        [
            np.linalg.lstsq(states, np.array([0.0, value]), rcond=None)[0]
            for value in wanted
        ]
    )
    objective = _skew_octupole_objective(orbit, wanted)
    volumes = np.asarray(ng.Integrate(1.0, mesh, element_wise=True))

    result = optimize_hdiv_mmm_magnet_from_third_order_taylor_map(
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
    np.testing.assert_allclose(result.realized_U, objective.target_U)
    assert result.taylor_map_max_band_ratio < 1.0e-7
    assert result.topology.valid
