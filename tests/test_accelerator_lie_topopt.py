import json
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from radia.accelerator_lie_topopt import (
    PlanarFourthOrderLieMapObjective,
    PlanarThirdOrderLieMapObjective,
    apply_dragt_finn_map,
    canonical_body_hamiltonian_jet,
    canonical_body_hamiltonian_rhs,
    canonical_poisson_matrix,
    canonical_vector_potential_hamiltonian_rhs,
    dragt_finn_factorize_fourth_order,
    dragt_finn_factorize_third_order,
    formal_fourth_order_symplectic_residual,
    fourth_order_lie_map_from_multipoles,
    optimize_hdiv_mmm_magnet_from_third_order_lie_map,
    planar_orbit_quartic_multipole_observations,
    run_fourth_order_lie_material_inverse_pipeline,
    run_third_order_lie_material_inverse_pipeline,
    third_order_lie_map_from_multipoles,
    track_canonical_hamiltonian_s,
)
from radia.accelerator_lie_topopt import (
    _canonical_vector_potential_hamiltonian_jet as _internal_hcurl_jet,
)
from radia.accelerator_lie_topopt import (
    _fourth_order_lie_map_from_vector_potential_polynomials as _internal_lie_from_jet,
)
from radia.accelerator_magnet_topopt import PlanarDesignOrbit


def _symbolic_reference():
    path = (
        Path(__file__).parents[1]
        / "validation_test"
        / "ffag_topopt"
        / "lie_map_symbolic_reference.json"
    )
    return json.loads(path.read_text(encoding="utf-8"))


def _straight_orbit(length=0.08, rigidity=3.0):
    return PlanarDesignOrbit(
        positions=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, length]]),
        tangents=np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]),
        magnetic_rigidity=rigidity,
        bend_axis=np.array([0.0, 1.0, 0.0]),
    )


def _skew_octupole_lie_objective(orbit, wanted):
    target = third_order_lie_map_from_multipoles(
        wanted,
        orbit.segment_lengths,
        orbit.magnetic_rigidity,
        maximum_step_m=0.002,
    )
    return PlanarThirdOrderLieMapObjective(
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
        maximum_step_m=0.002,
    )


def test_mathematica_fourth_degree_hamiltonian_golden_matches_python():
    reference = _symbolic_reference()
    parameters = reference["parameters"]
    jet = canonical_body_hamiltonian_jet(
        [
            parameters["curvature"],
            parameters["normal_quadrupole"],
            parameters["skew_quadrupole"],
            parameters["normal_sextupole"],
            parameters["skew_sextupole"],
            parameters["normal_octupole"],
            parameters["skew_octupole"],
            parameters["normal_decapole"],
            parameters["skew_decapole"],
        ],
        1.0,
        reference_beta=parameters["reference_beta"],
    )
    entries = reference["hamiltonian_tensor_entries"]
    actual = {
        "H2_xx": jet.H2[0, 0],
        "H2_xy": jet.H2[0, 2],
        "H2_x_delta": jet.H2[0, 5],
        "H2_delta_delta": jet.H2[5, 5],
        "H3_px_px_delta": jet.H3[1, 1, 5],
        "H3_x_px_px": jet.H3[0, 1, 1],
        "H3_x_x_x": jet.H3[0, 0, 0],
        "H3_x_x_y": jet.H3[0, 0, 2],
        "H4_px_px_px_px": jet.H4[1, 1, 1, 1],
        "H4_px_px_py_py": jet.H4[1, 1, 3, 3],
        "H4_x_px_px_delta": jet.H4[0, 1, 1, 5],
        "H4_x_x_x_x": jet.H4[0, 0, 0, 0],
        "H4_x_x_x_y": jet.H4[0, 0, 0, 2],
        "H4_delta_delta_delta_delta": jet.H4[5, 5, 5, 5],
        "H5_px_px_delta_delta_delta": jet.H5[1, 1, 5, 5, 5],
        "H5_px_px_px_px_delta": jet.H5[1, 1, 1, 1, 5],
        "H5_x_px_px_px_px": jet.H5[0, 1, 1, 1, 1],
        "H5_px_px_py_py_delta": jet.H5[1, 1, 3, 3, 5],
        "H5_x_x_x_x_x": jet.H5[0, 0, 0, 0, 0],
        "H5_x_x_x_x_y": jet.H5[0, 0, 0, 0, 2],
        "H5_delta_delta_delta_delta_delta": jet.H5[5, 5, 5, 5, 5],
    }
    for name, expected in entries.items():
        assert actual[name] == pytest.approx(expected, abs=3.0e-15)

    generator_entries = reference["linear_generator_entries"]
    actual_generator = {
        "A_x_px": jet.A[0, 1],
        "A_px_x": jet.A[1, 0],
        "A_px_y": jet.A[1, 2],
        "A_px_delta": jet.A[1, 5],
        "A_y_py": jet.A[2, 3],
        "A_py_x": jet.A[3, 0],
        "A_py_y": jet.A[3, 2],
        "A_ell_x": jet.A[4, 0],
        "A_ell_delta": jet.A[4, 5],
    }
    for name, expected in generator_entries.items():
        assert actual_generator[name] == pytest.approx(expected, abs=3.0e-15)


def test_design_orbit_curvature_is_independent_of_dipole_field_curvature():
    coefficients = np.array([0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    rigidity = 3.0
    field_curvature = coefficients[0] / rigidity
    design_curvature = 0.25
    jet = canonical_body_hamiltonian_jet(
        coefficients,
        rigidity,
        reference_curvature_per_m=design_curvature,
    )
    rhs = canonical_body_hamiltonian_rhs(
        np.zeros(6),
        coefficients,
        rigidity,
        reference_curvature_per_m=design_curvature,
    )

    assert jet.H2[0, 5] == pytest.approx(-design_curvature, abs=2.0e-15)
    assert jet.H2[0, 0] == pytest.approx(
        design_curvature * field_curvature, abs=2.0e-15
    )
    assert rhs[1] == pytest.approx(
        design_curvature - field_curvature, abs=2.0e-15
    )
    np.testing.assert_array_equal(rhs[[0, 2, 3, 4, 5]], 0.0)


def test_direct_As_Ay_hamiltonian_reduces_to_body_multipole_rhs():
    coefficients = np.array(
        [0.42, 1.3, -0.4, 7.0, 2.0, -15.0, 6.0, 25.0, -9.0]
    )
    rigidity = 2.8
    curvature = 0.17
    beta = 0.91
    state = np.array([2.0e-3, 1.5e-2, -1.2e-3, -8.0e-3, 0.0, 3.0e-2])
    x_value, y_value = state[0], state[2]
    coordinate = complex(x_value, y_value)
    field_curvature = coefficients[0] / rigidity
    potential_x = field_curvature + curvature * field_curvature * x_value
    potential_y = 0.0
    for order in range(1, 5):
        coefficient = complex(
            coefficients[2 * order - 1], coefficients[2 * order]
        ) / rigidity
        field = coefficient * coordinate**order
        potential_x += field.real
        potential_y -= field.imag
    gradient = np.zeros((3, 2))
    gradient[2] = [-potential_x, -potential_y]

    direct = canonical_vector_potential_hamiltonian_rhs(
        state,
        np.zeros(3),
        gradient,
        reference_curvature_per_m=curvature,
        reference_beta=beta,
    )
    body = canonical_body_hamiltonian_rhs(
        state,
        coefficients,
        rigidity,
        reference_curvature_per_m=curvature,
        reference_beta=beta,
    )

    np.testing.assert_allclose(direct, body, rtol=2.0e-15, atol=2.0e-15)


def test_direct_As_Ay_polynomial_jet_has_forward_ad_and_canonical_coupling():
    rigidity = 2.0
    Ay = np.zeros((3, 3))
    As = np.zeros((3, 3))
    Ay[1, 0] = 0.7
    normal_gradient = 1.6
    As[2, 0] = -0.5 * normal_gradient

    result = _internal_hcurl_jet(
        Ay,
        As,
        rigidity,
        reference_curvature_per_m=0.0,
        longitudinal_component="physical",
    )
    jet = result.jet
    normalized_Ay = Ay[1, 0] / rigidity

    assert result.parameter_names[:5] == (
        "Ay_x1_y0",
        "Ay_x0_y1",
        "Ay_x2_y0",
        "Ay_x1_y1",
        "Ay_x0_y2",
    )
    assert result.constant == pytest.approx(0.0, abs=2.0e-15)
    np.testing.assert_allclose(result.linear, 0.0, atol=2.0e-15)
    assert jet.H2[0, 3] == pytest.approx(-normalized_Ay, abs=2.0e-15)
    assert jet.H2[0, 0] == pytest.approx(
        normalized_Ay**2 + normal_gradient / rigidity,
        abs=2.0e-15,
    )
    assert jet.H2[1, 1] == pytest.approx(1.0, abs=2.0e-15)
    assert jet.H2[3, 3] == pytest.approx(1.0, abs=2.0e-15)
    assert jet.H2_jacobian[0, 0, 3] == pytest.approx(
        -1.0 / rigidity, abs=2.0e-15
    )
    assert jet.H2_jacobian[0, 0, 0] == pytest.approx(
        2.0 * normalized_Ay / rigidity, abs=2.0e-15
    )


def test_direct_As_Ay_polynomials_generate_fourth_order_lie_map_with_ad():
    rigidity = 2.0
    length = 0.06
    Ay = np.zeros((1, 3, 3))
    As = np.zeros((1, 3, 3))
    Ay[0, 1, 0] = 0.7
    As[0, 2, 0] = -0.8
    result = _internal_lie_from_jet(
        Ay,
        As,
        [length],
        rigidity,
        reference_curvature_per_m=[0.0],
        maximum_step_m=5.0e-4,
    )
    expected_jet = _internal_hcurl_jet(
        Ay[0],
        As[0],
        rigidity,
        reference_curvature_per_m=0.0,
    ).jet

    np.testing.assert_allclose(
        result.transfer.raw_R,
        expm(expected_jet.A * length),
        rtol=2.0e-12,
        atol=2.0e-13,
    )
    assert result.transfer.factorization.reconstructed_symplectic_residual.maximum < (
        2.0e-12
    )
    assert result.transfer.V.shape == (6, 6, 6, 6, 6)
    assert np.max(np.abs(result.transfer.T)) > 0.0
    assert np.max(np.abs(result.transfer.U)) > 0.0

    step = 2.0e-5
    plus = Ay.copy()
    minus = Ay.copy()
    plus[0, 1, 0] += step
    minus[0, 1, 0] -= step
    plus_map = _internal_lie_from_jet(
        plus,
        As,
        [length],
        rigidity,
        reference_curvature_per_m=[0.0],
        maximum_step_m=5.0e-4,
    )
    minus_map = _internal_lie_from_jet(
        minus,
        As,
        [length],
        rigidity,
        reference_curvature_per_m=[0.0],
        maximum_step_m=5.0e-4,
    )
    regression = (plus_map.transfer.R - minus_map.transfer.R) / (2.0 * step)
    np.testing.assert_allclose(
        result.transfer.R_jacobian[0],
        regression,
        rtol=1.0e-6,
        atol=5.0e-10,
    )


def test_direct_A_lie_map_rejects_a_nonreference_design_orbit():
    Ay = np.zeros((1, 3, 3))
    As = np.zeros_like(Ay)
    As[0, 1, 0] = 0.2

    with pytest.raises(ValueError, match="reference orbit is not"):
        _internal_lie_from_jet(
            Ay,
            As,
            [0.05],
            2.0,
            reference_curvature_per_m=[0.0],
            reference_orbit_tolerance=1.0e-10,
        )


def test_design_orbit_exposes_continuous_global_xyz_frame_and_h_of_s():
    radius = 2.4
    angles = np.linspace(0.0, 0.3, 7)
    positions = np.column_stack(
        (radius * np.sin(angles), np.zeros_like(angles), radius * np.cos(angles))
    )
    tangents = np.column_stack(
        (np.cos(angles), np.zeros_like(angles), -np.sin(angles))
    )
    stations = radius * angles
    orbit = PlanarDesignOrbit(
        positions,
        tangents,
        magnetic_rigidity=1.8,
        bend_axis=np.array([0.0, 1.0, 0.0]),
        path_length_stations=stations,
    )

    np.testing.assert_allclose(orbit.arc_length_stations, stations, atol=0.0)
    np.testing.assert_allclose(orbit.position_at(stations), positions, atol=2.0e-15)
    np.testing.assert_allclose(orbit.tangent_at(stations), tangents, atol=2.0e-15)
    np.testing.assert_allclose(
        orbit.signed_curvature_at(0.5 * (stations[:-1] + stations[1:])),
        1.0 / radius,
        atol=2.0e-15,
    )
    horizontal, vertical, evaluated_tangent = orbit.frame_at(stations)
    np.testing.assert_allclose(evaluated_tangent, tangents, atol=2.0e-15)
    np.testing.assert_allclose(
        vertical, np.broadcast_to(orbit.bend_axis, vertical.shape), atol=0.0
    )
    np.testing.assert_allclose(
        np.einsum("ij,ij->i", horizontal, evaluated_tangent), 0.0, atol=2.0e-15
    )
    np.testing.assert_allclose(
        orbit.local_to_global(stations, np.zeros_like(stations), np.zeros_like(stations)),
        positions,
        atol=2.0e-15,
    )


def test_s_runge_kutta_keeps_matched_uniform_field_design_orbit():
    radius = 3.2
    rigidity = 1.7
    angles = np.linspace(0.0, 0.25, 9)
    stations = radius * angles
    orbit = PlanarDesignOrbit(
        positions=np.column_stack(
            (
                radius * np.sin(angles),
                np.zeros_like(angles),
                radius * np.cos(angles),
            )
        ),
        tangents=np.column_stack(
            (np.cos(angles), np.zeros_like(angles), -np.sin(angles))
        ),
        magnetic_rigidity=rigidity,
        bend_axis=np.array([0.0, 1.0, 0.0]),
        path_length_stations=stations,
    )
    coefficients = np.zeros((9, len(orbit.segment_lengths)))
    coefficients[0] = rigidity / radius

    fixed = track_canonical_hamiltonian_s(
        orbit,
        coefficients.reshape(-1),
        np.zeros(6),
        integrator="RK4",
        maximum_step_m=0.004,
    )
    adaptive = track_canonical_hamiltonian_s(
        orbit,
        coefficients.reshape(-1),
        np.zeros(6),
        integrator="DOP853",
        maximum_step_m=0.004,
    )

    assert fixed.accepted_steps > len(orbit.segment_lengths)
    assert adaptive.accepted_steps > len(orbit.segment_lengths)
    np.testing.assert_allclose(fixed.canonical_states, 0.0, atol=2.0e-15)
    np.testing.assert_allclose(adaptive.canonical_states, 0.0, atol=2.0e-15)
    np.testing.assert_allclose(
        fixed.global_positions_m, fixed.reference_positions_m, atol=2.0e-15
    )
    np.testing.assert_allclose(
        fixed.reference_curvature_per_m, 1.0 / radius, atol=2.0e-15
    )


def test_dragt_finn_recovers_mathematica_self_cascade_and_quartic_kick():
    sample = _symbolic_reference()["lie_sample"]
    a = sample["f3_a"]
    b = sample["f4_b"]
    R = np.eye(6)
    T = np.zeros((6, 6, 6))
    U = np.zeros((6, 6, 6, 6))
    T[0, 0, 0] = a
    T[1, 0, 1] = T[1, 1, 0] = -a
    U[0, 0, 0, 0] = 1.5 * a * a
    U[1, 0, 0, 1] = U[1, 0, 1, 0] = U[1, 1, 0, 0] = 0.5 * a * a
    U[1, 0, 0, 0] = -6.0 * b

    factors = dragt_finn_factorize_third_order(R, T, U)

    assert factors.f3[0, 0, 1] == pytest.approx(a)
    assert factors.f4[0, 0, 0, 0] == pytest.approx(6.0 * b)
    assert factors.relative_reconstruction_error < 1.0e-15
    assert factors.reconstructed_symplectic_residual.maximum < 2.0e-16
    q_coefficients = sample["q_map_coefficients"]
    p_coefficients = sample["p_map_coefficients"]
    assert 0.5 * factors.T[0, 0, 0] == pytest.approx(q_coefficients["q2"])
    assert factors.T[1, 0, 1] == pytest.approx(p_coefficients["q_p"])
    assert factors.U[0, 0, 0, 0] / 6.0 == pytest.approx(q_coefficients["q3"])
    assert factors.U[1, 0, 0, 1] / 2.0 == pytest.approx(p_coefficients["q2_p"])
    assert factors.U[1, 0, 0, 0] / 6.0 == pytest.approx(p_coefficients["q3"])

    V = np.zeros((6, 6, 6, 6, 6))
    V[0, 0, 0, 0, 0] = 24.0 * q_coefficients["q4"]
    V[1, 0, 0, 0, 0] = 24.0 * p_coefficients["q4"]
    fourth = dragt_finn_factorize_fourth_order(R, T, U, V)
    assert fourth.f5[0, 0, 0, 0, 0] == pytest.approx(24.0 * sample["f5_c"])
    assert fourth.V[0, 0, 0, 0, 0] / 24.0 == pytest.approx(
        q_coefficients["q4"]
    )
    assert fourth.V[1, 0, 0, 0, 0] / 24.0 == pytest.approx(
        p_coefficients["q4"]
    )
    assert fourth.relative_reconstruction_error < 1.0e-14
    assert fourth.reconstructed_symplectic_residual.maximum < 3.0e-16


def test_factorization_reports_and_repairs_a_non_hamiltonian_quadratic_map():
    T = np.zeros((6, 6, 6))
    T[0, 0, 0] = 1.0
    factors = dragt_finn_factorize_third_order(np.eye(6), T, np.zeros((6, 6, 6, 6)))

    assert factors.f3_symmetry_defect > 0.1
    assert factors.relative_reconstruction_error > 0.1
    assert factors.raw_symplectic_residual.linear > 0.1
    assert factors.reconstructed_symplectic_residual.maximum < 2.0e-16


def test_finite_amplitude_dragt_finn_application_is_symplectic():
    raw = np.array([0.2, 1.0, 0.3, 2.0, 0.7, 4.0, -3.0])
    transfer = third_order_lie_map_from_multipoles(
        raw,
        [0.07],
        3.0,
        reference_beta=0.8,
        maximum_step_m=0.002,
    )
    state = np.array([1.0e-3, 2.0e-4, -7.0e-4, 1.0e-4, 2.0e-3, 1.0e-3])

    applied, jacobian = apply_dragt_finn_map(
        transfer.factorization, state, return_jacobian=True
    )
    polynomial = (
        transfer.R @ state
        + 0.5 * np.einsum("ijk,j,k->i", transfer.T, state, state)
        + np.einsum("ijkl,j,k,l->i", transfer.U, state, state, state) / 6.0
    )
    poisson = canonical_poisson_matrix()

    np.testing.assert_allclose(applied, polynomial, atol=1.2e-14, rtol=0.0)
    np.testing.assert_allclose(
        jacobian.T @ poisson @ jacobian, poisson, atol=3.0e-15, rtol=0.0
    )
    assert transfer.factorization.reconstructed_symplectic_residual.maximum < 3e-15


def test_lie_R_T_U_and_generators_forward_ad_match_centered_difference():
    raw = np.array([0.2, 1.0, 0.3, 2.0, 0.7, 4.0, -3.0])
    options = {
        "segment_lengths": [0.07],
        "magnetic_rigidity": 3.0,
        "reference_beta": 0.8,
        "maximum_step_m": 0.002,
    }
    differentiated = third_order_lie_map_from_multipoles(raw, **options)

    for name, jacobian in (
        ("R", differentiated.R_jacobian),
        ("T", differentiated.T_jacobian),
        ("U", differentiated.U_jacobian),
        ("f3", differentiated.f3_jacobian),
        ("f4", differentiated.f4_jacobian),
    ):
        for parameter in range(raw.size):
            step = 1.0e-6 * max(1.0, abs(raw[parameter]))
            direction = np.zeros_like(raw)
            direction[parameter] = step
            plus = third_order_lie_map_from_multipoles(raw + direction, **options)
            minus = third_order_lie_map_from_multipoles(raw - direction, **options)
            finite_difference = (getattr(plus, name) - getattr(minus, name)) / (
                2.0 * step
            )
            np.testing.assert_allclose(
                jacobian[parameter],
                finite_difference,
                rtol=8.0e-6,
                atol=8.0e-9,
            )


def test_truncated_lie_map_has_fourth_order_error_against_exact_hamiltonian():
    raw = np.array([0.2, 1.0, 0.3, 2.0, 0.7, 4.0, -3.0])
    length = 0.2
    transfer = third_order_lie_map_from_multipoles(
        raw,
        [length],
        3.0,
        reference_beta=0.8,
        maximum_step_m=0.001,
    )
    base = np.array([0.7, 0.2, -0.5, 0.1, 0.3, 0.15])
    errors = []
    for scale in (0.04, 0.02, 0.01):
        initial = scale * base
        exact = solve_ivp(
            lambda _, state: canonical_body_hamiltonian_rhs(
                state, raw, 3.0, reference_beta=0.8
            ),
            (0.0, length),
            initial,
            method="DOP853",
            rtol=2.0e-13,
            atol=2.0e-15,
        ).y[:, -1]
        polynomial = (
            transfer.R @ initial
            + 0.5 * np.einsum("ijk,j,k->i", transfer.T, initial, initial)
            + np.einsum("ijkl,j,k,l->i", transfer.U, initial, initial, initial) / 6.0
        )
        errors.append(float(np.max(np.abs(exact - polynomial))))

    assert errors[0] / errors[1] == pytest.approx(16.0, rel=0.03)
    assert errors[1] / errors[2] == pytest.approx(16.0, rel=0.03)


def test_lie_objective_rejects_non_symplectic_target():
    orbit = _straight_orbit()
    zero = third_order_lie_map_from_multipoles(
        np.zeros(7), orbit.segment_lengths, orbit.magnetic_rigidity
    )
    invalid_T = zero.T.copy()
    invalid_T[0, 0, 0] += 0.1

    with pytest.raises(ValueError, match="not formally symplectic"):
        PlanarThirdOrderLieMapObjective(
            orbit=orbit,
            target_R=zero.R,
            target_T=invalid_T,
            target_U=zero.U,
            R_band=1.0,
            T_band=1.0,
            U_band=1.0,
            normal_dipole_band=1.0,
            R_entries=((0, 0),),
            T_entries=((0, 0, 0),),
            U_entries=((0, 0, 0, 0),),
        )


def test_lie_material_pipeline_selects_skew_octupole_column():
    orbit = _straight_orbit()
    current = np.zeros(7)
    wanted = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0])
    objective = _skew_octupole_lie_objective(orbit, wanted)
    candidates = np.column_stack(
        (wanted, np.array([0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0]))
    )

    result = run_third_order_lie_material_inverse_pipeline(
        objective,
        current,
        candidate_elements=np.array([10, 11]),
        candidate_multipole_response_delta=candidates,
        candidate_volumes=np.ones(2),
        volume_budget=2.0,
        candidate_material_active=np.zeros(2, dtype=bool),
        maximum_changed_elements=2,
        field_inverse_relative_tolerance=1.0e-10,
        material_relative_tolerance=1.0e-10,
    )

    assert result.automatic_differentiation.backend == (
        "forward-mode-hamiltonian-lie-ad"
    )
    np.testing.assert_array_equal(result.material_selection.selected_elements, [10])
    assert result.proposed_exact_max_band_ratio < 1.0e-8


def test_whole_hex_topology_controls_canonical_skew_octupole_lie_term():
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
    rhs = np.asarray(mass @ np.random.default_rng(20260815).normal(size=fes.ndof))
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
    objective = _skew_octupole_lie_objective(orbit, wanted)
    volumes = np.asarray(ng.Integrate(1.0, mesh, element_wise=True))

    result = optimize_hdiv_mmm_magnet_from_third_order_lie_map(
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
    np.testing.assert_allclose(result.realized_multipole_response, wanted, atol=2e-13)
    assert result.taylor_map_max_band_ratio < 1.0e-7
    assert result.topology.valid


def test_fourth_order_lie_map_f5_and_forward_ad_match_centered_difference():
    raw = np.array([0.2, 1.0, 0.3, 2.0, 0.7, 4.0, -3.0, 5.0, -2.0])
    options = {
        "segment_lengths": [0.02],
        "magnetic_rigidity": 3.0,
        "reference_beta": 0.8,
        "maximum_step_m": 0.004,
    }
    differentiated = fourth_order_lie_map_from_multipoles(raw, **options)

    assert np.max(np.abs(differentiated.V)) > 1.0e-5
    assert np.max(np.abs(differentiated.f5)) > 1.0e-5
    assert differentiated.factorization.maximum_generator_symmetry_defect < 2e-10
    assert differentiated.factorization.relative_reconstruction_error < 2e-10
    assert (
        differentiated.factorization.reconstructed_symplectic_residual.maximum
        < 5e-15
    )
    for parameter in (0, 3, 6, 8):
        step = 1.0e-5 * max(1.0, abs(raw[parameter]))
        direction = np.zeros_like(raw)
        direction[parameter] = step
        plus = fourth_order_lie_map_from_multipoles(raw + direction, **options)
        minus = fourth_order_lie_map_from_multipoles(raw - direction, **options)
        for name in ("V", "f5"):
            finite_difference = (getattr(plus, name) - getattr(minus, name)) / (
                2.0 * step
            )
            np.testing.assert_allclose(
                getattr(differentiated, name + "_jacobian")[parameter],
                finite_difference,
                rtol=1.0e-6,
                atol=2.0e-10,
            )


def test_fourth_order_dragt_finn_application_and_formal_gate_are_symplectic():
    raw = np.array([0.2, 1.0, 0.3, 2.0, 0.7, 4.0, -3.0, 5.0, -2.0])
    transfer = fourth_order_lie_map_from_multipoles(
        raw, [0.03], 3.0, reference_beta=0.8, maximum_step_m=0.003
    )
    state = np.array([8e-4, 2e-4, -5e-4, 1e-4, 1e-3, 7e-4])
    applied, jacobian = apply_dragt_finn_map(
        transfer.factorization, state, return_jacobian=True
    )
    polynomial = (
        transfer.R @ state
        + 0.5 * np.einsum("ijk,j,k->i", transfer.T, state, state)
        + np.einsum("ijkl,j,k,l->i", transfer.U, state, state, state) / 6.0
        + np.einsum(
            "ijklm,j,k,l,m->i", transfer.V, state, state, state, state
        )
        / 24.0
    )
    poisson = canonical_poisson_matrix()

    np.testing.assert_allclose(applied, polynomial, atol=3.0e-14, rtol=0.0)
    np.testing.assert_allclose(
        jacobian.T @ poisson @ jacobian, poisson, atol=3.0e-15, rtol=0.0
    )
    damaged = transfer.V.copy()
    damaged[0, 0, 0, 0, 0] += 0.1
    assert formal_fourth_order_symplectic_residual(
        transfer.R, transfer.T, transfer.U, damaged
    ).cubic > 1.0e-3


def test_fourth_order_map_has_fifth_order_error_against_exact_hamiltonian():
    raw = np.array([0.2, 1.0, 0.3, 2.0, 0.7, 4.0, -3.0, 5.0, -2.0])
    length = 0.08
    transfer = fourth_order_lie_map_from_multipoles(
        raw,
        [length],
        3.0,
        reference_beta=0.8,
        maximum_step_m=0.002,
    )
    base = np.array([0.7, 0.2, -0.5, 0.1, 0.3, 0.15])
    errors = []
    for scale in (0.04, 0.02, 0.01):
        initial = scale * base
        exact = solve_ivp(
            lambda _, state: canonical_body_hamiltonian_rhs(
                state, raw, 3.0, reference_beta=0.8
            ),
            (0.0, length),
            initial,
            method="DOP853",
            rtol=2.0e-13,
            atol=2.0e-15,
        ).y[:, -1]
        polynomial = (
            transfer.R @ initial
            + 0.5 * np.einsum("ijk,j,k->i", transfer.T, initial, initial)
            + np.einsum("ijkl,j,k,l->i", transfer.U, initial, initial, initial)
            / 6.0
            + np.einsum(
                "ijklm,j,k,l,m->i",
                transfer.V,
                initial,
                initial,
                initial,
                initial,
            )
            / 24.0
        )
        errors.append(float(np.max(np.abs(exact - polynomial))))

    assert errors[0] / errors[1] == pytest.approx(32.0, rel=0.06)
    assert errors[1] / errors[2] == pytest.approx(32.0, rel=0.06)


def test_fourth_order_lie_objective_and_material_pipeline_control_v_term():
    orbit = _straight_orbit(length=0.02)
    current = np.zeros(9)
    wanted = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0])
    target = fourth_order_lie_map_from_multipoles(
        wanted, orbit.segment_lengths, orbit.magnetic_rigidity, maximum_step_m=0.004
    )
    objective = PlanarFourthOrderLieMapObjective(
        orbit=orbit,
        target_R=target.R,
        target_T=target.T,
        target_U=target.U,
        target_V=target.V,
        R_band=1.0e-6,
        T_band=1.0e-6,
        U_band=1.0e-6,
        V_band=1.0e-6,
        normal_dipole_band=1.0e-6,
        R_entries=((0, 0),),
        T_entries=((0, 1, 5),),
        U_entries=((1, 0, 0, 2),),
        V_entries=((1, 0, 0, 0, 2),),
        maximum_step_m=0.004,
    )
    candidates = np.column_stack(
        (
            wanted,
            np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0]),
        )
    )
    result = run_fourth_order_lie_material_inverse_pipeline(
        objective,
        current,
        candidate_elements=np.array([10, 11]),
        candidate_multipole_response_delta=candidates,
        candidate_volumes=np.ones(2),
        volume_budget=2.0,
        candidate_material_active=np.zeros(2, dtype=bool),
        maximum_changed_elements=2,
        field_inverse_relative_tolerance=1.0e-10,
        material_relative_tolerance=1.0e-10,
    )

    assert result.automatic_differentiation.backend == (
        "forward-mode-fourth-order-hamiltonian-lie-ad"
    )
    np.testing.assert_array_equal(result.material_selection.selected_elements, [10])
    assert result.proposed_exact_max_band_ratio < 1.0e-7
    np.testing.assert_allclose(result.proposed_V, target.V, atol=2.0e-12)


def test_quartic_orbit_observations_recover_normal_and_skew_decapole():
    orbit = _straight_orbit(length=0.02)
    points, weights = planar_orbit_quartic_multipole_observations(
        orbit, sample_radius=0.004
    )
    coordinate = points[:, 0] + 1.0j * points[:, 1]
    coefficient = 5.0 - 2.0j
    field_complex = coefficient * coordinate**4
    field = np.column_stack(
        (field_complex.imag, field_complex.real, np.zeros(points.shape[0]))
    )
    response = np.einsum("rpc,pc->r", weights, field)

    assert points.shape == (9, 3)
    assert response.shape == (9,)
    np.testing.assert_allclose(response[:7], 0.0, atol=2.0e-12)
    assert response[7] == pytest.approx(coefficient.real, abs=2.0e-12)
    assert response[8] == pytest.approx(coefficient.imag, abs=2.0e-12)
