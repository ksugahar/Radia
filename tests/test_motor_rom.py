from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from radia.motor_rom import (
    AnglePeriodicMotorROM,
    HysteresisEvaluation,
    LoadAnglePeriodicMotorROM,
    MotorPortContract,
    MotorROMFromHybridVIMSweep,
    MotorROMInput,
    PeriodicAngleTable,
    ProjectRigidRotationMotionalFluxGradient,
)
from radia.motor_rom_export import (
    MotorROMPortManifest,
    SaveMotorROMBundle,
    ValidateMotorROMBundle,
)


def _angles(count=33):
    return np.linspace(0.0, 2.0 * np.pi, count, endpoint=False)


def _motor(*, motion=True, skew=0.0, thermal=False):
    angles = _angles()
    L = np.array(
        [
            [
                [0.02 + 0.004 * np.cos(2.0 * a), 0.001 * np.sin(a)],
                [0.001 * np.sin(a), 0.006],
            ]
            for a in angles
        ]
    )
    R = np.repeat(np.diag([0.4, 2.0])[None, :, :], angles.size, axis=0)
    psi = np.array([[0.03 * np.cos(a), 0.01 * np.sin(a)] for a in angles])
    motion_values = np.array([[0.0, 0.003 * np.cos(a)] for a in angles])
    ports = MotorPortContract(("A",), ("eddy-bulk",))
    return AnglePeriodicMotorROM(
        ports,
        PeriodicAngleTable(angles, L),
        PeriodicAngleTable(angles, R),
        PeriodicAngleTable(angles, psi),
        motion_flux_gradient_Wb_per_rad=(
            PeriodicAngleTable(angles, motion_values) if motion else None
        ),
        inertia_kg_m2=0.02,
        viscous_friction_Nm_s=1.0e-3,
        skew_span_rad=skew,
        resistance_temperature_coefficient_per_K=(0.0039, 0.0039),
        thermal_capacity_J_per_K=100.0 if thermal else None,
        thermal_conductance_W_per_K=2.0 if thermal else 0.0,
    )


def test_periodic_fourier_value_derivative_and_continuous_skew():
    angles = _angles(33)
    values = 2.0 + 0.5 * np.cos(3.0 * angles) - 0.2 * np.sin(2.0 * angles)
    table = PeriodicAngleTable(angles, values)
    theta = 0.371
    expected = 2.0 + 0.5 * np.cos(3.0 * theta) - 0.2 * np.sin(2.0 * theta)
    derivative = -1.5 * np.sin(3.0 * theta) - 0.4 * np.cos(2.0 * theta)
    assert table.evaluate(theta) == pytest.approx(expected, abs=2.0e-14)
    assert table.evaluate(theta, derivative=1) == pytest.approx(derivative, abs=2.0e-13)

    span = 0.4
    skewed = (
        2.0
        + 0.5 * np.sinc(3.0 * span / (2.0 * np.pi)) * np.cos(3.0 * theta)
        - 0.2 * np.sinc(2.0 * span / (2.0 * np.pi)) * np.sin(2.0 * theta)
    )
    assert table.evaluate(theta, skew_span_rad=span) == pytest.approx(skewed, abs=2.0e-14)


def test_periodic_table_requires_odd_uniform_grid_without_duplicate_endpoint():
    with pytest.raises(ValueError, match="odd number"):
        PeriodicAngleTable(_angles(8), np.zeros(8))
    angles = np.linspace(0.0, 2.0 * np.pi, 9)
    with pytest.raises(ValueError, match="without a duplicate endpoint"):
        PeriodicAngleTable(angles, np.zeros(9))


def test_port_contract_zeroes_internal_eddy_voltage():
    ports = MotorPortContract(("A", "B", "C"), ("bulk", "bridge", "sibc"))
    voltage = ports.generalized_voltage((1.0, 2.0, 3.0))
    np.testing.assert_array_equal(voltage, (1.0, 2.0, 3.0, 0.0, 0.0, 0.0))
    assert ports.diagnostics()["n_eddy"] == 3


@dataclass
class _Basis:
    points: np.ndarray
    weights: np.ndarray
    modes: np.ndarray


def test_rigid_rotation_v_cross_b_projection_pairs_voltage_and_torque():
    basis = _Basis(
        points=np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]]),
        weights=np.array([0.5, 0.25]),
        modes=np.array(
            [
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]],
            ]
        ),
    )
    projected = ProjectRigidRotationMotionalFluxGradient(
        basis, (0.0, 0.0, 2.0)
    )
    # v/omega = ez x r; (v/omega) x B = (2*x, 2*y, 0).
    np.testing.assert_allclose(projected, (2.0, 0.0), atol=1.0e-15)
    q = np.array((3.0, -2.0))
    omega = 17.0
    assert q @ (omega * projected) == pytest.approx(omega * (q @ projected))


def test_motion_induces_eddy_current_and_implicit_step_closes_power_balance():
    motor = _motor(motion=True)
    state = motor.initial_state(
        rotor_angle_rad=0.2,
        rotor_speed_rad_s=40.0,
        phase_currents_A=(1.0,),
    )
    command = MotorROMInput(np.array((2.0,)), load_torque_Nm=0.1)
    worst_balance = 0.0
    for _ in range(1000):
        state, output = motor.step(state, command, 2.0e-6)
        worst_balance = max(worst_balance, abs(output.energy_balance_residual_W))
    assert abs(state.generalized_currents_A[1]) > 1.0e-3
    assert worst_balance < 1.0e-7
    assert output.resistive_loss_W >= 0.0
    assert motor.diagnostics()["passive"] is True
    assert motor.diagnostics()["has_motional_v_cross_b"] is True


def test_cogging_coenergy_drives_torque_and_is_in_power_balance():
    base = _motor(motion=False)
    angles = base.inductance_H.angles_rad
    cogging = PeriodicAngleTable(angles, 2.0e-3 * np.cos(6.0 * angles))
    motor = AnglePeriodicMotorROM(
        base.ports,
        base.inductance_H,
        base.resistance_ohm,
        base.pm_flux_linkage_Wb,
        inertia_kg_m2=base.inertia_kg_m2,
        viscous_friction_Nm_s=base.viscous_friction_Nm_s,
        cogging_coenergy_J=cogging,
    )
    theta = 0.23
    expected = -12.0e-3 * np.sin(6.0 * theta)
    assert motor.torque_components(theta, (0.0, 0.0))["cogging"] == pytest.approx(
        expected, abs=2.0e-13
    )
    assert motor.virtual_work_torque(theta, (0.0, 0.0)) == pytest.approx(
        expected, rel=2.0e-9
    )

    state = motor.initial_state(rotor_angle_rad=theta, rotor_speed_rad_s=10.0)
    command = MotorROMInput(np.array((0.0,)))
    worst_balance = 0.0
    for _ in range(100):
        state, output = motor.step(state, command, 1.0e-6)
        worst_balance = max(worst_balance, abs(output.energy_balance_residual_W))
    assert worst_balance < 1.0e-7
    assert motor.diagnostics()["has_cogging_coenergy"] is True


def test_temperature_scaling_thermal_state_and_end_effects_preserve_passivity():
    motor = _motor(motion=False, thermal=True)
    cold = motor.resistance(0.0, 293.15)
    hot = motor.resistance(0.0, 393.15)
    np.testing.assert_allclose(hot, 1.39 * cold, rtol=1.0e-14)
    state = motor.initial_state(phase_currents_A=(10.0,), temperature_K=293.15)
    next_state, output = motor.step(
        state,
        MotorROMInput(np.array((0.0,)), ambient_temperature_K=293.15),
        1.0e-3,
    )
    assert output.resistive_loss_W > 0.0
    assert next_state.temperature_K > state.temperature_K


def test_torque_audit_matches_fourier_derivative_virtual_work_and_maxwell_route():
    motor = _motor(motion=False)
    theta = 0.43
    q = np.array((2.0, -0.7))
    components = motor.torque_components(theta, q)
    maxwell = components["reluctance"] + components["permanent_magnet"]
    audit = motor.torque_audit(theta, q, maxwell_torque_Nm=maxwell)
    assert audit["relative_spread"] < 1.0e-9


class _FunctionalHysteresis:
    def __init__(self, n):
        self.n = n

    def evaluate(self, rotor_angle_rad, generalized_currents_A, committed_state):
        previous = 0 if committed_state is None else int(committed_state)
        q = np.asarray(generalized_currents_A)
        return HysteresisEvaluation(
            flux_linkage_Wb=1.0e-4 * q,
            torque_Nm=1.0e-3 * q[0],
            stored_energy_J=0.5e-4 * float(q @ q),
            dissipated_energy_increment_J=1.0e-8 * (previous + 1),
            state=previous + 1,
        )


def test_functional_hysteresis_state_is_carried_by_time_domain_rom():
    base = _motor(motion=False)
    motor = AnglePeriodicMotorROM(
        base.ports,
        base.inductance_H,
        base.resistance_ohm,
        base.pm_flux_linkage_Wb,
        inertia_kg_m2=base.inertia_kg_m2,
        hysteresis_port=_FunctionalHysteresis(base.ports.n_generalized),
    )
    state = motor.initial_state(phase_currents_A=(1.0,))
    assert state.hysteresis_state == 1
    next_state, output = motor.step(
        state, MotorROMInput(np.array((1.0,))), 1.0e-5
    )
    assert next_state.hysteresis_state > state.hysteresis_state
    assert output.torque_components_Nm["hysteresis"] != 0.0
    assert output.hysteresis_loss_W > 0.0


class _AcceptedStateHysteresis:
    def evaluate(self, rotor_angle_rad, generalized_currents_A, committed_state):
        q = np.asarray(generalized_currents_A, dtype=float).copy()
        return HysteresisEvaluation(
            flux_linkage_Wb=1.0e-4 * q,
            torque_Nm=0.0,
            stored_energy_J=0.5e-4 * float(q @ q),
            dissipated_energy_increment_J=0.0,
            state=(float(rotor_angle_rad), q),
        )


def test_hysteresis_state_is_committed_at_the_accepted_step_state():
    base = _motor(motion=False)
    motor = AnglePeriodicMotorROM(
        base.ports,
        base.inductance_H,
        base.resistance_ohm,
        base.pm_flux_linkage_Wb,
        inertia_kg_m2=base.inertia_kg_m2,
        hysteresis_port=_AcceptedStateHysteresis(),
    )
    state = motor.initial_state(phase_currents_A=(0.0,))
    next_state, _ = motor.step(
        state, MotorROMInput(np.array((100.0,))), 1.0e-3, tolerance=1.0
    )
    accepted_angle, accepted_currents = next_state.hysteresis_state
    assert accepted_angle == pytest.approx(next_state.rotor_angle_rad)
    np.testing.assert_allclose(
        accepted_currents, next_state.generalized_currents_A, atol=0.0, rtol=0.0
    )


def test_hysteresis_nonfinite_or_negative_dissipation_fails_loud():
    class BadHysteresis(_AcceptedStateHysteresis):
        def evaluate(self, rotor_angle_rad, generalized_currents_A, committed_state):
            result = super().evaluate(
                rotor_angle_rad, generalized_currents_A, committed_state
            )
            return HysteresisEvaluation(
                result.flux_linkage_Wb,
                result.torque_Nm,
                result.stored_energy_J,
                -1.0,
                result.state,
            )

    base = _motor(motion=False)
    motor = AnglePeriodicMotorROM(
        base.ports,
        base.inductance_H,
        base.resistance_ohm,
        base.pm_flux_linkage_Wb,
        inertia_kg_m2=base.inertia_kg_m2,
        hysteresis_port=BadHysteresis(),
    )
    with pytest.raises(ValueError, match="must be non-negative"):
        motor.initial_state(phase_currents_A=(0.0,))


@dataclass
class _Hybrid:
    resistance: np.ndarray
    inductance: np.ndarray
    surface_mass: np.ndarray
    blocks: dict[str, tuple[int, int]]

    @property
    def n_modes(self):
        return self.resistance.shape[0]


def test_hybrid_vim_factory_requires_time_domain_cln_when_sibc_is_active():
    angles = _angles()
    systems = tuple(
        _Hybrid(
            resistance=np.array([[2.0]]),
            inductance=np.array([[0.01]]),
            surface_mass=np.array([[1.0]]),
            blocks={"surface": (0, 1)},
        )
        for _ in angles
    )
    ports = MotorPortContract(("A",), ("sibc",))
    kwargs = dict(
        phase_inductance_H=np.array([[0.03]]),
        phase_resistance_ohm=np.array([[0.4]]),
        phase_eddy_mutual_H=np.zeros((1, 1)),
        pm_flux_linkage_Wb=np.zeros(2),
        inertia_kg_m2=0.01,
    )
    with pytest.raises(ValueError, match="positive-real time-domain CLN"):
        MotorROMFromHybridVIMSweep(angles, systems, ports, **kwargs)

    motor = MotorROMFromHybridVIMSweep(
        angles,
        systems,
        ports,
        time_domain_eddy_resistance_ohm=np.array([[2.5]]),
        time_domain_eddy_inductance_H=np.array([[0.015]]),
        **kwargs,
    )
    assert motor.diagnostics()["passive"] is True


def test_npz_export_preserves_platform_neutral_port_and_angle_arrays(tmp_path):
    motor = _motor()
    path = motor.save_npz(tmp_path / "motor_rom.npz")
    data = np.load(path)
    assert data["inductance_H"].shape == (33, 2, 2)
    assert data["phase_names"].tolist() == ["A"]
    assert data["eddy_names"].tolist() == ["eddy-bulk"]

    loaded = LoadAnglePeriodicMotorROM(path)
    for angle in (0.0, 0.37, 2.41):
        np.testing.assert_allclose(loaded.inductance(angle), motor.inductance(angle))
        np.testing.assert_allclose(
            loaded.resistance(angle, 310.0), motor.resistance(angle, 310.0)
        )
        np.testing.assert_allclose(loaded.pm_flux(angle), motor.pm_flux(angle))


def test_simulink_c_abi_and_fmi_boundary_bundle_is_synchronized(tmp_path):
    motor = _motor(thermal=True)
    paths = SaveMotorROMBundle(motor, tmp_path / "motor")
    assert all(Path(path).exists() for path in paths.values())

    manifest = MotorROMPortManifest(motor)
    assert manifest["c_abi"]["abi_version"] == 1
    assert manifest["fmi"]["version"] == "3.0.2"
    assert manifest["fmi"]["packaged_fmu"] is False
    assert manifest["simulink"]["mex_s_function"] == "radia_motor_rom_sfun"
    assert manifest["simulink"]["state_update"] == "internal-discrete-at-fixed-sample-time"
    assert manifest["generalized_current_order"] == ["A", "eddy-bulk"]
    from scipy.io import loadmat

    mat = loadmat(paths["mat"], squeeze_me=True)
    assert float(mat["n_phase"]) == motor.ports.n_phase
    assert float(mat["n_generalized"]) == motor.ports.n_generalized
    assert bool(mat["external_hysteresis_required"]) is False
    assert np.isfinite(float(mat["thermal_capacity_J_per_K"]))
    assert ValidateMotorROMBundle(tmp_path / "motor")["passed"]
