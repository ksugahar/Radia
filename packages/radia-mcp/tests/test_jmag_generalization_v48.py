from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.motor_semantic_identity_v48 import PWM, SKEW, validate_public_identity


PROMOTED_CASE_IDS = {
    "v48_public_skew_slice_weight_rotor_angle_torque_harmonic_phase_aggregation_mismatch",
    "v48_public_pwm_carrier_control_sample_switch_state_current_voltage_loss_owner_mismatch",
}


def _identity() -> dict[str, object]:
    skew_generation = "skew-aggregation-v48-901"
    pwm_generation = "pwm-timeline-v48-901"
    times = [0.0, 5.0e-6, 10.0e-6, 15.0e-6]
    states = [[1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 1, 1]]
    currents = [[10.0, -5.0, -5.0], [9.8, -4.7, -5.1], [9.5, -4.5, -5.0], [9.2, -4.2, -5.0]]
    voltages = [[300.0, -150.0, -150.0], [150.0, 150.0, -300.0], [-300.0, 150.0, 150.0], [-150.0, 300.0, -150.0]]
    losses = [2.0, 2.1, 1.9, 2.2]
    return {
        SKEW: {
            "generation": skew_generation,
            "slice_generation": skew_generation,
            "angle_generation": skew_generation,
            "harmonic_generation": skew_generation,
            "phase_generation": skew_generation,
            "result_generation": skew_generation,
            "slice_ids": ["slice:0", "slice:1", "slice:2"],
            "result_slice_ids": ["slice:0", "slice:1", "slice:2"],
            "slice_weights": [0.25, 0.50, 0.25],
            "result_slice_weights": [0.25, 0.50, 0.25],
            "rotor_angles_deg": [-5.0, 0.0, 5.0],
            "result_rotor_angles_deg": [-5.0, 0.0, 5.0],
            "torque_harmonic_phasors_nm": [[1.0, 0.0], [0.8, 0.1], [0.6, -0.1]],
            "result_torque_harmonic_phasors_nm": [[1.0, 0.0], [0.8, 0.1], [0.6, -0.1]],
            "phase_origins_deg": [0.0, 0.0, 0.0],
            "result_phase_origins_deg": [0.0, 0.0, 0.0],
            "machine_state_owner": "machine-state:skew-v48-901",
            "result_machine_state_owner": "machine-state:skew-v48-901",
            "result_sha256": "5" * 64,
            "accepted_result_sha256": "5" * 64,
        },
        PWM: {
            "generation": pwm_generation,
            "carrier_generation": pwm_generation,
            "control_generation": pwm_generation,
            "switch_generation": pwm_generation,
            "electrical_generation": pwm_generation,
            "loss_generation": pwm_generation,
            "result_generation": pwm_generation,
            "sample_times_s": times,
            "result_sample_times_s": times,
            "carrier_frequency_hz": 20000.0,
            "result_carrier_frequency_hz": 20000.0,
            "control_sample_divider": 2,
            "result_control_sample_divider": 2,
            "switch_states": states,
            "result_switch_states": states,
            "phase_current_a": currents,
            "result_phase_current_a": currents,
            "phase_voltage_v": voltages,
            "result_phase_voltage_v": voltages,
            "loss_w": losses,
            "result_loss_w": losses,
            "timeline_owner": "timeline:pwm-v48-901",
            "result_timeline_owner": "timeline:pwm-v48-901",
            "result_sha256": "6" * 64,
            "accepted_result_sha256": "6" * 64,
        },
    }


def test_v48_positive_skew_and_pwm_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v48_skew_permutation_and_phase_mutations_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[SKEW]["result_slice_weights"] = [0.50, 0.25, 0.25]
    identity[SKEW]["result_phase_origins_deg"] = [30.0, 0.0, 0.0]
    assert validate_public_identity(identity)["motor_v48_skew_slice_angle_harmonic_phase_owner"] is False


def test_v48_pwm_timeline_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[PWM]["result_sample_times_s"] = [0.0, 10.0e-6, 5.0e-6, 15.0e-6]
    identity[PWM]["result_timeline_owner"] = "timeline:pwm-v48-old"
    assert validate_public_identity(identity)["motor_v48_pwm_timeline_switch_electrical_loss_owner"] is False
