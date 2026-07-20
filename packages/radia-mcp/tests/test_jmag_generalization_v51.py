from copy import deepcopy

from radia_mcp.radia_ngsolve.motor_artifact_identity_v51 import TORQUE, WINDING, validate_public_identity


PROMOTED_CASE_IDS = {
    "v51_public_torque_ripple_rotor_angle_electrical_mechanical_period_fft_window_owner_mismatch",
    "v51_public_winding_temperature_resistance_endturn_length_fillfactor_copperloss_owner_mismatch",
}


def _identity() -> dict[str, object]:
    generation = "motor-public-v51"
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    mechanical = [index * 11.25 for index in range(9)]
    electrical = [value * 4.0 for value in mechanical]
    torque = [10.0, 10.4, 10.0, 9.6, 10.0, 10.4, 10.0, 9.6, 10.0]
    harmonics = [{"order": 0, "amplitude_nm": 10.0}, {"order": 2, "amplitude_nm": 0.4}]
    return {
        TORQUE: {
            "generation": generation, "rotor_generation": generation, "angle_generation": generation,
            "period_generation": generation, "window_generation": generation, "fft_generation": generation,
            "owner_generation": generation, "result_generation": generation, "pole_pairs": 4, "result_pole_pairs": 4,
            "rotor_angles_mechanical_deg": mechanical, "result_rotor_angles_mechanical_deg": mechanical,
            "rotor_angles_electrical_deg": electrical, "result_rotor_angles_electrical_deg": electrical,
            "mechanical_period_deg": 90.0, "result_mechanical_period_deg": 90.0, "electrical_period_deg": 360.0,
            "result_electrical_period_deg": 360.0, "sample_window_mechanical_deg": [0.0, 90.0],
            "result_sample_window_mechanical_deg": [0.0, 90.0], "torque_samples_nm": torque,
            "result_torque_samples_nm": torque, "fft_harmonics": harmonics, "result_fft_harmonics": harmonics,
            "torque_owner": "torque:motor-v51", "result_torque_owner": "torque:motor-v51", **result,
        },
        WINDING: {
            "generation": generation, "temperature_generation": generation, "resistance_generation": generation,
            "length_generation": generation, "fill_generation": generation, "loss_generation": generation,
            "owner_generation": generation, "result_generation": generation, "reference_temperature_c": 20.0,
            "winding_temperature_c": 120.0, "result_winding_temperature_c": 120.0,
            "copper_temperature_coefficient_per_k": 0.00393, "resistance_reference_ohm": 0.08,
            "resistance_at_temperature_ohm": 0.11144, "result_resistance_at_temperature_ohm": 0.11144,
            "active_length_m": 0.42, "end_turn_length_m": 0.18, "result_end_turn_length_m": 0.18,
            "slot_fill_factor": 0.48, "result_slot_fill_factor": 0.48, "current_rms_a": 20.0,
            "copper_loss_w": 44.576, "result_copper_loss_w": 44.576, "winding_owner": "winding:phase-u-v51",
            "result_winding_owner": "winding:phase-u-v51", **result,
        },
    }


def test_v51_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v51_frozen_counterfactuals_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[TORQUE].update({"result_electrical_period_deg": 90.0, "result_torque_owner": "torque:stale"})
    identity[WINDING].update({"result_resistance_at_temperature_ohm": 0.08, "result_winding_owner": "winding:stale"})
    assert not all(validate_public_identity(identity).values())


def test_v51_self_consistent_wrong_physics_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[TORQUE]["electrical_period_deg"] = identity[TORQUE]["result_electrical_period_deg"] = 90.0
    identity[WINDING]["resistance_at_temperature_ohm"] = identity[WINDING]["result_resistance_at_temperature_ohm"] = 0.08
    assert not all(validate_public_identity(identity).values())
