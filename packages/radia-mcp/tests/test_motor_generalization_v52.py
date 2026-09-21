from copy import deepcopy

from radia_mcp.radia_ngsolve.motor_artifact_identity_v52 import COGGING, IRON_LOSS, validate_public_identity


PROMOTED_CASE_IDS = {
    "v52_public_ironloss_fft_window_harmonic_rotation_frequency_coefficient_owner_mismatch",
    "v52_public_cogging_period_slotpole_sampling_phase_alignment_owner_mismatch",
}


def _generations(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _identity() -> dict[str, object]:
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    harmonics = [{"order": 1, "amplitude_t": 1.2}, {"order": 3, "amplitude_t": 0.08}]
    coefficients = {"hysteresis": 1.7, "classical_eddy": 0.012, "excess": 0.08}
    angles = [index * 0.5 for index in range(13)]
    torque = [0.12, 0.10, 0.04, -0.04, -0.10, -0.12, -0.10, -0.04, 0.04, 0.10, 0.12, 0.10, 0.12]
    return {
        IRON_LOSS: {
            **_generations("iron-loss-v52", ("window_generation", "harmonic_generation", "frequency_generation", "coefficient_generation", "owner_generation", "result_generation")),
            "fft_window": "hann_periodic", "result_fft_window": "hann_periodic",
            "sample_count": 1024, "result_sample_count": 1024,
            "harmonics": harmonics, "result_harmonics": harmonics,
            "rotation_frequency_hz": 50.0, "result_rotation_frequency_hz": 50.0,
            "pole_pairs": 4, "result_pole_pairs": 4,
            "electrical_frequency_hz": 200.0, "result_electrical_frequency_hz": 200.0,
            "loss_coefficients": coefficients, "result_loss_coefficients": coefficients,
            "waveform_owner": "waveform:iron-loss-v52", "result_waveform_owner": "waveform:iron-loss-v52",
            **result,
        },
        COGGING: {
            **_generations("cogging-v52", ("period_generation", "sampling_generation", "phase_generation", "owner_generation", "result_generation")),
            "slot_count": 12, "result_slot_count": 12, "pole_count": 10, "result_pole_count": 10,
            "cogging_period_mechanical_deg": 6.0, "result_cogging_period_mechanical_deg": 6.0,
            "sample_angles_mechanical_deg": angles, "result_sample_angles_mechanical_deg": angles,
            "torque_samples_nm": torque, "result_torque_samples_nm": torque,
            "phase_alignment": "slot_center_to_pole_center", "result_phase_alignment": "slot_center_to_pole_center",
            "torque_owner": "torque:cogging-v52", "result_torque_owner": "torque:cogging-v52", **result,
        },
    }


def test_v52_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v52_frozen_counterfactuals_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[IRON_LOSS]["result_fft_window"] = "rectangular"
    identity[COGGING]["result_cogging_period_mechanical_deg"] = 30.0
    assert not all(validate_public_identity(identity).values())


def test_v52_self_consistent_wrong_physics_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[IRON_LOSS]["fft_window"] = identity[IRON_LOSS]["result_fft_window"] = "rectangular"
    identity[COGGING]["cogging_period_mechanical_deg"] = identity[COGGING]["result_cogging_period_mechanical_deg"] = 30.0
    assert not all(validate_public_identity(identity).values())
