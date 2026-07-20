from copy import deepcopy

from radia_mcp.radia_ngsolve.motor_artifact_identity_v53 import DEMAG, SKEW, validate_public_identity


PROMOTED_CASE_IDS = {
    "v53_public_skew_slice_weight_angle_harmonic_torque_rotor_owner_mismatch",
    "v53_public_magnet_demag_operatingpoint_temperature_recoil_irreversible_owner_mismatch",
}


def _generations(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _identity():
    weights = [0.25, 0.5, 0.25]; angles = [-5.0, 0.0, 5.0]
    harmonics = [{"order": 1, "amplitude_nm": 12.0, "phase_deg": 0.0}, {"order": 6, "amplitude_nm": 0.35, "phase_deg": 25.0}]
    skew = {**_generations("skew-v53", ("slice_generation", "angle_generation", "harmonic_generation", "owner_generation", "result_generation")), "slice_weights": weights, "result_slice_weights": weights, "skew_angles_mechanical_deg": angles, "result_skew_angles_mechanical_deg": angles, "harmonic_torque": harmonics, "result_harmonic_torque": harmonics, "rotor_owner": "rotor:v53", "result_rotor_owner": "rotor:v53", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    operating = {"b_t": 0.62, "h_a_per_m": -420000.0}; recoil = {"relative_permeability": 1.05, "coercivity_a_per_m": 900000.0}
    demag = {**_generations("demag-v53", ("operating_generation", "temperature_generation", "recoil_generation", "irreversible_generation", "owner_generation", "result_generation")), "operating_point": operating, "result_operating_point": operating, "temperature_c": 140.0, "result_temperature_c": 140.0, "recoil_line": recoil, "result_recoil_line": recoil, "irreversible_demag_fraction": 0.03, "result_irreversible_demag_fraction": 0.03, "irreversible_state": "partially_demagnetized", "result_irreversible_state": "partially_demagnetized", "magnet_owner": "magnet:v53", "result_magnet_owner": "magnet:v53", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {SKEW: skew, DEMAG: demag}


def test_v53_positive_public_artifacts_are_accepted():
    assert all(validate_public_identity(_identity()).values())


def test_v53_frozen_counterfactuals_are_rejected():
    identity = deepcopy(_identity())
    identity[SKEW]["result_slice_weights"] = [1.0, 0.0, 0.0]
    identity[DEMAG]["result_temperature_c"] = 20.0
    assert not all(validate_public_identity(identity).values())


def test_v53_self_consistent_wrong_physics_is_rejected():
    identity = deepcopy(_identity())
    identity[SKEW]["slice_weights"] = identity[SKEW]["result_slice_weights"] = [0.5, 0.5, 0.5]
    identity[DEMAG]["irreversible_demag_fraction"] = identity[DEMAG]["result_irreversible_demag_fraction"] = 0.0
    assert not all(validate_public_identity(identity).values())
