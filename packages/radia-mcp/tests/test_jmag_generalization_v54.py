from copy import deepcopy

from radia_mcp.radia_ngsolve.motor_artifact_identity_v54 import DEMAG, TORQUE, validate_public_identity


CASE_IDS = {
    "v54_public_torqueripple_harmonic_mechanical_electrical_angle_polepair_owner_mismatch",
    "v54_public_demag_irreversible_knee_temperature_currentvector_recovery_owner_mismatch",
}


def _generations(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _payload():
    mechanical = [0.0, 5.0, 10.0]; electrical = [0.0, 20.0, 40.0]; harmonics = [{"order": 6, "amplitude_nm": 0.35, "phase_electrical_deg": 25.0}]
    torque = {**_generations("torque-v54", ("harmonic_generation", "mechanical_generation", "electrical_generation", "polepair_generation", "owner_generation", "result_generation")), "pole_pairs": 4, "result_pole_pairs": 4, "mechanical_angles_deg": mechanical, "result_mechanical_angles_deg": mechanical, "electrical_angles_deg": electrical, "result_electrical_angles_deg": electrical, "torque_harmonics": harmonics, "result_torque_harmonics": harmonics, "result_owner": "result:v54", "accepted_result_owner": "result:v54", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    knee = {"b_t": 0.24, "h_a_per_m": -640000.0, "criterion": "operating_point_below_knee"}; current = [120.0, -60.0, -60.0]
    demag = {**_generations("demag-v54", ("knee_generation", "temperature_generation", "current_generation", "recovery_generation", "owner_generation", "result_generation")), "knee_criterion": knee, "result_knee_criterion": knee, "temperature_c": 160.0, "result_temperature_c": 160.0, "current_vector_abc_a": current, "result_current_vector_abc_a": current, "irreversible_demag_fraction": 0.05, "result_irreversible_demag_fraction": 0.05, "post_recovery_remanence_fraction": 0.95, "result_post_recovery_remanence_fraction": 0.95, "recovery_state": "partially_demagnetized", "result_recovery_state": "partially_demagnetized", "magnet_owner": "magnet:v54", "result_magnet_owner": "magnet:v54", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {TORQUE: torque, DEMAG: demag}


def test_v54_positive_identities_are_accepted():
    assert all(validate_public_identity(_payload()).values())


def test_v54_frozen_mutations_are_rejected():
    payload = deepcopy(_payload())
    payload[TORQUE]["result_pole_pairs"] = 3
    payload[DEMAG]["result_temperature_c"] = 20.0
    assert not all(validate_public_identity(payload).values())


def test_v54_self_consistent_nonphysical_records_are_rejected():
    payload = deepcopy(_payload())
    payload[TORQUE]["electrical_angles_deg"] = payload[TORQUE]["result_electrical_angles_deg"] = [0.0, 15.0, 30.0]
    payload[DEMAG]["current_vector_abc_a"] = payload[DEMAG]["result_current_vector_abc_a"] = [120.0, -50.0, -50.0]
    assert not all(validate_public_identity(payload).values())


def test_v54_malformed_values_reject_without_raising():
    payload = deepcopy(_payload())
    payload[TORQUE]["torque_harmonics"] = [{"order": [6], "amplitude_nm": 0.35, "phase_electrical_deg": 25.0}]
    payload[DEMAG]["current_vector_abc_a"] = [[120.0], -60.0, -60.0]
    assert not all(validate_public_identity(payload).values())
