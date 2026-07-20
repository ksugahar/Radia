from copy import deepcopy

from radia_mcp.radia_ngsolve.femm_artifact_identity_v54 import FORCE, POWER, validate_public_identity


CASE_IDS = {
    "v54_public_harmonic_complexpower_active_reactive_loss_frequency_circuit_owner_mismatch",
    "v54_public_axisymmetric_weightedstress_force_radius_measure_selection_owner_mismatch",
}


def _generation(generation: str, names: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{name: generation for name in names}}


def _payload():
    power = {**_generation("power-v54", ("power_generation", "loss_generation", "frequency_generation", "circuit_generation", "owner_generation", "result_generation")), "complex_power_va": {"real": 120.0, "imag": 45.0}, "result_complex_power_va": {"real": 120.0, "imag": 45.0}, "active_power_w": 120.0, "result_active_power_w": 120.0, "reactive_power_var": 45.0, "result_reactive_power_var": 45.0, "loss_components_w": {"copper": 80.0, "core": 40.0}, "result_loss_components_w": {"copper": 80.0, "core": 40.0}, "frequency_hz": 400.0, "result_frequency_hz": 400.0, "circuit_id": "circuit:a", "result_circuit_id": "circuit:a", "circuit_owner": "circuit-owner:v54", "result_circuit_owner": "circuit-owner:v54", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    selection = ["block:armature", "interface:airgap"]
    force = {**_generation("force-v54", ("radius_generation", "stress_generation", "selection_generation", "direction_generation", "owner_generation", "result_generation")), "radius_weighting": "2*pi*r", "result_radius_weighting": "2*pi*r", "stress_measure": "weighted_stress_tensor", "result_stress_measure": "weighted_stress_tensor", "integration_selection": selection, "result_integration_selection": selection, "force_direction_rz": [1.0, 0.0], "result_force_direction_rz": [1.0, 0.0], "force_n": 12.5, "result_force_n": 12.5, "mesh_owner": "mesh:v54", "result_mesh_owner": "mesh:v54", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {POWER: power, FORCE: force}


def test_v54_positive_identities_are_accepted():
    assert all(validate_public_identity(_payload()).values())


def test_v54_frozen_mutations_are_rejected():
    payload = deepcopy(_payload())
    payload[POWER]["result_active_power_w"] = 80.0
    payload[FORCE]["result_radius_weighting"] = "planar"
    assert not all(validate_public_identity(payload).values())


def test_v54_self_consistent_nonphysical_records_are_rejected():
    payload = deepcopy(_payload())
    payload[POWER]["loss_components_w"] = payload[POWER]["result_loss_components_w"] = {"copper": 10.0}
    payload[FORCE]["force_direction_rz"] = payload[FORCE]["result_force_direction_rz"] = [2.0, 0.0]
    assert not all(validate_public_identity(payload).values())


def test_v54_malformed_values_reject_without_raising():
    payload = deepcopy(_payload())
    payload[POWER]["loss_components_w"] = {"copper": [80.0]}
    payload[FORCE]["integration_selection"] = [["block:armature"]]
    assert not all(validate_public_identity(payload).values())
