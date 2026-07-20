from copy import deepcopy
import math

from radia_mcp.radia_ngsolve.femm_artifact_identity_v55 import (
    CAPACITANCE,
    INDUCTION,
    validate_public_identity,
)


CASE_IDS = {
    "v55_public_induction_skin_depth_jouleloss_complexfield_frequency_conductor_owner_mismatch",
    "v55_public_electrostatic_capacitance_charge_voltage_energy_symmetry_owner_mismatch",
}


def _generation(generation: str, names: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{name: generation for name in names}}


def _payload():
    frequency = 1000.0; conductivity = 5.8e7; permeability = 4.0e-7 * math.pi
    depth = math.sqrt(2.0 / (2.0 * math.pi * frequency * permeability * conductivity))
    induction = {**_generation("induction-v55", ("skin_generation", "field_generation", "loss_generation", "frequency_generation", "material_generation", "owner_generation", "result_generation")), "frequency_hz": frequency, "result_frequency_hz": frequency, "conductivity_s_m": conductivity, "result_conductivity_s_m": conductivity, "permeability_h_m": permeability, "result_permeability_h_m": permeability, "skin_depth_m": depth, "result_skin_depth_m": depth, "complex_magnetic_field_a_m": {"real": 1200.0, "imag": -350.0}, "result_complex_magnetic_field_a_m": {"real": 1200.0, "imag": -350.0}, "joule_loss_w": 12.5, "result_joule_loss_w": 12.5, "conductor_owner": "conductor:v55", "result_conductor_owner": "conductor:v55", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    matrix = [[2.0e-12, -1.0e-12], [-1.0e-12, 2.0e-12]]
    capacitance = {**_generation("capacitance-v55", ("capacitance_generation", "charge_generation", "voltage_generation", "energy_generation", "symmetry_generation", "owner_generation", "result_generation")), "conductor_order": ["conductor:1", "conductor:2"], "result_conductor_order": ["conductor:1", "conductor:2"], "capacitance_matrix_f": matrix, "result_capacitance_matrix_f": matrix, "voltage_v": [1.0, 0.0], "result_voltage_v": [1.0, 0.0], "charge_c": [2.0e-12, -1.0e-12], "result_charge_c": [2.0e-12, -1.0e-12], "stored_energy_j": 1.0e-12, "result_stored_energy_j": 1.0e-12, "solution_owner": "solution:v55", "result_solution_owner": "solution:v55", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {INDUCTION: induction, CAPACITANCE: capacitance}


def test_v55_positive_identities_are_accepted():
    assert all(validate_public_identity(_payload()).values())


def test_v55_frozen_mutations_are_rejected():
    payload = deepcopy(_payload())
    payload[INDUCTION]["result_skin_depth_m"] = 0.1
    payload[CAPACITANCE]["result_conductor_order"] = ["conductor:2", "conductor:1"]
    assert not all(validate_public_identity(payload).values())


def test_v55_self_consistent_wrong_skin_depth_or_negative_loss_is_rejected():
    payload = deepcopy(_payload())
    payload[INDUCTION]["skin_depth_m"] = payload[INDUCTION]["result_skin_depth_m"] = 0.1
    payload[INDUCTION]["joule_loss_w"] = payload[INDUCTION]["result_joule_loss_w"] = -1.0
    assert not all(validate_public_identity(payload).values())


def test_v55_self_consistent_nonsymmetric_capacitance_or_bad_charge_is_rejected():
    payload = deepcopy(_payload())
    bad_matrix = [[2.0e-12, 1.0e-12], [-1.0e-12, 2.0e-12]]
    payload[CAPACITANCE]["capacitance_matrix_f"] = payload[CAPACITANCE]["result_capacitance_matrix_f"] = bad_matrix
    payload[CAPACITANCE]["charge_c"] = payload[CAPACITANCE]["result_charge_c"] = [9.0, 9.0]
    assert not all(validate_public_identity(payload).values())


def test_v55_numeric_sha256_values_are_rejected():
    payload = _payload()
    numeric_digest = int("9" * 64)
    for row in payload.values():
        row["result_sha256"] = numeric_digest
        row["accepted_result_sha256"] = numeric_digest
    assert not all(validate_public_identity(payload).values())
