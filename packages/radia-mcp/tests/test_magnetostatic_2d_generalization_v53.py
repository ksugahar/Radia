from copy import deepcopy
import math

from radia_mcp.radia_ngsolve.electromagnetic_artifact_identity_v53 import HYSTERESIS, VIRTUAL_WORK, validate_public_identity


PROMOTED_CASE_IDS = {
    "v53_public_hysteresis_complex_permeability_phasor_lossdensity_material_owner_mismatch",
    "v53_public_electrostatic_virtualwork_voltage_charge_constraint_force_owner_mismatch",
}


def _generations(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _identity():
    frequency = 400.0; mu = {"real": 220.0, "imag": -18.0}; h_rms = 120.0
    loss = 2.0 * math.pi * frequency * 4.0e-7 * math.pi * (-mu["imag"]) * h_rms**2
    hysteresis = {**_generations("hyst-v53", ("permeability_generation", "phasor_generation", "loss_generation", "material_generation", "owner_generation", "result_generation")), "relative_permeability": mu, "result_relative_permeability": mu, "phasor_convention": "exp(+j_omega_t)", "result_phasor_convention": "exp(+j_omega_t)", "frequency_hz": frequency, "result_frequency_hz": frequency, "h_rms_a_per_m": h_rms, "result_h_rms_a_per_m": h_rms, "loss_density_w_m3": loss, "result_loss_density_w_m3": loss, "material_id": "material:steel", "result_material_id": "material:steel", "material_owner": "material-owner:v53", "result_material_owner": "material-owner:v53", "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64}
    virtual_work = {**_generations("vw-v53", ("path_generation", "constraint_generation", "energy_generation", "force_generation", "owner_generation", "result_generation")), "virtual_work_path": "constant_voltage_coenergy", "result_virtual_work_path": "constant_voltage_coenergy", "constraint_mode": "fixed_voltage", "result_constraint_mode": "fixed_voltage", "voltage_v": 800.0, "result_voltage_v": 800.0, "charge_c": 2.0e-8, "result_charge_c": 2.0e-8, "virtual_displacement_m": 1.0e-6, "result_virtual_displacement_m": 1.0e-6, "coenergy_before_j": 0.01, "result_coenergy_before_j": 0.01, "coenergy_after_j": 0.010004, "result_coenergy_after_j": 0.010004, "force_n": 4.0, "result_force_n": 4.0, "force_owner": "force:v53", "result_force_owner": "force:v53", "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64}
    return {HYSTERESIS: hysteresis, VIRTUAL_WORK: virtual_work}


def test_v53_positive_public_artifacts_are_accepted():
    assert all(validate_public_identity(_identity()).values())


def test_v53_frozen_counterfactuals_are_rejected():
    identity = deepcopy(_identity())
    identity[HYSTERESIS]["result_phasor_convention"] = "exp(-j_omega_t)"
    identity[VIRTUAL_WORK]["result_constraint_mode"] = "fixed_charge"
    assert not all(validate_public_identity(identity).values())


def test_v53_self_consistent_wrong_physics_is_rejected():
    identity = deepcopy(_identity())
    identity[HYSTERESIS]["relative_permeability"]["imag"] = identity[HYSTERESIS]["result_relative_permeability"]["imag"] = 18.0
    identity[VIRTUAL_WORK]["force_n"] = identity[VIRTUAL_WORK]["result_force_n"] = -4.0
    assert not all(validate_public_identity(identity).values())
