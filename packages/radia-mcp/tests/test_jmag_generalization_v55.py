from copy import deepcopy

from radia_mcp.radia_ngsolve.motor_artifact_identity_v55 import DQ, IRON, validate_public_identity


CASE_IDS = {
    "v55_public_dq_fluxlinkage_current_angle_saliency_torque_owner_mismatch",
    "v55_public_ironloss_hysteresis_eddy_excess_frequency_fluxdensity_owner_mismatch",
}


def _payload():
    from math import atan2, degrees

    generation = lambda name, fields: {"generation": name, **{field: name for field in fields}}
    current = {"d": -20.0, "q": 80.0}; inductance = {"d": 1.0e-3, "q": 1.5e-3}; pm_flux = 0.08
    flux = {"d": pm_flux + inductance["d"] * current["d"], "q": inductance["q"] * current["q"]}; pole_pairs = 4
    angle = degrees(atan2(current["q"], current["d"])) % 360.0
    saliency = 1.5 * pole_pairs * (inductance["d"] - inductance["q"]) * current["d"] * current["q"]
    torque = 1.5 * pole_pairs * (flux["d"] * current["q"] - flux["q"] * current["d"])
    dq = {**generation("dq-v55", ("flux_generation", "current_generation", "angle_generation", "saliency_generation", "torque_generation", "owner_generation", "result_generation")), "pole_pairs": pole_pairs, "result_pole_pairs": pole_pairs, "current_dq_a": current, "result_current_dq_a": current, "inductance_dq_h": inductance, "result_inductance_dq_h": inductance, "pm_flux_linkage_wb": pm_flux, "result_pm_flux_linkage_wb": pm_flux, "flux_linkage_dq_wb": flux, "result_flux_linkage_dq_wb": flux, "current_electrical_angle_deg": angle, "result_current_electrical_angle_deg": angle, "saliency_torque_nm": saliency, "result_saliency_torque_nm": saliency, "torque_nm": torque, "result_torque_nm": torque, "result_owner": "result:dq-v55", "accepted_result_owner": "result:dq-v55", "result_sha256": "5" * 64, "accepted_result_sha256": "5" * 64}
    components = {"hysteresis_w": 12.0, "eddy_w": 8.0, "excess_w": 2.0}; waveform = [0.0, 1.4, 0.0, -1.4, 0.0]
    iron = {**generation("iron-v55", ("component_generation", "frequency_generation", "waveform_generation", "material_generation", "owner_generation", "result_generation")), "loss_components": components, "result_loss_components": components, "total_iron_loss_w": 22.0, "result_total_iron_loss_w": 22.0, "frequency_hz": 400.0, "result_frequency_hz": 400.0, "flux_density_waveform_t": waveform, "result_flux_density_waveform_t": waveform, "peak_flux_density_t": 1.4, "result_peak_flux_density_t": 1.4, "material_revision": "steel-v55-r7", "result_material_revision": "steel-v55-r7", "material_owner": "material:steel-v55", "result_material_owner": "material:steel-v55", "result_sha256": "6" * 64, "accepted_result_sha256": "6" * 64}
    return {DQ: dq, IRON: iron}


def test_v55_positive_identities_are_accepted():
    assert all(validate_public_identity(_payload()).values())


def test_v55_frozen_mutations_are_rejected():
    payload = deepcopy(_payload()); payload[DQ]["accepted_result_owner"] = "result:stale"; payload[IRON]["result_material_owner"] = "material:stale"
    assert not all(validate_public_identity(payload).values())


def test_v55_self_consistent_nonphysical_records_are_rejected():
    payload = deepcopy(_payload())
    payload[DQ]["flux_linkage_dq_wb"] = payload[DQ]["result_flux_linkage_dq_wb"] = {"d": 0.4, "q": 0.4}
    payload[IRON]["loss_components"] = payload[IRON]["result_loss_components"] = {"hysteresis_w": -1.0, "eddy_w": 8.0, "excess_w": 2.0}
    payload[IRON]["total_iron_loss_w"] = payload[IRON]["result_total_iron_loss_w"] = 9.0
    assert not all(validate_public_identity(payload).values())


def test_v55_malformed_values_reject_without_raising():
    payload = deepcopy(_payload()); payload[DQ]["current_dq_a"] = {"d": [1.0], "q": 2.0}; payload[IRON]["flux_density_waveform_t"] = [0.0, [1.0], 0.0]
    assert not all(validate_public_identity(payload).values())
