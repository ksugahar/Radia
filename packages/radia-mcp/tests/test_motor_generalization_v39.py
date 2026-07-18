from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import pwm_controlled_motor_loss_gate
from test_motor_generalization_v38 import _payload_v38


_SKEW = "skewed_rotor_slice_angle_phase_weight_torque_ripple_power_model_owner_result_identity"
_DEMAG = "pm_irreversible_demag_temperature_recoil_knee_operating_flux_torque_mesh_owner_result_identity"
_PROMOTED_CASE_IDS = (
    "v39_public_skewed_rotor_slice_angle_phase_weight_torque_ripple_power_mismatch",
    "v39_public_pm_irreversible_demag_recoil_temperature_operatingpoint_flux_torque_mismatch",
)


def _payload_v39():
    payload = _payload_v38()
    identity = payload["artifact_identity"]
    generation = "skewed-rotor-271"
    angles = [-5.0, 0.0, 5.0]
    pole_pairs = 2
    phases = [pole_pairs * angle for angle in angles]
    weights = [0.25, 0.5, 0.25]
    torque = [10.0, 10.4, 10.0]
    phasors = [[math.cos(math.radians(phase)), math.sin(math.radians(phase))] for phase in phases]
    mean_torque = sum(weight * item for weight, item in zip(weights, torque))
    ripple = math.hypot(
        sum(weight * pair[0] for weight, pair in zip(weights, phasors)),
        sum(weight * pair[1] for weight, pair in zip(weights, phasors)),
    )
    identity[_SKEW] = {
        "skew_generation": generation,
        **{key: generation for key in ("slice_generation", "phase_generation", "weight_generation", "torque_generation", "ripple_generation", "power_generation", "owner_generation", "result_generation")},
        "pole_pairs": pole_pairs,
        "result_pole_pairs": pole_pairs,
        "slice_angles_mechanical_deg": angles,
        "result_slice_angles_mechanical_deg": angles,
        "slice_phase_offsets_electrical_deg": phases,
        "result_slice_phase_offsets_electrical_deg": phases,
        "axial_weights": weights,
        "result_axial_weights": weights,
        "slice_mean_torque_nm": torque,
        "result_slice_mean_torque_nm": torque,
        "weighted_mean_torque_nm": mean_torque,
        "result_weighted_mean_torque_nm": mean_torque,
        "slice_ripple_phasor": phasors,
        "result_slice_ripple_phasor": phasors,
        "weighted_ripple_residual": ripple,
        "result_weighted_ripple_residual": ripple,
        "mechanical_speed_rad_s": 100.0,
        "result_mechanical_speed_rad_s": 100.0,
        "mechanical_power_w": mean_torque * 100.0,
        "result_mechanical_power_w": mean_torque * 100.0,
        "model_owner": "motor:skewed-rotor-271",
        "accepted_model_owner": "motor:skewed-rotor-271",
        "skew_result_sha256": "1" * 64,
        "accepted_skew_result_sha256": "1" * 64,
    }

    generation = "pm-demag-271"
    temperature = 120.0
    br_reference = 1.2
    coefficient = -1.1e-3
    br_temperature = br_reference * (1.0 + coefficient * (temperature - 20.0))
    recoil_mu = 1.05
    operating_h = -8.0e5
    operating_b = br_temperature + 4.0e-7 * math.pi * recoil_mu * operating_h
    loss = 0.08
    identity[_DEMAG] = {
        "demag_generation": generation,
        **{key: generation for key in ("temperature_generation", "recoil_generation", "knee_generation", "operating_generation", "remanence_generation", "flux_generation", "torque_generation", "mesh_generation", "owner_generation", "result_generation")},
        "reference_temperature_c": 20.0,
        "result_reference_temperature_c": 20.0,
        "magnet_temperature_c": temperature,
        "result_magnet_temperature_c": temperature,
        "remanence_reference_t": br_reference,
        "result_remanence_reference_t": br_reference,
        "remanence_temperature_coefficient_per_k": coefficient,
        "result_remanence_temperature_coefficient_per_k": coefficient,
        "temperature_adjusted_remanence_t": br_temperature,
        "result_temperature_adjusted_remanence_t": br_temperature,
        "recoil_relative_permeability": recoil_mu,
        "result_recoil_relative_permeability": recoil_mu,
        "operating_h_a_per_m": operating_h,
        "result_operating_h_a_per_m": operating_h,
        "operating_b_t": operating_b,
        "result_operating_b_t": operating_b,
        "knee_h_a_per_m": -7.0e5,
        "result_knee_h_a_per_m": -7.0e5,
        "irreversible_region": True,
        "result_irreversible_region": True,
        "remanence_loss_fraction": loss,
        "result_remanence_loss_fraction": loss,
        "airgap_flux_before_wb": 1.0e-2,
        "result_airgap_flux_before_wb": 1.0e-2,
        "airgap_flux_after_wb": 1.0e-2 * (1.0 - loss),
        "result_airgap_flux_after_wb": 1.0e-2 * (1.0 - loss),
        "torque_before_nm": 10.0,
        "result_torque_before_nm": 10.0,
        "torque_after_nm": 10.0 * (1.0 - loss),
        "result_torque_after_nm": 10.0 * (1.0 - loss),
        "mesh_owner": "mesh:pm-demag-271",
        "accepted_mesh_owner": "mesh:pm-demag-271",
        "demag_result_sha256": "2" * 64,
        "accepted_demag_result_sha256": "2" * 64,
    }
    return payload


def test_v39_public_positive_skew_and_irreversible_demag_closure():
    assert pwm_controlled_motor_loss_gate(_payload_v39())["status"] == "ok"


def test_v39_public_skewed_rotor_slice_angle_phase_weight_torque_ripple_power_mismatch():
    payload = _payload_v39()
    payload["artifact_identity"][_SKEW].update({"phase_generation": "skewed-rotor-270", "power_generation": "skewed-rotor-269", "result_generation": "skewed-rotor-268", "result_slice_phase_offsets_electrical_deg": [5.0, 0.0, -5.0], "result_axial_weights": [0.5, 0.5, 0.5], "result_slice_mean_torque_nm": [9.0, 9.0, 9.0], "result_weighted_mean_torque_nm": -1.0, "result_weighted_ripple_residual": -1.0, "result_mechanical_power_w": -100.0, "accepted_model_owner": "stale:motor", "accepted_skew_result_sha256": "a" * 64})
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v39_public_pm_irreversible_demag_recoil_temperature_operatingpoint_flux_torque_mismatch():
    payload = _payload_v39()
    payload["artifact_identity"][_DEMAG].update({"temperature_generation": "pm-demag-270", "operating_generation": "pm-demag-269", "result_generation": "pm-demag-268", "result_temperature_adjusted_remanence_t": -1.0, "result_operating_h_a_per_m": 8.0e5, "result_operating_b_t": -1.0, "result_irreversible_region": False, "result_remanence_loss_fraction": -0.1, "result_airgap_flux_after_wb": 2.0e-2, "result_torque_after_nm": 20.0, "accepted_mesh_owner": "stale:mesh", "accepted_demag_result_sha256": "b" * 64})
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v39_public_rejects_self_consistent_wrong_skew_phase_conversion():
    payload = _payload_v39()
    row = payload["artifact_identity"][_SKEW]
    phases = row["slice_angles_mechanical_deg"]
    row["slice_phase_offsets_electrical_deg"] = phases
    row["result_slice_phase_offsets_electrical_deg"] = phases
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v39_public_rejects_self_consistent_reversible_demag_claim():
    payload = _payload_v39()
    row = payload["artifact_identity"][_DEMAG]
    row["irreversible_region"] = False
    row["result_irreversible_region"] = False
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"
