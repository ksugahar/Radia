from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_motor_generalization_v37 import _payload_v37


_PROMOTED_CASE_IDS = (
    "v38_public_wound_field_synchronous_excitation_flux_torque_powerfactor_fieldloss_energy_mismatch",
    "v38_public_flux_switching_pm_slot_pole_polarity_harmonic_backemf_torque_ripple_mismatch",
)


def _payload_v38():
    payload = _payload_v37()
    identity = payload["artifact_identity"]
    generation = "wound-field-258"
    field_current = 5.0
    excitation_inductance = 0.04
    excitation_flux = field_current * excitation_inductance
    torque_angle = math.pi / 6.0
    pole_pairs = 2
    stator_current = 10.0
    torque = (
        1.5
        * pole_pairs
        * excitation_flux
        * math.sqrt(2.0)
        * stator_current
        * math.sin(torque_angle)
    )
    speed = 100.0
    mechanical_power = torque * speed
    field_resistance = 1.0
    stator_resistance = 0.2
    field_loss = field_current**2 * field_resistance
    stator_loss = 3.0 * stator_current**2 * stator_resistance
    active_power = mechanical_power + field_loss + stator_loss
    line_voltage = 40.0
    apparent_power = math.sqrt(3.0) * line_voltage * stator_current
    identity[
        "wound_field_synchronous_excitation_flux_torque_angle_powerfactor_field_stator_loss_mechanical_energy_mesh_owner_result_identity"
    ] = {
        "wound_field_generation": generation,
        **{
            key: generation
            for key in (
                "excitation_generation", "flux_generation",
                "torque_generation", "powerfactor_generation",
                "field_loss_generation", "stator_loss_generation",
                "mechanical_generation", "energy_generation",
                "mesh_generation", "owner_generation", "result_generation",
            )
        },
        "field_current_a": field_current,
        "result_field_current_a": field_current,
        "excitation_inductance_h": excitation_inductance,
        "result_excitation_inductance_h": excitation_inductance,
        "excitation_flux_linkage_wb_turn": excitation_flux,
        "result_excitation_flux_linkage_wb_turn": excitation_flux,
        "torque_angle_rad": torque_angle,
        "result_torque_angle_rad": torque_angle,
        "pole_pairs": pole_pairs,
        "result_pole_pairs": pole_pairs,
        "stator_current_rms_a": stator_current,
        "result_stator_current_rms_a": stator_current,
        "electromagnetic_torque_nm": torque,
        "result_electromagnetic_torque_nm": torque,
        "mechanical_speed_rad_s": speed,
        "result_mechanical_speed_rad_s": speed,
        "mechanical_power_w": mechanical_power,
        "result_mechanical_power_w": mechanical_power,
        "field_resistance_ohm": field_resistance,
        "result_field_resistance_ohm": field_resistance,
        "field_copper_loss_w": field_loss,
        "result_field_copper_loss_w": field_loss,
        "stator_phase_resistance_ohm": stator_resistance,
        "result_stator_phase_resistance_ohm": stator_resistance,
        "stator_copper_loss_w": stator_loss,
        "result_stator_copper_loss_w": stator_loss,
        "line_voltage_rms_v": line_voltage,
        "result_line_voltage_rms_v": line_voltage,
        "apparent_power_va": apparent_power,
        "result_apparent_power_va": apparent_power,
        "active_input_power_w": active_power,
        "result_active_input_power_w": active_power,
        "power_factor": active_power / apparent_power,
        "result_power_factor": active_power / apparent_power,
        "energy_balance_residual_w": 0.0,
        "result_energy_balance_residual_w": 0.0,
        "energy_tolerance_w": 1.0e-9,
        "result_energy_tolerance_w": 1.0e-9,
        "mesh_owner": "mesh:wound-field-258",
        "accepted_mesh_owner": "mesh:wound-field-258",
        "motor_result_sha256": "1" * 64,
        "accepted_motor_result_sha256": "1" * 64,
    }

    generation = "flux-switching-258"
    torque_samples = [10.0, 11.0, 9.0, 10.0]
    average_torque = sum(torque_samples) / len(torque_samples)
    torque_ripple = (max(torque_samples) - min(torque_samples)) / average_torque
    phase_angles = [0.0, -2.0 * math.pi / 3.0, 2.0 * math.pi / 3.0]
    identity[
        "flux_switching_pm_slot_pole_polarity_phase_harmonic_backemf_torque_ripple_periodicity_mesh_owner_result_identity"
    ] = {
        "flux_switching_generation": generation,
        **{
            key: generation
            for key in (
                "slot_pole_generation", "polarity_generation",
                "phase_generation", "harmonic_generation",
                "backemf_generation", "torque_generation",
                "periodicity_generation", "mesh_generation",
                "owner_generation", "result_generation",
            )
        },
        "slot_count": 12,
        "result_slot_count": 12,
        "pole_count": 10,
        "result_pole_count": 10,
        "magnet_polarity_sequence": [1, -1] * 5,
        "result_magnet_polarity_sequence": [1, -1] * 5,
        "phase_sequence": "ABC",
        "result_phase_sequence": "ABC",
        "working_harmonic_order": 5,
        "result_working_harmonic_order": 5,
        "backemf_phase_angles_rad": phase_angles,
        "result_backemf_phase_angles_rad": phase_angles,
        "torque_samples_nm": torque_samples,
        "result_torque_samples_nm": torque_samples,
        "average_torque_nm": average_torque,
        "result_average_torque_nm": average_torque,
        "torque_ripple_ratio": torque_ripple,
        "result_torque_ripple_ratio": torque_ripple,
        "periodic_multiplier": 2,
        "result_periodic_multiplier": 2,
        "sector_slot_count": 6,
        "result_sector_slot_count": 6,
        "sector_pole_count": 5,
        "result_sector_pole_count": 5,
        "mesh_owner": "mesh:flux-switching-258",
        "accepted_mesh_owner": "mesh:flux-switching-258",
        "flux_switching_result_sha256": "2" * 64,
        "accepted_flux_switching_result_sha256": "2" * 64,
    }
    return payload


def test_v38_public_positive_wound_field_and_flux_switching_closure():
    assert pwm_controlled_motor_loss_gate(_payload_v38())["status"] == "ok"


def test_v38_public_wound_field_synchronous_excitation_flux_torque_powerfactor_fieldloss_energy_mismatch():
    payload = _payload_v38()
    row = payload["artifact_identity"][
        "wound_field_synchronous_excitation_flux_torque_angle_powerfactor_field_stator_loss_mechanical_energy_mesh_owner_result_identity"
    ]
    row.update(
        {
            "flux_generation": "wound-field-257",
            "result_field_current_a": -5.0,
            "result_excitation_flux_linkage_wb_turn": -0.2,
            "result_torque_angle_rad": math.pi,
            "result_electromagnetic_torque_nm": -1.0,
            "result_power_factor": 1.5,
            "result_field_copper_loss_w": -25.0,
            "result_energy_balance_residual_w": 100.0,
            "accepted_mesh_owner": "stale:wound-field",
        }
    )
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v38_public_flux_switching_pm_slot_pole_polarity_harmonic_backemf_torque_ripple_mismatch():
    payload = _payload_v38()
    row = payload["artifact_identity"][
        "flux_switching_pm_slot_pole_polarity_phase_harmonic_backemf_torque_ripple_periodicity_mesh_owner_result_identity"
    ]
    row.update(
        {
            "polarity_generation": "flux-switching-257",
            "result_slot_count": 10,
            "result_pole_count": 12,
            "result_magnet_polarity_sequence": [1] * 10,
            "result_phase_sequence": "ACB",
            "result_working_harmonic_order": 3,
            "result_backemf_phase_angles_rad": [0.0, 0.0, 0.0],
            "result_torque_ripple_ratio": -1.0,
            "result_periodic_multiplier": 1,
            "accepted_mesh_owner": "stale:flux-switching",
        }
    )
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_wrong_wound_field_torque():
    payload = _payload_v38()
    row = payload["artifact_identity"][
        "wound_field_synchronous_excitation_flux_torque_angle_powerfactor_field_stator_loss_mechanical_energy_mesh_owner_result_identity"
    ]
    row["electromagnetic_torque_nm"] *= 2.0
    row["result_electromagnetic_torque_nm"] = row["electromagnetic_torque_nm"]
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_wrong_flux_switching_periodicity():
    payload = _payload_v38()
    row = payload["artifact_identity"][
        "flux_switching_pm_slot_pole_polarity_phase_harmonic_backemf_torque_ripple_periodicity_mesh_owner_result_identity"
    ]
    row["periodic_multiplier"] = 1
    row["result_periodic_multiplier"] = 1
    row["sector_slot_count"] = 12
    row["result_sector_slot_count"] = 12
    row["sector_pole_count"] = 10
    row["result_sector_pole_count"] = 10
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"
