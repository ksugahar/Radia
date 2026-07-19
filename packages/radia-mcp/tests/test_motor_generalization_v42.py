from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_motor_generalization_v41 import _payload_v41


_IPM = (
    "ipm_dq_current_flux_torque_voltage_powerfactor_mtpv_energy_mesh_"
    "result_generation_identity"
)
_INDUCTION = (
    "inductionmotor_slip_rotorfrequency_copperloss_torque_airgappower_"
    "mechanicalpower_efficiency_result_generation_identity"
)
_PROMOTED_CASE_IDS = (
    "v42_public_ipm_dq_fluxmap_torque_voltage_powerfactor_mtpv_energy_mismatch",
    "v42_public_inductionmotor_slip_rotorfrequency_copperloss_torque_airgap_power_efficiency_mismatch",
)


def _payload_v42():
    payload = _payload_v41()
    identity = payload["artifact_identity"]
    generation = "ipm-dq-842"
    pole_pairs = 4
    current_d = -20.0
    current_q = 50.0
    pm_flux = 0.1
    inductance_d = 1.0e-3
    inductance_q = 2.0e-3
    flux_d = pm_flux + inductance_d * current_d
    flux_q = inductance_q * current_q
    torque = 1.5 * pole_pairs * (flux_d * current_q - flux_q * current_d)
    resistance = 0.05
    electrical_speed = 1000.0
    voltage_d = resistance * current_d - electrical_speed * flux_q
    voltage_q = resistance * current_q + electrical_speed * flux_d
    voltage_magnitude = math.hypot(voltage_d, voltage_q)
    voltage_limit = 150.0
    active_power = 1.5 * (voltage_d * current_d + voltage_q * current_q)
    apparent_power = 1.5 * voltage_magnitude * math.hypot(current_d, current_q)
    power_factor = active_power / apparent_power
    field_energy = 0.5 * (inductance_d * current_d**2 + inductance_q * current_q**2)
    coenergy = pm_flux * current_d + field_energy
    identity[_IPM] = {
        "dq_generation": generation,
        **{
            key: generation
            for key in (
                "current_generation", "flux_generation", "torque_generation",
                "voltage_generation", "powerfactor_generation", "mtpv_generation",
                "energy_generation", "mesh_generation", "result_generation",
            )
        },
        "pole_pairs": pole_pairs,
        "result_pole_pairs": pole_pairs,
        "current_d_a": current_d,
        "result_current_d_a": current_d,
        "current_q_a": current_q,
        "result_current_q_a": current_q,
        "pm_flux_linkage_wb": pm_flux,
        "result_pm_flux_linkage_wb": pm_flux,
        "inductance_d_h": inductance_d,
        "result_inductance_d_h": inductance_d,
        "inductance_q_h": inductance_q,
        "result_inductance_q_h": inductance_q,
        "flux_d_wb": flux_d,
        "result_flux_d_wb": flux_d,
        "flux_q_wb": flux_q,
        "result_flux_q_wb": flux_q,
        "torque_nm": torque,
        "result_torque_nm": torque,
        "phase_resistance_ohm": resistance,
        "result_phase_resistance_ohm": resistance,
        "electrical_speed_rad_s": electrical_speed,
        "result_electrical_speed_rad_s": electrical_speed,
        "voltage_d_v": voltage_d,
        "result_voltage_d_v": voltage_d,
        "voltage_q_v": voltage_q,
        "result_voltage_q_v": voltage_q,
        "voltage_magnitude_v": voltage_magnitude,
        "result_voltage_magnitude_v": voltage_magnitude,
        "voltage_limit_v": voltage_limit,
        "result_voltage_limit_v": voltage_limit,
        "active_power_w": active_power,
        "result_active_power_w": active_power,
        "apparent_power_va": apparent_power,
        "result_apparent_power_va": apparent_power,
        "power_factor": power_factor,
        "result_power_factor": power_factor,
        "mtpv_branch": "negative_id_high_speed",
        "result_mtpv_branch": "negative_id_high_speed",
        "mtpv_voltage_margin_v": voltage_limit - voltage_magnitude,
        "result_mtpv_voltage_margin_v": voltage_limit - voltage_magnitude,
        "field_energy_j": field_energy,
        "result_field_energy_j": field_energy,
        "coenergy_j": coenergy,
        "result_coenergy_j": coenergy,
        "mesh_owner": "mesh:ipm-dq-842",
        "accepted_mesh_owner": "mesh:ipm-dq-842",
        "ipm_result_sha256": "1" * 64,
        "accepted_ipm_result_sha256": "1" * 64,
    }

    generation = "induction-power-842"
    pole_pairs = 2
    supply_frequency = 50.0
    synchronous_speed = 2.0 * math.pi * supply_frequency / pole_pairs
    rotor_speed = 150.0
    slip = (synchronous_speed - rotor_speed) / synchronous_speed
    rotor_frequency = slip * supply_frequency
    torque = 20.0
    airgap_power = torque * synchronous_speed
    rotor_loss = slip * airgap_power
    converted = (1.0 - slip) * airgap_power
    mechanical_loss = 30.0
    mechanical_power = converted - mechanical_loss
    input_power = airgap_power + 100.0 + 50.0
    efficiency = mechanical_power / input_power
    identity[_INDUCTION] = {
        "induction_generation": generation,
        **{
            key: generation
            for key in (
                "slip_generation", "frequency_generation", "loss_generation",
                "torque_generation", "airgap_power_generation",
                "mechanical_power_generation", "efficiency_generation", "result_generation",
            )
        },
        "pole_pairs": pole_pairs,
        "result_pole_pairs": pole_pairs,
        "supply_frequency_hz": supply_frequency,
        "result_supply_frequency_hz": supply_frequency,
        "synchronous_speed_rad_s": synchronous_speed,
        "result_synchronous_speed_rad_s": synchronous_speed,
        "rotor_speed_rad_s": rotor_speed,
        "result_rotor_speed_rad_s": rotor_speed,
        "slip": slip,
        "result_slip": slip,
        "rotor_electrical_frequency_hz": rotor_frequency,
        "result_rotor_electrical_frequency_hz": rotor_frequency,
        "torque_nm": torque,
        "result_torque_nm": torque,
        "airgap_power_w": airgap_power,
        "result_airgap_power_w": airgap_power,
        "stator_copper_loss_w": 100.0,
        "result_stator_copper_loss_w": 100.0,
        "rotor_copper_loss_w": rotor_loss,
        "result_rotor_copper_loss_w": rotor_loss,
        "core_loss_w": 50.0,
        "result_core_loss_w": 50.0,
        "mechanical_loss_w": mechanical_loss,
        "result_mechanical_loss_w": mechanical_loss,
        "converted_power_w": converted,
        "result_converted_power_w": converted,
        "mechanical_power_w": mechanical_power,
        "result_mechanical_power_w": mechanical_power,
        "input_power_w": input_power,
        "result_input_power_w": input_power,
        "efficiency": efficiency,
        "result_efficiency": efficiency,
        "motor_result_sha256": "2" * 64,
        "accepted_motor_result_sha256": "2" * 64,
    }
    return payload


def test_v42_public_positive_ipm_and_induction_power_closure():
    assert pwm_controlled_motor_loss_gate(_payload_v42())["status"] == "ok"


def test_v42_public_ipm_dq_mismatch():
    payload = _payload_v42()
    payload["artifact_identity"][_IPM].update(
        {
            "flux_generation": "ipm-dq-841",
            "energy_generation": "ipm-dq-840",
            "result_generation": "ipm-dq-839",
            "result_current_d_a": 20.0,
            "result_flux_d_wb": -0.08,
            "result_torque_nm": -36.0,
            "result_voltage_d_v": 101.0,
            "result_voltage_magnitude_v": 200.0,
            "result_power_factor": -0.5,
            "result_mtpv_branch": "stale_branch",
            "result_mtpv_voltage_margin_v": -50.0,
            "result_field_energy_j": -2.7,
            "result_coenergy_j": -0.7,
            "accepted_mesh_owner": "stale:mesh",
            "accepted_ipm_result_sha256": "9" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "ipm_dq_maps_close_currents_flux_torque_voltage_powerfactor_mtpv_energy_mesh_and_result"
    ]


def test_v42_public_induction_motor_mismatch():
    payload = _payload_v42()
    payload["artifact_identity"][_INDUCTION].update(
        {
            "slip_generation": "induction-power-841",
            "loss_generation": "induction-power-840",
            "result_generation": "induction-power-839",
            "result_slip": -0.1,
            "result_rotor_electrical_frequency_hz": -5.0,
            "result_torque_nm": -20.0,
            "result_airgap_power_w": -1.0,
            "result_rotor_copper_loss_w": -1.0,
            "result_converted_power_w": -1.0,
            "result_mechanical_power_w": -1.0,
            "result_efficiency": 1.5,
            "accepted_motor_result_sha256": "a" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "induction_motors_close_slip_rotor_frequency_losses_torque_airgap_mechanical_power_efficiency_and_result"
    ]


def test_v42_public_rejects_self_consistent_wrong_ipm_torque():
    payload = _payload_v42()
    record = payload["artifact_identity"][_IPM]
    record["torque_nm"] = 30.0
    record["result_torque_nm"] = 30.0
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v42_public_rejects_self_consistent_wrong_rotor_loss():
    payload = _payload_v42()
    record = payload["artifact_identity"][_INDUCTION]
    record["rotor_copper_loss_w"] *= 2.0
    record["result_rotor_copper_loss_w"] = record["rotor_copper_loss_w"]
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"
