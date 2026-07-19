from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import pwm_controlled_motor_loss_gate
from test_motor_generalization_v39 import _payload_v39


_INDUCTION = "induction_cage_slip_rotorbar_current_loss_endring_airgap_torque_power_model_mesh_result_identity"
_FIELDWEAKENING = "ipm_fieldweakening_dq_flux_voltage_limit_current_angle_speed_torque_power_model_result_identity"
_PROMOTED_CASE_IDS = (
    "v40_public_induction_cage_rotorbar_slip_current_loss_torque_endring_power_mismatch",
    "v40_public_ipm_fieldweakening_dq_flux_voltage_limit_current_angle_torque_power_mismatch",
)


def _payload_v40():
    payload = _payload_v39()
    identity = payload["artifact_identity"]
    generation = "induction-cage-280"
    synchronous_speed = 157.07963267948966
    rotor_speed = 141.3716694115407
    slip = (synchronous_speed - rotor_speed) / synchronous_speed
    bar_count = 24
    bar_current = 80.0
    bar_resistance = 1.5e-4
    ring_current = 260.0
    ring_resistance = 3.0e-5
    bar_loss = bar_count * bar_current**2 * bar_resistance
    ring_loss = 2.0 * ring_current**2 * ring_resistance
    rotor_loss = bar_loss + ring_loss
    airgap_power = rotor_loss / slip
    torque = airgap_power / synchronous_speed
    mechanical_power = torque * rotor_speed
    identity[_INDUCTION] = {
        "induction_generation": generation,
        **{key: generation for key in ("slip_generation", "rotorbar_generation", "endring_generation", "loss_generation", "airgap_generation", "torque_generation", "power_generation", "model_generation", "mesh_generation", "result_generation")},
        "synchronous_mechanical_speed_rad_s": synchronous_speed,
        "result_synchronous_mechanical_speed_rad_s": synchronous_speed,
        "rotor_mechanical_speed_rad_s": rotor_speed,
        "result_rotor_mechanical_speed_rad_s": rotor_speed,
        "slip": slip,
        "result_slip": slip,
        "rotor_bar_count": bar_count,
        "result_rotor_bar_count": bar_count,
        "rotor_bar_current_rms_a": bar_current,
        "result_rotor_bar_current_rms_a": bar_current,
        "rotor_bar_resistance_ohm": bar_resistance,
        "result_rotor_bar_resistance_ohm": bar_resistance,
        "endring_current_rms_a": ring_current,
        "result_endring_current_rms_a": ring_current,
        "endring_segment_resistance_ohm": ring_resistance,
        "result_endring_segment_resistance_ohm": ring_resistance,
        "rotor_bar_loss_w": bar_loss,
        "result_rotor_bar_loss_w": bar_loss,
        "endring_loss_w": ring_loss,
        "result_endring_loss_w": ring_loss,
        "rotor_copper_loss_w": rotor_loss,
        "result_rotor_copper_loss_w": rotor_loss,
        "airgap_power_w": airgap_power,
        "result_airgap_power_w": airgap_power,
        "electromagnetic_torque_nm": torque,
        "result_electromagnetic_torque_nm": torque,
        "mechanical_power_w": mechanical_power,
        "result_mechanical_power_w": mechanical_power,
        "model_owner": "motor:induction-cage-280",
        "accepted_model_owner": "motor:induction-cage-280",
        "mesh_owner": "mesh:induction-cage-280",
        "accepted_mesh_owner": "mesh:induction-cage-280",
        "induction_result_sha256": "5" * 64,
        "accepted_induction_result_sha256": "5" * 64,
    }

    generation = "ipm-fieldweakening-280"
    pole_pairs = 4
    resistance = 0.05
    ld = 5.0e-3
    lq = 9.0e-3
    magnet_flux = 0.1
    current_d = -8.0
    current_q = 12.0
    electrical_speed = 400.0
    flux_d = magnet_flux + ld * current_d
    flux_q = lq * current_q
    voltage_d = resistance * current_d - electrical_speed * flux_q
    voltage_q = resistance * current_q + electrical_speed * flux_d
    current_magnitude = math.hypot(current_d, current_q)
    voltage_magnitude = math.hypot(voltage_d, voltage_q)
    current_angle = math.degrees(math.atan2(-current_d, current_q))
    torque = 1.5 * pole_pairs * (flux_d * current_q - flux_q * current_d)
    mechanical_speed = electrical_speed / pole_pairs
    copper_loss = 1.5 * resistance * current_magnitude**2
    electrical_power = 1.5 * (voltage_d * current_d + voltage_q * current_q)
    mechanical_power = torque * mechanical_speed
    values = {
        "phase_resistance_ohm": resistance, "ld_h": ld, "lq_h": lq,
        "magnet_flux_wb": magnet_flux, "current_d_a": current_d,
        "current_q_a": current_q, "flux_d_wb": flux_d, "flux_q_wb": flux_q,
        "electrical_speed_rad_s": electrical_speed,
        "mechanical_speed_rad_s": mechanical_speed, "voltage_d_v": voltage_d,
        "voltage_q_v": voltage_q, "current_magnitude_a": current_magnitude,
        "current_limit_a": 20.0, "voltage_magnitude_v": voltage_magnitude,
        "voltage_limit_v": 60.0, "current_angle_deg": current_angle,
        "electromagnetic_torque_nm": torque, "copper_loss_w": copper_loss,
        "electrical_power_w": electrical_power, "mechanical_power_w": mechanical_power,
    }
    identity[_FIELDWEAKENING] = {
        "fieldweakening_generation": generation,
        **{key: generation for key in ("flux_generation", "voltage_generation", "current_generation", "angle_generation", "speed_generation", "torque_generation", "power_generation", "model_generation", "result_generation")},
        "pole_pairs": pole_pairs,
        "result_pole_pairs": pole_pairs,
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "model_owner": "motor:ipm-fieldweakening-280",
        "accepted_model_owner": "motor:ipm-fieldweakening-280",
        "fieldweakening_result_sha256": "6" * 64,
        "accepted_fieldweakening_result_sha256": "6" * 64,
    }
    return payload


def test_v40_public_positive_induction_and_fieldweakening_closure():
    assert pwm_controlled_motor_loss_gate(_payload_v40())["status"] == "ok"


def test_v40_public_induction_cage_rotorbar_slip_current_loss_torque_endring_power_mismatch():
    payload = _payload_v40()
    payload["artifact_identity"][_INDUCTION].update({"slip_generation": "induction-cage-279", "result_slip": 0.2, "result_airgap_power_w": -1.0, "accepted_mesh_owner": "stale:mesh"})
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v40_public_ipm_fieldweakening_dq_flux_voltage_limit_current_angle_torque_power_mismatch():
    payload = _payload_v40()
    payload["artifact_identity"][_FIELDWEAKENING].update({"flux_generation": "ipm-fieldweakening-279", "result_flux_d_wb": -1.0, "result_voltage_d_v": 200.0, "result_electrical_power_w": -1.0, "accepted_model_owner": "stale:motor"})
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_wrong_induction_slip():
    payload = _payload_v40()
    row = payload["artifact_identity"][_INDUCTION]
    row["slip"] = row["result_slip"] = 0.2
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_fieldweakening_power_gap():
    payload = _payload_v40()
    row = payload["artifact_identity"][_FIELDWEAKENING]
    row["electrical_power_w"] = row["result_electrical_power_w"] = row["mechanical_power_w"]
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"
