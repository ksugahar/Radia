from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_motor_generalization_v36 import _payload_v36


_PROMOTED_CASE_IDS = (
    "v37_public_induction_motor_slip_airgap_power_torque_rotor_loss_efficiency_owner_mismatch",
    "v37_public_axial_flux_motor_periodicity_airgap_flux_torque_ripple_backemf_owner_mismatch",
)


def _payload_v37():
    payload = _payload_v36()
    identity = payload["artifact_identity"]
    generation = "induction-power-246"
    frequency, pole_pairs, slip = 50.0, 2, 0.05
    synchronous_speed = 2.0 * math.pi * frequency / pole_pairs
    mechanical_speed = (1.0 - slip) * synchronous_speed
    airgap_power, input_power = 1000.0, 1200.0
    torque = airgap_power / synchronous_speed
    rotor_loss = slip * airgap_power
    mechanical_output = (1.0 - slip) * airgap_power
    identity[
        "induction_motor_slip_synchronous_speed_airgap_power_torque_rotor_loss_mechanical_output_efficiency_owner_result_identity"
    ] = {
        "induction_generation": generation,
        **{key: generation for key in (
            "frequency_generation", "speed_generation", "slip_generation",
            "airgap_generation", "torque_generation", "rotor_loss_generation",
            "mechanical_generation", "efficiency_generation", "owner_generation",
            "result_generation")},
        "electrical_frequency_hz": frequency, "result_electrical_frequency_hz": frequency,
        "pole_pairs": pole_pairs, "result_pole_pairs": pole_pairs,
        "slip": slip, "result_slip": slip,
        "synchronous_speed_rad_s": synchronous_speed,
        "result_synchronous_speed_rad_s": synchronous_speed,
        "mechanical_speed_rad_s": mechanical_speed,
        "result_mechanical_speed_rad_s": mechanical_speed,
        "airgap_power_w": airgap_power, "result_airgap_power_w": airgap_power,
        "electromagnetic_torque_nm": torque,
        "result_electromagnetic_torque_nm": torque,
        "rotor_copper_loss_w": rotor_loss,
        "result_rotor_copper_loss_w": rotor_loss,
        "mechanical_output_w": mechanical_output,
        "result_mechanical_output_w": mechanical_output,
        "input_power_w": input_power, "result_input_power_w": input_power,
        "efficiency": mechanical_output / input_power,
        "result_efficiency": mechanical_output / input_power,
        "motor_owner": "motor/induction-246",
        "accepted_motor_owner": "motor/induction-246",
        "motor_result_sha256": "1" * 64,
        "accepted_motor_result_sha256": "1" * 64,
    }
    generation = "axial-flux-246"
    sector_factor = 12
    sector_torque = [2.0, 2.1, 1.9]
    full_torque = [sector_factor * item for item in sector_torque]
    average_torque = sum(full_torque) / len(full_torque)
    ripple = (max(full_torque) - min(full_torque)) / average_torque
    phases = [0.0, -2.0 * math.pi / 3.0, 2.0 * math.pi / 3.0]
    identity[
        "axial_flux_motor_sector_periodicity_dual_airgap_axial_flux_torque_ripple_backemf_frame_mesh_owner_result_identity"
    ] = {
        "axial_flux_generation": generation,
        **{key: generation for key in (
            "sector_generation", "airgap_generation", "flux_generation",
            "torque_generation", "ripple_generation", "backemf_generation",
            "frame_generation", "mesh_generation", "owner_generation",
            "result_generation")},
        "sector_factor": sector_factor, "result_sector_factor": sector_factor,
        "sector_angle_rad": 2.0 * math.pi / sector_factor,
        "result_sector_angle_rad": 2.0 * math.pi / sector_factor,
        "dual_airgap_m": [0.001, 0.001],
        "result_dual_airgap_m": [0.001, 0.001],
        "sector_axial_flux_per_gap_wb": [0.002, 0.002],
        "result_sector_axial_flux_per_gap_wb": [0.002, 0.002],
        "sector_torque_samples_nm": sector_torque,
        "result_sector_torque_samples_nm": sector_torque,
        "full_machine_torque_samples_nm": full_torque,
        "result_full_machine_torque_samples_nm": full_torque,
        "average_torque_nm": average_torque,
        "result_average_torque_nm": average_torque,
        "torque_ripple_ratio": ripple, "result_torque_ripple_ratio": ripple,
        "backemf_phase_angles_rad": phases,
        "result_backemf_phase_angles_rad": phases,
        "coordinate_frame": "cylindrical_z_axial",
        "result_coordinate_frame": "cylindrical_z_axial",
        "mesh_owner": "mesh/axial-flux-246",
        "accepted_mesh_owner": "mesh/axial-flux-246",
        "axial_flux_result_sha256": "2" * 64,
        "accepted_axial_flux_result_sha256": "2" * 64,
    }
    return payload


def test_v37_public_positive_induction_and_axial_flux_closure():
    assert pwm_controlled_motor_loss_gate(_payload_v37())["status"] == "ok"


def test_v37_public_induction_motor_slip_airgap_power_torque_rotor_loss_efficiency_owner_mismatch():
    payload = _payload_v37()
    row = payload["artifact_identity"][
        "induction_motor_slip_synchronous_speed_airgap_power_torque_rotor_loss_mechanical_output_efficiency_owner_result_identity"
    ]
    row.update({
        "slip_generation": "induction-power-245", "torque_generation": "induction-power-244",
        "result_generation": "induction-power-243", "result_slip": 0.2,
        "result_synchronous_speed_rad_s": -1.0, "result_mechanical_speed_rad_s": 200.0,
        "result_airgap_power_w": -1000.0, "result_electromagnetic_torque_nm": -5.0,
        "result_rotor_copper_loss_w": 500.0, "result_mechanical_output_w": 100.0,
        "result_input_power_w": 100.0, "result_efficiency": 2.0,
        "accepted_motor_owner": "stale/motor", "accepted_motor_result_sha256": "a" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "induction_motor_closes_slip_speed_airgap_power_torque_rotor_loss_output_efficiency_owner_and_result"
    ]


def test_v37_public_axial_flux_motor_periodicity_airgap_flux_torque_ripple_backemf_owner_mismatch():
    payload = _payload_v37()
    row = payload["artifact_identity"][
        "axial_flux_motor_sector_periodicity_dual_airgap_axial_flux_torque_ripple_backemf_frame_mesh_owner_result_identity"
    ]
    row.update({
        "sector_generation": "axial-flux-245", "backemf_generation": "axial-flux-244",
        "result_generation": "axial-flux-243", "result_sector_factor": 6,
        "result_sector_angle_rad": math.pi, "result_dual_airgap_m": [0.001, -0.001],
        "result_sector_axial_flux_per_gap_wb": [0.002, -0.002],
        "result_sector_torque_samples_nm": [2.0, -2.1, 1.9],
        "result_full_machine_torque_samples_nm": [1.0, 2.0, 3.0],
        "result_average_torque_nm": -1.0, "result_torque_ripple_ratio": -0.5,
        "result_backemf_phase_angles_rad": [0.0, 0.0, 0.0],
        "result_coordinate_frame": "cartesian_x_axial",
        "accepted_mesh_owner": "stale/mesh",
        "accepted_axial_flux_result_sha256": "b" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axial_flux_motor_closes_sector_dual_airgap_flux_torque_ripple_backemf_frame_mesh_owner_and_result"
    ]


def test_v37_public_rejects_self_consistent_induction_power_imbalance():
    payload = _payload_v37()
    row = payload["artifact_identity"][
        "induction_motor_slip_synchronous_speed_airgap_power_torque_rotor_loss_mechanical_output_efficiency_owner_result_identity"
    ]
    row["rotor_copper_loss_w"] = row["result_rotor_copper_loss_w"] = 500.0
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v37_public_rejects_self_consistent_axial_flux_phase_collapse():
    payload = _payload_v37()
    row = payload["artifact_identity"][
        "axial_flux_motor_sector_periodicity_dual_airgap_axial_flux_torque_ripple_backemf_frame_mesh_owner_result_identity"
    ]
    row["backemf_phase_angles_rad"] = row["result_backemf_phase_angles_rad"] = [0.0, 0.0, 0.0]
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"
