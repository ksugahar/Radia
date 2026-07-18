from __future__ import annotations

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_motor_generalization_v31 import _payload_v31


_PROMOTED_CASE_IDS = (
    "v32_public_ipm_demagnetization_knee_temperature_current_angle_region_fraction_mismatch",
    "v32_public_synrm_dq_map_angle_saturation_cross_coupling_mtpa_torque_derivative_mismatch",
)


def _payload_v32():
    payload = _payload_v31()
    identity = payload["artifact_identity"]
    generation = "ipm-demag-191"
    identity[
        "ipm_demagnetization_knee_temperature_current_angle_region_fraction_mesh_result_identity"
    ] = {
        "demag_generation": generation,
        **{
            key: generation
            for key in (
                "knee_demag_generation",
                "temperature_demag_generation",
                "current_demag_generation",
                "angle_demag_generation",
                "region_demag_generation",
                "fraction_demag_generation",
                "mesh_demag_generation",
                "result_demag_generation",
            )
        },
        "knee_criterion": "b_parallel_below_temperature_knee",
        "result_knee_criterion": "b_parallel_below_temperature_knee",
        "magnet_temperature_c": 120.0,
        "result_magnet_temperature_c": 120.0,
        "phase_current_rms_a": 200.0,
        "result_phase_current_rms_a": 200.0,
        "current_angle_electrical_deg": 135.0,
        "result_current_angle_electrical_deg": 135.0,
        "irreversible_region_labels": ["magnet_1/edge", "magnet_2/edge"],
        "result_irreversible_region_labels": ["magnet_1/edge", "magnet_2/edge"],
        "demagnetized_fraction": 0.015,
        "result_demagnetized_fraction": 0.015,
        "operating_point_owner": "ipm/case-191/120C/200A/135deg",
        "result_operating_point_owner": "ipm/case-191/120C/200A/135deg",
        "mesh_sha256": "1" * 64,
        "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "synrm-dq-191"
    rows = [
        {"id_a": -100.0, "iq_a": 100.0, "psi_d_wb": 0.18, "psi_q_wb": 0.07},
        {"id_a": -50.0, "iq_a": 150.0, "psi_d_wb": 0.21, "psi_q_wb": 0.09},
        {"id_a": 0.0, "iq_a": 180.0, "psi_d_wb": 0.24, "psi_q_wb": 0.11},
    ]
    identity[
        "synrm_dq_map_angle_saturation_cross_coupling_mtpa_torque_mesh_result_identity"
    ] = {
        "map_generation": generation,
        **{
            key: generation
            for key in (
                "angle_map_generation",
                "saturation_map_generation",
                "cross_coupling_map_generation",
                "mtpa_map_generation",
                "torque_map_generation",
                "mesh_map_generation",
                "result_map_generation",
            )
        },
        "electrical_angle_deg": [0.0, 30.0, 60.0, 90.0],
        "result_electrical_angle_deg": [0.0, 30.0, 60.0, 90.0],
        "saturation_branch": "nonlinear_forward",
        "result_saturation_branch": "nonlinear_forward",
        "flux_map_rows": rows,
        "result_flux_map_rows": [dict(row) for row in rows],
        "dpsi_d_diq_h": [1.5e-4, 1.7e-4, 1.9e-4],
        "result_dpsi_d_diq_h": [1.5e-4, 1.7e-4, 1.9e-4],
        "dpsi_q_did_h": [1.5e-4, 1.7e-4, 1.9e-4],
        "result_dpsi_q_did_h": [1.5e-4, 1.7e-4, 1.9e-4],
        "mtpa_row_indices": [0, 1, 2],
        "result_mtpa_row_indices": [0, 1, 2],
        "pole_pairs": 2,
        "result_pole_pairs": 2,
        "torque_reconstruction": "1.5*p*(psi_d*iq-psi_q*id)",
        "result_torque_reconstruction": "1.5*p*(psi_d*iq-psi_q*id)",
        "torque_nm": [75.0, 81.0, 86.4],
        "result_torque_nm": [75.0, 81.0, 86.4],
        "mesh_sha256": "3" * 64,
        "result_mesh_sha256": "3" * 64,
        "map_sha256": "4" * 64,
        "accepted_map_sha256": "4" * 64,
    }
    return payload


def test_v32_public_positive_ipm_demagnetization_and_synrm_dq_map_closure():
    assert pwm_controlled_motor_loss_gate(_payload_v32())["status"] == "ok"


def test_v32_public_ipm_demagnetization_knee_temperature_current_angle_region_fraction_mismatch():
    payload = _payload_v32()
    record = payload["artifact_identity"][
        "ipm_demagnetization_knee_temperature_current_angle_region_fraction_mesh_result_identity"
    ]
    record.update(
        {
            "knee_demag_generation": "ipm-demag-190",
            "temperature_demag_generation": "ipm-demag-189",
            "result_demag_generation": "ipm-demag-188",
            "result_knee_criterion": "b_magnitude_below_room_temperature_knee",
            "result_magnet_temperature_c": 20.0,
            "result_phase_current_rms_a": 120.0,
            "result_current_angle_electrical_deg": 90.0,
            "result_irreversible_region_labels": ["magnet_3/center"],
            "result_demagnetized_fraction": 0.15,
            "result_operating_point_owner": "ipm/old-case",
            "result_mesh_sha256": "8" * 64,
            "accepted_result_sha256": "9" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "ipm_demagnetization_uses_current_knee_temperature_current_angle_regions_fraction_mesh_owner_and_result"
    ]


def test_v32_public_synrm_dq_map_angle_saturation_cross_coupling_mtpa_torque_derivative_mismatch():
    payload = _payload_v32()
    record = payload["artifact_identity"][
        "synrm_dq_map_angle_saturation_cross_coupling_mtpa_torque_mesh_result_identity"
    ]
    record.update(
        {
            "angle_map_generation": "synrm-dq-190",
            "cross_coupling_map_generation": "synrm-dq-189",
            "result_map_generation": "synrm-dq-188",
            "result_electrical_angle_deg": [0.0, 15.0, 30.0, 45.0],
            "result_saturation_branch": "linearized",
            "result_flux_map_rows": [],
            "result_dpsi_d_diq_h": [0.0, 0.0, 0.0],
            "result_dpsi_q_did_h": [-1.0e-3, -1.0e-3, -1.0e-3],
            "result_mtpa_row_indices": [2, 1, 0],
            "result_pole_pairs": 3,
            "result_torque_reconstruction": "1.5*(psi_d*iq-psi_q*id)",
            "result_torque_nm": [10.0, 20.0, 30.0],
            "result_mesh_sha256": "a" * 64,
            "accepted_map_sha256": "b" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "synrm_dq_map_uses_current_angles_saturation_cross_coupling_mtpa_torque_mesh_and_result"
    ]
