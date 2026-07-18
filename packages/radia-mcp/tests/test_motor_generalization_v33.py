from __future__ import annotations

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_motor_generalization_v32 import _payload_v32


_PROMOTED_CASE_IDS = (
    "v33_public_srm_commutation_current_chop_dwell_overlap_coenergy_torque_loss_mismatch",
    "v33_public_axial_flux_pm_sector_airgap_end_effect_torque_axial_force_surface_mismatch",
)


def _payload_v33():
    payload = _payload_v32()
    identity = payload["artifact_identity"]
    generation = "srm-commutation-201"
    identity[
        "srm_commutation_phase_dwell_chop_overlap_coenergy_torque_loss_angle_mesh_result_identity"
    ] = {
        "srm_generation": generation,
        **{
            key: generation
            for key in (
                "phase_generation",
                "dwell_generation",
                "chop_generation",
                "overlap_generation",
                "coenergy_generation",
                "torque_generation",
                "loss_generation",
                "angle_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "phase_sequence": ["A", "B", "C"],
        "result_phase_sequence": ["A", "B", "C"],
        "turn_on_deg": [0.0, 30.0, 60.0],
        "result_turn_on_deg": [0.0, 30.0, 60.0],
        "turn_off_deg": [20.0, 50.0, 80.0],
        "result_turn_off_deg": [20.0, 50.0, 80.0],
        "current_chop_a": 100.0,
        "result_current_chop_a": 100.0,
        "overlap_deg": 5.0,
        "result_overlap_deg": 5.0,
        "angle_grid_rad": [0.0, 0.1, 0.2],
        "result_angle_grid_rad": [0.0, 0.1, 0.2],
        "coenergy_j": [0.0, 0.5, 1.0],
        "result_coenergy_j": [0.0, 0.5, 1.0],
        "torque_nm": [5.0, 5.0, 5.0],
        "result_torque_nm": [5.0, 5.0, 5.0],
        "copper_loss_w": 120.0,
        "result_copper_loss_w": 120.0,
        "iron_loss_w": 30.0,
        "result_iron_loss_w": 30.0,
        "total_loss_w": 150.0,
        "result_total_loss_w": 150.0,
        "mesh_sha256": "1" * 64,
        "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "axial-flux-pm-201"
    identity[
        "axial_flux_pm_sector_airgap_end_effect_torque_force_surface_direction_frame_mesh_result_identity"
    ] = {
        "axial_generation": generation,
        **{
            key: generation
            for key in (
                "sector_generation",
                "airgap_generation",
                "end_effect_generation",
                "torque_generation",
                "force_generation",
                "direction_generation",
                "frame_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "sector_multiplier": 12,
        "result_sector_multiplier": 12,
        "air_gaps_m": [0.001, 0.001],
        "result_air_gaps_m": [0.001, 0.001],
        "end_effect_factor": 0.96,
        "result_end_effect_factor": 0.96,
        "surface_coordinates": [[0.0, 0.0], [10.0, 100.0], [20.0, 200.0]],
        "result_surface_coordinates": [
            [0.0, 0.0],
            [10.0, 100.0],
            [20.0, 200.0],
        ],
        "torque_surface_nm": [40.0, 55.0, 60.0],
        "result_torque_surface_nm": [40.0, 55.0, 60.0],
        "axial_force_surface_n": [0.0, 12.0, 18.0],
        "result_axial_force_surface_n": [0.0, 12.0, 18.0],
        "force_direction": "+z",
        "result_force_direction": "+z",
        "axial_frame": "rotor_global_z",
        "result_axial_frame": "rotor_global_z",
        "mesh_sha256": "3" * 64,
        "result_mesh_sha256": "3" * 64,
        "result_lineage_sha256": "4" * 64,
        "accepted_result_lineage_sha256": "4" * 64,
    }
    return payload


def test_v33_public_positive_srm_commutation_and_axial_flux_pm_closure():
    assert pwm_controlled_motor_loss_gate(_payload_v33())["status"] == "ok"


def test_v33_public_srm_commutation_current_chop_dwell_overlap_coenergy_torque_loss_mismatch():
    payload = _payload_v33()
    record = payload["artifact_identity"][
        "srm_commutation_phase_dwell_chop_overlap_coenergy_torque_loss_angle_mesh_result_identity"
    ]
    record.update(
        {
            "phase_generation": "srm-commutation-200",
            "coenergy_generation": "srm-commutation-199",
            "result_generation": "srm-commutation-198",
            "result_phase_sequence": ["C", "B", "A"],
            "result_turn_on_deg": [5.0, 35.0, 65.0],
            "result_turn_off_deg": [10.0, 40.0, 70.0],
            "result_current_chop_a": 50.0,
            "result_overlap_deg": -5.0,
            "result_angle_grid_rad": [0.0, 0.2, 0.1],
            "result_coenergy_j": [0.0, 1.0, 0.5],
            "result_torque_nm": [-5.0, 0.0, 5.0],
            "result_copper_loss_w": 10.0,
            "result_iron_loss_w": 300.0,
            "result_total_loss_w": 100.0,
            "result_mesh_sha256": "8" * 64,
            "accepted_result_sha256": "9" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "srm_commutation_uses_current_phases_dwell_chop_overlap_coenergy_torque_loss_mesh_and_result"
    ]


def test_v33_public_axial_flux_pm_sector_airgap_end_effect_torque_axial_force_surface_mismatch():
    payload = _payload_v33()
    record = payload["artifact_identity"][
        "axial_flux_pm_sector_airgap_end_effect_torque_force_surface_direction_frame_mesh_result_identity"
    ]
    record.update(
        {
            "sector_generation": "axial-flux-pm-200",
            "force_generation": "axial-flux-pm-199",
            "result_generation": "axial-flux-pm-198",
            "result_sector_multiplier": 6,
            "result_air_gaps_m": [0.001, 0.002],
            "result_end_effect_factor": 1.2,
            "result_surface_coordinates": [[0.0, 0.0]],
            "result_torque_surface_nm": [5.0],
            "result_axial_force_surface_n": [-18.0],
            "result_force_direction": "-x",
            "result_axial_frame": "stale_local_y",
            "result_mesh_sha256": "a" * 64,
            "accepted_result_lineage_sha256": "b" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axial_flux_pm_uses_current_sector_airgaps_end_effect_torque_force_surface_frame_mesh_and_result"
    ]
