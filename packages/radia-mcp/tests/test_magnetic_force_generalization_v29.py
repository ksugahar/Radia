from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.force_coenergy_gate import force_coenergy_displacement_gate
from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from test_force_coenergy_gate import _artifact_identity, _quadratic_case
from test_magnetic_force_generalization_v28 import _summary_v28


_PROMOTED_CASE_IDS = (
    "v29_public_airgap_stress_harmonic_sector_periodicity_sampling_alias_torque_mismatch",
    "v29_public_laminated_core_hysteresis_eddy_excess_loss_frequency_flux_volume_mismatch",
    "v29_public_magnet_demag_recoil_knee_temperature_local_field_volume_fraction_mismatch",
    "v29_public_linear_motor_end_effect_phase_sequence_travel_wave_force_ripple_pitch_mismatch",
)


def _summary_v29():
    summary = _summary_v28()
    identity = summary["artifact_identity"]
    identity["magnet_demag_recoil_knee_field_volume_generation_identity"] = {
        "demag_generation": "demag-state-331",
        "recoil_demag_generation": "demag-state-331",
        "knee_demag_generation": "demag-state-331",
        "temperature_demag_generation": "demag-state-331",
        "field_demag_generation": "demag-state-331",
        "mask_demag_generation": "demag-state-331",
        "volume_demag_generation": "demag-state-331",
        "mesh_demag_generation": "demag-state-331",
        "result_demag_generation": "demag-state-331",
        "recoil_relative_permeability": 1.05,
        "result_recoil_relative_permeability": 1.05,
        "knee_field_a_m": -700000.0,
        "result_knee_field_a_m": -700000.0,
        "temperature_c": 120.0,
        "result_temperature_c": 120.0,
        "element_ids": [1, 2, 3],
        "result_element_ids": [1, 2, 3],
        "local_recoil_axis_field_a_m": [-600000.0, -800000.0, -500000.0],
        "result_local_recoil_axis_field_a_m": [-600000.0, -800000.0, -500000.0],
        "irreversible_mask": [False, True, False],
        "result_irreversible_mask": [False, True, False],
        "element_volumes_m3": [1.0e-6, 2.0e-6, 1.0e-6],
        "result_element_volumes_m3": [1.0e-6, 2.0e-6, 1.0e-6],
        "magnet_volume_m3": 4.0e-6,
        "result_magnet_volume_m3": 4.0e-6,
        "irreversible_volume_fraction": 0.5,
        "result_irreversible_volume_fraction": 0.5,
        "material_state_sha256": "1" * 64,
        "result_material_state_sha256": "1" * 64,
        "mesh_sha256": "2" * 64,
        "result_mesh_sha256": "2" * 64,
        "result_sha256": "3" * 64,
        "accepted_result_sha256": "3" * 64,
    }
    identity["linear_motor_end_phase_wave_pitch_force_generation_identity"] = {
        "linear_motor_generation": "linear-force-331",
        "end_effect_linear_motor_generation": "linear-force-331",
        "phase_linear_motor_generation": "linear-force-331",
        "wave_linear_motor_generation": "linear-force-331",
        "pitch_linear_motor_generation": "linear-force-331",
        "position_linear_motor_generation": "linear-force-331",
        "force_linear_motor_generation": "linear-force-331",
        "ripple_linear_motor_generation": "linear-force-331",
        "result_linear_motor_generation": "linear-force-331",
        "phase_sequence": ["U", "V", "W"],
        "result_phase_sequence": ["U", "V", "W"],
        "traveling_wave_direction": "global-x-positive",
        "result_traveling_wave_direction": "global-x-positive",
        "pole_pitch_m": 0.06,
        "result_pole_pitch_m": 0.06,
        "position_m": [0.0, 0.015, 0.03, 0.045, 0.06],
        "result_position_m": [0.0, 0.015, 0.03, 0.045, 0.06],
        "end_effect_factor": [0.9, 1.0, 1.0, 1.0, 0.9],
        "result_end_effect_factor": [0.9, 1.0, 1.0, 1.0, 0.9],
        "force_n": [100.0, 110.0, 100.0, 90.0, 100.0],
        "result_force_n": [100.0, 110.0, 100.0, 90.0, 100.0],
        "mean_force_n": 100.0,
        "result_mean_force_n": 100.0,
        "force_ripple_peak_to_peak_n": 20.0,
        "result_force_ripple_peak_to_peak_n": 20.0,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def test_v29_public_positive_demag_and_linear_motor_identities():
    assert magnetic_force_method_profile_gate(_summary_v29())["status"] == "ok"


def test_v29_public_magnet_demag_recoil_knee_temperature_local_field_volume_fraction_mismatch():
    summary = _summary_v29()
    identity = summary["artifact_identity"][
        "magnet_demag_recoil_knee_field_volume_generation_identity"
    ]
    identity.update({
        "recoil_demag_generation": "demag-state-330",
        "field_demag_generation": "demag-state-329",
        "result_recoil_relative_permeability": 1.2,
        "result_knee_field_a_m": -500000.0,
        "result_temperature_c": 20.0,
        "result_element_ids": [3, 2, 1],
        "result_local_recoil_axis_field_a_m": [-400000.0, -450000.0, -300000.0],
        "result_irreversible_mask": [False, False, False],
        "result_element_volumes_m3": [1.0e-6, 1.0e-6, 1.0e-6],
        "result_magnet_volume_m3": 3.0e-6,
        "result_irreversible_volume_fraction": 0.0,
        "result_material_state_sha256": "9" * 64,
        "result_mesh_sha256": "a" * 64,
        "accepted_result_sha256": "b" * 64,
    })
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "magnet_demag_uses_current_recoil_knee_temperature_local_field_mask_volume_and_result"
    ]


def test_v29_public_linear_motor_end_effect_phase_sequence_travel_wave_force_ripple_pitch_mismatch():
    summary = _summary_v29()
    identity = summary["artifact_identity"][
        "linear_motor_end_phase_wave_pitch_force_generation_identity"
    ]
    identity.update({
        "phase_linear_motor_generation": "linear-force-330",
        "end_effect_linear_motor_generation": "linear-force-329",
        "result_phase_sequence": ["U", "W", "V"],
        "result_traveling_wave_direction": "global-x-negative",
        "result_pole_pitch_m": 0.08,
        "result_position_m": [0.0, 0.02, 0.04, 0.08],
        "result_end_effect_factor": [1.0, 1.0, 1.0, 1.0],
        "result_force_n": [80.0, 120.0, 70.0, 130.0],
        "result_mean_force_n": 90.0,
        "result_force_ripple_peak_to_peak_n": 60.0,
        "accepted_result_sha256": "c" * 64,
    })
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "linear_motor_force_uses_current_end_effect_phase_wave_pitch_positions_ripple_and_result"
    ]


def _force_summary_v29():
    positions, _, _ = _quadratic_case()
    identity = _artifact_identity(len(positions))
    generation = "airgap-stress-161"
    identity["airgap_stress_harmonic_sector_periodicity_origin_sampling_alias_radius_torque_generation_identity"] = {
        "stress_generation": generation, "harmonic_stress_generation": generation,
        "sector_stress_generation": generation, "sampling_stress_generation": generation,
        "alias_stress_generation": generation, "geometry_stress_generation": generation,
        "mesh_stress_generation": generation, "result_stress_generation": generation,
        "sector_pitch_deg": 30.0, "result_sector_pitch_deg": 30.0,
        "sector_count": 12, "result_sector_count": 12,
        "angular_origin_deg": 0.0, "result_angular_origin_deg": 0.0,
        "angular_sample_count": 720, "result_angular_sample_count": 720,
        "sector_sample_count": 60, "result_sector_sample_count": 60,
        "harmonic_orders": [0, 6, 12, 18], "result_harmonic_orders": [0, 6, 12, 18],
        "torque_harmonics_nm": [5.0, 0.2, 0.05, 0.01],
        "result_torque_harmonics_nm": [5.0, 0.2, 0.05, 0.01],
        "alias_filter": "truncate_below_nyquist", "result_alias_filter": "truncate_below_nyquist",
        "alias_cutoff_order": 24, "result_alias_cutoff_order": 24,
        "airgap_radius_m": 0.05, "result_airgap_radius_m": 0.05,
        "axial_length_m": 0.1, "result_axial_length_m": 0.1,
        "torque_nm": 5.26, "result_torque_nm": 5.26,
        "airgap_mesh_sha256": "1" * 64, "result_airgap_mesh_sha256": "1" * 64,
        "torque_result_sha256": "2" * 64, "accepted_torque_result_sha256": "2" * 64,
    }
    frequency, flux, thickness, volume = 400.0, 1.2, 0.00035, 0.001
    kh, alpha, ke, kex = 50.0, 1.6, 0.02, 0.5
    hysteresis = kh * frequency * flux**alpha * volume
    eddy = ke * frequency**2 * flux**2 * thickness**2 * volume
    excess = kex * frequency**1.5 * flux**1.5 * volume
    generation = "laminated-loss-161"
    identity["laminated_core_hysteresis_eddy_excess_frequency_flux_lamination_volume_result_generation_identity"] = {
        "loss_generation": generation, "hysteresis_loss_generation": generation,
        "eddy_loss_generation": generation, "excess_loss_generation": generation,
        "frequency_loss_generation": generation, "flux_loss_generation": generation,
        "lamination_loss_generation": generation, "volume_loss_generation": generation,
        "result_loss_generation": generation,
        "frequency_hz": frequency, "result_frequency_hz": frequency,
        "peak_flux_density_t": flux, "result_peak_flux_density_t": flux,
        "lamination_thickness_m": thickness, "result_lamination_thickness_m": thickness,
        "magnetic_volume_m3": volume, "result_magnetic_volume_m3": volume,
        "hysteresis_coefficient": kh, "result_hysteresis_coefficient": kh,
        "hysteresis_exponent": alpha, "result_hysteresis_exponent": alpha,
        "eddy_coefficient": ke, "result_eddy_coefficient": ke,
        "excess_coefficient": kex, "result_excess_coefficient": kex,
        "hysteresis_loss_w": hysteresis, "result_hysteresis_loss_w": hysteresis,
        "eddy_loss_w": eddy, "result_eddy_loss_w": eddy,
        "excess_loss_w": excess, "result_excess_loss_w": excess,
        "total_core_loss_w": hysteresis + eddy + excess,
        "result_total_core_loss_w": hysteresis + eddy + excess,
        "material_sha256": "3" * 64, "result_material_sha256": "3" * 64,
        "loss_result_sha256": "4" * 64, "accepted_loss_result_sha256": "4" * 64,
    }
    return identity


def _force_gate(identity):
    positions, coenergy, forces = _quadratic_case()
    return force_coenergy_displacement_gate(
        positions, coenergy, forces, artifact_identity=identity
    )


def test_v29_public_positive_airgap_harmonics_and_laminated_loss():
    assert _force_gate(_force_summary_v29())["status"] == "ok"


def test_v29_public_airgap_stress_harmonic_sector_periodicity_sampling_alias_torque_mismatch():
    summary = _force_summary_v29()
    summary["airgap_stress_harmonic_sector_periodicity_origin_sampling_alias_radius_torque_generation_identity"].update({
        "harmonic_stress_generation": "airgap-stress-160", "sampling_stress_generation": "airgap-stress-159",
        "result_stress_generation": "airgap-stress-158", "result_sector_pitch_deg": 45.0,
        "result_sector_count": 10, "result_angular_origin_deg": 15.0,
        "result_angular_sample_count": 30, "result_sector_sample_count": 3,
        "result_harmonic_orders": [0, 7, 13, 25], "result_torque_harmonics_nm": [5.0, -1.0],
        "result_alias_filter": "none", "result_alias_cutoff_order": 400,
        "result_airgap_radius_m": 0.06, "result_axial_length_m": 0.08,
        "result_torque_nm": 8.0, "result_airgap_mesh_sha256": "b" * 64,
        "accepted_torque_result_sha256": "c" * 64,
    })
    result = _force_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["airgap_torque_uses_current_sector_sampling_alias_harmonics_geometry_mesh_and_result"]


def test_v29_public_laminated_core_hysteresis_eddy_excess_loss_frequency_flux_volume_mismatch():
    summary = _force_summary_v29()
    summary["laminated_core_hysteresis_eddy_excess_frequency_flux_lamination_volume_result_generation_identity"].update({
        "hysteresis_loss_generation": "laminated-loss-160", "frequency_loss_generation": "laminated-loss-159",
        "result_loss_generation": "laminated-loss-158", "result_frequency_hz": 50.0,
        "result_peak_flux_density_t": 0.8, "result_lamination_thickness_m": 0.0005,
        "result_magnetic_volume_m3": 0.01, "result_hysteresis_coefficient": 30.0,
        "result_hysteresis_exponent": 2.0, "result_eddy_coefficient": 0.2,
        "result_excess_coefficient": 2.0, "result_hysteresis_loss_w": 1.0,
        "result_eddy_loss_w": 2.0, "result_excess_loss_w": -1.0,
        "result_total_core_loss_w": 99.0, "result_material_sha256": "d" * 64,
        "accepted_loss_result_sha256": "e" * 64,
    })
    result = _force_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["laminated_core_loss_uses_current_frequency_flux_lamination_volume_components_and_result"]
