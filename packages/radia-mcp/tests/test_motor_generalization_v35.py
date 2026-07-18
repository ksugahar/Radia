from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import pwm_controlled_motor_loss_gate
from test_motor_generalization_v34 import _payload_v34


_PROMOTED_CASE_IDS = (
    "v35_public_skew_slice_torque_harmonic_axial_weight_phase_periodicity_mismatch",
    "v35_public_ironloss_hysteresis_eddy_excess_waveform_frequency_volume_temperature_mismatch",
)


def _payload_v35():
    payload = _payload_v34()
    identity = payload["artifact_identity"]
    generation = "skew-slice-torque-235"
    angles = [-0.1, 0.0, 0.1]
    phases = [[harmonic, [harmonic * 4 * angle for angle in angles]] for harmonic in [1, 6]]
    identity["skew_slice_torque_angle_axial_weight_harmonic_phase_pole_periodicity_mean_ripple_mesh_owner_result_identity"] = {
        "skew_generation": generation,
        **{key: generation for key in (
            "slice_generation", "weight_generation", "phase_generation", "periodicity_generation",
            "torque_generation", "ripple_generation", "mesh_generation", "owner_generation",
            "result_generation")},
        "slice_angles_rad": angles, "result_slice_angles_rad": angles,
        "axial_weights": [0.25, 0.5, 0.25], "result_axial_weights": [0.25, 0.5, 0.25],
        "harmonic_phase_shifts_rad": phases, "result_harmonic_phase_shifts_rad": phases,
        "pole_pairs": 4, "result_pole_pairs": 4,
        "pole_periodicity_angle_rad": math.pi / 2.0, "result_pole_periodicity_angle_rad": math.pi / 2.0,
        "slice_mean_torque_nm": [47.0, 50.0, 49.0], "result_slice_mean_torque_nm": [47.0, 50.0, 49.0],
        "weighted_mean_torque_nm": 49.0, "result_weighted_mean_torque_nm": 49.0,
        "torque_ripple_spectrum_nm": [[0, 49.0], [6, 2.0]],
        "result_torque_ripple_spectrum_nm": [[0, 49.0], [6, 2.0]],
        "skew_mesh_sha256": "1" * 64, "result_skew_mesh_sha256": "1" * 64,
        "skew_result_owner": "motor/skew-235", "accepted_skew_result_owner": "motor/skew-235",
        "skew_result_sha256": "2" * 64, "accepted_skew_result_sha256": "2" * 64,
    }
    generation = "iron-loss-separation-235"
    frequency, b_peak, waveform_factor, temperature_factor = 100.0, 1.2, 1.1, 1.1
    components = [
        2.0 * frequency * b_peak**2 * waveform_factor * temperature_factor,
        0.1 * frequency**2 * b_peak**2 * waveform_factor,
        0.05 * frequency**1.5 * b_peak**1.5 * waveform_factor,
    ]
    identity["ironloss_hysteresis_eddy_excess_waveform_frequency_coeff_volume_temperature_total_owner_result_identity"] = {
        "ironloss_generation": generation,
        **{key: generation for key in (
            "component_generation", "waveform_generation", "frequency_generation", "coefficient_generation",
            "volume_generation", "temperature_generation", "total_generation", "owner_generation",
            "result_generation")},
        "b_waveform_peak_t": b_peak, "result_b_waveform_peak_t": b_peak,
        "waveform_factor": waveform_factor, "result_waveform_factor": waveform_factor,
        "frequency_hz": frequency, "result_frequency_hz": frequency,
        "material_coefficients": [2.0, 0.1, 0.05], "result_material_coefficients": [2.0, 0.1, 0.05],
        "active_volume_m3": 0.001, "result_active_volume_m3": 0.001,
        "temperature_c": 80.0, "result_temperature_c": 80.0,
        "temperature_factor": temperature_factor, "result_temperature_factor": temperature_factor,
        "loss_components_w_m3": components, "result_loss_components_w_m3": components,
        "total_iron_loss_w": sum(components) * 0.001, "result_total_iron_loss_w": sum(components) * 0.001,
        "waveform_sha256": "3" * 64, "result_waveform_sha256": "3" * 64,
        "ironloss_owner": "motor/iron-loss-235", "accepted_ironloss_owner": "motor/iron-loss-235",
        "ironloss_result_sha256": "4" * 64, "accepted_ironloss_result_sha256": "4" * 64,
    }
    return payload


def test_v35_public_positive_skew_and_ironloss_closure():
    assert pwm_controlled_motor_loss_gate(_payload_v35())["status"] == "ok"


def test_v35_public_skew_slice_torque_harmonic_axial_weight_phase_periodicity_mismatch():
    payload = _payload_v35()
    row = payload["artifact_identity"]["skew_slice_torque_angle_axial_weight_harmonic_phase_pole_periodicity_mean_ripple_mesh_owner_result_identity"]
    row.update({"slice_generation": "skew-slice-torque-234", "phase_generation": "skew-slice-torque-233",
                "result_generation": "skew-slice-torque-232", "result_slice_angles_rad": [0.1, 0.0, -0.1],
                "result_axial_weights": [0.8, 0.8, -0.6],
                "result_harmonic_phase_shifts_rad": [[1, [0.0, 0.0, 0.0]], [6, [1.0, 1.0, 1.0]]],
                "result_pole_pairs": 2, "result_pole_periodicity_angle_rad": math.pi,
                "result_slice_mean_torque_nm": [10.0, -20.0, 100.0], "result_weighted_mean_torque_nm": -49.0,
                "result_torque_ripple_spectrum_nm": [[0, -49.0], [5, 20.0]],
                "result_skew_mesh_sha256": "9" * 64, "accepted_skew_result_owner": "stale/skew",
                "accepted_skew_result_sha256": "a" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["skew_slice_torque_closes_angles_axial_weights_harmonic_phases_pole_periodicity_mean_ripple_mesh_owner_and_result"]


def test_v35_public_ironloss_hysteresis_eddy_excess_waveform_frequency_volume_temperature_mismatch():
    payload = _payload_v35()
    row = payload["artifact_identity"]["ironloss_hysteresis_eddy_excess_waveform_frequency_coeff_volume_temperature_total_owner_result_identity"]
    row.update({"component_generation": "iron-loss-separation-234", "temperature_generation": "iron-loss-separation-233",
                "result_generation": "iron-loss-separation-232", "result_b_waveform_peak_t": -1.2,
                "result_waveform_factor": 0.0, "result_frequency_hz": -100.0,
                "result_material_coefficients": [2.0, -0.1, 0.0], "result_active_volume_m3": -0.001,
                "result_temperature_c": 20.0, "result_temperature_factor": -1.0,
                "result_loss_components_w_m3": [100.0, -200.0, 0.0], "result_total_iron_loss_w": 99.0,
                "result_waveform_sha256": "b" * 64, "accepted_ironloss_owner": "stale/loss",
                "accepted_ironloss_result_sha256": "c" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["iron_loss_closes_hysteresis_eddy_excess_waveform_frequency_coefficients_volume_temperature_total_owner_and_result"]


def test_v35_rejects_self_consistent_skew_weight_sum_error():
    payload = _payload_v35()
    row = payload["artifact_identity"]["skew_slice_torque_angle_axial_weight_harmonic_phase_pole_periodicity_mean_ripple_mesh_owner_result_identity"]
    row["axial_weights"] = row["result_axial_weights"] = [0.5, 0.5, 0.5]
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v35_rejects_self_consistent_ironloss_total_error():
    payload = _payload_v35()
    row = payload["artifact_identity"]["ironloss_hysteresis_eddy_excess_waveform_frequency_coeff_volume_temperature_total_owner_result_identity"]
    row["total_iron_loss_w"] = row["result_total_iron_loss_w"] = 99.0
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"
