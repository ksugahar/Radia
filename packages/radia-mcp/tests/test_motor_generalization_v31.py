from __future__ import annotations

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import pwm_controlled_motor_loss_gate
from test_motor_generalization_v30 import _payload_v30


_PROMOTED_CASE_IDS = (
    "v31_public_pwm_iron_loss_sampling_carrier_sideband_angle_alias_energy_balance_mismatch",
    "v31_public_skew_slice_torque_weight_phase_offset_periodicity_ripple_mismatch",
)


def _payload_v31():
    payload = _payload_v30(); identity = payload["artifact_identity"]
    generation = "pwm-loss-181"
    identity["pwm_iron_loss_sampling_sideband_angle_alias_volume_energy_result_identity"] = {
        "loss_generation": generation,
        **{key: generation for key in ("sampling_loss_generation", "sideband_loss_generation", "angle_loss_generation", "alias_loss_generation", "volume_loss_generation", "energy_loss_generation", "mesh_loss_generation", "result_loss_generation")},
        "sample_period_s": 2.5e-6, "result_sample_period_s": 2.5e-6,
        "samples_per_fundamental_cycle": 4000, "result_samples_per_fundamental_cycle": 4000,
        "carrier_frequency_hz": 10000.0, "result_carrier_frequency_hz": 10000.0,
        "fundamental_frequency_hz": 100.0, "result_fundamental_frequency_hz": 100.0,
        "carrier_sidebands_hz": [9900.0, 10100.0], "result_carrier_sidebands_hz": [9900.0, 10100.0],
        "pole_pairs": 4, "result_pole_pairs": 4,
        "mechanical_angle_deg": [0.0, 5.0, 10.0], "result_mechanical_angle_deg": [0.0, 5.0, 10.0],
        "electrical_angle_deg": [0.0, 20.0, 40.0], "result_electrical_angle_deg": [0.0, 20.0, 40.0],
        "alias_filter": "nyquist_guard_and_sideband_keep", "result_alias_filter": "nyquist_guard_and_sideband_keep",
        "gross_volume_m3": 0.001, "result_gross_volume_m3": 0.001,
        "active_volume_m3": 0.00095, "result_active_volume_m3": 0.00095,
        "stacking_factor": 0.95, "result_stacking_factor": 0.95,
        "mean_iron_loss_w": 50.0, "result_mean_iron_loss_w": 50.0,
        "cycle_energy_j": 0.5, "result_cycle_energy_j": 0.5,
        "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    generation = "skew-torque-181"
    identity["skew_slice_torque_weight_axial_phase_periodicity_ripple_mesh_result_identity"] = {
        "skew_generation": generation,
        **{key: generation for key in ("weight_skew_generation", "axial_skew_generation", "phase_skew_generation", "periodicity_skew_generation", "torque_skew_generation", "ripple_skew_generation", "mesh_skew_generation", "result_skew_generation")},
        "axial_locations_m": [-0.01, 0.0, 0.01], "result_axial_locations_m": [-0.01, 0.0, 0.01],
        "slice_weights": [0.25, 0.5, 0.25], "result_slice_weights": [0.25, 0.5, 0.25],
        "electrical_phase_offsets_deg": [-2.0, 0.0, 2.0], "result_electrical_phase_offsets_deg": [-2.0, 0.0, 2.0],
        "periodic_wrap_electrical_deg": 360.0, "result_periodic_wrap_electrical_deg": 360.0,
        "slice_mean_torque_nm": [10.0, 10.2, 10.0], "result_slice_mean_torque_nm": [10.0, 10.2, 10.0],
        "skew_mean_torque_nm": 10.1, "result_skew_mean_torque_nm": 10.1,
        "slice_ripple_harmonics_nm": [{"6": 0.3}, {"6": 0.2}, {"6": 0.3}],
        "result_slice_ripple_harmonics_nm": [{"6": 0.3}, {"6": 0.2}, {"6": 0.3}],
        "skew_ripple_harmonics_nm": {"6": 0.25}, "result_skew_ripple_harmonics_nm": {"6": 0.25},
        "mesh_sha256": "3" * 64, "result_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return payload


def test_v31_public_positive_pwm_sampling_and_skew_torque():
    assert pwm_controlled_motor_loss_gate(_payload_v31())["status"] == "ok"


def test_v31_public_pwm_iron_loss_sampling_carrier_sideband_angle_alias_energy_balance_mismatch():
    payload = _payload_v31(); record = payload["artifact_identity"]["pwm_iron_loss_sampling_sideband_angle_alias_volume_energy_result_identity"]
    record.update({"sampling_loss_generation": "pwm-loss-180", "result_sample_period_s": 2.5e-5, "result_carrier_sidebands_hz": [9800.0, 10200.0], "result_pole_pairs": 3, "result_alias_filter": "none", "result_active_volume_m3": 0.001, "result_cycle_energy_j": 0.2, "accepted_result_sha256": "b" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["pwm_iron_loss_uses_current_sampling_sidebands_angles_alias_volume_energy_mesh_and_result"]


def test_v31_public_skew_slice_torque_weight_phase_offset_periodicity_ripple_mismatch():
    payload = _payload_v31(); record = payload["artifact_identity"]["skew_slice_torque_weight_axial_phase_periodicity_ripple_mesh_result_identity"]
    record.update({"weight_skew_generation": "skew-torque-180", "result_slice_weights": [1.0, 1.0, 1.0], "result_electrical_phase_offsets_deg": [2.0, 0.0, -2.0], "result_periodic_wrap_electrical_deg": 180.0, "result_skew_mean_torque_nm": 30.2, "result_skew_ripple_harmonics_nm": {"6": 0.8}, "accepted_result_sha256": "c" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["skew_torque_uses_current_slice_weights_axial_phase_periodicity_ripple_mesh_and_result"]
