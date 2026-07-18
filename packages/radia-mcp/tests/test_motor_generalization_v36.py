from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import pwm_controlled_motor_loss_gate
from test_motor_generalization_v35 import _payload_v35


_PROMOTED_CASE_IDS = (
    "v36_public_dq_flux_inductance_torque_mtpa_current_angle_speed_convention_mismatch",
    "v36_public_iron_loss_hysteresis_eddy_excess_frequency_flux_energy_balance_mismatch",
)


def _payload_v36():
    payload = _payload_v35()
    identity = payload["artifact_identity"]
    generation = "dq-mtpa-236"
    pole_pairs = 4
    i_d, i_q = -50.0, 100.0
    l_d, l_q, flux_pm = 0.001, 0.0015, 0.075
    flux = [flux_pm + l_d * i_d, l_q * i_q]
    torque = 1.5 * pole_pairs * (flux[0] * i_q - flux[1] * i_d)
    identity["dq_flux_inductance_torque_mtpa_current_angle_speed_convention_owner_result_identity"] = {
        "dq_generation": generation,
        **{key: generation for key in (
            "park_generation", "flux_generation", "inductance_generation", "torque_generation",
            "mtpa_generation", "current_generation", "angle_generation", "speed_generation",
            "owner_generation", "result_generation")},
        "park_convention": "power_invariant_q_leads_d", "result_park_convention": "power_invariant_q_leads_d",
        "pole_pairs": pole_pairs, "result_pole_pairs": pole_pairs,
        "current_dq_a": [i_d, i_q], "result_current_dq_a": [i_d, i_q],
        "current_magnitude_a": math.hypot(i_d, i_q), "result_current_magnitude_a": math.hypot(i_d, i_q),
        "current_angle_rad": math.atan2(i_q, i_d), "result_current_angle_rad": math.atan2(i_q, i_d),
        "pm_flux_linkage_wb_turn": flux_pm, "result_pm_flux_linkage_wb_turn": flux_pm,
        "flux_linkage_dq_wb_turn": flux, "result_flux_linkage_dq_wb_turn": flux,
        "differential_inductance_dq_h": [l_d, l_q], "result_differential_inductance_dq_h": [l_d, l_q],
        "torque_nm": torque, "result_torque_nm": torque,
        "mechanical_speed_rad_s": 100.0, "result_mechanical_speed_rad_s": 100.0,
        "electrical_speed_rad_s": 400.0, "result_electrical_speed_rad_s": 400.0,
        "dq_owner": "motor/dq-236", "accepted_dq_owner": "motor/dq-236",
        "dq_result_sha256": "1" * 64, "accepted_dq_result_sha256": "1" * 64,
    }
    generation = "iron-loss-energy-236"
    frequency, flux_peak, temperature_factor = 100.0, 1.2, 1.1
    coefficients = [2.0, 0.1, 0.05]
    components = [
        coefficients[0] * frequency * flux_peak**2 * temperature_factor,
        coefficients[1] * frequency**2 * flux_peak**2,
        coefficients[2] * frequency**1.5 * flux_peak**1.5,
    ]
    regions = [["stator", 0.0007], ["rotor", 0.0003]]
    total_power = sum(components) * sum(row[1] for row in regions)
    identity["iron_loss_component_frequency_flux_region_thermal_energy_balance_owner_result_identity"] = {
        "iron_loss_generation": generation,
        **{key: generation for key in (
            "component_generation", "frequency_generation", "flux_generation", "region_generation",
            "thermal_generation", "power_generation", "energy_generation", "owner_generation",
            "result_generation")},
        "frequency_hz": frequency, "result_frequency_hz": frequency,
        "flux_peak_t": flux_peak, "result_flux_peak_t": flux_peak,
        "loss_coefficients": coefficients, "result_loss_coefficients": coefficients,
        "loss_components_w_m3": components, "result_loss_components_w_m3": components,
        "regional_volumes_m3": regions, "result_regional_volumes_m3": regions,
        "temperature_c": 80.0, "result_temperature_c": 80.0,
        "temperature_factor": temperature_factor, "result_temperature_factor": temperature_factor,
        "total_iron_loss_w": total_power, "result_total_iron_loss_w": total_power,
        "integration_duration_s": 0.2, "result_integration_duration_s": 0.2,
        "loss_energy_j": total_power * 0.2, "result_loss_energy_j": total_power * 0.2,
        "iron_loss_owner": "motor/iron-loss-236", "accepted_iron_loss_owner": "motor/iron-loss-236",
        "iron_loss_result_sha256": "2" * 64, "accepted_iron_loss_result_sha256": "2" * 64,
    }
    return payload


def test_v36_public_positive_dq_mtpa_and_iron_loss_energy_closure():
    assert pwm_controlled_motor_loss_gate(_payload_v36())["status"] == "ok"


def test_v36_public_dq_flux_inductance_torque_mtpa_current_angle_speed_convention_mismatch():
    payload = _payload_v36()
    row = payload["artifact_identity"]["dq_flux_inductance_torque_mtpa_current_angle_speed_convention_owner_result_identity"]
    row.update({"park_generation": "dq-mtpa-235", "speed_generation": "dq-mtpa-234",
                "result_generation": "dq-mtpa-233", "result_park_convention": "amplitude_invariant_d_leads_q",
                "result_pole_pairs": 2, "result_current_dq_a": [100.0, -50.0],
                "result_current_magnitude_a": -1.0, "result_current_angle_rad": -2.0,
                "result_pm_flux_linkage_wb_turn": -0.075, "result_flux_linkage_dq_wb_turn": [0.15, 0.025],
                "result_differential_inductance_dq_h": [-0.001, 0.0], "result_torque_nm": -60.0,
                "result_mechanical_speed_rad_s": -100.0, "result_electrical_speed_rad_s": 100.0,
                "accepted_dq_owner": "stale/dq", "accepted_dq_result_sha256": "a" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["dq_map_closes_park_flux_inductance_torque_mtpa_current_angle_speed_owner_and_result"]


def test_v36_public_iron_loss_hysteresis_eddy_excess_frequency_flux_energy_balance_mismatch():
    payload = _payload_v36()
    row = payload["artifact_identity"]["iron_loss_component_frequency_flux_region_thermal_energy_balance_owner_result_identity"]
    row.update({"component_generation": "iron-loss-energy-235", "thermal_generation": "iron-loss-energy-234",
                "result_generation": "iron-loss-energy-233", "result_frequency_hz": -100.0,
                "result_flux_peak_t": -1.2, "result_loss_coefficients": [2.0, -0.1, 0.0],
                "result_loss_components_w_m3": [1.0, -2.0, 3.0],
                "result_regional_volumes_m3": [["stator", -0.0007], ["old", 0.0]],
                "result_temperature_c": 20.0, "result_temperature_factor": -1.0,
                "result_total_iron_loss_w": -5.0, "result_integration_duration_s": -0.2,
                "result_loss_energy_j": 99.0, "accepted_iron_loss_owner": "stale/loss",
                "accepted_iron_loss_result_sha256": "b" * 64})
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["iron_loss_closes_components_frequency_flux_regions_thermal_power_energy_owner_and_result"]


def test_v36_rejects_self_consistent_dq_torque_error():
    payload = _payload_v36()
    row = payload["artifact_identity"]["dq_flux_inductance_torque_mtpa_current_angle_speed_convention_owner_result_identity"]
    row["torque_nm"] = row["result_torque_nm"] = -60.0
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v36_rejects_self_consistent_iron_loss_energy_error():
    payload = _payload_v36()
    row = payload["artifact_identity"]["iron_loss_component_frequency_flux_region_thermal_energy_balance_owner_result_identity"]
    row["loss_energy_j"] = row["result_loss_energy_j"] = 99.0
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"
