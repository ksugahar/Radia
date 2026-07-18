from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.pwm_controlled_motor_loss_gate import (
    pwm_controlled_motor_loss_gate,
)
from test_motor_generalization_v33 import _payload_v33


_PROMOTED_CASE_IDS = (
    "v34_public_demagnetization_temperature_recoil_loadline_operating_point_margin_mismatch",
    "v34_public_eccentricity_unbalanced_magnetic_pull_harmonic_frame_force_torque_mismatch",
)


def _payload_v34():
    payload = _payload_v33()
    identity = payload["artifact_identity"]

    generation = "pm-demag-operating-point-211"
    identity[
        "pm_demagnetization_temperature_recoil_loadline_operating_point_knee_margin_angle_mesh_owner_result_identity"
    ] = {
        "demag_generation": generation,
        **{
            key: generation
            for key in (
                "temperature_generation",
                "recoil_generation",
                "loadline_generation",
                "operating_point_generation",
                "knee_generation",
                "margin_generation",
                "angle_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "reference_temperature_c": 20.0,
        "result_reference_temperature_c": 20.0,
        "operating_temperature_c": 120.0,
        "result_operating_temperature_c": 120.0,
        "remanence_reference_t": 1.2,
        "result_remanence_reference_t": 1.2,
        "remanence_temperature_coefficient_per_c": -0.001,
        "result_remanence_temperature_coefficient_per_c": -0.001,
        "remanence_operating_t": 1.08,
        "result_remanence_operating_t": 1.08,
        "coercivity_reference_a_m": 900000.0,
        "result_coercivity_reference_a_m": 900000.0,
        "coercivity_temperature_coefficient_per_c": -0.002,
        "result_coercivity_temperature_coefficient_per_c": -0.002,
        "coercivity_operating_a_m": 720000.0,
        "result_coercivity_operating_a_m": 720000.0,
        "recoil_permeability_relative": 1.05,
        "result_recoil_permeability_relative": 1.05,
        "loadline_slope_t_per_a_m": 1.0e-6,
        "result_loadline_slope_t_per_a_m": 1.0e-6,
        "operating_field_a_m": -600000.0,
        "result_operating_field_a_m": -600000.0,
        "operating_flux_density_t": 0.48,
        "result_operating_flux_density_t": 0.48,
        "knee_field_a_m": -700000.0,
        "result_knee_field_a_m": -700000.0,
        "irreversible_margin_a_m": 100000.0,
        "result_irreversible_margin_a_m": 100000.0,
        "rotor_angle_rad": 0.3,
        "result_rotor_angle_rad": 0.3,
        "demag_mesh_sha256": "1" * 64,
        "result_demag_mesh_sha256": "1" * 64,
        "demag_result_owner": "motor/pm-demag-211",
        "accepted_demag_result_owner": "motor/pm-demag-211",
        "demag_result_sha256": "2" * 64,
        "accepted_demag_result_sha256": "2" * 64,
    }

    generation = "eccentricity-ump-211"
    angle_grid = [index * math.pi / 8.0 for index in range(5)]
    identity[
        "eccentricity_static_dynamic_frame_radial_force_harmonic_ump_torque_pole_periodicity_angle_owner_result_identity"
    ] = {
        "eccentricity_generation": generation,
        **{
            key: generation
            for key in (
                "static_generation",
                "dynamic_generation",
                "frame_generation",
                "harmonic_generation",
                "force_generation",
                "torque_generation",
                "periodicity_generation",
                "angle_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "static_eccentricity_m": [1.0e-4, 0.0],
        "result_static_eccentricity_m": [1.0e-4, 0.0],
        "dynamic_eccentricity_amplitude_m": 5.0e-5,
        "result_dynamic_eccentricity_amplitude_m": 5.0e-5,
        "mechanical_frame": "stator_global_xy",
        "result_mechanical_frame": "stator_global_xy",
        "radial_force_harmonics_n": [[0, 0.0, 0.0], [1, 100.0, 0.0], [2, 0.0, 0.0]],
        "result_radial_force_harmonics_n": [[0, 0.0, 0.0], [1, 100.0, 0.0], [2, 0.0, 0.0]],
        "unbalanced_magnetic_pull_n": [100.0, 0.0],
        "result_unbalanced_magnetic_pull_n": [100.0, 0.0],
        "torque_nm": 20.0,
        "result_torque_nm": 20.0,
        "pole_pairs": 4,
        "result_pole_pairs": 4,
        "periodicity_angle_rad": math.pi / 2.0,
        "result_periodicity_angle_rad": math.pi / 2.0,
        "angle_grid_rad": angle_grid,
        "result_angle_grid_rad": angle_grid,
        "eccentricity_mesh_sha256": "3" * 64,
        "result_eccentricity_mesh_sha256": "3" * 64,
        "eccentricity_result_owner": "motor/eccentricity-211",
        "accepted_eccentricity_result_owner": "motor/eccentricity-211",
        "eccentricity_result_sha256": "4" * 64,
        "accepted_eccentricity_result_sha256": "4" * 64,
    }
    return payload


def test_v34_public_positive_demagnetization_and_eccentricity_closure():
    assert pwm_controlled_motor_loss_gate(_payload_v34())["status"] == "ok"


def test_v34_public_demagnetization_temperature_recoil_loadline_operating_point_margin_mismatch():
    payload = _payload_v34()
    record = payload["artifact_identity"][
        "pm_demagnetization_temperature_recoil_loadline_operating_point_knee_margin_angle_mesh_owner_result_identity"
    ]
    record.update(
        {
            "temperature_generation": "pm-demag-operating-point-210",
            "margin_generation": "pm-demag-operating-point-209",
            "result_generation": "pm-demag-operating-point-208",
            "result_operating_temperature_c": 20.0,
            "result_remanence_operating_t": 1.2,
            "result_coercivity_operating_a_m": 900000.0,
            "result_recoil_permeability_relative": -1.05,
            "result_loadline_slope_t_per_a_m": -1.0e-6,
            "result_operating_field_a_m": -800000.0,
            "result_operating_flux_density_t": -0.2,
            "result_knee_field_a_m": -600000.0,
            "result_irreversible_margin_a_m": -200000.0,
            "result_rotor_angle_rad": -0.3,
            "result_demag_mesh_sha256": "9" * 64,
            "accepted_demag_result_owner": "stale/demag",
            "accepted_demag_result_sha256": "a" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "pm_demagnetization_uses_current_temperature_recoil_loadline_knee_margin_angle_mesh_owner_and_result"
    ]


def test_v34_public_eccentricity_unbalanced_magnetic_pull_harmonic_frame_force_torque_mismatch():
    payload = _payload_v34()
    record = payload["artifact_identity"][
        "eccentricity_static_dynamic_frame_radial_force_harmonic_ump_torque_pole_periodicity_angle_owner_result_identity"
    ]
    record.update(
        {
            "static_generation": "eccentricity-ump-210",
            "harmonic_generation": "eccentricity-ump-209",
            "result_generation": "eccentricity-ump-208",
            "result_static_eccentricity_m": [-1.0e-4, 0.0],
            "result_dynamic_eccentricity_amplitude_m": -5.0e-5,
            "result_mechanical_frame": "rotor_local_yz",
            "result_radial_force_harmonics_n": [[0, 0.0, 0.0], [1, -100.0, 50.0]],
            "result_unbalanced_magnetic_pull_n": [0.0, -100.0],
            "result_torque_nm": -20.0,
            "result_pole_pairs": 2,
            "result_periodicity_angle_rad": math.pi,
            "result_angle_grid_rad": [0.0, 1.0, 0.5],
            "result_eccentricity_mesh_sha256": "b" * 64,
            "accepted_eccentricity_result_owner": "stale/eccentricity",
            "accepted_eccentricity_result_sha256": "c" * 64,
        }
    )
    result = pwm_controlled_motor_loss_gate(payload)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "eccentricity_ump_uses_current_static_dynamic_frame_harmonics_force_torque_periodicity_angles_owner_and_result"
    ]


def test_v34_public_rejects_self_consistent_but_temperature_wrong_remanence():
    payload = _payload_v34()
    record = payload["artifact_identity"][
        "pm_demagnetization_temperature_recoil_loadline_operating_point_knee_margin_angle_mesh_owner_result_identity"
    ]
    record["remanence_operating_t"] = record["result_remanence_operating_t"] = 1.1
    record["operating_flux_density_t"] = record["result_operating_flux_density_t"] = 0.5
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"


def test_v34_public_rejects_self_consistent_ump_opposite_static_eccentricity():
    payload = _payload_v34()
    record = payload["artifact_identity"][
        "eccentricity_static_dynamic_frame_radial_force_harmonic_ump_torque_pole_periodicity_angle_owner_result_identity"
    ]
    harmonics = [[0, 0.0, 0.0], [1, -100.0, 0.0], [2, 0.0, 0.0]]
    record["radial_force_harmonics_n"] = harmonics
    record["result_radial_force_harmonics_n"] = harmonics
    record["unbalanced_magnetic_pull_n"] = [-100.0, 0.0]
    record["result_unbalanced_magnetic_pull_n"] = [-100.0, 0.0]
    assert pwm_controlled_motor_loss_gate(payload)["status"] == "needs_attention"
