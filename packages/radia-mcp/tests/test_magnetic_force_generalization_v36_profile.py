from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from test_magnetic_force_generalization_v35_profile import _summary_v35


_PROMOTED_CASE_IDS = (
    "v36_public_magnetic_bearing_stiffness_force_displacement_bias_current_linearization_mismatch",
    "v36_public_pm_demag_recoil_knee_loadline_temperature_irreversible_loss_mismatch",
)


def _summary_v36():
    summary = _summary_v35()
    identity = summary["artifact_identity"]
    generation = "bearing-linearization-402"
    identity[
        "magnetic_bearing_bias_displacement_force_stiffness_crosscoupling_frame_owner_result_identity"
    ] = {
        "bearing_generation": generation,
        **{
            key: generation
            for key in (
                "bias_generation", "displacement_generation", "force_generation",
                "stiffness_generation", "crosscoupling_generation", "frame_generation",
                "owner_generation", "result_generation",
            )
        },
        "bias_current_a": 5.0, "result_bias_current_a": 5.0,
        "displacement_samples_m": [-0.001, 0.0, 0.001],
        "result_displacement_samples_m": [-0.001, 0.0, 0.001],
        "force_x_samples_n": [10.0, 0.0, -10.0],
        "result_force_x_samples_n": [10.0, 0.0, -10.0],
        "stiffness_matrix_n_m": [[10000.0, 100.0], [100.0, 9000.0]],
        "result_stiffness_matrix_n_m": [[10000.0, 100.0], [100.0, 9000.0]],
        "coordinate_frame": "global_xyz_right_handed",
        "result_coordinate_frame": "global_xyz_right_handed",
        "bearing_owner": "bearing/case-402",
        "accepted_bearing_owner": "bearing/case-402",
        "bearing_result_sha256": "1" * 64,
        "accepted_bearing_result_sha256": "1" * 64,
    }
    generation = "pm-demag-402"
    reference_temperature = 20.0
    operating_temperature = 100.0
    remanence_reference = 1.2
    coefficient = -0.001
    remanence_temperature = remanence_reference * (
        1.0 + coefficient * (operating_temperature - reference_temperature)
    )
    recoil_mu = 1.05
    h_points = [-900000.0, -800000.0, -600000.0]
    b_points = [
        remanence_temperature + 4.0e-7 * math.pi * recoil_mu * field
        for field in h_points
    ]
    identity[
        "pm_demag_recoil_knee_loadline_temperature_irreversible_orientation_mesh_owner_result_identity"
    ] = {
        "demag_generation": generation,
        **{
            key: generation
            for key in (
                "recoil_generation", "knee_generation", "loadline_generation",
                "temperature_generation", "irreversible_generation", "orientation_generation",
                "mesh_generation", "owner_generation", "result_generation",
            )
        },
        "reference_temperature_c": reference_temperature,
        "result_reference_temperature_c": reference_temperature,
        "operating_temperature_c": operating_temperature,
        "result_operating_temperature_c": operating_temperature,
        "remanence_reference_t": remanence_reference,
        "result_remanence_reference_t": remanence_reference,
        "remanence_temperature_coefficient_per_c": coefficient,
        "result_remanence_temperature_coefficient_per_c": coefficient,
        "temperature_adjusted_remanence_t": remanence_temperature,
        "result_temperature_adjusted_remanence_t": remanence_temperature,
        "recoil_relative_permeability": recoil_mu,
        "result_recoil_relative_permeability": recoil_mu,
        "knee_field_a_m": -800000.0, "result_knee_field_a_m": -800000.0,
        "loadline_h_a_m": h_points, "result_loadline_h_a_m": h_points,
        "loadline_b_t": b_points, "result_loadline_b_t": b_points,
        "knee_crossed": True, "result_knee_crossed": True,
        "irreversible_flux_loss_fraction": 0.05,
        "result_irreversible_flux_loss_fraction": 0.05,
        "field_orientation": "magnetization_antiparallel_h",
        "result_field_orientation": "magnetization_antiparallel_h",
        "mesh_sha256": "2" * 64, "result_mesh_sha256": "2" * 64,
        "demag_owner": "pm/demag-402", "accepted_demag_owner": "pm/demag-402",
        "demag_result_sha256": "3" * 64,
        "accepted_demag_result_sha256": "3" * 64,
    }
    return summary


def test_v36_public_positive_bearing_and_pm_demag_closure():
    assert magnetic_force_method_profile_gate(_summary_v36())["status"] == "ok"


def test_v36_public_magnetic_bearing_stiffness_force_displacement_bias_current_linearization_mismatch():
    summary = _summary_v36()
    record = summary["artifact_identity"][
        "magnetic_bearing_bias_displacement_force_stiffness_crosscoupling_frame_owner_result_identity"
    ]
    record.update(
        {
            "bias_generation": "bearing-linearization-401",
            "stiffness_generation": "bearing-linearization-400",
            "result_generation": "bearing-linearization-399",
            "result_bias_current_a": -5.0,
            "result_displacement_samples_m": [0.0, 0.001],
            "result_force_x_samples_n": [10.0, 10.0],
            "result_stiffness_matrix_n_m": [[-10000.0, 5000.0], [-100.0, -9000.0]],
            "result_coordinate_frame": "rotor_left_handed",
            "accepted_bearing_owner": "bearing/old",
            "accepted_bearing_result_sha256": "a" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "magnetic_bearing_bias_sweep_uses_current_bias_symmetric_force_derivative_crosscoupling_frame_owner_and_result"
    ]


def test_v36_public_pm_demag_recoil_knee_loadline_temperature_irreversible_loss_mismatch():
    summary = _summary_v36()
    record = summary["artifact_identity"][
        "pm_demag_recoil_knee_loadline_temperature_irreversible_orientation_mesh_owner_result_identity"
    ]
    record.update(
        {
            "recoil_generation": "pm-demag-401",
            "temperature_generation": "pm-demag-400",
            "result_generation": "pm-demag-399",
            "result_temperature_adjusted_remanence_t": -1.0,
            "result_recoil_relative_permeability": -1.05,
            "result_knee_crossed": False,
            "result_irreversible_flux_loss_fraction": -0.2,
            "result_field_orientation": "parallel_h",
            "result_mesh_sha256": "b" * 64,
            "accepted_demag_owner": "pm/old",
            "accepted_demag_result_sha256": "c" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "pm_demag_uses_temperature_adjusted_recoil_knee_loadline_irreversible_loss_orientation_mesh_owner_and_result"
    ]


def test_v36_public_rejects_self_consistent_bearing_derivative_mismatch():
    summary = _summary_v36()
    record = summary["artifact_identity"][
        "magnetic_bearing_bias_displacement_force_stiffness_crosscoupling_frame_owner_result_identity"
    ]
    stiffness = [[8000.0, 100.0], [100.0, 9000.0]]
    record["stiffness_matrix_n_m"] = stiffness
    record["result_stiffness_matrix_n_m"] = stiffness
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v36_public_rejects_self_consistent_demag_loadline_mismatch():
    summary = _summary_v36()
    record = summary["artifact_identity"][
        "pm_demag_recoil_knee_loadline_temperature_irreversible_orientation_mesh_owner_result_identity"
    ]
    loadline = [0.5, 0.5, 0.5]
    record["loadline_b_t"] = loadline
    record["result_loadline_b_t"] = loadline
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"
