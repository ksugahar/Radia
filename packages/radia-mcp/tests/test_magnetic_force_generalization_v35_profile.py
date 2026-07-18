from __future__ import annotations

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from test_magnetic_force_generalization_v34_profile import _summary_v34


_PROMOTED_CASE_IDS = (
    "v35_public_magnetic_gear_polepair_harmonic_torque_phase_power_balance_mismatch",
    "v35_public_demag_bem_surface_charge_neutrality_normal_farfield_energy_mesh_mismatch",
)


def _summary_v35():
    summary = _summary_v34()
    identity = summary["artifact_identity"]
    generation = "magnetic-gear-401"
    identity[
        "magnetic_gear_pole_harmonic_torque_phase_power_frame_mesh_owner_result_identity"
    ] = {
        "magnetic_gear_generation": generation,
        **{
            key: generation
            for key in (
                "pole_generation", "harmonic_generation", "torque_generation",
                "phase_generation", "power_generation", "frame_generation",
                "mesh_generation", "owner_generation", "result_generation",
            )
        },
        "high_speed_pole_pairs": 4, "result_high_speed_pole_pairs": 4,
        "low_speed_pole_pairs": 22, "result_low_speed_pole_pairs": 22,
        "modulator_pole_count": 26, "result_modulator_pole_count": 26,
        "transmitted_harmonic_order": 22, "result_transmitted_harmonic_order": 22,
        "high_speed_torque_nm": -10.0, "result_high_speed_torque_nm": -10.0,
        "low_speed_torque_nm": 55.0, "result_low_speed_torque_nm": 55.0,
        "high_speed_angular_velocity_rad_s": 110.0,
        "result_high_speed_angular_velocity_rad_s": 110.0,
        "low_speed_angular_velocity_rad_s": 20.0,
        "result_low_speed_angular_velocity_rad_s": 20.0,
        "high_speed_harmonic_phase_rad": 0.5,
        "result_high_speed_harmonic_phase_rad": 0.5,
        "low_speed_harmonic_phase_rad": -0.1,
        "result_low_speed_harmonic_phase_rad": -0.1,
        "modulator_phase_rad": 0.2, "result_modulator_phase_rad": 0.2,
        "transmitted_phase_rad": 0.8, "result_transmitted_phase_rad": 0.8,
        "coordinate_frame": "global_xyz_right_handed",
        "result_coordinate_frame": "global_xyz_right_handed",
        "gear_mesh_sha256": "1" * 64, "result_gear_mesh_sha256": "1" * 64,
        "gear_result_owner": "magnetic-gear/case-401",
        "accepted_gear_result_owner": "magnetic-gear/case-401",
        "gear_result_sha256": "2" * 64,
        "accepted_gear_result_sha256": "2" * 64,
    }
    generation = "demag-bem-401"
    identity[
        "demag_bem_surface_charge_normal_jump_farfield_energy_mesh_owner_solution_identity"
    ] = {
        "demag_bem_generation": generation,
        **{
            key: generation
            for key in (
                "charge_generation", "normal_generation", "jump_generation",
                "farfield_generation", "energy_generation", "mesh_generation",
                "owner_generation", "solution_generation",
            )
        },
        "panel_areas_m2": [1.0, 1.0, 2.0, 2.0],
        "result_panel_areas_m2": [1.0, 1.0, 2.0, 2.0],
        "surface_charge_density_a_m": [2.0, -2.0, 1.0, -1.0],
        "result_surface_charge_density_a_m": [2.0, -2.0, 1.0, -1.0],
        "surface_charge_integral_a_m": 0.0, "result_surface_charge_integral_a_m": 0.0,
        "outward_normals": [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]],
        "result_outward_normals": [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0]],
        "outward_orientation_verified": True, "result_outward_orientation_verified": True,
        "normal_field_jump_a_m": [2.0, -2.0, 1.0, -1.0],
        "result_normal_field_jump_a_m": [2.0, -2.0, 1.0, -1.0],
        "farfield_radius_m": [2.0, 4.0, 8.0], "result_farfield_radius_m": [2.0, 4.0, 8.0],
        "farfield_potential_a": [0.25, 0.0625, 0.015625],
        "result_farfield_potential_a": [0.25, 0.0625, 0.015625],
        "farfield_field_a_m": [0.125, 0.015625, 0.001953125],
        "result_farfield_field_a_m": [0.125, 0.015625, 0.001953125],
        "magnetic_energy_j": 0.75, "result_magnetic_energy_j": 0.75,
        "boundary_mesh_sha256": "3" * 64, "result_boundary_mesh_sha256": "3" * 64,
        "demag_solution_owner": "demag-bem/case-401",
        "accepted_demag_solution_owner": "demag-bem/case-401",
        "demag_solution_sha256": "4" * 64,
        "accepted_demag_solution_sha256": "4" * 64,
    }
    return summary


def test_v35_public_positive_magnetic_gear_and_demag_bem_closure():
    assert magnetic_force_method_profile_gate(_summary_v35())["status"] == "ok"


def test_v35_public_magnetic_gear_polepair_harmonic_torque_phase_power_balance_mismatch():
    summary = _summary_v35()
    record = summary["artifact_identity"][
        "magnetic_gear_pole_harmonic_torque_phase_power_frame_mesh_owner_result_identity"
    ]
    record.update(
        {
            "pole_generation": "magnetic-gear-400",
            "power_generation": "magnetic-gear-399",
            "result_generation": "magnetic-gear-398",
            "result_high_speed_pole_pairs": 5,
            "result_low_speed_pole_pairs": 20,
            "result_modulator_pole_count": 24,
            "result_transmitted_harmonic_order": 21,
            "result_high_speed_torque_nm": 10.0,
            "result_low_speed_torque_nm": 20.0,
            "result_high_speed_angular_velocity_rad_s": 100.0,
            "result_low_speed_angular_velocity_rad_s": 40.0,
            "result_high_speed_harmonic_phase_rad": -0.5,
            "result_low_speed_harmonic_phase_rad": 0.7,
            "result_modulator_phase_rad": -0.2,
            "result_transmitted_phase_rad": 2.5,
            "result_coordinate_frame": "rotor_left_handed",
            "result_gear_mesh_sha256": "9" * 64,
            "accepted_gear_result_owner": "magnetic-gear/old",
            "accepted_gear_result_sha256": "a" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "magnetic_gear_uses_current_poles_harmonic_torque_phase_power_frame_mesh_owner_and_result"
    ]


def test_v35_public_demag_bem_surface_charge_neutrality_normal_farfield_energy_mesh_mismatch():
    summary = _summary_v35()
    record = summary["artifact_identity"][
        "demag_bem_surface_charge_normal_jump_farfield_energy_mesh_owner_solution_identity"
    ]
    record.update(
        {
            "charge_generation": "demag-bem-400",
            "farfield_generation": "demag-bem-399",
            "solution_generation": "demag-bem-398",
            "result_panel_areas_m2": [1.0, -1.0],
            "result_surface_charge_density_a_m": [2.0, 2.0],
            "result_surface_charge_integral_a_m": 4.0,
            "result_outward_normals": [[0.0, 0.0, 0.0]],
            "result_outward_orientation_verified": False,
            "result_normal_field_jump_a_m": [-2.0, 2.0, -1.0, 1.0],
            "result_farfield_radius_m": [8.0, 4.0, 2.0],
            "result_farfield_potential_a": [0.25, 0.25, 0.25],
            "result_farfield_field_a_m": [0.125, 0.125, 0.125],
            "result_magnetic_energy_j": -0.75,
            "result_boundary_mesh_sha256": "b" * 64,
            "accepted_demag_solution_owner": "demag-bem/old",
            "accepted_demag_solution_sha256": "c" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "demag_bem_uses_neutral_surface_charge_outward_normals_jump_farfield_energy_mesh_owner_and_solution"
    ]


def test_v35_public_rejects_self_consistent_non_neutral_demag_charge():
    summary = _summary_v35()
    record = summary["artifact_identity"][
        "demag_bem_surface_charge_normal_jump_farfield_energy_mesh_owner_solution_identity"
    ]
    charges = [2.0, 2.0, 1.0, 1.0]
    record["surface_charge_density_a_m"] = charges
    record["result_surface_charge_density_a_m"] = charges
    record["normal_field_jump_a_m"] = charges
    record["result_normal_field_jump_a_m"] = charges
    record["surface_charge_integral_a_m"] = 8.0
    record["result_surface_charge_integral_a_m"] = 8.0
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v35_public_rejects_self_consistent_magnetic_gear_power_imbalance():
    summary = _summary_v35()
    record = summary["artifact_identity"][
        "magnetic_gear_pole_harmonic_torque_phase_power_frame_mesh_owner_result_identity"
    ]
    record["low_speed_torque_nm"] = 50.0
    record["result_low_speed_torque_nm"] = 50.0
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"
