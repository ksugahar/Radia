from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from test_magnetic_force_generalization_v36_profile import _summary_v36


_PROMOTED_CASE_IDS = (
    "v37_public_maglev_dynamic_stiffness_frequency_damping_phase_bias_equilibrium_owner_mismatch",
    "v37_public_bem_demag_reciprocity_energy_field_magnetization_surface_owner_mismatch",
)


def _summary_v37():
    summary = _summary_v36()
    identity = summary["artifact_identity"]
    generation = "maglev-dynamic-246"
    frequency, stiffness, damping = 100.0, 10000.0, 10.0
    omega = 2.0 * math.pi * frequency
    displacement = [1.0e-5, 0.0]
    dynamic_stiffness = [stiffness, omega * damping]
    force = [dynamic_stiffness[0] * displacement[0], dynamic_stiffness[1] * displacement[0]]
    identity["maglev_bias_equilibrium_frequency_complex_stiffness_damping_force_displacement_phase_frame_owner_result_identity"] = {
        "maglev_generation": generation,
        **{key: generation for key in (
            "bias_generation", "equilibrium_generation", "frequency_generation",
            "stiffness_generation", "damping_generation", "force_generation",
            "displacement_generation", "phase_generation", "frame_generation",
            "owner_generation", "result_generation")},
        "bias_current_a": 5.0, "result_bias_current_a": 5.0,
        "equilibrium_gap_m": 0.005, "result_equilibrium_gap_m": 0.005,
        "equilibrium_force_n": 100.0, "result_equilibrium_force_n": 100.0,
        "supported_load_n": 100.0, "result_supported_load_n": 100.0,
        "excitation_frequency_hz": frequency, "result_excitation_frequency_hz": frequency,
        "complex_stiffness_n_m": dynamic_stiffness, "result_complex_stiffness_n_m": dynamic_stiffness,
        "viscous_damping_n_s_m": damping, "result_viscous_damping_n_s_m": damping,
        "displacement_phasor_m": displacement, "result_displacement_phasor_m": displacement,
        "force_phasor_n": force, "result_force_phasor_n": force,
        "force_displacement_phase_rad": math.atan2(force[1], force[0]),
        "result_force_displacement_phase_rad": math.atan2(force[1], force[0]),
        "coordinate_frame": "global_z_up_force_positive",
        "result_coordinate_frame": "global_z_up_force_positive",
        "maglev_owner": "maglev/dynamic-246", "accepted_maglev_owner": "maglev/dynamic-246",
        "maglev_result_sha256": "1" * 64, "accepted_maglev_result_sha256": "1" * 64,
    }
    generation = "bem-demag-246"
    identity["bem_demag_reciprocity_interaction_energy_field_magnetization_surface_volume_mesh_solution_result_identity"] = {
        "demag_generation": generation,
        **{key: generation for key in (
            "reciprocity_generation", "energy_generation", "field_generation",
            "magnetization_generation", "surface_generation", "volume_generation",
            "mesh_generation", "solution_generation", "result_generation")},
        "interaction_energy_12_j": -0.01, "result_interaction_energy_12_j": -0.01,
        "interaction_energy_21_j": -0.01, "result_interaction_energy_21_j": -0.01,
        "field_1_due_2_a_m": [-1000.0, 0.0, 0.0], "result_field_1_due_2_a_m": [-1000.0, 0.0, 0.0],
        "field_2_due_1_a_m": [0.0, -1000.0, 0.0], "result_field_2_due_1_a_m": [0.0, -1000.0, 0.0],
        "magnetization_1_a_m": [800000.0, 0.0, 0.0], "result_magnetization_1_a_m": [800000.0, 0.0, 0.0],
        "magnetization_2_a_m": [0.0, 800000.0, 0.0], "result_magnetization_2_a_m": [0.0, 800000.0, 0.0],
        "surface_orientation": "outward_right_handed", "result_surface_orientation": "outward_right_handed",
        "region_volumes_m3": [1.0e-5, 1.0e-5], "result_region_volumes_m3": [1.0e-5, 1.0e-5],
        "mesh_owner": "mesh/bem-demag-246", "accepted_mesh_owner": "mesh/bem-demag-246",
        "solution_owner": "solution/bem-demag-246", "accepted_solution_owner": "solution/bem-demag-246",
        "demag_result_sha256": "2" * 64, "accepted_demag_result_sha256": "2" * 64,
    }
    return summary


def test_v37_public_positive_maglev_dynamic_and_bem_demag_closure():
    assert magnetic_force_method_profile_gate(_summary_v37())["status"] == "ok"


def test_v37_public_maglev_dynamic_stiffness_frequency_damping_phase_bias_equilibrium_owner_mismatch():
    summary = _summary_v37()
    row = summary["artifact_identity"]["maglev_bias_equilibrium_frequency_complex_stiffness_damping_force_displacement_phase_frame_owner_result_identity"]
    row.update({"bias_generation": "maglev-dynamic-245", "phase_generation": "maglev-dynamic-244",
                "result_generation": "maglev-dynamic-243", "result_bias_current_a": -5.0,
                "result_equilibrium_gap_m": -0.005, "result_equilibrium_force_n": -100.0,
                "result_supported_load_n": 50.0, "result_excitation_frequency_hz": -100.0,
                "result_complex_stiffness_n_m": [-10000.0, -1.0], "result_viscous_damping_n_s_m": -10.0,
                "result_displacement_phasor_m": [0.0, 1.0e-5], "result_force_phasor_n": [-0.1, -0.1],
                "result_force_displacement_phase_rad": -2.0, "result_coordinate_frame": "rotor_down_left_handed",
                "accepted_maglev_owner": "stale/maglev", "accepted_maglev_result_sha256": "a" * 64})
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["maglev_dynamics_close_bias_equilibrium_frequency_stiffness_damping_phase_frame_owner_and_result"]


def test_v37_public_bem_demag_reciprocity_energy_field_magnetization_surface_owner_mismatch():
    summary = _summary_v37()
    row = summary["artifact_identity"]["bem_demag_reciprocity_interaction_energy_field_magnetization_surface_volume_mesh_solution_result_identity"]
    row.update({"reciprocity_generation": "bem-demag-245", "surface_generation": "bem-demag-244",
                "result_generation": "bem-demag-243", "result_interaction_energy_21_j": 0.02,
                "result_field_1_due_2_a_m": [1000.0, 0.0, 0.0], "result_field_2_due_1_a_m": [0.0, 1000.0, 0.0],
                "result_surface_orientation": "inward_left_handed", "result_region_volumes_m3": [1.0e-5, -1.0e-5],
                "accepted_mesh_owner": "stale/mesh", "accepted_solution_owner": "stale/solution",
                "accepted_demag_result_sha256": "b" * 64})
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["bem_demag_closes_reciprocal_energy_field_magnetization_surface_volume_mesh_solution_and_result"]


def test_v37_public_rejects_self_consistent_negative_dynamic_damping():
    summary = _summary_v37()
    row = summary["artifact_identity"]["maglev_bias_equilibrium_frequency_complex_stiffness_damping_force_displacement_phase_frame_owner_result_identity"]
    row["viscous_damping_n_s_m"] = row["result_viscous_damping_n_s_m"] = -10.0
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v37_public_rejects_self_consistent_bem_nonreciprocity():
    summary = _summary_v37()
    row = summary["artifact_identity"]["bem_demag_reciprocity_interaction_energy_field_magnetization_surface_volume_mesh_solution_result_identity"]
    row["interaction_energy_21_j"] = row["result_interaction_energy_21_j"] = 0.02
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"
