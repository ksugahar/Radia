"""Payload builders for ``test_magnetic_force_demagnetization_generalization.py``.

Not collected by pytest (leading underscore).  Both builders feed
``magnetic_force_v44_identity.validate_public_identity`` as independent
roots (magnetic-bearing dynamics / demagnetization minor loop, and
partial-solve force / demagnetization branch restart).
"""

from __future__ import annotations


def _identity_v45():
    generation = "test-845"
    return {
        "v45_public_magnetic_bearing_dynamic_stiffness_phase_damping_force_power_stability_owner_mismatch": {
            "generation": generation, **{key: generation for key in ("stiffness_generation", "phase_generation", "damping_generation", "force_generation", "power_generation", "stability_generation", "mesh_generation", "result_generation")},
            "frequency_hz": [100.0, 200.0], "result_frequency_hz": [100.0, 200.0], "dynamic_stiffness_n_per_m": [100.0, 110.0], "result_dynamic_stiffness_n_per_m": [100.0, 110.0], "phase_deg": [0.0, 10.0], "result_phase_deg": [0.0, 10.0], "damping_n_s_per_m": [2.0, 2.5], "result_damping_n_s_per_m": [2.0, 2.5], "force_n": [10.0, 11.0], "result_force_n": [10.0, 11.0], "power_w": [1.0, 1.2], "result_power_w": [1.0, 1.2], "stability_sign": "stable", "result_stability_sign": "stable", "mesh_owner": "mesh:test", "result_mesh_owner": "mesh:test", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
        },
        "v45_public_demagnetization_minor_loop_field_path_remanence_coercivity_loss_temperature_mismatch": {
            "generation": generation, **{key: generation for key in ("field_path_generation", "remanence_generation", "coercivity_generation", "loss_generation", "temperature_generation", "material_generation", "result_generation")}, "field_path_a_per_m": [-1.0, 0.5, -0.5, 1.0], "result_field_path_a_per_m": [-1.0, 0.5, -0.5, 1.0], "remanence_a_per_m": 0.8, "result_remanence_a_per_m": 0.8, "coercivity_a_per_m": 0.4, "result_coercivity_a_per_m": 0.4, "loss_energy_j": 0.16, "result_loss_energy_j": 0.16, "temperature_k": 293.15, "result_temperature_k": 293.15, "material_owner": "material:test", "result_material_owner": "material:test", "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64,
        },
    }


def _identity_v46():
    generation = "test-846"
    return {
        "v46_public_magnetic_force_partial_solve_unit_scale_coordinate_frame_nan_mismatch": {
            "generation": generation,
            **{key: generation for key in ("solve_generation", "unit_scale_generation", "frame_generation", "force_generation", "finite_generation", "result_generation")},
            "solve_completion": "complete", "result_solve_completion": "complete", "unit_scale_to_si": 1.0, "result_unit_scale_to_si": 1.0,
            "coordinate_frame": "global_cartesian", "result_coordinate_frame": "global_cartesian", "force_n": [1.0, 2.0, 3.0], "result_force_n": [1.0, 2.0, 3.0],
            "nonfinite_value_count": 0, "result_nonfinite_value_count": 0, "finite_values": True, "result_finite_values": True,
            "mesh_owner": "mesh:test", "result_mesh_owner": "mesh:test", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
        },
        "v46_public_demagnetization_curve_branch_restart_temperature_window_mismatch": {
            "generation": generation,
            **{key: generation for key in ("branch_generation", "restart_generation", "temperature_generation", "path_generation", "completion_generation", "result_generation")},
            "branch_mode": "continuous", "result_branch_mode": "continuous", "result_restart_generation": generation,
            "temperature_window_k": [293.15, 353.15], "result_temperature_window_k": [293.15, 353.15], "field_path_a_per_m": [-1.0, 0.0, 1.0], "result_field_path_a_per_m": [-1.0, 0.0, 1.0],
            "partial_path_status": "complete", "result_partial_path_status": "complete", "material_owner": "material:test", "result_material_owner": "material:test", "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64,
        },
    }
