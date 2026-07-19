from __future__ import annotations

from radia_mcp.radia_ngsolve.magnetic_force_v44_identity import validate_public_identity


def _identity():
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


def test_v46_public_identity_accepts_closed_artifacts():
    checks = validate_public_identity(_identity())
    assert checks and all(checks.values())


def test_v46_public_identity_rejects_partial_nan_and_restart_mutations():
    identity = _identity()
    identity["v46_public_magnetic_force_partial_solve_unit_scale_coordinate_frame_nan_mismatch"]["result_finite_values"] = False
    identity["v46_public_demagnetization_curve_branch_restart_temperature_window_mismatch"]["result_partial_path_status"] = "partial"
    checks = validate_public_identity(identity)
    assert checks and not all(checks.values())
