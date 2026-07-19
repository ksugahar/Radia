from __future__ import annotations

from radia_mcp.radia_ngsolve.comsol_v46_identity import validate_public_identity


PROMOTED_CASE_IDS = {
    "v46_public_time_adaptive_partial_solution_nan_inf_restart_window_mismatch",
    "v46_public_unit_scale_coordinate_frame_complex_field_vector_mismatch",
}


def _time_record(**updates: object) -> dict[str, object]:
    generation = "td-partial-restart-v46-901"
    record: dict[str, object] = {
        "generation": generation,
        "adaptive_step_generation": generation,
        "restart_window_generation": generation,
        "finite_value_generation": generation,
        "result_generation": generation,
        "restart_window_s": [0.4, 0.8],
        "result_restart_window_s": [0.4, 0.8],
        "accepted_step_count": 8,
        "result_accepted_step_count": 8,
        "adaptive_step_policy": "bdf_adaptive_accepted_only",
        "result_adaptive_step_policy": "bdf_adaptive_accepted_only",
        "finite_value_status": "finite",
        "result_finite_value_status": "finite",
        "restart_state": "partial_restartable",
        "result_restart_state": "partial_restartable",
        "owner": "model/td-partial-v46-901",
        "accepted_owner": "model/td-partial-v46-901",
        "result_sha256": "5" * 64,
        "accepted_result_sha256": "5" * 64,
    }
    record.update(updates)
    return record


def _field_record(**updates: object) -> dict[str, object]:
    generation = "field-frame-v46-901"
    record: dict[str, object] = {
        "generation": generation,
        "unit_scale_generation": generation,
        "coordinate_frame_generation": generation,
        "complex_vector_generation": generation,
        "result_generation": generation,
        "unit_name": "V_per_m",
        "result_unit_name": "V_per_m",
        "unit_scale_to_si": 1.0,
        "result_unit_scale_to_si": 1.0,
        "coordinate_frame": "global_cartesian",
        "result_coordinate_frame": "global_cartesian",
        "complex_vector_convention": "phasor_real_imag",
        "result_complex_vector_convention": "phasor_real_imag",
        "owner": "model/field-frame-v46-901",
        "accepted_owner": "model/field-frame-v46-901",
        "result_sha256": "6" * 64,
        "accepted_result_sha256": "6" * 64,
    }
    record.update(updates)
    return record


def test_v46_positive_replays_are_accepted() -> None:
    result = validate_public_identity(
        {
            "time_adaptive_partial_solution_nan_inf_restart_window_identity": _time_record(),
            "unit_scale_coordinate_frame_complex_field_vector_identity": _field_record(),
        }
    )
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_v46_partial_restart_mutation_is_rejected() -> None:
    result = validate_public_identity(
        {
            "time_adaptive_partial_solution_nan_inf_restart_window_identity": _time_record(
                result_restart_window_s=[0.0, 0.4],
                result_accepted_step_count=0,
                result_finite_value_status="contains_nan",
                accepted_owner="model/old",
                accepted_result_sha256="a" * 64,
            )
        }
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["v46_time_partial_restart_identity"] is False


def test_v46_unit_frame_mutation_is_rejected() -> None:
    result = validate_public_identity(
        {
            "unit_scale_coordinate_frame_complex_field_vector_identity": _field_record(
                result_unit_scale_to_si=1000.0,
                result_coordinate_frame="global_cylindrical",
                result_complex_vector_convention="real_imag_components",
                accepted_owner="model/old",
                accepted_result_sha256="b" * 64,
            )
        }
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["v46_field_unit_frame_identity"] is False
