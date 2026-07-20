from __future__ import annotations

import math
from collections.abc import Mapping


def _sha(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _closed(row: Mapping[str, object], fields: tuple[str, ...]) -> bool:
    generation = str(row.get("generation", "")).strip()
    return bool(generation) and all(row.get(field) == generation for field in fields)


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    torque = identity.get("v46_public_transient_torque_sampling_periodic_window_partial_solution_mismatch")
    if isinstance(torque, Mapping):
        checks["motor_v46_torque_generation_closure"] = _closed(torque, ("torque_generation", "sampling_generation", "periodic_window_generation", "partial_solution_generation", "result_generation"))
        samples = torque.get("sample_times_s")
        checks["motor_v46_torque_sampling_window_state"] = (
            isinstance(samples, list)
            and samples == torque.get("result_sample_times_s")
            and len(samples) >= 3
            and all(math.isfinite(float(value)) for value in samples)
            and torque.get("periodic_window_deg") == torque.get("result_periodic_window_deg") == 360.0
            and torque.get("partial_transient_status") == torque.get("result_partial_transient_status") == "complete"
        )
        checks["motor_v46_torque_owner_digest"] = (
            str(torque.get("mesh_owner", "")).startswith("mesh:")
            and torque.get("result_mesh_owner") == torque.get("mesh_owner")
            and _sha(torque.get("result_sha256"))
            and torque.get("accepted_result_sha256") == torque.get("result_sha256")
        )
    thermal = identity.get("v46_public_thermal_coupling_unit_scale_temperature_coordinate_frame_mismatch")
    if isinstance(thermal, Mapping):
        checks["motor_v46_thermal_generation_closure"] = _closed(thermal, ("thermal_generation", "unit_generation", "temperature_generation", "coordinate_frame_generation", "result_generation"))
        checks["motor_v46_thermal_units_frame"] = (
            thermal.get("thermal_unit") == thermal.get("result_thermal_unit") == "kelvin"
            and thermal.get("temperature_frame") == thermal.get("result_temperature_frame") == "absolute_kelvin"
            and thermal.get("coordinate_frame") == thermal.get("result_coordinate_frame") == "global_cartesian"
            and math.isfinite(float(thermal.get("temperature_k")))
            and thermal.get("result_temperature_k") == thermal.get("temperature_k")
        )
        checks["motor_v46_thermal_owner_digest"] = (
            str(thermal.get("study_owner", "")).startswith("study:")
            and thermal.get("result_study_owner") == thermal.get("study_owner")
            and _sha(thermal.get("result_sha256"))
            and thermal.get("accepted_result_sha256") == thermal.get("result_sha256")
        )
    return checks


def validate_source_identity(replay_identities: object) -> dict[str, bool]:
    if not isinstance(replay_identities, list):
        return {}
    checks: dict[str, bool] = {}
    studies = [row.get("v46_source_tool_study_restart_parameter_scope_partial_result_tree_mismatch") for row in replay_identities if isinstance(row, Mapping) and isinstance(row.get("v46_source_tool_study_restart_parameter_scope_partial_result_tree_mismatch"), Mapping)]
    if studies:
        checks["motor_v46_source_study_generation_closure"] = len(studies) == len(replay_identities) and all(_closed(row, ("study_generation", "restart_generation", "parameter_scope_generation", "partial_result_tree_generation", "result_generation")) for row in studies)
        checks["motor_v46_source_study_replay_state"] = all(row.get("restart_point") == row.get("replayed_restart_point") == "step_12" and row.get("parameter_scope") == row.get("replayed_parameter_scope") == "global" and row.get("partial_result_tree") == row.get("replayed_partial_result_tree") == "complete" for row in studies)
        checks["motor_v46_source_study_owner_digest"] = all(str(row.get("study_owner", "")).startswith("study:") and row.get("replayed_study_owner") == row.get("study_owner") and _sha(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256") for row in studies)
    materials = [row.get("v46_source_tool_material_curve_interpolation_nan_temperature_step_mismatch") for row in replay_identities if isinstance(row, Mapping) and isinstance(row.get("v46_source_tool_material_curve_interpolation_nan_temperature_step_mismatch"), Mapping)]
    if materials:
        checks["motor_v46_source_material_generation_closure"] = len(materials) == len(replay_identities) and all(_closed(row, ("material_curve_generation", "interpolation_generation", "finite_sample_generation", "temperature_step_generation", "restart_generation", "result_generation")) for row in materials)
        checks["motor_v46_source_material_replay_state"] = all(row.get("interpolation") == row.get("replayed_interpolation") == "linear" and row.get("finite_sample_status") == row.get("replayed_finite_sample_status") == "none" and row.get("temperature_step_k") == row.get("replayed_temperature_step_k") == 5.0 and row.get("restart_state") == row.get("replayed_restart_state") == "cold_start" for row in materials)
        checks["motor_v46_source_material_owner_digest"] = all(str(row.get("study_owner", "")).startswith("study:") and row.get("replayed_study_owner") == row.get("study_owner") and _sha(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256") for row in materials)
    return checks
