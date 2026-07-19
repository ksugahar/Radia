from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


_FORCE = "v46_public_magnetic_force_partial_solve_unit_scale_coordinate_frame_nan_mismatch"
_DEMAG = "v46_public_demagnetization_curve_branch_restart_temperature_window_mismatch"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _finite_sequence(value: object, *, minimum: int = 1) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and len(value) >= minimum
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def _linked(row: Mapping[str, object], fields: tuple[str, ...]) -> bool:
    generation = str(row.get("generation", "")).strip()
    present = [field for field in fields if field in row]
    return bool(generation) and all(row.get(field) == generation for field in present)


def _force_ok(row: Mapping[str, object]) -> bool:
    return (
        _linked(row, ("solve_generation", "unit_scale_generation", "frame_generation", "force_generation", "finite_generation", "result_generation"))
        and row.get("solve_completion") == row.get("result_solve_completion") == "complete"
        and row.get("unit_scale_to_si") == row.get("result_unit_scale_to_si") == 1.0
        and row.get("coordinate_frame") == row.get("result_coordinate_frame") == "global_cartesian"
        and _finite_sequence(row.get("force_n"))
        and row.get("force_n") == row.get("result_force_n")
        and row.get("nonfinite_value_count") == row.get("result_nonfinite_value_count") == 0
        and row.get("finite_values") == row.get("result_finite_values") is True
        and str(row.get("mesh_owner", "")).startswith("mesh:")
        and row.get("result_mesh_owner") == row.get("mesh_owner")
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def _demag_ok(row: Mapping[str, object]) -> bool:
    temperature = row.get("temperature_window_k")
    return (
        _linked(row, ("branch_generation", "restart_generation", "temperature_generation", "path_generation", "completion_generation", "result_generation"))
        and row.get("branch_mode") == row.get("result_branch_mode") == "continuous"
        and row.get("restart_generation") == row.get("result_restart_generation")
        and _finite_sequence(temperature, minimum=2)
        and temperature == row.get("result_temperature_window_k")
        and float(temperature[0]) < float(temperature[-1])
        and _finite_sequence(row.get("field_path_a_per_m"), minimum=3)
        and row.get("field_path_a_per_m") == row.get("result_field_path_a_per_m")
        and row.get("partial_path_status") == row.get("result_partial_path_status") == "complete"
        and str(row.get("material_owner", "")).startswith("material:")
        and row.get("result_material_owner") == row.get("material_owner")
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def validate_public_v46_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    if _FORCE in identity:
        row = identity[_FORCE]
        checks["elf_v46_partial_force_identity"] = isinstance(row, Mapping) and _force_ok(row)
    if _DEMAG in identity:
        row = identity[_DEMAG]
        checks["elf_v46_demagnetization_restart_identity"] = isinstance(row, Mapping) and _demag_ok(row)
    return checks
