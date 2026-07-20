from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _finite(value: object) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and bool(value) and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)


def _closed(row: Mapping[str, object], fields: tuple[str, ...]) -> bool:
    generation = str(row.get("generation", "")).strip()
    return bool(generation) and all(row.get(field) == generation for field in fields)


def validate_public_v45_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    dynamic = identity.get("v45_public_magnetic_bearing_dynamic_stiffness_phase_damping_force_power_stability_owner_mismatch")
    if isinstance(dynamic, Mapping):
        checks["elf_v45_dynamic_generation"] = _closed(dynamic, ("stiffness_generation", "phase_generation", "damping_generation", "force_generation", "power_generation", "stability_generation", "mesh_generation", "result_generation"))
        checks["elf_v45_dynamic_values"] = all(_finite(dynamic.get(key)) and dynamic.get(f"result_{key}") == dynamic.get(key) for key in ("frequency_hz", "dynamic_stiffness_n_per_m", "phase_deg", "damping_n_s_per_m", "force_n", "power_w")) and all(float(item) >= 0.0 for item in dynamic.get("damping_n_s_per_m", [])) and all(float(item) >= 0.0 for item in dynamic.get("power_w", [])) and dynamic.get("stability_sign") == dynamic.get("result_stability_sign") == "stable"
        checks["elf_v45_dynamic_owner"] = str(dynamic.get("mesh_owner", "")).startswith("mesh:") and dynamic.get("result_mesh_owner") == dynamic.get("mesh_owner") and _digest(dynamic.get("result_sha256")) and dynamic.get("accepted_result_sha256") == dynamic.get("result_sha256")
    demag = identity.get("v45_public_demagnetization_minor_loop_field_path_remanence_coercivity_loss_temperature_mismatch")
    if isinstance(demag, Mapping):
        checks["elf_v45_demag_generation"] = _closed(demag, ("field_path_generation", "remanence_generation", "coercivity_generation", "loss_generation", "temperature_generation", "material_generation", "result_generation"))
        checks["elf_v45_demag_values"] = demag.get("field_path_a_per_m") == demag.get("result_field_path_a_per_m") and float(demag.get("remanence_a_per_m")) >= 0.0 and demag.get("result_remanence_a_per_m") == demag.get("remanence_a_per_m") and float(demag.get("coercivity_a_per_m")) >= 0.0 and demag.get("result_coercivity_a_per_m") == demag.get("coercivity_a_per_m") and float(demag.get("loss_energy_j")) >= 0.0 and demag.get("result_loss_energy_j") == demag.get("loss_energy_j") and float(demag.get("temperature_k")) > 0.0 and demag.get("result_temperature_k") == demag.get("temperature_k")
        checks["elf_v45_demag_owner"] = str(demag.get("material_owner", "")).startswith("material:") and demag.get("result_material_owner") == demag.get("material_owner") and _digest(demag.get("result_sha256")) and demag.get("accepted_result_sha256") == demag.get("result_sha256")
    return checks
