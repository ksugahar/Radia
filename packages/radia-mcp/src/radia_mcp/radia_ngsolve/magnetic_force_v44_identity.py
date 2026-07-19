"""Neutral identity checks for dynamic magnetic-force result artifacts."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .elf_v45_identity import validate_public_v45_identity
from .elf_v46_identity import validate_public_v46_identity
from .magnetic_force_artifact_lineage_v47 import validate_public_v47_identity
from .bem_hysteresis_identity_v48 import validate_public_identity as validate_public_v48_identity
from .demag_virtual_work_identity_v49 import validate_public_identity as validate_public_v49_identity


_DYNAMIC = "magneticbearing_dynamicstiffness_phase_damping_force_power_stability_mesh_result_identity"
_DEMAG = "demag_minorloop_fieldpath_remanence_loss_energy_temperature_material_mesh_result_identity"


def _same(row: Mapping[str, object], left: str, right: str) -> bool:
    return row.get(left) == row.get(right)


def _finite_sequence(value: object, *, minimum: int = 1) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and len(value) >= minimum
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _dynamic_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("bearing_dynamic_generation", "")).strip()
    linked = (
        "dynamic_stiffness_generation", "phase_generation", "damping_generation",
        "force_generation", "power_generation", "stability_generation",
        "mesh_generation", "result_generation",
    )
    arrays = (
        ("frequency_hz", "result_frequency_hz"),
        ("dynamic_stiffness_n_per_m", "result_dynamic_stiffness_n_per_m"),
        ("phase_deg", "result_phase_deg"),
        ("damping_n_s_per_m", "result_damping_n_s_per_m"),
        ("force_n", "result_force_n"),
        ("power_w", "result_power_w"),
    )
    return (
        bool(generation)
        and all(row.get(key) == generation for key in linked)
        and all(_finite_sequence(row.get(left)) and _same(row, left, right) for left, right in arrays)
        and all(float(item) >= 0.0 for item in row.get("damping_n_s_per_m", []))
        and all(float(item) >= 0.0 for item in row.get("power_w", []))
        and row.get("stability_sign") == "stable"
        and row.get("result_stability_sign") == row.get("stability_sign")
        and str(row.get("mesh_owner", "")).startswith("mesh:")
        and row.get("result_mesh_owner") == row.get("mesh_owner")
        and _digest(row.get("bearing_dynamic_result_sha256"))
        and row.get("accepted_bearing_dynamic_result_sha256") == row.get("bearing_dynamic_result_sha256")
    )


def _demag_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("demag_minor_generation", "")).strip()
    linked = (
        "fieldpath_generation", "remanence_generation", "branch_generation",
        "loss_generation", "energy_generation", "temperature_generation",
        "material_generation", "mesh_generation", "result_generation",
    )
    arrays = (
        ("field_path_a_per_m", "result_field_path_a_per_m"),
        ("magnetization_a_per_m", "result_magnetization_a_per_m"),
    )
    scalars = (
        ("remanence_a_per_m", "result_remanence_a_per_m"),
        ("coercivity_a_per_m", "result_coercivity_a_per_m"),
        ("loss_energy_j", "result_loss_energy_j"),
        ("temperature_k", "result_temperature_k"),
    )
    return (
        bool(generation)
        and all(row.get(key) == generation for key in linked)
        and all(_finite_sequence(row.get(left)) and _same(row, left, right) for left, right in arrays)
        and all(
            isinstance(row.get(left), (int, float))
            and math.isfinite(float(row.get(left)))
            and row.get(left) >= 0.0
            and _same(row, left, right)
            for left, right in scalars[:3]
        )
        and isinstance(row.get("temperature_k"), (int, float))
        and float(row["temperature_k"]) > 0.0
        and _same(row, "temperature_k", "result_temperature_k")
        and str(row.get("material_owner", "")).startswith("material:")
        and row.get("result_material_owner") == row.get("material_owner")
        and str(row.get("mesh_owner", "")).startswith("mesh:")
        and row.get("result_mesh_owner") == row.get("mesh_owner")
        and _digest(row.get("demag_minor_result_sha256"))
        and row.get("accepted_demag_minor_result_sha256") == row.get("demag_minor_result_sha256")
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    """Return only checks for v44 records present in a neutral artifact."""
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    if _DYNAMIC in identity:
        row = identity[_DYNAMIC]
        checks["magnetic_force_v44_dynamic_bearing_identity"] = isinstance(row, Mapping) and _dynamic_ok(row)
    if _DEMAG in identity:
        row = identity[_DEMAG]
        checks["magnetic_force_v44_demag_minorloop_identity"] = isinstance(row, Mapping) and _demag_ok(row)
    checks.update(validate_public_v45_identity(identity))
    checks.update(validate_public_v46_identity(identity))
    checks.update(validate_public_v47_identity(identity))
    checks.update(validate_public_v48_identity(identity))
    checks.update(validate_public_v49_identity(identity))
    return checks
