"""Assembly, sketch, STEP, and Boolean identity checks for v49 CAD artifacts."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


ASSEMBLY = "assembly_mass_density_material_occurrence_transform_suppression_owner_identity"
SKETCH = "sketch_constraint_dof_plane_unit_profile_wire_owner_identity"
STEP = "step_schema_assembly_color_layer_unit_tolerance_owner_identity"
BOOLEAN = "boolean_deleted_subshape_selector_adjacency_mass_cache_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_vector(value: object, length: int | None = None) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and (length is None or len(value) == length)
        and bool(value)
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def _assembly_ok(row: Mapping[str, object]) -> bool:
    materials = row.get("material_by_occurrence")
    densities = row.get("density_kg_m3_by_occurrence")
    transforms = row.get("occurrence_transforms")
    suppressed = row.get("suppressed_occurrences")
    mass = row.get("assembly_mass_kg")
    occurrences = set(materials) if isinstance(materials, Mapping) else set()
    return (
        _generations(row, "density_generation", "material_generation", "occurrence_generation", "transform_generation", "suppression_generation", "mass_generation", "result_generation")
        and bool(occurrences)
        and all(isinstance(value, str) and value for value in materials.values())
        and row.get("result_material_by_occurrence") == materials
        and isinstance(densities, Mapping)
        and set(densities) == occurrences
        and all(isinstance(value, (int, float)) and math.isfinite(float(value)) and float(value) > 0.0 for value in densities.values())
        and row.get("result_density_kg_m3_by_occurrence") == densities
        and isinstance(transforms, Mapping)
        and set(transforms) == occurrences
        and all(_finite_vector(value, 4) for value in transforms.values())
        and row.get("result_occurrence_transforms") == transforms
        and isinstance(suppressed, list)
        and len(set(suppressed)) == len(suppressed)
        and set(suppressed).issubset(occurrences)
        and row.get("result_suppressed_occurrences") == suppressed
        and isinstance(mass, (int, float))
        and math.isfinite(float(mass))
        and float(mass) >= 0.0
        and row.get("result_assembly_mass_kg") == mass
        and str(row.get("assembly_owner") or "").startswith("assembly:")
        and row.get("result_assembly_owner") == row.get("assembly_owner")
        and _result(row)
    )


def _sketch_ok(row: Mapping[str, object]) -> bool:
    constraints = row.get("constraints")
    dof = row.get("remaining_dof")
    wires = row.get("profile_wires")
    return (
        _generations(row, "constraint_generation", "dof_generation", "plane_generation", "unit_generation", "wire_generation", "result_generation")
        and isinstance(constraints, list)
        and bool(constraints)
        and len(set(constraints)) == len(constraints)
        and all(isinstance(value, str) and value for value in constraints)
        and row.get("result_constraints") == constraints
        and isinstance(dof, int)
        and dof >= 0
        and row.get("result_remaining_dof") == dof
        and row.get("work_plane") in {"Plane.XY", "Plane.XZ", "Plane.YZ"}
        and row.get("result_work_plane") == row.get("work_plane")
        and row.get("length_unit") in {"m", "mm", "cm", "in"}
        and row.get("result_length_unit") == row.get("length_unit")
        and isinstance(wires, list)
        and bool(wires)
        and len(set(wires)) == len(wires)
        and all(isinstance(value, str) and value.startswith("wire:") for value in wires)
        and row.get("result_profile_wires") == wires
        and str(row.get("sketch_owner") or "").startswith("sketch:")
        and row.get("result_sketch_owner") == row.get("sketch_owner")
        and _result(row)
    )


def _step_ok(row: Mapping[str, object]) -> bool:
    structure = row.get("assembly_structure")
    colors = row.get("color_map")
    layers = row.get("layer_map")
    tolerance = row.get("model_tolerance_m")
    parts = {part for values in structure.values() for part in values} if isinstance(structure, Mapping) else set()
    return (
        _generations(row, "schema_generation", "assembly_generation", "metadata_generation", "unit_generation", "tolerance_generation", "result_generation")
        and row.get("ap_schema") in {"AP203", "AP214", "AP242"}
        and row.get("result_ap_schema") == row.get("ap_schema")
        and isinstance(structure, Mapping)
        and bool(structure)
        and all(isinstance(values, list) and bool(values) and len(set(values)) == len(values) for values in structure.values())
        and row.get("result_assembly_structure") == structure
        and isinstance(colors, Mapping)
        and set(colors) == parts
        and all(_finite_vector(value, 3) and all(0.0 <= float(channel) <= 1.0 for channel in value) for value in colors.values())
        and row.get("result_color_map") == colors
        and isinstance(layers, Mapping)
        and set(layers) == parts
        and all(isinstance(value, str) and value for value in layers.values())
        and row.get("result_layer_map") == layers
        and row.get("length_unit") in {"m", "mm", "cm", "in"}
        and row.get("result_length_unit") == row.get("length_unit")
        and isinstance(tolerance, (int, float))
        and math.isfinite(float(tolerance))
        and float(tolerance) > 0.0
        and row.get("result_model_tolerance_m") == tolerance
        and str(row.get("import_owner") or "").startswith("import:")
        and row.get("result_import_owner") == row.get("import_owner")
        and _result(row)
    )


def _entity_map(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and bool(value)
        and all(isinstance(rows, list) and bool(rows) and len(set(rows)) == len(rows) and all(isinstance(item, str) and item for item in rows) for rows in value.values())
    )


def _boolean_ok(row: Mapping[str, object]) -> bool:
    deleted = row.get("deleted_subshapes")
    selectors = row.get("selector_results")
    adjacency = row.get("adjacency_map")
    mass = row.get("mass_kg")
    return (
        _generations(row, "deletion_generation", "selector_generation", "adjacency_generation", "mass_generation", "cache_generation", "result_generation")
        and isinstance(deleted, list)
        and bool(deleted)
        and len(set(deleted)) == len(deleted)
        and row.get("result_deleted_subshapes") == deleted
        and _entity_map(selectors)
        and row.get("result_selector_results") == selectors
        and _entity_map(adjacency)
        and row.get("result_adjacency_map") == adjacency
        and isinstance(mass, (int, float))
        and math.isfinite(float(mass))
        and float(mass) >= 0.0
        and row.get("cached_mass_kg") == mass
        and row.get("result_mass_kg") == mass
        and _digest(row.get("cache_shape_sha256"))
        and row.get("result_cache_shape_sha256") == row.get("cache_shape_sha256")
        and str(row.get("history_owner") or "").startswith("history:")
        and row.get("result_history_owner") == row.get("history_owner")
        and _result(row)
    )


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {"policy": policy, "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, ok in checks.items() if not ok]}


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    rows: list[Mapping[str, object]] = []
    if isinstance(payload.get("reference"), list):
        rows.extend(row for row in payload["reference"] if isinstance(row, Mapping))
    measured = payload.get("measured")
    if isinstance(measured, Mapping):
        for values in measured.values():
            if isinstance(values, list):
                rows.extend(row for row in values if isinstance(row, Mapping))
    checks: dict[str, bool] = {}
    assemblies = [row.get(ASSEMBLY) for row in rows if ASSEMBLY in row]
    sketches = [row.get(SKETCH) for row in rows if SKETCH in row]
    if assemblies:
        checks["v49_assembly_material_density_transform_suppression_mass_owner"] = len(assemblies) == len(rows) and all(isinstance(item, Mapping) and _assembly_ok(item) for item in assemblies)
    if sketches:
        checks["v49_sketch_constraints_dof_plane_unit_wires_owner"] = len(sketches) == len(rows) and all(isinstance(item, Mapping) and _sketch_ok(item) for item in sketches)
    return _report("build123d_v49_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping) or not isinstance(payload.get("replay_identity"), Mapping):
        return {}
    identity = payload["replay_identity"]
    checks: dict[str, bool] = {}
    step = identity.get(STEP)
    boolean = identity.get(BOOLEAN)
    if step is not None:
        checks["v49_step_schema_assembly_metadata_unit_tolerance_owner"] = isinstance(step, Mapping) and _step_ok(step)
    if boolean is not None:
        checks["v49_boolean_deleted_selector_adjacency_mass_cache_owner"] = isinstance(boolean, Mapping) and _boolean_ok(boolean)
    return _report("build123d_v49_source_identity_v1", checks) if checks else {}

