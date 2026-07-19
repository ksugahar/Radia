"""Neutral v46 build123d CAD replay identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _public_mixed(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation", "")).strip()
    placement = row.get("placement_m")
    inertia = row.get("inertia_frame")
    return (
        bool(generation)
        and row.get("placement_generation") == row.get("result_placement_generation") == generation
        and row.get("unit_scale_generation") == row.get("result_unit_scale_generation") == generation
        and row.get("inertia_frame_generation") == row.get("result_inertia_frame_generation") == generation
        and row.get("mass_property_generation") == row.get("result_mass_property_generation") == generation
        and isinstance(placement, list)
        and len(placement) == 3
        and all(isinstance(x, (int, float)) and math.isfinite(float(x)) for x in placement)
        and row.get("result_placement_m") == placement
        and isinstance(row.get("unit_scale_to_si"), (int, float))
        and math.isfinite(float(row["unit_scale_to_si"]))
        and float(row["unit_scale_to_si"]) > 0.0
        and row.get("result_unit_scale_to_si") == row.get("unit_scale_to_si")
        and isinstance(inertia, str)
        and inertia == row.get("result_inertia_frame")
        and row.get("finite_mass_property_status") == row.get("result_finite_mass_property_status") == "finite"
        and bool(str(row.get("shape_owner") or ""))
        and row.get("accepted_shape_owner") == row.get("shape_owner")
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def _public_boolean(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation", "")).strip()
    topology = row.get("topology_signature")
    return (
        bool(generation)
        and row.get("boolean_generation") == row.get("result_boolean_generation") == generation
        and row.get("fillet_generation") == row.get("result_fillet_generation") == generation
        and row.get("recovery_generation") == row.get("result_recovery_generation") == generation
        and row.get("topology_generation") == row.get("result_topology_generation") == generation
        and row.get("boolean_status") == row.get("result_boolean_status") == "recovered"
        and isinstance(topology, dict)
        and topology == row.get("result_topology_signature")
        and row.get("partial_shape_status") == row.get("result_partial_shape_status") == "complete"
        and bool(str(row.get("shape_owner") or ""))
        and row.get("accepted_shape_owner") == row.get("shape_owner")
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    rows = []
    if isinstance(payload.get("reference"), list):
        rows.extend(payload["reference"])
    measured = payload.get("measured")
    if isinstance(measured, Mapping):
        for values in measured.values():
            if isinstance(values, list):
                rows.extend(values)
    mixed = [row.get("placement_unit_scale_inertia_coordinate_frame_nonfinite_identity") for row in rows if isinstance(row, Mapping) and "placement_unit_scale_inertia_coordinate_frame_nonfinite_identity" in row]
    boolean = [row.get("boolean_partial_shape_fillet_recovery_topology_identity") for row in rows if isinstance(row, Mapping) and "boolean_partial_shape_fillet_recovery_topology_identity" in row]
    checks = {}
    if mixed:
        checks["v46_placement_unit_inertia_frame_identity"] = all(isinstance(row, Mapping) and _public_mixed(row) for row in mixed) and len(mixed) == len(rows)
    if boolean:
        checks["v46_boolean_partial_topology_identity"] = all(isinstance(row, Mapping) and _public_boolean(row) for row in boolean) and len(boolean) == len(rows)
    if not checks:
        return {}
    return {"policy": "build123d_v46_public_identity_v1", "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, ok in checks.items() if not ok]}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    identity = payload.get("replay_identity")
    if not isinstance(identity, Mapping):
        return {}
    sketch = identity.get("sketch_solver_partial_constraint_warning_plane_unit_identity")
    step = identity.get("step_import_partial_face_tolerance_coordinate_frame_checksum_identity")
    checks = {}
    if sketch is not None:
        generation = str(sketch.get("generation", "")).strip() if isinstance(sketch, Mapping) else ""
        checks["v46_sketch_partial_constraint_identity"] = isinstance(sketch, Mapping) and bool(generation) and sketch.get("solver_warning") == sketch.get("result_solver_warning") == "none" and sketch.get("constraint_state") == sketch.get("result_constraint_state") == "fully_constrained" and sketch.get("plane") == sketch.get("result_plane") == "XY" and sketch.get("unit_scale_to_si") == sketch.get("result_unit_scale_to_si") == 1.0 and sketch.get("result_generation") == generation and _digest(sketch.get("result_sha256")) and sketch.get("accepted_result_sha256") == sketch.get("result_sha256")
    if step is not None:
        generation = str(step.get("generation", "")).strip() if isinstance(step, Mapping) else ""
        checksum = step.get("checksum_sha256") if isinstance(step, Mapping) else None
        checks["v46_step_partial_face_identity"] = isinstance(step, Mapping) and bool(generation) and step.get("partial_face_count") == step.get("result_partial_face_count") == 0 and step.get("tolerance_m") == step.get("result_tolerance_m") and step.get("coordinate_frame") == step.get("result_coordinate_frame") == "global_cartesian" and _digest(checksum) and step.get("result_checksum_sha256") == checksum and step.get("result_generation") == generation and _digest(step.get("result_sha256")) and step.get("accepted_result_sha256") == step.get("result_sha256")
    if not checks:
        return {}
    return {"policy": "build123d_v46_source_identity_v1", "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, ok in checks.items() if not ok]}
