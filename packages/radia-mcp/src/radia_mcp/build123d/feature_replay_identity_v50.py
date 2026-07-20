"""Assembly-mate, feature-history, healing, and tessellation checks for v50."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


MATE = "assembly_mate_constraint_dof_frame_occurrence_transform_owner_identity"
FEATURE = "fillet_chamfer_edge_selector_radius_topology_history_owner_identity"
HEALING = "occt_tolerance_healing_sewing_shell_solid_orientation_owner_identity"
STL = "stl_tessellation_linear_angular_deflection_triangle_normal_unit_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _unique_strings(value: object, prefix: str | None = None) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and bool(value)
        and all(isinstance(item, str) and item and (prefix is None or item.startswith(prefix)) for item in value)
        and len(set(value)) == len(value)
    )


def _rigid_frames(value: object, prefix: str | None = None) -> bool:
    if not isinstance(value, Mapping) or not value:
        return False
    for name, frame in value.items():
        if not isinstance(name, str) or not name or (prefix is not None and not name.startswith(prefix)):
            return False
        if not isinstance(frame, Sequence) or isinstance(frame, (str, bytes)) or len(frame) != 7:
            return False
        try:
            values = [float(item) for item in frame]
        except (TypeError, ValueError):
            return False
        if not all(math.isfinite(item) for item in values):
            return False
        quaternion_norm = math.sqrt(sum(item * item for item in values[3:]))
        if not math.isclose(quaternion_norm, 1.0, rel_tol=1e-9, abs_tol=1e-12):
            return False
    return True


def _mate_ok(row: Mapping[str, object]) -> bool:
    constraints = row.get("mate_constraints")
    dof = row.get("remaining_dof")
    frames = row.get("mate_frames")
    transforms = row.get("occurrence_transforms")
    frame_names = set(frames) if isinstance(frames, Mapping) else set()
    occurrence_names = (
        {name.split(":", 1)[1] for name in transforms if isinstance(name, str) and name.startswith("occurrence:")}
        if isinstance(transforms, Mapping)
        else set()
    )
    return (
        _generations(
            row,
            "constraint_generation",
            "dof_generation",
            "frame_generation",
            "occurrence_generation",
            "transform_generation",
            "owner_generation",
            "result_generation",
        )
        and _unique_strings(constraints)
        and row.get("result_mate_constraints") == constraints
        and isinstance(dof, int)
        and not isinstance(dof, bool)
        and dof >= 0
        and row.get("result_remaining_dof") == dof
        and _rigid_frames(frames)
        and row.get("result_mate_frames") == frames
        and _rigid_frames(transforms, "occurrence:")
        and frame_names == occurrence_names
        and row.get("result_occurrence_transforms") == transforms
        and str(row.get("assembly_owner") or "").startswith("assembly:")
        and row.get("result_assembly_owner") == row.get("assembly_owner")
        and _result(row)
    )


def _entity_map(value: object, key_prefix: str | None = None) -> bool:
    return (
        isinstance(value, Mapping)
        and bool(value)
        and all(
            isinstance(name, str)
            and name
            and (key_prefix is None or name.startswith(key_prefix))
            and _unique_strings(items)
            for name, items in value.items()
        )
    )


def _positive_finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) and float(value) > 0.0


def _feature_ok(row: Mapping[str, object]) -> bool:
    selectors = row.get("edge_selectors")
    history = row.get("topology_history")
    fillet = row.get("fillet_radius_m")
    chamfer = row.get("chamfer_distance_m")
    if not _entity_map(selectors) or set(selectors) != {"fillet", "chamfer"}:
        return False
    selected_edges = [edge for edges in selectors.values() for edge in edges]
    return (
        _generations(
            row,
            "selector_generation",
            "radius_generation",
            "topology_generation",
            "history_generation",
            "owner_generation",
            "result_generation",
        )
        and len(selected_edges) == len(set(selected_edges))
        and all(edge.startswith("edge:") for edge in selected_edges)
        and row.get("result_edge_selectors") == selectors
        and _positive_finite(fillet)
        and row.get("result_fillet_radius_m") == fillet
        and _positive_finite(chamfer)
        and row.get("result_chamfer_distance_m") == chamfer
        and _entity_map(history, "edge:")
        and set(history) == set(selected_edges)
        and row.get("result_topology_history") == history
        and str(row.get("shape_owner") or "").startswith("shape:")
        and row.get("result_shape_owner") == row.get("shape_owner")
        and _result(row)
    )


def _healing_ok(row: Mapping[str, object]) -> bool:
    input_tolerance = row.get("input_tolerance_m")
    sewing_tolerance = row.get("sewing_tolerance_m")
    shells = row.get("shell_count")
    solids = row.get("solid_count")
    orientations = row.get("shell_orientation_signs")
    return (
        _generations(
            row,
            "tolerance_generation",
            "healing_generation",
            "sewing_generation",
            "shell_generation",
            "solid_generation",
            "orientation_generation",
            "owner_generation",
            "result_generation",
        )
        and _positive_finite(input_tolerance)
        and row.get("result_input_tolerance_m") == input_tolerance
        and row.get("healing_applied") is True
        and row.get("result_healing_applied") is True
        and _positive_finite(sewing_tolerance)
        and float(sewing_tolerance) >= float(input_tolerance)
        and row.get("result_sewing_tolerance_m") == sewing_tolerance
        and isinstance(shells, int)
        and not isinstance(shells, bool)
        and shells > 0
        and row.get("result_shell_count") == shells
        and isinstance(solids, int)
        and not isinstance(solids, bool)
        and solids > 0
        and row.get("result_solid_count") == solids
        and isinstance(orientations, Mapping)
        and len(orientations) == shells
        and all(isinstance(name, str) and name.startswith("shell:") and sign in {-1, 1} for name, sign in orientations.items())
        and row.get("result_shell_orientation_signs") == orientations
        and str(row.get("shape_owner") or "").startswith("shape:")
        and row.get("result_shape_owner") == row.get("shape_owner")
        and _result(row)
    )


def _stl_ok(row: Mapping[str, object]) -> bool:
    linear = row.get("linear_deflection_m")
    angular = row.get("angular_deflection_rad")
    triangles = row.get("triangle_count")
    normals = row.get("normal_counts")
    return (
        _generations(
            row,
            "linear_generation",
            "angular_generation",
            "triangle_generation",
            "normal_generation",
            "unit_generation",
            "owner_generation",
            "result_generation",
        )
        and _positive_finite(linear)
        and row.get("result_linear_deflection_m") == linear
        and _positive_finite(angular)
        and float(angular) <= math.pi
        and row.get("result_angular_deflection_rad") == angular
        and isinstance(triangles, int)
        and not isinstance(triangles, bool)
        and triangles > 0
        and row.get("result_triangle_count") == triangles
        and isinstance(normals, Mapping)
        and set(normals) == {"outward", "inward", "degenerate"}
        and all(isinstance(value, int) and not isinstance(value, bool) and value >= 0 for value in normals.values())
        and sum(normals.values()) == triangles
        and normals["outward"] == triangles
        and normals["inward"] == 0
        and normals["degenerate"] == 0
        and row.get("result_normal_counts") == normals
        and row.get("length_unit") == "m"
        and row.get("result_length_unit") == "m"
        and str(row.get("mesh_owner") or "").startswith("mesh:")
        and row.get("result_mesh_owner") == row.get("mesh_owner")
        and _result(row)
    )


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {
        "policy": policy,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
    }


def _public_rows(payload: Mapping[str, object]) -> list[Mapping[str, object]]:
    rows: list[Mapping[str, object]] = []
    if isinstance(payload.get("reference"), list):
        rows.extend(row for row in payload["reference"] if isinstance(row, Mapping))
    measured = payload.get("measured")
    if isinstance(measured, Mapping):
        for values in measured.values():
            if isinstance(values, list):
                rows.extend(row for row in values if isinstance(row, Mapping))
    return rows


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    rows = _public_rows(payload)
    checks: dict[str, bool] = {}
    mates = [row.get(MATE) for row in rows if MATE in row]
    features = [row.get(FEATURE) for row in rows if FEATURE in row]
    if mates:
        checks["v50_assembly_mates_dof_frames_transforms_owner"] = (
            len(mates) == len(rows) and all(isinstance(item, Mapping) and _mate_ok(item) for item in mates)
        )
    if features:
        checks["v50_fillet_chamfer_selectors_radii_topology_owner"] = (
            len(features) == len(rows) and all(isinstance(item, Mapping) and _feature_ok(item) for item in features)
        )
    return _report("build123d_v50_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping) or not isinstance(payload.get("replay_identity"), Mapping):
        return {}
    identity = payload["replay_identity"]
    checks: dict[str, bool] = {}
    healing = identity.get(HEALING)
    stl = identity.get(STL)
    if healing is not None:
        checks["v50_occt_tolerance_healing_sewing_closure_orientation_owner"] = (
            isinstance(healing, Mapping) and _healing_ok(healing)
        )
    if stl is not None:
        checks["v50_stl_deflection_triangles_normals_units_owner"] = (
            isinstance(stl, Mapping) and _stl_ok(stl)
        )
    return _report("build123d_v50_source_identity_v1", checks) if checks else {}
