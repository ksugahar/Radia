"""Selector, workplane, BREP, and glTF identity checks for v52."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


SELECTOR = "selector_query_order_topology_label_owner_identity"
WORKPLANE = "workplane_local_coordinate_pending_edge_wire_owner_identity"
BREP = "brep_occversion_location_tshape_serialization_owner_identity"
GLTF = "gltf_axis_scale_material_instance_scene_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], *names: str) -> bool:
    value = str(row.get("generation") or "")
    return bool(value) and all(row.get(name) == value for name in names)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_vector(value: object, length: int) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and len(value) == length
        and all(isinstance(item, (int, float)) and not isinstance(item, bool) and math.isfinite(float(item)) for item in value)
    )


def _selector_ok(row: Mapping[str, object]) -> bool:
    selected = row.get("selected_topology_order")
    labels = row.get("topology_labels")
    return (
        _generation(row, "query_generation", "order_generation", "topology_generation", "label_generation", "owner_generation", "result_generation")
        and "sort_by(" in str(row.get("selector_query") or "")
        and row.get("result_selector_query") == row.get("selector_query")
        and isinstance(selected, Sequence)
        and not isinstance(selected, (str, bytes))
        and bool(selected)
        and len(selected) == len(set(selected))
        and all(isinstance(item, str) and item.startswith(("face:", "edge:", "vertex:")) for item in selected)
        and row.get("result_selected_topology_order") == selected
        and str(row.get("topology_revision") or "").startswith("topology:")
        and row.get("result_topology_revision") == row.get("topology_revision")
        and isinstance(labels, Mapping)
        and set(labels) == set(selected)
        and all(isinstance(label, str) and label for label in labels.values())
        and row.get("result_topology_labels") == labels
        and str(row.get("shape_owner") or "").startswith("shape:")
        and row.get("result_shape_owner") == row.get("shape_owner")
        and _result(row)
    )


def _workplane_ok(row: Mapping[str, object]) -> bool:
    frame = row.get("local_frame")
    edges = row.get("pending_edge_order")
    if not isinstance(frame, Mapping) or set(frame) != {"origin", "x_dir", "z_dir"}:
        return False
    if not all(_finite_vector(frame[name], 3) for name in frame):
        return False
    x_dir = [float(item) for item in frame["x_dir"]]
    z_dir = [float(item) for item in frame["z_dir"]]
    frame_ok = math.isclose(sum(item * item for item in x_dir), 1.0, abs_tol=1.0e-12) and math.isclose(sum(item * item for item in z_dir), 1.0, abs_tol=1.0e-12) and math.isclose(sum(a * b for a, b in zip(x_dir, z_dir)), 0.0, abs_tol=1.0e-12)
    return (
        _generation(row, "coordinate_generation", "edge_generation", "closure_generation", "owner_generation", "result_generation")
        and frame_ok
        and row.get("result_local_frame") == frame
        and isinstance(edges, Sequence)
        and not isinstance(edges, (str, bytes))
        and len(edges) >= 3
        and len(edges) == len(set(edges))
        and all(isinstance(edge, str) and edge.startswith("edge:") for edge in edges)
        and row.get("result_pending_edge_order") == edges
        and row.get("wire_closed") is True
        and row.get("result_wire_closed") is True
        and str(row.get("builder_owner") or "").startswith("builder:")
        and row.get("result_builder_owner") == row.get("builder_owner")
        and _result(row)
    )


def _brep_ok(row: Mapping[str, object]) -> bool:
    version = str(row.get("occt_version") or "")
    return (
        _generation(row, "version_generation", "location_generation", "tshape_generation", "serialization_generation", "owner_generation", "result_generation")
        and len(version.split(".")) == 3
        and all(part.isdigit() for part in version.split("."))
        and row.get("replayed_occt_version") == version
        and _finite_vector(row.get("shape_location"), 7)
        and row.get("replayed_shape_location") == row.get("shape_location")
        and _digest(row.get("tshape_sha256"))
        and row.get("replayed_tshape_sha256") == row.get("tshape_sha256")
        and _digest(row.get("serialization_sha256"))
        and row.get("replayed_serialization_sha256") == row.get("serialization_sha256")
        and str(row.get("shape_owner") or "").startswith("shape:")
        and row.get("replayed_shape_owner") == row.get("shape_owner")
        and _result(row)
    )


def _gltf_ok(row: Mapping[str, object]) -> bool:
    materials = row.get("part_materials")
    instances = row.get("instance_transforms")
    scale = row.get("length_scale_to_m")
    return (
        _generation(row, "axis_generation", "scale_generation", "material_generation", "instance_generation", "owner_generation", "result_generation")
        and row.get("axis_convention") == "Y_up_right_handed"
        and row.get("replayed_axis_convention") == row.get("axis_convention")
        and isinstance(scale, (int, float))
        and not isinstance(scale, bool)
        and math.isfinite(float(scale))
        and float(scale) > 0.0
        and row.get("replayed_length_scale_to_m") == scale
        and isinstance(materials, Mapping)
        and bool(materials)
        and all(isinstance(part, str) and part.startswith("part:") and isinstance(material, str) and material.startswith("material:") for part, material in materials.items())
        and row.get("replayed_part_materials") == materials
        and isinstance(instances, Mapping)
        and bool(instances)
        and all(isinstance(name, str) and name.startswith("instance:") and _finite_vector(transform, 7) for name, transform in instances.items())
        and row.get("replayed_instance_transforms") == instances
        and str(row.get("scene_owner") or "").startswith("scene:")
        and row.get("replayed_scene_owner") == row.get("scene_owner")
        and _result(row)
    )


def _public_rows(payload: Mapping[str, object]) -> list[Mapping[str, object]]:
    rows: list[Mapping[str, object]] = []
    reference = payload.get("reference")
    if isinstance(reference, Sequence) and not isinstance(reference, (str, bytes)):
        rows.extend(item for item in reference if isinstance(item, Mapping))
    measured = payload.get("measured")
    if isinstance(measured, Mapping):
        for family in measured.values():
            if isinstance(family, Sequence) and not isinstance(family, (str, bytes)):
                rows.extend(item for item in family if isinstance(item, Mapping))
    return rows


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {"policy": policy, "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, accepted in checks.items() if not accepted]}


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    rows = _public_rows(payload)
    checks: dict[str, bool] = {}
    selectors = [row.get(SELECTOR) for row in rows if SELECTOR in row]
    workplanes = [row.get(WORKPLANE) for row in rows if WORKPLANE in row]
    if selectors:
        checks["v52_selector_query_order_topology_label_owner"] = len(selectors) == len(rows) and all(isinstance(item, Mapping) and _selector_ok(item) for item in selectors)
    if workplanes:
        checks["v52_workplane_coordinate_edge_closure_owner"] = len(workplanes) == len(rows) and all(isinstance(item, Mapping) and _workplane_ok(item) for item in workplanes)
    return _report("build123d_v52_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping) or not isinstance(payload.get("replay_identity"), Mapping):
        return {}
    identity = payload["replay_identity"]
    checks: dict[str, bool] = {}
    if identity.get(BREP) is not None:
        checks["v52_brep_version_location_tshape_serialization_owner"] = isinstance(identity[BREP], Mapping) and _brep_ok(identity[BREP])
    if identity.get(GLTF) is not None:
        checks["v52_gltf_axis_scale_material_instance_owner"] = isinstance(identity[GLTF], Mapping) and _gltf_ok(identity[GLTF])
    return _report("build123d_v52_source_identity_v1", checks) if checks else {}
