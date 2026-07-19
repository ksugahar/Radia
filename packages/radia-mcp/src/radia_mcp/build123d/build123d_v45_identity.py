"""Identity checks for held-out CAD result and source replay records."""

from __future__ import annotations

from collections.abc import Mapping, Sequence


def _sha(value: object) -> bool:
    text = str(value or "")
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text.lower())


def _same(identity: Mapping[str, object], *names: str) -> bool:
    return all(identity.get(f"result_{name}") == identity.get(name) for name in names)


def boolean_v45_ok(value: object) -> bool:
    if not isinstance(value, Mapping):
        return False
    try:
        thickness = float(value["shell_thickness_m"])
        radius = float(value["fillet_radius_m"])
        volume = float(value["volume_m3"])
        area = float(value["surface_area_m2"])
    except (KeyError, TypeError, ValueError):
        return False
    return (
        bool(str(value.get("generation") or ""))
        and thickness > 0.0 and radius > 0.0 and volume > 0.0 and area > 0.0
        and _same(value, "operation", "shell_thickness_m", "fillet_radius_m", "center_of_mass_m", "volume_m3", "surface_area_m2", "topology_signature", "shape_owner")
        and bool(str(value.get("shape_owner") or "").startswith("part:"))
        and _sha(value.get("boolean_brep_sha256"))
        and value.get("accepted_boolean_brep_sha256") == value.get("boolean_brep_sha256")
    )


def loft_v45_ok(value: object) -> bool:
    if not isinstance(value, Mapping):
        return False
    try:
        count = int(value["section_count"])
        volume = float(value["volume_m3"])
        inertia = value["inertia_tensor_kg_m2"]
    except (KeyError, TypeError, ValueError):
        return False
    return (
        bool(str(value.get("generation") or ""))
        and count >= 2 and volume > 0.0
        and value.get("section_orientation") == value.get("result_section_orientation") == "consistent_ccw"
        and value.get("tangent_continuity") == value.get("result_tangent_continuity") == "C1"
        and isinstance(inertia, Sequence) and inertia == value.get("result_inertia_tensor_kg_m2")
        and _same(value, "section_count", "section_area_m2", "volume_m3", "shape_owner")
        and bool(str(value.get("shape_owner") or "").startswith("part:"))
        and _sha(value.get("loft_brep_sha256"))
        and value.get("accepted_loft_brep_sha256") == value.get("loft_brep_sha256")
    )


def sketch_v45_ok(value: object) -> bool:
    if not isinstance(value, Mapping):
        return False
    order = value.get("constraint_order")
    return (
        bool(str(value.get("generation") or ""))
        and isinstance(order, Sequence) and list(order) == list(value.get("replayed_constraint_order") or [])
        and value.get("plane_frame") == value.get("replayed_plane_frame")
        and value.get("parameter_cache_key") == value.get("replayed_parameter_cache_key")
        and value.get("solver_status") == value.get("replayed_solver_status") == "solved"
        and value.get("shape_generation_id") == value.get("replayed_shape_generation_id")
        and value.get("shape_owner") == value.get("replayed_shape_owner")
        and bool(str(value.get("shape_owner") or "").startswith("headless:"))
        and _sha(value.get("result_sha256"))
        and value.get("accepted_result_sha256") == value.get("result_sha256")
    )


def step_v45_ok(value: object) -> bool:
    if not isinstance(value, Mapping):
        return False
    topology = value.get("face_topology")
    return (
        bool(str(value.get("generation") or ""))
        and value.get("length_unit") == value.get("replayed_length_unit") == "mm"
        and float(value.get("tessellation_tolerance_m", 0.0)) > 0.0
        and value.get("tessellation_tolerance_m") == value.get("replayed_tessellation_tolerance_m")
        and isinstance(topology, Mapping) and topology == value.get("replayed_face_topology")
        and value.get("brep_generation_id") == value.get("replayed_brep_generation_id")
        and value.get("export_owner") == value.get("replayed_export_owner")
        and bool(str(value.get("export_owner") or "").startswith("headless:"))
        and _sha(value.get("step_digest_sha256"))
        and value.get("replayed_step_digest_sha256") == value.get("step_digest_sha256")
        and value.get("accepted_result_sha256") == value.get("step_digest_sha256")
    )
