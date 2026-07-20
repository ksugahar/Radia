"""Assembly, STEP import, history, and mass-property identity checks for v53."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .assembly_tessellation_identity_v54 import (
    validate_public_identity as validate_public_v54_identity,
    validate_source_identity as validate_source_v54_identity,
)


MATE = "assembly_mate_frame_handedness_axis_offset_component_owner_identity"
STEP = "step_import_unit_tolerance_color_hierarchy_owner_identity"
HISTORY = "boolean_face_ancestry_fillet_chamfer_shape_owner_identity"
MASS = "massproperty_inertia_frame_centroid_density_owner_identity"
_LENGTH_UNITS = {"m", "mm", "cm", "inch"}


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], *names: str) -> bool:
    value = str(row.get("generation") or "")
    return bool(value) and all(row.get(name) == value for name in names)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _vector(value: object, length: int = 3) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and len(value) == length and all(_finite(item) for item in value)


def _unit(value: object) -> bool:
    return _vector(value) and math.isclose(sum(float(item) ** 2 for item in value), 1.0, rel_tol=0.0, abs_tol=1.0e-12)


def _frame(value: object) -> bool:
    if not isinstance(value, Mapping) or set(value) != {"origin", "x_dir", "z_dir"}:
        return False
    if not _vector(value["origin"]) or not _unit(value["x_dir"]) or not _unit(value["z_dir"]):
        return False
    return math.isclose(sum(float(a) * float(b) for a, b in zip(value["x_dir"], value["z_dir"])), 0.0, rel_tol=0.0, abs_tol=1.0e-12)


def _mate_ok(row: Mapping[str, object]) -> bool:
    pair = row.get("component_pair")
    return (
        _generation(row, "frame_generation", "handedness_generation", "axis_generation", "offset_generation", "component_generation", "owner_generation", "result_generation")
        and _frame(row.get("mate_frame"))
        and row.get("result_mate_frame") == row.get("mate_frame")
        and row.get("handedness") == "right"
        and row.get("result_handedness") == "right"
        and _unit(row.get("mate_axis"))
        and row.get("result_mate_axis") == row.get("mate_axis")
        and _finite(row.get("offset_m"))
        and row.get("result_offset_m") == row.get("offset_m")
        and isinstance(pair, Sequence)
        and not isinstance(pair, (str, bytes))
        and len(pair) == 2
        and len(set(pair)) == 2
        and all(isinstance(item, str) and item.startswith("component:") for item in pair)
        and row.get("result_component_pair") == pair
        and str(row.get("assembly_owner") or "").startswith("assembly:")
        and row.get("result_assembly_owner") == row.get("assembly_owner")
        and _result(row)
    )


def _colors(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and bool(value)
        and all(
            isinstance(component, str)
            and component.startswith("component:")
            and _vector(color)
            and all(0.0 <= float(channel) <= 1.0 for channel in color)
            for component, color in value.items()
        )
    )


def _hierarchy(value: object, components: set[str]) -> bool:
    if not isinstance(value, Mapping) or not value:
        return False
    leaves: list[str] = []
    for assembly, children in value.items():
        if not isinstance(assembly, str) or not assembly.startswith("assembly:"):
            return False
        if not isinstance(children, Sequence) or isinstance(children, (str, bytes)) or not children:
            return False
        if len(children) != len(set(children)) or not all(isinstance(child, str) and child.startswith("component:") for child in children):
            return False
        leaves.extend(children)
    return len(leaves) == len(set(leaves)) and set(leaves) == components


def _step_ok(row: Mapping[str, object]) -> bool:
    colors = row.get("component_colors")
    hierarchy = row.get("component_hierarchy")
    return (
        _generation(row, "unit_generation", "tolerance_generation", "color_generation", "hierarchy_generation", "owner_generation", "result_generation")
        and row.get("length_unit") in _LENGTH_UNITS
        and row.get("result_length_unit") == row.get("length_unit")
        and _finite(row.get("linear_tolerance_m"))
        and float(row["linear_tolerance_m"]) > 0.0
        and row.get("result_linear_tolerance_m") == row.get("linear_tolerance_m")
        and _colors(colors)
        and row.get("result_component_colors") == colors
        and _hierarchy(hierarchy, set(colors))
        and row.get("result_component_hierarchy") == hierarchy
        and str(row.get("document_owner") or "").startswith("document:")
        and row.get("result_document_owner") == row.get("document_owner")
        and _result(row)
    )


def _history_ok(row: Mapping[str, object]) -> bool:
    ancestry = row.get("face_ancestry")
    ancestry_ok = isinstance(ancestry, Mapping) and bool(ancestry)
    if ancestry_ok:
        ancestry_ok = all(
            isinstance(face, str)
            and face.startswith("face:")
            and isinstance(parents, Sequence)
            and not isinstance(parents, (str, bytes))
            and bool(parents)
            and len(parents) == len(set(parents))
            and all(isinstance(parent, str) and parent.startswith("face:") for parent in parents)
            for face, parents in ancestry.items()
        )
    fillets = row.get("fillet_edges")
    chamfers = row.get("chamfer_edges")
    edge_lists_ok = all(
        isinstance(items, Sequence)
        and not isinstance(items, (str, bytes))
        and len(items) == len(set(items))
        and all(isinstance(edge, str) and edge.startswith("edge:") for edge in items)
        for items in (fillets, chamfers)
    )
    return (
        _generation(row, "boolean_generation", "ancestry_generation", "fillet_generation", "chamfer_generation", "owner_generation", "result_generation")
        and row.get("boolean_operation") in {"cut", "fuse", "intersect"}
        and row.get("replayed_boolean_operation") == row.get("boolean_operation")
        and ancestry_ok
        and row.get("replayed_face_ancestry") == ancestry
        and edge_lists_ok
        and set(fillets).isdisjoint(chamfers)
        and row.get("replayed_fillet_edges") == fillets
        and row.get("replayed_chamfer_edges") == chamfers
        and str(row.get("shape_owner") or "").startswith("shape:")
        and row.get("replayed_shape_owner") == row.get("shape_owner")
        and _result(row)
    )


def _inertia_tensor(value: object) -> bool:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 3:
        return False
    if not all(_vector(row) for row in value):
        return False
    matrix = [[float(item) for item in row] for row in value]
    if not all(math.isclose(matrix[i][j], matrix[j][i], rel_tol=0.0, abs_tol=1.0e-12) for i in range(3) for j in range(3)):
        return False
    a, b, c = matrix[0]
    _, d, e = matrix[1]
    _, _, f = matrix[2]
    determinant = a * (d * f - e * e) - b * (b * f - c * e) + c * (b * e - c * d)
    positive_definite = a > 0.0 and a * d - b * b > 0.0 and determinant > 0.0
    triangle = a + d >= f - 1.0e-12 and a + f >= d - 1.0e-12 and d + f >= a - 1.0e-12
    return positive_definite and triangle


def _mass_ok(row: Mapping[str, object]) -> bool:
    return (
        _generation(row, "tensor_generation", "frame_generation", "centroid_generation", "density_generation", "owner_generation", "result_generation")
        and _inertia_tensor(row.get("inertia_tensor_kg_m2"))
        and row.get("replayed_inertia_tensor_kg_m2") == row.get("inertia_tensor_kg_m2")
        and _frame(row.get("reference_frame"))
        and row.get("replayed_reference_frame") == row.get("reference_frame")
        and _vector(row.get("centroid_m"))
        and row.get("replayed_centroid_m") == row.get("centroid_m")
        and _finite(row.get("density_kg_m3"))
        and float(row["density_kg_m3"]) > 0.0
        and row.get("replayed_density_kg_m3") == row.get("density_kg_m3")
        and str(row.get("solid_owner") or "").startswith("solid:")
        and row.get("replayed_solid_owner") == row.get("solid_owner")
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
    v54 = validate_public_v54_identity(payload)
    if v54:
        checks.update(v54["checks"])
    mates = [row.get(MATE) for row in rows if MATE in row]
    steps = [row.get(STEP) for row in rows if STEP in row]
    if mates:
        checks["v53_mate_frame_handedness_axis_offset_component_owner"] = len(mates) == len(rows) and all(isinstance(item, Mapping) and _mate_ok(item) for item in mates)
    if steps:
        checks["v53_step_unit_tolerance_color_hierarchy_owner"] = len(steps) == len(rows) and all(isinstance(item, Mapping) and _step_ok(item) for item in steps)
    return _report("build123d_v53_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping) or not isinstance(payload.get("replay_identity"), Mapping):
        return {}
    identity = payload["replay_identity"]
    checks: dict[str, bool] = {}
    v54 = validate_source_v54_identity(payload)
    if v54:
        checks.update(v54["checks"])
    if identity.get(HISTORY) is not None:
        checks["v53_boolean_face_ancestry_fillet_chamfer_owner"] = isinstance(identity[HISTORY], Mapping) and _history_ok(identity[HISTORY])
    if identity.get(MASS) is not None:
        checks["v53_mass_inertia_frame_centroid_density_owner"] = isinstance(identity[MASS], Mapping) and _mass_ok(identity[MASS])
    return _report("build123d_v53_source_identity_v1", checks) if checks else {}
