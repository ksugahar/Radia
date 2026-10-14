"""Solver-neutral labelled-pole and CoilBuilder exchange contracts."""

from __future__ import annotations

import hashlib
import json
import math
from typing import Any


SCHEMA = "cae-ai-lab.shared-pole-coil.v1"
_COIL_KIND = "rounded_rectangle_racetrack"
_BODY_KIND = "axis_aligned_box"


def _canonical_number(value: Any, name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite")
    if abs(number) < 1.0e-15:
        return 0.0
    return float(f"{number:.15g}")


def _canonical_sha256(payload: dict[str, Any]) -> str:
    encoded = json.dumps(
        payload, ensure_ascii=True, sort_keys=True, separators=(",", ":")
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _nonempty_text(value: Any, name: str) -> str:
    text = str(value or "").strip()
    if not text:
        raise ValueError(f"{name} must be nonempty")
    return text


def _vector3(value: Any, name: str) -> list[float]:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"{name} must contain three coordinates")
    return [_canonical_number(item, f"{name}[{index}]") for index, item in enumerate(value)]


def _right_handed_frame(value: Any) -> list[list[float]]:
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError("orientation_rows must be a 3x3 matrix")
    rows = [_vector3(row, f"orientation_rows[{index}]") for index, row in enumerate(value)]
    for index, row in enumerate(rows):
        norm = math.sqrt(sum(item * item for item in row))
        if abs(norm - 1.0) > 1.0e-10:
            raise ValueError(f"orientation row {index} must have unit norm")
    for first, second in ((0, 1), (0, 2), (1, 2)):
        dot = sum(rows[first][index] * rows[second][index] for index in range(3))
        if abs(dot) > 1.0e-10:
            raise ValueError("orientation_rows must be orthogonal")
    determinant = (
        rows[0][0] * (rows[1][1] * rows[2][2] - rows[1][2] * rows[2][1])
        - rows[0][1] * (rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0])
        + rows[0][2] * (rows[1][0] * rows[2][1] - rows[1][1] * rows[2][0])
    )
    if determinant <= 0.0 or abs(determinant - 1.0) > 1.0e-10:
        raise ValueError("orientation_rows must form a right-handed orthonormal frame")
    return rows


def _normalize_body(raw: Any, index: int) -> dict[str, Any]:
    if not isinstance(raw, dict):
        raise ValueError(f"bodies[{index}] must be a mapping")
    if raw.get("kind") != _BODY_KIND:
        raise ValueError(f"unsupported body kind: {raw.get('kind')!r}")
    if raw.get("role") != "magnetic_pole":
        raise ValueError("shared bodies must be explicitly labelled magnetic_pole")
    lower = _vector3(raw.get("lower_m"), f"bodies[{index}].lower_m")
    upper = _vector3(raw.get("upper_m"), f"bodies[{index}].upper_m")
    if any(high <= low for low, high in zip(lower, upper)):
        raise ValueError("axis-aligned box upper_m must exceed lower_m on every axis")
    return {
        "name": _nonempty_text(raw.get("name"), f"bodies[{index}].name"),
        "role": "magnetic_pole",
        "kind": _BODY_KIND,
        "lower_m": lower,
        "upper_m": upper,
        "material_label": _nonempty_text(
            raw.get("material_label"), f"bodies[{index}].material_label"
        ),
    }


def _normalize_coil(raw: Any, index: int) -> dict[str, Any]:
    if not isinstance(raw, dict):
        raise ValueError(f"coils[{index}] must be a mapping")
    if raw.get("kind") != _COIL_KIND:
        raise ValueError(f"unsupported coil kind: {raw.get('kind')!r}")
    positive_names = (
        "corner_centerline_radius_m",
        "cross_section_width_m",
        "cross_section_height_m",
    )
    values = {
        name: _canonical_number(raw.get(name), f"coils[{index}].{name}")
        for name in positive_names
    }
    for name, value in values.items():
        if value <= 0.0:
            raise ValueError(f"{name} must be positive")
    straight_x = _canonical_number(raw.get("straight_x_m"), "straight_x_m")
    straight_y = _canonical_number(raw.get("straight_y_m"), "straight_y_m")
    if straight_x < 0.0 or straight_y < 0.0 or (straight_x == 0.0 and straight_y == 0.0):
        raise ValueError("racetrack straight lengths must be nonnegative and not both zero")
    current = _canonical_number(raw.get("current_A"), f"coils[{index}].current_A")
    turns_value = _canonical_number(raw.get("turns"), f"coils[{index}].turns")
    turns = int(turns_value)
    if turns < 1 or turns_value != turns:
        raise ValueError("turns must be a positive integer")
    conductor_model = _nonempty_text(raw.get("conductor_model"), "conductor_model")
    if conductor_model != "solid_uniform_current_density":
        raise ValueError(f"unsupported conductor_model: {conductor_model!r}")
    density = raw.get("current_density_A_per_m2")
    density_value = None
    if density is not None:
        density_value = _canonical_number(density, "current_density_A_per_m2")
        expected_current = density_value * values["cross_section_width_m"] * values["cross_section_height_m"]
        scale = max(abs(current), abs(expected_current), 1.0)
        if abs(current - expected_current) > 1.0e-12 * scale:
            raise ValueError("current_A is inconsistent with current density and cross-section")
    return {
        "name": _nonempty_text(raw.get("name"), f"coils[{index}].name"),
        "kind": _COIL_KIND,
        "corner_centerline_radius_m": values["corner_centerline_radius_m"],
        "straight_x_m": straight_x,
        "straight_y_m": straight_y,
        "cross_section_width_m": values["cross_section_width_m"],
        "cross_section_height_m": values["cross_section_height_m"],
        "centre_m": _vector3(raw.get("centre_m"), f"coils[{index}].centre_m"),
        "orientation_rows": _right_handed_frame(raw.get("orientation_rows")),
        "current_A": current,
        "turns": turns,
        "current_direction": "segment_path_start_to_end",
        "conductor_model": conductor_model,
        "region_label": _nonempty_text(raw.get("region_label"), "region_label"),
        **(
            {"current_density_A_per_m2": density_value}
            if density_value is not None
            else {}
        ),
    }


def normalize_shared_pole_coil_contract(payload: dict[str, Any]) -> dict[str, Any]:
    """Return the canonical SI subset used by both source and Radia adapters."""

    if not isinstance(payload, dict):
        raise ValueError("payload must be a mapping")
    if payload.get("schema") != SCHEMA:
        raise ValueError(f"schema must be {SCHEMA!r}")
    if payload.get("coordinate_system") != "right-handed Cartesian":
        raise ValueError("coordinate_system must be right-handed Cartesian")
    if payload.get("length_unit") != "m" or payload.get("current_unit") != "A":
        raise ValueError("shared exchange requires SI length_unit=m and current_unit=A")
    raw_bodies = payload.get("bodies")
    raw_coils = payload.get("coils")
    if not isinstance(raw_bodies, list) or not raw_bodies:
        raise ValueError("bodies must contain at least one labelled magnetic pole")
    if not isinstance(raw_coils, list) or not raw_coils:
        raise ValueError("coils must contain at least one supported CoilBuilder path")
    bodies = [_normalize_body(raw, index) for index, raw in enumerate(raw_bodies)]
    coils = [_normalize_coil(raw, index) for index, raw in enumerate(raw_coils)]
    names = [item["name"] for item in bodies + coils]
    if len(set(names)) != len(names):
        raise ValueError("body and coil names must be unique")
    return {
        "schema": SCHEMA,
        "component_label": _nonempty_text(payload.get("component_label"), "component_label"),
        "coordinate_system": "right-handed Cartesian",
        "length_unit": "m",
        "current_unit": "A",
        "bodies": bodies,
        "coils": coils,
    }


def _coilbuilder_racetrack_identity(coil: dict[str, Any]) -> dict[str, Any]:
    """Reproduce the CoilBuilder racetrack identity without a native import.

    FastMCP executes synchronous tools in a worker thread. Importing the full
    numerical/CAD stack there can stall some Windows runtimes, so the MCP
    adapter emits the exact straight/quarter-arc manifest with scalar math.
    Solver code remains responsible for constructing the actual CoilBuilder.
    """

    radius = coil["corner_centerline_radius_m"]
    straight_x = coil["straight_x_m"]
    straight_y = coil["straight_y_m"]
    width = coil["cross_section_width_m"]
    height = coil["cross_section_height_m"]
    centre = coil["centre_m"]
    frame = [list(row) for row in coil["orientation_rows"]]

    def add(*vectors: list[float]) -> list[float]:
        return [sum(vector[index] for vector in vectors) for index in range(3)]

    def scaled(scale: float, vector: list[float]) -> list[float]:
        return [scale * value for value in vector]

    def vector(values: list[float]) -> list[float]:
        return [_canonical_number(value, "coil identity coordinate") for value in values]

    position = add(
        centre,
        scaled(0.5 * straight_x + radius, frame[0]),
        scaled(-0.5 * straight_y, frame[1]),
    )
    first_position = list(position)
    segments: list[dict[str, Any]] = []
    operations = (
        ("straight", straight_y),
        ("arc", radius),
        ("straight", straight_x),
        ("arc", radius),
        ("straight", straight_y),
        ("arc", radius),
        ("straight", straight_x),
        ("arc", radius),
    )
    for index, (kind, value) in enumerate(operations):
        start = list(position)
        orientation = [list(row) for row in frame]
        if kind == "straight":
            position = add(position, scaled(value, frame[1]))
            entry = {
                "index": index,
                "kind": "straight",
                "start_m": vector(start),
                "end_m": vector(position),
                "orientation": [vector(row) for row in orientation],
                "cross_section_width_m": width,
                "cross_section_height_m": height,
                "length_m": value,
            }
        else:
            arc_center = add(position, scaled(-radius, frame[0]))
            # CoilBuilder.add_arc(..., 90 deg) rotates the local row frame.
            frame = [list(frame[1]), scaled(-1.0, frame[0]), list(frame[2])]
            position = add(arc_center, scaled(radius, frame[0]))
            entry = {
                "index": index,
                "kind": "arc",
                "start_m": vector(start),
                "end_m": vector(position),
                "orientation": [vector(row) for row in orientation],
                "cross_section_width_m": width,
                "cross_section_height_m": height,
                "corner_centerline_radius_m": radius,
                "arc_angle_deg": 90.0,
                "arc_center_m": vector(arc_center),
            }
        segments.append(entry)

    closure_error = math.sqrt(
        sum((end - start) ** 2 for start, end in zip(first_position, position))
    )
    payload = {
        "schema": "radia.coil-geometry-excitation-identity.v1",
        "coordinate_system": "right-handed Cartesian",
        "length_unit": "m",
        "current_unit": "A",
        "current_A": coil["current_A"],
        "turns": coil["turns"],
        "ampere_turns_A": _canonical_number(
            coil["current_A"] * coil["turns"], "ampere_turns_A"
        ),
        "current_direction": "segment_path_start_to_end",
        "conductor_model": coil["conductor_model"],
        "region_label": coil["region_label"],
        "closed": closure_error <= 1.0e-12,
        "closure_error_m": _canonical_number(closure_error, "closure_error_m"),
        "segments": segments,
    }
    return {**payload, "identity_sha256": _canonical_sha256(payload)}


def build_shared_pole_coil_exchange(payload: dict[str, Any]) -> dict[str, Any]:
    """Build CoilBuilder-compatible identities from one neutral exchange payload."""

    try:
        contract = normalize_shared_pole_coil_contract(payload)
        identity = _canonical_sha256(contract)
        claimed = payload.get("identity_sha256")
        if claimed is not None and str(claimed) != identity:
            raise ValueError("identity_sha256 does not match the canonical shared contract")

        adapters = []
        source_rows = []
        coil_geometry_rows = []
        for coil in contract["coils"]:
            coil_identity = _coilbuilder_racetrack_identity(coil)
            if not coil_identity["closed"]:
                raise ValueError(f"coil {coil['name']!r} is not closed")
            first = coil_identity["segments"][0]
            delta = [end - start for start, end in zip(first["start_m"], first["end_m"])]
            norm = math.sqrt(sum(value * value for value in delta))
            if norm <= 0.0:
                raise ValueError(f"coil {coil['name']!r} has a zero-length first segment")
            direction = [_canonical_number(value / norm, "source direction") for value in delta]
            adapters.append(
                {
                    "name": coil["name"],
                    "factory": "CoilBuilder.rounded_rectangle_racetrack",
                    "factory_arguments": {
                        "current": coil["current_A"],
                        "corner_centerline_radius": coil["corner_centerline_radius_m"],
                        "straight_x": coil["straight_x_m"],
                        "straight_y": coil["straight_y_m"],
                        "width": coil["cross_section_width_m"],
                        "height": coil["cross_section_height_m"],
                        "centre": coil["centre_m"],
                        "orientation": coil["orientation_rows"],
                    },
                    "coil_identity": coil_identity,
                }
            )
            source_rows.append(
                {
                    "name": coil["name"],
                    "kind": "coil",
                    "domain": coil["region_label"],
                    "n_turns": coil["turns"],
                    "current": coil["current_A"],
                    "direction": direction,
                    "magnitude": 0.0,
                }
            )
            coil_geometry_rows.append(
                {
                    key: value
                    for key, value in coil.items()
                    if key not in {
                        "current_A",
                        "turns",
                        "current_direction",
                        "conductor_model",
                        "current_density_A_per_m2",
                    }
                }
            )

        return {
            "status": "accepted",
            "accepted": True,
            "policy": "shared_pole_coil_exchange_v1",
            "identity_sha256": identity,
            "shared_contract": {**contract, "identity_sha256": identity},
            "radia_adapters": adapters,
            "modelir_fragment": {
                "name": contract["component_label"],
                "formulation": "magnetostatic_A",
                "length_unit": "m",
                "geometry": {
                    "schema": "modelir.labelled-pole-coil-geometry.v1",
                    "bodies": contract["bodies"],
                    "coil_paths": coil_geometry_rows,
                },
                "sources": source_rows,
            },
            "constitutive_policy": {
                "owner": "Radia/NGSolve",
                "changed_by_exchange": False,
            },
            "issues": [],
            "claim_boundary": (
                "This binds SI geometry, pole labels, coil path, frame, cross-section, "
                "and excitation. It does not claim mesh or field-result equivalence."
            ),
        }
    except (TypeError, ValueError) as exc:
        return {
            "status": "rejected",
            "accepted": False,
            "policy": "shared_pole_coil_exchange_v1",
            "issues": [str(exc)],
        }


def build_shared_pole_coil_exchange_json(payload_json: str) -> dict[str, Any]:
    """Parse JSON and build one shared pole/coil exchange contract."""

    try:
        payload = json.loads(payload_json)
    except (json.JSONDecodeError, TypeError) as exc:
        return {
            "status": "invalid_input",
            "accepted": False,
            "policy": "shared_pole_coil_exchange_v1",
            "issues": [str(exc)],
        }
    return build_shared_pole_coil_exchange(payload)
