"""Gear, pipe-sweep, STEP-occurrence, and mesh-exchange identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


GEAR = "involutegear_module_toothcount_pressureangle_pitchdiameter_volume_owner_identity"
PIPE = "pipesweep_pathlength_frame_twist_selfintersection_volume_owner_identity"
STEP = "step_occurrence_transform_unit_product_assemblyframe_owner_identity"
MESH = "meshformat_watertight_manifold_unit_signedvolume_owner_identity"
_UNIT_SCALE = {"m": 1.0, "mm": 1.0e-3, "cm": 1.0e-2, "inch": 0.0254}


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _number(value: object, *, positive: bool = False) -> bool:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return False
    number = float(value)
    return math.isfinite(number) and (not positive or number > 0.0)


def _close(left: object, right: object) -> bool:
    return _number(left) and _number(right) and math.isclose(
        float(left), float(right), rel_tol=1.0e-10, abs_tol=1.0e-12
    )


def _generation(row: Mapping[str, object], *names: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get(
        "accepted_result_sha256"
    ) == row.get("result_sha256")


def _gear_ok(row: Mapping[str, object]) -> bool:
    module = row.get("module_m")
    teeth = row.get("tooth_count")
    pressure = row.get("pressure_angle_rad")
    pitch = row.get("pitch_diameter_m")
    names = (
        "module_m",
        "tooth_count",
        "pressure_angle_rad",
        "pitch_diameter_m",
        "volume_m3",
        "shape_owner",
    )
    return (
        _generation(
            row,
            "module_generation",
            "tooth_generation",
            "pressure_generation",
            "pitch_generation",
            "volume_generation",
            "owner_generation",
            "result_generation",
        )
        and _number(module, positive=True)
        and isinstance(teeth, int)
        and not isinstance(teeth, bool)
        and teeth >= 6
        and _number(pressure, positive=True)
        and float(pressure) < math.pi / 2.0
        and _number(pitch, positive=True)
        and _close(pitch, float(module) * teeth)
        and _number(row.get("volume_m3"), positive=True)
        and str(row.get("shape_owner") or "").startswith("shape:")
        and all(row.get("result_" + name) == row.get(name) for name in names)
        and _result(row)
    )


def _frame_ok(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and set(value) == {"method", "samples", "closed"}
        and value.get("method") in {"parallel_transport", "frenet"}
        and isinstance(value.get("samples"), int)
        and not isinstance(value.get("samples"), bool)
        and int(value["samples"]) >= 3
        and isinstance(value.get("closed"), bool)
    )


def _pipe_ok(row: Mapping[str, object]) -> bool:
    names = (
        "path_length_m",
        "transport_frame",
        "net_twist_rad",
        "self_intersection",
        "volume_m3",
        "shape_owner",
    )
    return (
        _generation(
            row,
            "path_generation",
            "frame_generation",
            "twist_generation",
            "intersection_generation",
            "volume_generation",
            "owner_generation",
            "result_generation",
        )
        and _number(row.get("path_length_m"), positive=True)
        and _frame_ok(row.get("transport_frame"))
        and _number(row.get("net_twist_rad"))
        and isinstance(row.get("self_intersection"), bool)
        and row.get("self_intersection") is False
        and _number(row.get("volume_m3"), positive=True)
        and str(row.get("shape_owner") or "").startswith("shape:")
        and all(row.get("result_" + name) == row.get(name) for name in names)
        and _result(row)
    )


def _matrix4(value: object) -> list[list[float]] | None:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 4:
        return None
    result: list[list[float]] = []
    for row in value:
        if (
            not isinstance(row, Sequence)
            or isinstance(row, (str, bytes))
            or len(row) != 4
            or not all(_number(item) for item in row)
        ):
            return None
        result.append([float(item) for item in row])
    return result


def _rigid_transform(value: object) -> bool:
    matrix = _matrix4(value)
    if matrix is None or matrix[3] != [0.0, 0.0, 0.0, 1.0]:
        return False
    rotation = [row[:3] for row in matrix[:3]]
    orthonormal = all(
        math.isclose(
            sum(rotation[row][index] * rotation[column][index] for index in range(3)),
            1.0 if row == column else 0.0,
            abs_tol=1.0e-12,
        )
        for row in range(3)
        for column in range(3)
    )
    determinant = (
        rotation[0][0] * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
        - rotation[0][1] * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
        + rotation[0][2] * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0])
    )
    return orthonormal and math.isclose(determinant, 1.0, abs_tol=1.0e-12)


def _step_ok(row: Mapping[str, object]) -> bool:
    occurrence = row.get("product_occurrence")
    frame = row.get("assembly_frame")
    unit = row.get("length_unit")
    semantic_ok = (
        isinstance(occurrence, Mapping)
        and set(occurrence) == {"id", "product"}
        and str(occurrence.get("id") or "").startswith("occurrence:")
        and str(occurrence.get("product") or "").startswith("product:")
        and isinstance(frame, Mapping)
        and set(frame) == {"parent", "child"}
        and str(frame.get("parent") or "").startswith("assembly:")
        and frame.get("child") == occurrence.get("id")
    )
    names = (
        "occurrence_transform_4x4",
        "length_unit",
        "product_occurrence",
        "assembly_frame",
        "document_owner",
    )
    return (
        _generation(
            row,
            "transform_generation",
            "unit_generation",
            "product_generation",
            "frame_generation",
            "owner_generation",
            "result_generation",
        )
        and _rigid_transform(row.get("occurrence_transform_4x4"))
        and isinstance(unit, str)
        and unit in _UNIT_SCALE
        and semantic_ok
        and str(row.get("document_owner") or "").startswith("document:")
        and all(row.get("replayed_" + name) == row.get(name) for name in names)
        and _result(row)
    )


def _mesh_ok(row: Mapping[str, object]) -> bool:
    unit = row.get("length_unit")
    mesh_format = row.get("mesh_format")
    names = (
        "mesh_format",
        "watertight",
        "manifold",
        "length_unit",
        "unit_scale_to_m",
        "signed_volume_m3",
        "mesh_owner",
    )
    return (
        _generation(
            row,
            "format_generation",
            "watertight_generation",
            "manifold_generation",
            "unit_generation",
            "volume_generation",
            "owner_generation",
            "result_generation",
        )
        and isinstance(mesh_format, str)
        and mesh_format in {"stl", "3mf"}
        and row.get("watertight") is True
        and row.get("manifold") is True
        and isinstance(unit, str)
        and unit in _UNIT_SCALE
        and _close(row.get("unit_scale_to_m"), _UNIT_SCALE[unit])
        and _number(row.get("signed_volume_m3"), positive=True)
        and str(row.get("mesh_owner") or "").startswith("mesh:")
        and all(row.get("replayed_" + name) == row.get(name) for name in names)
        and _result(row)
    )


def _public_rows(payload: Mapping[str, object]) -> list[Mapping[str, object]] | None:
    rows: list[Mapping[str, object]] = []
    reference = payload.get("reference")
    if reference is not None:
        if not isinstance(reference, Sequence) or isinstance(reference, (str, bytes)):
            return None
        for item in reference:
            if not isinstance(item, Mapping):
                return None
            rows.append(item)
    measured = payload.get("measured")
    if measured is not None:
        if not isinstance(measured, Mapping):
            return None
        for family in measured.values():
            if not isinstance(family, Sequence) or isinstance(family, (str, bytes)):
                return None
            for item in family:
                if not isinstance(item, Mapping):
                    return None
                rows.append(item)
    return rows


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {
        "policy": policy,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
    }


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    rows = _public_rows(payload)
    if rows is None:
        return _report("build123d_v56_public_identity_v1", {"v56_public_payload_rows_are_mappings": False})
    checks: dict[str, bool] = {}
    gears = [row.get(GEAR) for row in rows if GEAR in row]
    pipes = [row.get(PIPE) for row in rows if PIPE in row]
    if gears:
        checks["v56_involute_gear_pitch_volume_owner"] = len(gears) == len(rows) and all(
            isinstance(item, Mapping) and _gear_ok(item) for item in gears
        )
    if pipes:
        checks["v56_pipe_sweep_frame_intersection_volume_owner"] = len(pipes) == len(rows) and all(
            isinstance(item, Mapping) and _pipe_ok(item) for item in pipes
        )
    return _report("build123d_v56_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping) or not isinstance(payload.get("replay_identity"), Mapping):
        return {}
    identity = payload["replay_identity"]
    checks: dict[str, bool] = {}
    if identity.get(STEP) is not None:
        checks["v56_step_occurrence_transform_unit_frame_owner"] = isinstance(
            identity[STEP], Mapping
        ) and _step_ok(identity[STEP])
    if identity.get(MESH) is not None:
        checks["v56_mesh_format_manifold_unit_volume_owner"] = isinstance(
            identity[MESH], Mapping
        ) and _mesh_ok(identity[MESH])
    return _report("build123d_v56_source_identity_v1", checks) if checks else {}
