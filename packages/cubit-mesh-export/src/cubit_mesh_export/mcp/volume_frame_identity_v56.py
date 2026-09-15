"""Headless hex-volume, curved-face, batch, and Exodus-frame checks."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


VOLUME = "hexmesh_blockvolume_cadvolume_orientation_owner_identity"
CURVED = "curvedhex_faceconformity_nodeorder_jacobian_geometry_owner_identity"
BATCH = "batchjournal_errorstatus_rollback_outputdatabase_owner_identity"
EXODUS = "exodus_coordinateframe_qarecord_timestep_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _number(value: object, *, positive: bool = False, nonnegative: bool = False) -> bool:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return False
    number = float(value)
    return math.isfinite(number) and (not positive or number > 0.0) and (not nonnegative or number >= 0.0)


def _close(left: object, right: object) -> bool:
    return _number(left) and _number(right) and math.isclose(float(left), float(right), rel_tol=1.0e-10, abs_tol=1.0e-12)


def _generation(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _face_nodes_conform(side_a: object, side_b: object) -> bool:
    sides = (side_a, side_b)
    if not all(
        isinstance(side, Sequence)
        and not isinstance(side, (str, bytes))
        and len(side) == 9
        and all(isinstance(node, int) and not isinstance(node, bool) and node > 0 for node in side)
        for side in sides
    ):
        return False
    a = list(side_a)
    b = list(side_b)
    if len(set(a)) != 9 or len(set(b)) != 9:
        return False
    a_corners, b_corners = a[:4], b[:4]
    a_edges, b_edges = a[4:8], b[4:8]
    if set(a_corners) != set(b_corners) or set(a_edges) != set(b_edges) or a[8] != b[8]:
        return False
    edge_by_vertices = {
        frozenset((a_corners[index], a_corners[(index + 1) % 4])): a_edges[index]
        for index in range(4)
    }
    return len(edge_by_vertices) == 4 and all(
        edge_by_vertices.get(frozenset((b_corners[index], b_corners[(index + 1) % 4]))) == b_edges[index]
        for index in range(4)
    )


def _volume_ok(row: Mapping[str, object]) -> bool:
    elements = row.get("signed_element_volume_m3")
    block_sums = row.get("block_volume_sum_m3")
    if not isinstance(elements, Mapping) or not elements or not isinstance(block_sums, Mapping) or set(block_sums) != set(elements):
        return False
    computed: dict[str, float] = {}
    for block, values in elements.items():
        if not isinstance(block, str) or not block.startswith("block:") or not isinstance(values, Sequence) or isinstance(values, (str, bytes)) or not values or not all(_number(value, positive=True) for value in values):
            return False
        computed[block] = sum(float(value) for value in values)
    sums_ok = all(_number(block_sums[block], positive=True) and _close(block_sums[block], computed[block]) for block in computed)
    cad_volume = row.get("cad_volume_m3")
    return (
        _generation(row, "element_generation", "block_generation", "cad_generation", "orientation_generation", "owner_generation", "result_generation")
        and sums_ok and _number(cad_volume, positive=True) and _close(cad_volume, sum(computed.values()))
        and row.get("orientation") == "positive"
        and all(row.get("result_" + field) == row.get(field) for field in ("signed_element_volume_m3", "block_volume_sum_m3", "cad_volume_m3", "orientation", "mesh_owner"))
        and str(row.get("mesh_owner") or "").startswith("headless:") and _result(row)
    )


def _curved_ok(row: Mapping[str, object]) -> bool:
    faces = row.get("shared_face_nodes")
    jacobians = row.get("high_order_jacobian_samples")
    faces_ok = isinstance(faces, Sequence) and not isinstance(faces, (str, bytes)) and bool(faces)
    names: list[str] = []
    if faces_ok:
        for face in faces:
            if not isinstance(face, Mapping) or set(face) != {"face", "side_a_nodes", "side_b_nodes"}:
                faces_ok = False
                break
            side_a = face["side_a_nodes"]
            side_b = face["side_b_nodes"]
            if not isinstance(face["face"], str) or not face["face"].startswith("face:") or not _face_nodes_conform(side_a, side_b):
                faces_ok = False
                break
            names.append(face["face"])
        faces_ok = faces_ok and len(names) == len(set(names))
    return (
        _generation(row, "face_generation", "order_generation", "jacobian_generation", "geometry_generation", "owner_generation", "result_generation")
        and faces_ok and row.get("element_order") == "HEX27"
        and isinstance(jacobians, Sequence) and not isinstance(jacobians, (str, bytes)) and bool(jacobians) and all(_number(value, positive=True) for value in jacobians)
        and _digest(row.get("geometry_revision_sha256"))
        and all(row.get("result_" + field) == row.get(field) for field in ("shared_face_nodes", "element_order", "high_order_jacobian_samples", "geometry_revision_sha256", "mesh_owner"))
        and str(row.get("mesh_owner") or "").startswith("headless:") and _result(row)
    )


def _batch_ok(row: Mapping[str, object]) -> bool:
    return (
        _generation(row, "status_generation", "rollback_generation", "database_generation", "owner_generation", "result_generation")
        and row.get("error_status") == "none"
        and isinstance(row.get("rollback_applied"), bool) and not row.get("rollback_applied")
        and str(row.get("rollback_checkpoint") or "").startswith("checkpoint:")
        and str(row.get("output_database_revision") or "").startswith("database:post")
        and all(row.get("replayed_" + field) == row.get(field) for field in ("error_status", "rollback_checkpoint", "rollback_applied", "output_database_revision", "session_owner"))
        and str(row.get("session_owner") or "").startswith("headless:") and _result(row)
    )


def _matrix4(value: object) -> list[list[float]] | None:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 4:
        return None
    matrix: list[list[float]] = []
    for line in value:
        if not isinstance(line, Sequence) or isinstance(line, (str, bytes)) or len(line) != 4 or not all(_number(item) for item in line):
            return None
        matrix.append([float(item) for item in line])
    return matrix


def _frame_ok(matrix: list[list[float]] | None) -> bool:
    if matrix is None or matrix[3] != [0.0, 0.0, 0.0, 1.0]:
        return False
    rotation = [line[:3] for line in matrix[:3]]
    orthonormal = all(math.isclose(sum(rotation[row][index] * rotation[column][index] for index in range(3)), 1.0 if row == column else 0.0, abs_tol=1.0e-12) for row in range(3) for column in range(3))
    determinant = (
        rotation[0][0] * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
        - rotation[0][1] * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
        + rotation[0][2] * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0])
    )
    return orthonormal and math.isclose(determinant, 1.0, abs_tol=1.0e-12)


def _exodus_ok(row: Mapping[str, object]) -> bool:
    transform = row.get("coordinate_transform_4x4")
    qa = row.get("qa_records")
    steps = row.get("time_steps")
    qa_ok = isinstance(qa, Sequence) and not isinstance(qa, (str, bytes)) and bool(qa) and all(isinstance(record, Mapping) and set(record) == {"program", "version", "date", "time"} and all(isinstance(record[field], str) and record[field] for field in record) for record in qa)
    steps_ok = isinstance(steps, Sequence) and not isinstance(steps, (str, bytes)) and len(steps) >= 2
    times: list[float] = []
    if steps_ok:
        for index, step in enumerate(steps):
            step_index = step.get("index") if isinstance(step, Mapping) else None
            if not isinstance(step, Mapping) or set(step) != {"index", "time_s"} or not isinstance(step_index, int) or isinstance(step_index, bool) or step_index != index or not _number(step.get("time_s"), nonnegative=True):
                steps_ok = False
                break
            times.append(float(step["time_s"]))
        steps_ok = steps_ok and times[0] == 0.0 and all(left < right for left, right in zip(times, times[1:]))
    return (
        _generation(row, "frame_generation", "qa_generation", "timestep_generation", "owner_generation", "result_generation")
        and _frame_ok(_matrix4(transform)) and qa_ok and steps_ok
        and all(row.get("replayed_" + field) == row.get(field) for field in ("coordinate_transform_4x4", "qa_records", "time_steps", "file_owner"))
        and str(row.get("file_owner") or "").startswith("headless:") and _result(row)
    )


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {"policy": policy, "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, accepted in checks.items() if not accepted]}


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    if payload.get(VOLUME) is not None:
        checks["v56_hex_block_cad_volume_orientation_owner"] = isinstance(payload[VOLUME], Mapping) and _volume_ok(payload[VOLUME])
    if payload.get(CURVED) is not None:
        checks["v56_curved_hex_face_order_jacobian_geometry_owner"] = isinstance(payload[CURVED], Mapping) and _curved_ok(payload[CURVED])
    return _report("cubit_v56_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    if payload.get(BATCH) is not None:
        checks["v56_batch_status_rollback_database_owner"] = isinstance(payload[BATCH], Mapping) and _batch_ok(payload[BATCH])
    if payload.get(EXODUS) is not None:
        checks["v56_exodus_frame_qa_timestep_owner"] = isinstance(payload[EXODUS], Mapping) and _exodus_ok(payload[EXODUS])
    return _report("cubit_v56_source_identity_v1", checks) if checks else {}
