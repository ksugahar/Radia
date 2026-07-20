"""Hex-quality, periodic, Sculpt, and Exodus identity checks for v51."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


_QUALITY = "hex_quality_metric_reference_jacobian_samples_order_mesh_owner_identity"
_PERIODIC = "periodic_highorder_node_edge_face_parametric_rotation_owner_identity"
_SCULPT = "parallel_sculpt_partition_ghost_seed_merge_revision_owner_identity"
_EXODUS = "exodus_int64_idmap_names_qa_time_global_order_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *names: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_vector(value: object, size: int) -> list[float] | None:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != size:
        return None
    try:
        values = [float(item) for item in value]
    except (TypeError, ValueError):
        return None
    return values if all(math.isfinite(item) for item in values) else None


def _quality_ok(row: Mapping[str, object]) -> bool:
    points = row.get("sample_points")
    parsed = (
        [_finite_vector(point, 3) for point in points]
        if isinstance(points, Sequence) and not isinstance(points, (str, bytes)) and points
        else []
    )
    order = row.get("element_order")
    return (
        _generations(row, "metric_generation", "jacobian_generation", "sample_generation", "order_generation", "owner_generation", "result_generation")
        and row.get("quality_metric") == "scaled_jacobian"
        and row.get("result_quality_metric") == row.get("quality_metric")
        and _digest(row.get("reference_jacobian_sha256"))
        and row.get("result_reference_jacobian_sha256") == row.get("reference_jacobian_sha256")
        and bool(parsed)
        and all(point is not None and all(-1.0 <= coordinate <= 1.0 for coordinate in point) for point in parsed)
        and len({tuple(point or ()) for point in parsed}) == len(parsed)
        and row.get("result_sample_points") == points
        and isinstance(order, int)
        and not isinstance(order, bool)
        and 1 <= order <= 3
        and row.get("result_element_order") == order
        and str(row.get("mesh_owner") or "").startswith("headless:")
        and row.get("result_mesh_owner") == row.get("mesh_owner")
        and _result(row)
    )


def _pair_list(value: object) -> list[list[int]] | None:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        return None
    pairs: list[list[int]] = []
    for pair in value:
        if not isinstance(pair, Sequence) or isinstance(pair, (str, bytes)) or len(pair) != 2:
            return None
        if not all(isinstance(item, int) and not isinstance(item, bool) and item > 0 for item in pair):
            return None
        pairs.append(list(pair))
    return pairs


def _parametric_map(value: object, dimension: int) -> bool:
    if not isinstance(value, Mapping) or not value:
        return False
    for name, coordinates in value.items():
        if not isinstance(name, str) or not name:
            return False
        if not isinstance(coordinates, Sequence) or isinstance(coordinates, (str, bytes)) or not coordinates:
            return False
        rows = coordinates if dimension > 1 else [[coordinate] for coordinate in coordinates]
        parsed = [_finite_vector(row, dimension) for row in rows]
        if any(row is None or any(not 0.0 <= item <= 1.0 for item in row) for row in parsed):
            return False
    return True


def _proper_rotation(value: object) -> bool:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 3:
        return False
    rows = [_finite_vector(row, 3) for row in value]
    if any(row is None for row in rows):
        return False
    matrix = rows  # type: ignore[assignment]
    dot = lambda left, right: sum(a * b for a, b in zip(left, right))
    orthogonal = all(math.isclose(dot(matrix[i], matrix[j]), 1.0 if i == j else 0.0, abs_tol=1e-12) for i in range(3) for j in range(3))
    determinant = (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )
    return orthogonal and math.isclose(determinant, 1.0, abs_tol=1e-12)


def _periodic_ok(row: Mapping[str, object]) -> bool:
    pairs = _pair_list(row.get("node_pairs"))
    return (
        _generations(row, "node_generation", "edge_generation", "face_generation", "rotation_generation", "owner_generation", "result_generation")
        and pairs is not None
        and len({pair[0] for pair in pairs}) == len(pairs)
        and len({pair[1] for pair in pairs}) == len(pairs)
        and row.get("result_node_pairs") == row.get("node_pairs")
        and _parametric_map(row.get("edge_parametric_coordinates"), 1)
        and row.get("result_edge_parametric_coordinates") == row.get("edge_parametric_coordinates")
        and _parametric_map(row.get("face_parametric_coordinates"), 2)
        and row.get("result_face_parametric_coordinates") == row.get("face_parametric_coordinates")
        and _proper_rotation(row.get("rotation_map"))
        and row.get("result_rotation_map") == row.get("rotation_map")
        and str(row.get("periodic_owner") or "").startswith("headless:")
        and row.get("result_periodic_owner") == row.get("periodic_owner")
        and _result(row)
    )


def _partition_map(value: object) -> dict[str, set[int]] | None:
    if not isinstance(value, Mapping) or not value:
        return None
    parsed: dict[str, set[int]] = {}
    for name, elements in value.items():
        if not isinstance(name, str) or not name.startswith("rank:"):
            return None
        if not isinstance(elements, Sequence) or isinstance(elements, (str, bytes)) or not elements:
            return None
        if not all(isinstance(item, int) and not isinstance(item, bool) and item > 0 for item in elements):
            return None
        parsed[name] = set(elements)
    all_ids = [item for values in parsed.values() for item in values]
    return parsed if len(all_ids) == len(set(all_ids)) else None


def _sculpt_ok(row: Mapping[str, object]) -> bool:
    partitions = _partition_map(row.get("partitions"))
    overlaps = _pair_list(row.get("ghost_overlap"))
    merge_order = row.get("merge_order")
    seed = row.get("random_seed")
    overlap_ok = partitions is not None and overlaps is not None and all(
        any(pair[0] in left and pair[1] in right for left_name, left in partitions.items() for right_name, right in partitions.items() if left_name != right_name)
        for pair in overlaps
    )
    return (
        _generations(row, "partition_generation", "ghost_generation", "seed_generation", "merge_generation", "revision_generation", "owner_generation", "result_generation")
        and partitions is not None
        and row.get("result_partitions") == row.get("partitions")
        and overlap_ok
        and row.get("result_ghost_overlap") == row.get("ghost_overlap")
        and isinstance(seed, int)
        and not isinstance(seed, bool)
        and seed >= 0
        and row.get("result_random_seed") == seed
        and isinstance(merge_order, Sequence)
        and not isinstance(merge_order, (str, bytes))
        and list(merge_order) == list(partitions)
        and row.get("result_merge_order") == merge_order
        and str(row.get("mesh_revision") or "").startswith("mesh-revision:")
        and row.get("result_mesh_revision") == row.get("mesh_revision")
        and str(row.get("job_owner") or "").startswith("headless:")
        and row.get("result_job_owner") == row.get("job_owner")
        and _result(row)
    )


def _exodus_ok(row: Mapping[str, object]) -> bool:
    ids = row.get("int64_ids")
    names = row.get("entity_names")
    qa = row.get("qa_records")
    times = row.get("time_values")
    globals_ = row.get("global_variable_order")
    try:
        parsed_times = [float(value) for value in times] if isinstance(times, Sequence) and not isinstance(times, (str, bytes)) else []
    except (TypeError, ValueError):
        parsed_times = []
    return (
        _generations(row, "id_generation", "name_generation", "qa_generation", "time_generation", "global_generation", "owner_generation", "result_generation")
        and isinstance(ids, Mapping)
        and bool(ids)
        and all(isinstance(name, str) and name and isinstance(value, int) and not isinstance(value, bool) and value > 2**31 for name, value in ids.items())
        and len(set(ids.values())) == len(ids)
        and row.get("result_int64_ids") == ids
        and isinstance(names, Mapping)
        and bool(names)
        and all(isinstance(name, str) and name and isinstance(value, str) and value for name, value in names.items())
        and row.get("result_entity_names") == names
        and isinstance(qa, Sequence)
        and not isinstance(qa, (str, bytes))
        and bool(qa)
        and all(isinstance(record, Sequence) and not isinstance(record, (str, bytes)) and len(record) == 4 and all(isinstance(item, str) and item for item in record) for record in qa)
        and row.get("result_qa_records") == qa
        and bool(parsed_times)
        and all(math.isfinite(value) and value >= 0.0 for value in parsed_times)
        and all(left < right for left, right in zip(parsed_times, parsed_times[1:]))
        and row.get("result_time_values") == times
        and isinstance(globals_, Sequence)
        and not isinstance(globals_, (str, bytes))
        and bool(globals_)
        and all(isinstance(name, str) and name for name in globals_)
        and len(set(globals_)) == len(globals_)
        and row.get("result_global_variable_order") == globals_
        and str(row.get("database_owner") or "").startswith("headless:")
        and row.get("result_database_owner") == row.get("database_owner")
        and _result(row)
    )


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {"policy": policy, "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, accepted in checks.items() if not accepted]}


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    quality = payload.get(_QUALITY)
    periodic = payload.get(_PERIODIC)
    if quality is not None:
        checks["v51_hex_quality_metric_jacobian_samples_order_owner"] = isinstance(quality, Mapping) and _quality_ok(quality)
    if periodic is not None:
        checks["v51_periodic_highorder_nodes_parametric_rotation_owner"] = isinstance(periodic, Mapping) and _periodic_ok(periodic)
    return _report("cubit_v51_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    sculpt = payload.get(_SCULPT)
    exodus = payload.get(_EXODUS)
    if sculpt is not None:
        checks["v51_parallel_sculpt_partition_ghost_seed_merge_revision_owner"] = isinstance(sculpt, Mapping) and _sculpt_ok(sculpt)
    if exodus is not None:
        checks["v51_exodus_int64_names_qa_time_global_order_owner"] = isinstance(exodus, Mapping) and _exodus_ok(exodus)
    return _report("cubit_v51_source_identity_v1", checks) if checks else {}
