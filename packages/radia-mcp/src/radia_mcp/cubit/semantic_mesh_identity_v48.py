"""Semantic identity checks for periodic, curved, replay, and Exodus artifacts."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _finite_sequence(value: object, *, length: int | None = None) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and (length is None or len(value) == length)
        and bool(value)
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def _result_ok(row: Mapping[str, object]) -> bool:
    owner = str(row.get("owner") or "")
    return (
        owner.startswith("headless:")
        and row.get("accepted_owner") == owner
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def _generations_ok(row: Mapping[str, object], fields: tuple[str, ...]) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _periodic_ok(row: Mapping[str, object]) -> bool:
    pairs = row.get("node_pairs")
    orientation = row.get("face_orientation")
    translation = row.get("translation")
    tolerance = row.get("pair_tolerance")
    error = row.get("maximum_pair_error")
    partition_owner = str(row.get("partition_owner") or "")
    valid_pairs = (
        isinstance(pairs, list)
        and bool(pairs)
        and all(
            isinstance(pair, list)
            and len(pair) == 2
            and all(isinstance(node, int) and not isinstance(node, bool) and node > 0 for node in pair)
            and pair[0] != pair[1]
            for pair in pairs
        )
        and len({pair[0] for pair in pairs}) == len(pairs)
        and len({pair[1] for pair in pairs}) == len(pairs)
    )
    return (
        _generations_ok(
            row,
            (
                "pair_generation",
                "orientation_generation",
                "translation_generation",
                "partition_generation",
                "result_generation",
            ),
        )
        and valid_pairs
        and row.get("result_node_pairs") == pairs
        and isinstance(orientation, list)
        and len(orientation) == 2
        and set(orientation) == {-1, 1}
        and row.get("result_face_orientation") == orientation
        and _finite_sequence(translation, length=3)
        and any(abs(float(value)) > 0.0 for value in translation)
        and row.get("result_translation") == translation
        and isinstance(tolerance, (int, float))
        and math.isfinite(float(tolerance))
        and float(tolerance) > 0.0
        and isinstance(error, (int, float))
        and math.isfinite(float(error))
        and 0.0 <= float(error) <= float(tolerance)
        and row.get("result_maximum_pair_error") == error
        and partition_owner.startswith("headless:")
        and row.get("result_partition_owner") == partition_owner
        and _result_ok(row)
    )


def _control_points(value: object) -> bool:
    if not isinstance(value, Mapping) or not value:
        return False
    identifiers: list[int] = []
    for entity, points in value.items():
        if not str(entity).strip() or not isinstance(points, list) or not points:
            return False
        if not all(isinstance(point, int) and not isinstance(point, bool) and point > 0 for point in points):
            return False
        identifiers.extend(points)
    return len(set(identifiers)) == len(identifiers)


def _curved_ok(row: Mapping[str, object]) -> bool:
    order = row.get("element_order")
    edges = row.get("edge_control_point_order")
    faces = row.get("face_control_point_order")
    error = row.get("maximum_projection_error")
    tolerance = row.get("projection_tolerance")
    return (
        _generations_ok(
            row,
            (
                "edge_generation",
                "face_generation",
                "cad_generation",
                "mesh_generation",
                "result_generation",
            ),
        )
        and isinstance(order, int)
        and not isinstance(order, bool)
        and order >= 2
        and row.get("result_element_order") == order
        and _control_points(edges)
        and row.get("result_edge_control_point_order") == edges
        and _control_points(faces)
        and row.get("result_face_control_point_order") == faces
        and _digest(row.get("cad_geometry_signature"))
        and row.get("result_cad_geometry_signature") == row.get("cad_geometry_signature")
        and isinstance(error, (int, float))
        and isinstance(tolerance, (int, float))
        and math.isfinite(float(error))
        and math.isfinite(float(tolerance))
        and 0.0 <= float(error) <= float(tolerance)
        and float(tolerance) > 0.0
        and row.get("result_maximum_projection_error") == error
        and _result_ok(row)
    )


def _scheme_ok(row: Mapping[str, object]) -> bool:
    fallback = row.get("scheme_fallback")
    sweep_map = row.get("sweep_source_target_map")
    count = row.get("fallback_count")
    valid_fallback = (
        isinstance(fallback, Mapping)
        and set(fallback) == {"requested", "applied", "reason"}
        and all(bool(str(fallback.get(field) or "").strip()) for field in fallback)
        and fallback.get("requested") != fallback.get("applied")
    )
    valid_map = (
        isinstance(sweep_map, Mapping)
        and bool(sweep_map)
        and all(
            bool(str(volume).strip())
            and isinstance(surfaces, list)
            and bool(surfaces)
            and len(set(surfaces)) == len(surfaces)
            and all(bool(str(surface).strip()) for surface in surfaces)
            for volume, surfaces in sweep_map.items()
        )
    )
    return (
        _generations_ok(
            row,
            (
                "scheme_generation",
                "fallback_generation",
                "sweep_map_generation",
                "volume_generation",
                "result_generation",
            ),
        )
        and valid_fallback
        and row.get("result_scheme_fallback") == fallback
        and valid_map
        and row.get("result_sweep_source_target_map") == sweep_map
        and isinstance(count, int)
        and not isinstance(count, bool)
        and count == 1
        and row.get("result_fallback_count") == count
        and _result_ok(row)
    )


def _exodus_ok(row: Mapping[str, object]) -> bool:
    times = row.get("timesteps")
    names = row.get("global_variable_names")
    rows = row.get("global_variable_rows")
    qa = row.get("qa_records")
    revision = str(row.get("mesh_revision") or "")
    valid_times = (
        _finite_sequence(times)
        and all(float(times[index]) < float(times[index + 1]) for index in range(len(times) - 1))
    )
    valid_names = (
        isinstance(names, list)
        and bool(names)
        and len(set(names)) == len(names)
        and all(bool(str(name).strip()) for name in names)
    )
    valid_rows = (
        isinstance(rows, list)
        and isinstance(times, Sequence)
        and len(rows) == len(times)
        and valid_names
        and all(_finite_sequence(values, length=len(names)) for values in rows)
    )
    valid_qa = (
        isinstance(qa, list)
        and bool(qa)
        and all(
            isinstance(record, list)
            and len(record) == 4
            and all(bool(str(value).strip()) for value in record)
            for record in qa
        )
    )
    return (
        _generations_ok(
            row,
            (
                "timestep_generation",
                "global_variable_generation",
                "qa_generation",
                "mesh_generation",
                "result_generation",
            ),
        )
        and valid_times
        and row.get("result_timesteps") == times
        and valid_names
        and row.get("result_global_variable_names") == names
        and valid_rows
        and row.get("result_global_variable_rows") == rows
        and valid_qa
        and row.get("result_qa_records") == qa
        and bool(revision)
        and row.get("result_mesh_revision") == revision
        and _result_ok(row)
    )


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {
        "policy": policy,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
    }


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    periodic = payload.get("periodic_node_pair_orientation_translation_tolerance_partition_owner_identity")
    curved = payload.get("high_order_curved_edge_face_control_point_cad_signature_identity")
    if periodic is not None:
        checks["v48_periodic_node_pair_transform_partition_owner"] = (
            isinstance(periodic, Mapping) and _periodic_ok(periodic)
        )
    if curved is not None:
        checks["v48_curved_control_point_cad_signature"] = isinstance(curved, Mapping) and _curved_ok(curved)
    return _report("cubit_v48_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    scheme = payload.get("mesh_scheme_fallback_provenance_sweep_source_target_map_owner_identity")
    exodus = payload.get("exodus_timestep_global_variable_qa_record_mesh_revision_identity")
    if scheme is not None:
        checks["v48_mesh_scheme_fallback_sweep_map_owner"] = isinstance(scheme, Mapping) and _scheme_ok(scheme)
    if exodus is not None:
        checks["v48_exodus_timestep_global_qa_mesh_revision"] = isinstance(exodus, Mapping) and _exodus_ok(exodus)
    return _report("cubit_v48_source_identity_v1", checks) if checks else {}
