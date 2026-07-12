"""Solver-neutral evidence gate for heterogeneous part-mesh replay drift."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


def _finite(value: object, name: str, *, positive: bool = False) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{name} must be finite")
    if positive and parsed <= 0.0:
        raise ValueError(f"{name} must be positive")
    return parsed


def _count(value: object, name: str) -> int:
    parsed = _finite(value, name, positive=True)
    integer = int(parsed)
    if parsed != integer:
        raise ValueError(f"{name} must be an integer")
    return integer


def _rows(value: object, name: str, minimum: int) -> list[Mapping[str, object]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    rows = list(value)
    if len(rows) < minimum or any(not isinstance(row, Mapping) for row in rows):
        raise ValueError(f"{name} must contain at least {minimum} objects")
    return rows


def _relative(actual: float, reference: float) -> float:
    return abs(actual - reference) / max(abs(reference), 1.0e-300)


def heterogeneous_part_mesh_replay_gate(
    summary: Mapping[str, object],
) -> dict[str, Any]:
    """Classify deterministic remesh drift without promoting it to solver-ready."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    references = _rows(summary.get("reference_evidence"), "reference_evidence", 2)
    routes = _rows(summary.get("mesh_routes"), "mesh_routes", 3)
    replays = _rows(summary.get("live_replays"), "live_replays", 2)
    warning = summary.get("warning")
    if not isinstance(warning, Mapping):
        raise ValueError("warning must be an object")
    maximum_drift = _finite(
        summary.get("maximum_reference_live_relative_error"),
        "maximum_reference_live_relative_error",
        positive=True,
    )

    reference_ids = [str(row.get("evidence_id") or "").strip() for row in references]
    reference_counts = [
        (_count(row.get("nodes"), f"reference_evidence[{index}].nodes"),
         _count(row.get("elements"), f"reference_evidence[{index}].elements"))
        for index, row in enumerate(references)
    ]
    replay_ids = [str(row.get("replay_id") or "").strip() for row in replays]
    replay_counts = [
        (_count(row.get("nodes"), f"live_replays[{index}].nodes"),
         _count(row.get("elements"), f"live_replays[{index}].elements"))
        for index, row in enumerate(replays)
    ]
    route_names = [str(row.get("route") or "").strip() for row in routes]
    reference_nodes, reference_elements = reference_counts[0]
    live_nodes, live_elements = replay_counts[0]
    node_drift = _relative(live_nodes, reference_nodes)
    element_drift = _relative(live_elements, reference_elements)
    replay_node_span = max(
        _relative(nodes, live_nodes) for nodes, _ in replay_counts
    )
    replay_element_span = max(
        _relative(elements, live_elements) for _, elements in replay_counts
    )
    changed_artifacts_ok = all(
        isinstance(row.get("changed_artifacts"), Sequence)
        and not isinstance(row.get("changed_artifacts"), (str, bytes))
        and len(row.get("changed_artifacts")) > 0
        for row in replays
    )
    checks = {
        "independent_reference_records_are_identified": len(set(reference_ids))
        == len(reference_ids)
        and all(reference_ids)
        and all(row.get("independent") is True for row in references),
        "independent_reference_mesh_counts_agree": all(
            count == reference_counts[0] for count in reference_counts[1:]
        ),
        "heterogeneous_mesh_routes_are_explicit": len(set(route_names))
        == len(route_names)
        and all(route_names)
        and all(row.get("observed") is True for row in routes)
        and "external_part_mesh" in route_names,
        "source_warning_is_classified_not_ignored": _count(
            warning.get("count"), "warning.count"
        )
        >= 1
        and bool(str(warning.get("code") or "").strip())
        and warning.get("observed_in_report") is True
        and warning.get("disposition") == "source_configuration_requires_review",
        "two_or_more_fresh_replays_are_identified": len(set(replay_ids)) == len(replay_ids)
        and all(replay_ids),
        "fresh_replay_mesh_counts_are_deterministic": replay_node_span <= 1.0e-12
        and replay_element_span <= 1.0e-12,
        "fresh_replays_preserve_source_and_process_ownership": all(
            row.get("source_preserved") is True
            and row.get("temporary_work_copy") is True
            and row.get("pass_marker") is True
            and int(row.get("owned_processes_after", -1)) == 0
            for row in replays
        ),
        "fresh_replays_create_mesh_without_physics_result": all(
            row.get("has_mesh_before") is False
            and row.get("has_mesh_after") is True
            and row.get("has_mesh_any_part_after") is True
            and row.get("has_result_after") is False
            for row in replays
        ),
        "fresh_replays_record_changed_mesh_artifacts": changed_artifacts_ok,
        "reference_to_live_mesh_drift_is_material": max(node_drift, element_drift)
        > maximum_drift,
        "material_drift_is_not_promoted_to_solver_ready": summary.get("classification")
        == "reproducible_remesh_drift_not_solver_ready"
        and summary.get("solver_ready") is False,
    }
    return {
        "policy": "heterogeneous_part_mesh_replay_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "classification": (
            "reproducible_remesh_drift_not_solver_ready"
            if all(checks.values())
            else "incomplete_or_contradictory_mesh_evidence"
        ),
        "solver_ready": False,
        "metrics": {
            "reference_nodes": reference_nodes,
            "reference_elements": reference_elements,
            "live_nodes": live_nodes,
            "live_elements": live_elements,
            "live_node_retention": live_nodes / reference_nodes,
            "live_element_retention": live_elements / reference_elements,
            "node_relative_drift": node_drift,
            "element_relative_drift": element_drift,
            "replay_node_relative_span": replay_node_span,
            "replay_element_relative_span": replay_element_span,
        },
        "lesson": (
            "A deterministic remesh is not equivalent to the archived heterogeneous mesh. "
            "Require two independent reference records, repeated fresh generation, explicit route and warning provenance, "
            "and reject solver handoff when the counts materially drift instead of hiding the change with a loose tolerance."
        ),
    }
