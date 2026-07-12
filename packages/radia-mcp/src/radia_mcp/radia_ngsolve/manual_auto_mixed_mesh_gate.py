"""Solver-neutral preservation gate for manual-plus-automatic 2D meshes."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


def _count(value: object, name: str) -> int:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed) or parsed <= 0.0 or int(parsed) != parsed:
        raise ValueError(f"{name} must be a positive integer")
    return int(parsed)


def _finite(value: object, name: str, *, nonnegative: bool = False) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed) or (nonnegative and parsed < 0.0):
        raise ValueError(f"{name} must be finite")
    return parsed


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    return value


def _rows(value: object, name: str) -> list[Mapping[str, object]]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be an array")
    rows = list(value)
    if len(rows) < 2 or any(not isinstance(row, Mapping) for row in rows):
        raise ValueError(f"{name} must contain at least two objects")
    return rows


def _relative(actual: float, reference: float) -> float:
    return abs(actual - reference) / max(abs(reference), 1.0e-300)


def _parse_mesh(
    row: Mapping[str, object], name: str, *, require_vertices: bool = True
) -> dict[str, object]:
    regions = _mapping(row.get("region_elements"), f"{name}.region_elements")
    sizes = _mapping(row.get("auto_mesh_size_m"), f"{name}.auto_mesh_size_m")
    return {
        "total_elements": _count(row.get("total_elements"), f"{name}.total_elements"),
        "vertices": (
            _count(row.get("vertices"), f"{name}.vertices")
            if require_vertices
            else None
        ),
        "manual_region_elements": _count(
            regions.get("manual_region"), f"{name}.region_elements.manual_region"
        ),
        "automatic_region_elements": _count(
            regions.get("automatic_region"),
            f"{name}.region_elements.automatic_region",
        ),
        "manual_region_auto_size": _finite(
            sizes.get("manual_region"),
            f"{name}.auto_mesh_size_m.manual_region",
            nonnegative=True,
        ),
        "automatic_region_auto_size": _finite(
            sizes.get("automatic_region"),
            f"{name}.auto_mesh_size_m.automatic_region",
            nonnegative=True,
        ),
        "keep_existing_mesh": row.get("keep_existing_mesh") is True,
        "partial_mesh_before": row.get("partial_mesh_before") is True,
        "physics_result_after": row.get("physics_result_after") is True,
        "source_preserved": row.get("source_preserved") is True,
        "temporary_work_copy": row.get("temporary_work_copy") is True,
        "pass_marker": row.get("pass_marker") is True,
        "owned_processes_after": int(row.get("owned_processes_after", -1)),
    }


def manual_auto_mixed_mesh_preservation_gate(
    summary: Mapping[str, object],
    *,
    maximum_automatic_region_relative_drift: float = 0.10,
    maximum_total_element_relative_drift: float = 0.10,
) -> dict[str, Any]:
    """Gate exact manual-region preservation and bounded automatic remeshing.

    A manual element region must survive completion exactly.  The automatic
    region may change across mesher versions, but its drift must be explicit,
    bounded, and deterministic across two fresh replays.  This evidence is a
    mesh-preservation diagnosis; it is not solver-ready proof without an
    element-family/interface audit.
    """
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    contract = _mapping(summary.get("model_contract"), "model_contract")
    archived = _mapping(summary.get("archived_reference"), "archived_reference")
    fresh = _mapping(summary.get("fresh_replays"), "fresh_replays")
    archived_auto = _parse_mesh(
        _mapping(archived.get("automatic_only"), "archived_reference.automatic_only"),
        "archived_reference.automatic_only",
        require_vertices=False,
    )
    archived_mixed = _parse_mesh(
        _mapping(archived.get("manual_plus_automatic"), "archived_reference.manual_plus_automatic"),
        "archived_reference.manual_plus_automatic",
        require_vertices=False,
    )
    fresh_auto = [
        _parse_mesh(row, f"fresh_replays.automatic_only[{index}]")
        for index, row in enumerate(_rows(fresh.get("automatic_only"), "fresh_replays.automatic_only"))
    ]
    fresh_mixed = [
        _parse_mesh(row, f"fresh_replays.manual_plus_automatic[{index}]")
        for index, row in enumerate(
            _rows(fresh.get("manual_plus_automatic"), "fresh_replays.manual_plus_automatic")
        )
    ]
    automatic_limit = _finite(
        maximum_automatic_region_relative_drift,
        "maximum_automatic_region_relative_drift",
        nonnegative=True,
    )
    total_limit = _finite(
        maximum_total_element_relative_drift,
        "maximum_total_element_relative_drift",
        nonnegative=True,
    )

    auto_replay = fresh_auto[0]
    mixed_replay = fresh_mixed[0]
    auto_region_drift = _relative(
        int(auto_replay["automatic_region_elements"]),
        int(archived_auto["automatic_region_elements"]),
    )
    mixed_auto_region_drift = _relative(
        int(mixed_replay["automatic_region_elements"]),
        int(archived_mixed["automatic_region_elements"]),
    )
    auto_total_drift = _relative(
        int(auto_replay["total_elements"]), int(archived_auto["total_elements"])
    )
    mixed_total_drift = _relative(
        int(mixed_replay["total_elements"]), int(archived_mixed["total_elements"])
    )
    metrics = {
        "manual_region_reference_elements": archived_mixed["manual_region_elements"],
        "manual_region_live_elements": mixed_replay["manual_region_elements"],
        "automatic_only_region_relative_drift": auto_region_drift,
        "mixed_automatic_region_relative_drift": mixed_auto_region_drift,
        "automatic_only_total_relative_drift": auto_total_drift,
        "mixed_total_relative_drift": mixed_total_drift,
        "automatic_only_live_elements": auto_replay["total_elements"],
        "mixed_live_elements": mixed_replay["total_elements"],
    }
    checks = {
        "two_dimensional_manual_auto_element_contract": contract
        == {
            "dimension": 2,
            "same_two_part_geometry": True,
            "manual_region_element_family": "quadrilateral",
            "automatic_region_element_family": "triangle",
            "manual_automatic_interface": "conformal",
        },
        "two_deterministic_automatic_replays": all(
            row["total_elements"] == auto_replay["total_elements"]
            and row["vertices"] == auto_replay["vertices"]
            and row["manual_region_elements"] == auto_replay["manual_region_elements"]
            and row["automatic_region_elements"] == auto_replay["automatic_region_elements"]
            for row in fresh_auto[1:]
        ),
        "two_deterministic_mixed_replays": all(
            row["total_elements"] == mixed_replay["total_elements"]
            and row["vertices"] == mixed_replay["vertices"]
            and row["manual_region_elements"] == mixed_replay["manual_region_elements"]
            and row["automatic_region_elements"] == mixed_replay["automatic_region_elements"]
            for row in fresh_mixed[1:]
        ),
        "manual_region_is_preserved_exactly": archived_mixed[
            "manual_region_elements"
        ]
        == mixed_replay["manual_region_elements"]
        and all(
            row["manual_region_elements"] == archived_mixed["manual_region_elements"]
            for row in fresh_mixed
        ),
        "automatic_mesher_skips_manual_region": archived_mixed[
            "manual_region_auto_size"
        ]
        == 0.0
        and all(row["manual_region_auto_size"] == 0.0 for row in fresh_mixed)
        and archived_auto["manual_region_auto_size"] > 0.0
        and all(row["manual_region_auto_size"] > 0.0 for row in fresh_auto),
        "automatic_region_mesher_remains_active": archived_mixed[
            "automatic_region_auto_size"
        ]
        > 0.0
        and all(row["automatic_region_auto_size"] > 0.0 for row in fresh_mixed),
        "keep_existing_mesh_and_partial_mesh_are_explicit": all(
            row["keep_existing_mesh"] and row["partial_mesh_before"]
            for row in fresh_mixed
        )
        and all(not row["keep_existing_mesh"] for row in fresh_auto),
        "automatic_region_version_drift_is_bounded": max(
            auto_region_drift, mixed_auto_region_drift
        )
        <= automatic_limit,
        "total_element_version_drift_is_bounded": max(
            auto_total_drift, mixed_total_drift
        )
        <= total_limit,
        "automatic_and_mixed_routes_are_distinct": auto_replay["total_elements"]
        != mixed_replay["total_elements"]
        and auto_replay["manual_region_elements"]
        != mixed_replay["manual_region_elements"],
        "fresh_replays_preserve_source_and_process_ownership": all(
            row["source_preserved"]
            and row["temporary_work_copy"]
            and row["pass_marker"]
            and row["owned_processes_after"] == 0
            for row in fresh_auto + fresh_mixed
        ),
        "mesh_only_evidence_not_physics_result": all(
            not row["physics_result_after"] for row in fresh_auto + fresh_mixed
        ),
        "not_promoted_without_element_interface_audit": summary.get("solver_ready")
        is False
        and summary.get("classification")
        == "manual_region_preserved_auto_region_version_drift_recorded",
    }
    return {
        "policy": "manual_auto_mixed_mesh_preservation_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": metrics,
        "tolerances": {
            "maximum_automatic_region_relative_drift": automatic_limit,
            "maximum_total_element_relative_drift": total_limit,
        },
        "solver_ready": False,
        "lesson": (
            "Enable existing-mesh preservation before completing a mixed mesh. "
            "The manual region must retain its exact element inventory, while "
            "the automatic region may show bounded deterministic mesher-version "
            "drift that must be recorded separately."
        ),
    }
