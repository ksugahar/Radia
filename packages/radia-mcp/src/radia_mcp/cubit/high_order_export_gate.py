"""Replay gates for headless high-order Netgen exports from Cubit."""

from __future__ import annotations

import math
from typing import Any, Iterable, Mapping


def cubit_loft_high_order_vol_series_gate(
    rows: Iterable[Mapping[str, Any]],
    *,
    expected_orders: Iterable[int] = (1, 2, 3, 4, 5),
    min_quality: float = 0.2,
) -> dict[str, Any]:
    """Gate topology, curved payload, sidecars, and quality across export orders."""

    records = [dict(row) for row in rows]
    orders_expected = [int(value) for value in expected_orders]
    quality_limit = float(min_quality)
    if not records:
        raise ValueError("rows must not be empty")
    if not orders_expected or len(set(orders_expected)) != len(orders_expected):
        raise ValueError("expected_orders must be nonempty and unique")
    if quality_limit <= 0.0 or not math.isfinite(quality_limit):
        raise ValueError("min_quality must be finite and positive")

    normalized = []
    for index, row in enumerate(records):
        inventory = row.get("inventory")
        sidecar = row.get("sidecar")
        if not isinstance(inventory, Mapping) or not isinstance(sidecar, Mapping):
            raise ValueError(f"row {index} needs inventory and sidecar mappings")
        try:
            order = int(row["order"])
            curved_node_count = int(row["curved_node_count"])
            quality_minimum = float(row["quality_minimum"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"row {index} has malformed order/curving/quality") from exc
        if curved_node_count < 0 or not math.isfinite(quality_minimum):
            raise ValueError(f"row {index} has invalid curving or quality evidence")
        normalized.append(
            {
                "order": order,
                "curved_node_count": curved_node_count,
                "quality_minimum": quality_minimum,
                "inventory": dict(inventory),
                "sidecar": dict(sidecar),
            }
        )

    orders = [row["order"] for row in normalized]
    first_inventory = normalized[0]["inventory"]
    topology_keys = (
        "volume_elements",
        "surface_elements",
        "points",
        "volume_kind_counts",
        "surface_kind_counts",
        "materials",
        "boundary_names",
        "surface_section",
        "routing_hint",
    )
    topology_invariant = all(
        all(row["inventory"].get(key) == first_inventory.get(key) for key in topology_keys)
        for row in normalized
    )
    curved_counts = [row["curved_node_count"] for row in normalized]
    expected_curve_presence = all(
        bool(row["inventory"].get("curvedelements_present")) == (row["order"] > 1)
        for row in normalized
    )
    sidecar_matches = all(
        int(row["sidecar"].get("order", -1)) == row["order"]
        and int(row["sidecar"].get("n_elements", -1)) == int(row["inventory"].get("volume_elements", -2))
        and int(row["sidecar"].get("n_points", -1)) == int(row["inventory"].get("points", -2))
        for row in normalized
    )
    sidecar_names_invariant = all(
        set(row["sidecar"].get("materials", {})) == set(normalized[0]["sidecar"].get("materials", {}))
        and set(row["sidecar"].get("boundaries", {})) == set(normalized[0]["sidecar"].get("boundaries", {}))
        for row in normalized
    )
    checks = {
        "orders_match_expected_series": orders == orders_expected,
        "topology_and_labels_invariant": topology_invariant,
        "all_hex_volume_and_quad_boundary": (
            first_inventory.get("volume_kind_counts") == {"hex": 24}
            and first_inventory.get("surface_kind_counts") == {"quad": 40}
        ),
        "surfaceelementsuv_counted": first_inventory.get("surface_section") == "surfaceelementsuv",
        "cubit_hex_route_retained": first_inventory.get("routing_hint") == "cubit_hex_or_mixed_path",
        "curved_section_matches_order": expected_curve_presence,
        "curved_nodes_start_zero_then_strictly_increase": (
            curved_counts[0] == 0
            and all(right > left for left, right in zip(curved_counts[1:], curved_counts[2:]))
            and curved_counts[1] > 0
        ),
        "sidecar_counts_and_orders_match_inventory": sidecar_matches,
        "sidecar_material_and_boundary_names_invariant": sidecar_names_invariant,
        "mesh_quality_meets_floor": all(row["quality_minimum"] >= quality_limit for row in normalized),
    }
    return {
        "policy": "cubit_loft_high_order_vol_series_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "orders": orders,
        "curved_node_counts": curved_counts,
        "volume_elements": first_inventory.get("volume_elements"),
        "surface_elements": first_inventory.get("surface_elements"),
        "base_points": first_inventory.get("points"),
        "minimum_quality": min(row["quality_minimum"] for row in normalized),
        "lesson": (
            "High-order .vol exports keep the first-order hex/quad inventory and base point count; "
            "validate order through curvedelements evidence and the exporter sidecar instead of expecting "
            "the points or element rows to grow."
        ),
    }


def cubit_headless_netgen_export_gate(
    summary: Mapping[str, Any],
    *,
    expected_orders: Iterable[int] = (1, 2, 3, 4, 5),
    min_quality: float = 0.2,
) -> dict[str, Any]:
    """Gate GUI-free migration from a legacy plugin command to native Netgen export."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    expected = [int(value) for value in expected_orders]
    produced = [int(value) for value in summary.get("produced_orders", [])]
    quality = float(summary.get("quality_minimum", math.nan))
    exit_code = int(summary.get("exit_code", -999))
    known_startup_diagnostics = bool(summary.get("known_startup_plugin_path_diagnostics"))
    complete = bool(summary.get("artifact_set_complete"))
    checks = {
        "source_native_journal_recorded": bool(str(summary.get("source_journal") or "").strip()),
        "headless_batch_recorded": summary.get("headless") is True,
        "no_persistent_gui_started": summary.get("persistent_gui_started") is False,
        "legacy_gui_plugin_command_identified": str(summary.get("legacy_command") or "").strip().lower() == "radia_export",
        "native_export_netgen_command_used": str(summary.get("replay_command") or "").strip().lower() == "export netgen",
        "legacy_command_unavailable_headless": summary.get("legacy_command_available_headless") is False,
        "native_command_available_headless": summary.get("native_command_available_headless") is True,
        "all_expected_orders_and_sidecars_created": (
            produced == expected
            and int(summary.get("vol_file_count", -1)) == len(expected)
            and int(summary.get("sidecar_file_count", -1)) == len(expected)
            and complete
        ),
        "all_hex_mesh_completed": int(summary.get("hex_count", 0)) > 0,
        "quality_meets_floor": math.isfinite(quality) and quality >= float(min_quality),
        "exit_code_explained_by_known_startup_diagnostics": (
            exit_code == 0 or (exit_code == 2 and known_startup_diagnostics and complete)
        ),
        "no_model_or_export_errors": summary.get("model_or_export_errors") is False,
    }
    return {
        "policy": "cubit_headless_netgen_export_command_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "produced_orders": produced,
        "exit_code": exit_code,
        "quality_minimum": quality,
        "lesson": (
            "A GUI-installed export command is not a headless contract. Replay source journals with Cubit's "
            "native export netgen command, then judge completion from the full .vol+sidecar set, mesh quality, "
            "and classified startup diagnostics rather than the process exit code alone."
        ),
    }
