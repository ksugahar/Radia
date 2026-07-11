"""Topology and mass-property validation for a Platonic-solid CAD family."""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


EXPECTED = {
    "tetrahedron": (4, 6, 4),
    "cube": (8, 12, 6),
    "octahedron": (6, 12, 8),
    "dodecahedron": (20, 30, 12),
    "icosahedron": (12, 30, 20),
}


def platonic_solid_family_gate(
    summary: Mapping[str, Any],
    *,
    analytic_volume_rtol: float = 1.0e-12,
    self_roundtrip_rtol: float = 5.0e-12,
    external_volume_rtol: float = 5.0e-12,
) -> dict[str, object]:
    """Gate all five Platonic solids through topology and two CAD kernels."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    rows_raw = summary.get("rows")
    if not isinstance(rows_raw, Sequence) or isinstance(rows_raw, (str, bytes)):
        raise ValueError("rows must be a sequence")
    tolerances = (analytic_volume_rtol, self_roundtrip_rtol, external_volume_rtol)
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")
    rows = []
    for raw in rows_raw:
        if not isinstance(raw, Mapping):
            raise ValueError("each row must be a mapping")
        name = str(raw.get("name") or "").strip().lower()
        try:
            row = {
                "name": name,
                "vertices": int(raw["vertices"]),
                "edges": int(raw["edges"]),
                "faces": int(raw["faces"]),
                "analytic_relative_error": float(raw["analytic_relative_error"]),
                "self_roundtrip_relative_error": float(raw["self_roundtrip_relative_error"]),
                "external_relative_error": float(raw["external_relative_error"]),
                "external_volume_count": int(raw["external_volume_count"]),
                "valid": raw.get("valid") is True,
            }
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("row topology and error metrics are required") from exc
        rows.append(row)
    by_name = {row["name"]: row for row in rows}
    finite_errors = all(
        math.isfinite(row[key]) and row[key] >= 0.0
        for row in rows
        for key in (
            "analytic_relative_error",
            "self_roundtrip_relative_error",
            "external_relative_error",
        )
    )
    checks = {
        "upstream_identity_recorded": bool(str(summary.get("upstream_commit") or "").strip())
        and bool(str(summary.get("source_sha256") or "").strip()),
        "exact_source_execution_recorded": summary.get("source_execution_mode")
        == "exact_source_with_display_stub",
        "installed_api_contract_recorded": summary.get("shape_valid_access") == "property",
        "five_unique_named_solids": len(rows) == len(by_name) == 5 and set(by_name) == set(EXPECTED),
        "known_vertex_edge_face_counts": set(by_name) == set(EXPECTED)
        and all(
            (by_name[name]["vertices"], by_name[name]["edges"], by_name[name]["faces"])
            == counts
            for name, counts in EXPECTED.items()
        ),
        "euler_characteristic_two": bool(rows)
        and all(row["vertices"] - row["edges"] + row["faces"] == 2 for row in rows),
        "all_shapes_valid": bool(rows) and all(row["valid"] for row in rows),
        "error_metrics_finite": finite_errors,
        "analytic_volumes_match": finite_errors
        and all(row["analytic_relative_error"] <= analytic_volume_rtol for row in rows),
        "same_kernel_roundtrips_match": finite_errors
        and all(row["self_roundtrip_relative_error"] <= self_roundtrip_rtol for row in rows),
        "external_kernel_volumes_match": finite_errors
        and all(row["external_relative_error"] <= external_volume_rtol for row in rows),
        "external_imports_are_single_volumes": bool(rows)
        and all(row["external_volume_count"] == 1 for row in rows),
    }
    return {
        "policy": "build123d_platonic_solid_family_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "solid_count": len(rows),
            "max_analytic_relative_error": max(
                (row["analytic_relative_error"] for row in rows), default=math.inf
            ),
            "max_self_roundtrip_relative_error": max(
                (row["self_roundtrip_relative_error"] for row in rows), default=math.inf
            ),
            "max_external_relative_error": max(
                (row["external_relative_error"] for row in rows), default=math.inf
            ),
        },
        "lesson": (
            "Validate a custom CAD object family with exact V-E+F topology, "
            "closed-form volume, same-kernel STEP replay, and an independent CAD "
            "kernel. Matching volume alone does not prove the expected polyhedron."
        ),
    }
