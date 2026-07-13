"""Validation gate for first-order sphere tri/tet geometry refinement."""
from __future__ import annotations

import math


def _relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(abs(actual), abs(expected), 1.0e-300)


def linear_sphere_geometry_convergence_gate(
    rows,
    *,
    analytic_volume: float,
    analytic_surface_area: float,
    replay,
    max_reader_relative_error: float = 1.0e-12,
    max_surface_radius_error: float = 1.0e-12,
    max_final_geometry_relative_error: float = 3.0e-3,
    min_asymptotic_order: float = 1.8,
):
    if not isinstance(rows, list) or len(rows) < 4:
        raise ValueError("rows must contain at least four refinement levels")
    if not isinstance(replay, dict):
        raise ValueError("replay must be an object")
    analytic_volume = float(analytic_volume)
    analytic_surface_area = float(analytic_surface_area)
    tolerances = (
        max_reader_relative_error,
        max_surface_radius_error,
        max_final_geometry_relative_error,
        min_asymptotic_order,
    )
    if (
        not math.isfinite(analytic_volume)
        or not math.isfinite(analytic_surface_area)
        or analytic_volume <= 0.0
        or analytic_surface_area <= 0.0
    ):
        raise ValueError("analytic volume and surface area must be finite and positive")
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")

    required = {
        "level",
        "points",
        "triangles",
        "tets",
        "volume",
        "surface_area",
        "boundary_orientation",
        "maximum_surface_radius_error",
        "volume_reader_relative_error",
        "surface_reader_relative_error",
    }
    normalized = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError(f"row {index} must be an object")
        missing = sorted(required - row.keys())
        if missing:
            raise ValueError(f"row {index} is missing: {', '.join(missing)}")
        item = {
            "level": int(row["level"]),
            "points": int(row["points"]),
            "triangles": int(row["triangles"]),
            "tets": int(row["tets"]),
            "volume": float(row["volume"]),
            "surface_area": float(row["surface_area"]),
            "boundary_orientation": str(row["boundary_orientation"]).strip().lower(),
            "maximum_surface_radius_error": float(row["maximum_surface_radius_error"]),
            "volume_reader_relative_error": float(row["volume_reader_relative_error"]),
            "surface_reader_relative_error": float(row["surface_reader_relative_error"]),
        }
        numeric = [
            value
            for key, value in item.items()
            if key != "boundary_orientation"
        ]
        if not all(math.isfinite(value) for value in numeric):
            raise ValueError(f"row {index} contains a non-finite value")
        if min(item["points"], item["triangles"], item["tets"]) <= 0:
            raise ValueError(f"row {index} mesh counts must be positive")
        if item["volume"] <= 0.0 or item["surface_area"] <= 0.0:
            raise ValueError(f"row {index} geometry measures must be positive")
        if min(
            item["maximum_surface_radius_error"],
            item["volume_reader_relative_error"],
            item["surface_reader_relative_error"],
        ) < 0.0:
            raise ValueError(f"row {index} error measures must be nonnegative")
        item["volume_relative_error"] = _relative_error(
            item["volume"], analytic_volume
        )
        item["surface_area_relative_error"] = _relative_error(
            item["surface_area"], analytic_surface_area
        )
        normalized.append(item)

    levels = [row["level"] for row in normalized]
    volume_errors = [row["volume_relative_error"] for row in normalized]
    area_errors = [row["surface_area_relative_error"] for row in normalized]
    volume_orders = []
    area_orders = []
    for left, right, left_error, right_error, left_area, right_area in zip(
        normalized,
        normalized[1:],
        volume_errors,
        volume_errors[1:],
        area_errors,
        area_errors[1:],
    ):
        delta = right["level"] - left["level"]
        if delta <= 0 or min(left_error, right_error, left_area, right_area) <= 0.0:
            volume_orders.append(-math.inf)
            area_orders.append(-math.inf)
        else:
            volume_orders.append(math.log(left_error / right_error) / (delta * math.log(2.0)))
            area_orders.append(math.log(left_area / right_area) / (delta * math.log(2.0)))

    replay_level = int(replay.get("level", -1))
    replay_reference = next(
        (row for row in normalized if row["level"] == replay_level), None
    )
    replay_error = math.inf
    replay_counts_match = False
    if replay_reference is not None:
        replay_counts_match = all(
            int(replay.get(key, -1)) == replay_reference[key]
            for key in ("points", "triangles", "tets")
        )
        replay_error = max(
            _relative_error(float(replay.get("volume", math.nan)), replay_reference["volume"]),
            _relative_error(
                float(replay.get("surface_area", math.nan)),
                replay_reference["surface_area"],
            ),
        )

    count_scaling = all(
        right["triangles"] == left["triangles"] * 4 ** (right["level"] - left["level"])
        and right["tets"] == left["tets"] * 4 ** (right["level"] - left["level"])
        for left, right in zip(normalized, normalized[1:])
        if right["level"] > left["level"]
    )
    checks = {
        "levels_strictly_increase": all(
            right > left for left, right in zip(levels, levels[1:])
        ),
        "mesh_counts_strictly_increase": all(
            right["points"] > left["points"]
            and right["triangles"] > left["triangles"]
            and right["tets"] > left["tets"]
            for left, right in zip(normalized, normalized[1:])
        ),
        "triangle_and_tet_counts_follow_four_way_refinement": count_scaling,
        "all_boundaries_are_outward": all(
            row["boundary_orientation"] == "outward" for row in normalized
        ),
        "surface_vertices_stay_on_analytic_sphere": max(
            row["maximum_surface_radius_error"] for row in normalized
        )
        <= max_surface_radius_error,
        "independent_readers_agree": max(
            max(
                row["volume_reader_relative_error"],
                row["surface_reader_relative_error"],
            )
            for row in normalized
        )
        <= max_reader_relative_error,
        "volume_error_strictly_decreases": all(
            right < left for left, right in zip(volume_errors, volume_errors[1:])
        ),
        "surface_area_error_strictly_decreases": all(
            right < left for left, right in zip(area_errors, area_errors[1:])
        ),
        "finest_geometry_errors_within_tolerance": max(
            volume_errors[-1], area_errors[-1]
        )
        <= max_final_geometry_relative_error,
        "asymptotic_order_near_two": min(volume_orders[-1], area_orders[-1])
        >= min_asymptotic_order,
        "independent_replay_is_exact": replay_counts_match
        and replay_error <= max_reader_relative_error,
    }
    return {
        "policy": "linear_sphere_geometry_convergence_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "rows": normalized,
        "metrics": {
            "volume_relative_errors": volume_errors,
            "surface_area_relative_errors": area_errors,
            "volume_observed_orders": volume_orders,
            "surface_area_observed_orders": area_orders,
            "maximum_reader_relative_error": max(
                max(
                    row["volume_reader_relative_error"],
                    row["surface_reader_relative_error"],
                )
                for row in normalized
            ),
            "replay_relative_error": replay_error,
        },
        "lesson": (
            "For a first-order faceted sphere, require outward tri/tet topology, "
            "independent-reader agreement, monotone volume/area error, near-second-order "
            "asymptotics, and an exact replay before promoting the mesh family."
        ),
    }
