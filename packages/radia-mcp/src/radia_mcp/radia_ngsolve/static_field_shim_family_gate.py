"""Solver-neutral quality gate for a static-field source and shim family."""

from __future__ import annotations

import math
from typing import Any


def static_field_shim_family_gate(
    cases: list[dict[str, Any]],
    *,
    interaction_axis: str = "z",
    min_paired_source_field_ratio: float = 1.2,
    max_paired_source_uniformity_ratio: float = 0.5,
    min_shim_center_field_delta_relative: float = 0.01,
    max_center_transverse_relative: float = 1.0e-4,
    max_central_divergence_relative: float = 0.05,
) -> dict[str, Any]:
    """Gate field scale, ROI uniformity, shim sensitivity, and map quality."""

    if not isinstance(cases, list) or len(cases) != 4:
        raise ValueError("cases must contain exactly four role records")
    axis_index = {"x": 0, "y": 1, "z": 2}.get(interaction_axis)
    if axis_index is None:
        raise ValueError("interaction_axis must be x, y, or z")
    if (
        not math.isfinite(float(min_paired_source_field_ratio))
        or float(min_paired_source_field_ratio) <= 1.0
    ):
        raise ValueError("min_paired_source_field_ratio must be finite and greater than one")
    bounded = {
        "max_paired_source_uniformity_ratio": max_paired_source_uniformity_ratio,
        "min_shim_center_field_delta_relative": min_shim_center_field_delta_relative,
        "max_center_transverse_relative": max_center_transverse_relative,
        "max_central_divergence_relative": max_central_divergence_relative,
    }
    if any(
        not math.isfinite(float(value)) or not 0.0 <= float(value) <= 1.0
        for value in bounded.values()
    ):
        raise ValueError("relative tolerances must be finite and between zero and one")

    expected_roles = {"single_source", "paired_source", "balanced_shim", "offset_shim"}
    parsed: dict[str, dict[str, Any]] = {}
    metadata_ok = True
    for index, case in enumerate(cases):
        if not isinstance(case, dict):
            raise ValueError(f"case {index} must be an object")
        role = str(case.get("role") or "").strip()
        if role in parsed:
            raise ValueError(f"duplicate role: {role}")
        try:
            center = [float(value) for value in case["center_b_t"]]
            p2p = float(case["central_axial_peak_to_peak_ppm"])
            rms = float(case["central_axial_rms_ppm"])
            transverse = float(case["center_transverse_relative"])
            divergence = float(case["central_divergence_max_relative"])
            shape = [int(value) for value in case["grid_shape"]]
            row_count = int(case["row_count"])
            sample_count = int(case["central_sample_count"])
            roi_half_width = float(case["roi_half_width_mm"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"case {index} has invalid field-map metadata") from exc
        numeric = [*center, p2p, rms, transverse, divergence, roi_half_width]
        if len(center) != 3 or not all(math.isfinite(value) for value in numeric):
            raise ValueError(f"case {index} field metrics must be finite")
        if len(shape) != 3 or any(value < 2 for value in shape):
            raise ValueError(f"case {index}.grid_shape must contain three dimensions >= 2")
        parsed[role] = {
            "center": center,
            "p2p": p2p,
            "rms": rms,
            "transverse": transverse,
            "divergence": divergence,
            "shape": shape,
            "row_count": row_count,
            "sample_count": sample_count,
            "roi_half_width": roi_half_width,
        }
        metadata_ok = metadata_ok and (
            case.get("coordinate_unit") == "mm"
            and case.get("field_unit") == "T"
            and row_count == math.prod(shape)
            and sample_count >= 27
            and roi_half_width > 0.0
            and p2p >= 0.0
            and rms >= 0.0
            and rms <= p2p
            and transverse >= 0.0
            and divergence >= 0.0
        )
    if set(parsed) != expected_roles:
        raise ValueError(f"roles must be exactly {sorted(expected_roles)}")

    axial = {role: row["center"][axis_index] for role, row in parsed.items()}
    nonzero_axial = all(abs(value) > 0.0 for value in axial.values())
    same_sign = nonzero_axial and all(
        value * axial["single_source"] > 0.0 for value in axial.values()
    )
    paired_field_ratio = abs(axial["paired_source"]) / max(
        abs(axial["single_source"]), 1.0e-30
    )
    paired_uniformity_ratio = parsed["paired_source"]["p2p"] / max(
        parsed["single_source"]["p2p"], 1.0e-30
    )
    shim_center_delta_relative = abs(
        abs(axial["balanced_shim"]) - abs(axial["offset_shim"])
    ) / max(abs(axial["balanced_shim"]), abs(axial["offset_shim"]), 1.0e-30)

    checks = {
        "field_map_units_shape_and_roi_recorded": metadata_ok,
        "common_nonzero_axial_direction": same_sign,
        "paired_source_increases_center_field": paired_field_ratio
        >= min_paired_source_field_ratio,
        "paired_source_improves_roi_uniformity": paired_uniformity_ratio
        <= max_paired_source_uniformity_ratio,
        "both_shim_candidates_improve_over_single_source": max(
            parsed["balanced_shim"]["p2p"], parsed["offset_shim"]["p2p"]
        )
        < parsed["single_source"]["p2p"],
        "shim_variants_are_resolved": shim_center_delta_relative
        >= min_shim_center_field_delta_relative,
        "center_transverse_field_is_bounded": max(
            row["transverse"] for row in parsed.values()
        )
        <= max_center_transverse_relative,
        "coarse_map_divergence_is_bounded": max(
            row["divergence"] for row in parsed.values()
        )
        <= max_central_divergence_relative,
    }
    return {
        "policy": "static_field_shim_family_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "interaction_axis": interaction_axis,
        "metrics": {
            "center_axial_field_by_role_t": axial,
            "paired_source_field_ratio": paired_field_ratio,
            "paired_source_uniformity_ratio": paired_uniformity_ratio,
            "shim_center_field_delta_relative": shim_center_delta_relative,
            "max_center_transverse_relative": max(
                row["transverse"] for row in parsed.values()
            ),
            "max_central_divergence_relative": max(
                row["divergence"] for row in parsed.values()
            ),
        },
        "lesson": (
            "Evaluate a static-field design as a family: establish the single-source "
            "baseline, require the paired source to increase field and improve ROI "
            "uniformity, then prove that shim variants create a resolved change. Record "
            "the actual field-grid shape and keep transverse-field and div(B) diagnostics."
        ),
    }
