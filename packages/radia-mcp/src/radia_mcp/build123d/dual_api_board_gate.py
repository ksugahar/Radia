"""Validation gate for equivalent Builder and Algebra CAD implementations."""

from __future__ import annotations

import math
from typing import Any, Mapping


def dual_api_perforated_board_gate(
    summary: Mapping[str, Any],
    *,
    max_analytic_relative_error: float = 1.0e-12,
    max_self_roundtrip_relative_error: float = 5.0e-12,
    max_external_volume_relative_error: float = 2.0e-6,
    max_external_area_relative_error: float = 5.0e-12,
    max_external_bbox_absolute_error: float = 1.0e-10,
) -> dict[str, Any]:
    """Gate two source APIs and two independent STEP import modes together."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    tolerances = (
        max_analytic_relative_error,
        max_self_roundtrip_relative_error,
        max_external_volume_relative_error,
        max_external_area_relative_error,
        max_external_bbox_absolute_error,
    )
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")

    api_rows = summary.get("api_rows") or []
    external_rows = summary.get("external_rows") or []
    dimensions = summary.get("board_dimensions") or {}
    holes = summary.get("hole_contract") or {}
    if not isinstance(api_rows, list) or not isinstance(external_rows, list):
        raise ValueError("api_rows and external_rows must be lists")
    if not isinstance(dimensions, Mapping) or not isinstance(holes, Mapping):
        raise ValueError("board_dimensions and hole_contract must be mappings")

    api_by_name = {
        str(row.get("implementation") or ""): row
        for row in api_rows
        if isinstance(row, Mapping)
    }
    external_by_key = {
        (
            str(row.get("import_mode") or ""),
            str(row.get("implementation") or ""),
        ): row
        for row in external_rows
        if isinstance(row, Mapping)
    }
    expected_apis = {"builder", "algebra"}
    expected_modes = {"noheal", "default_heal"}
    length = float(dimensions.get("length", math.nan))
    width = float(dimensions.get("width", math.nan))
    height = float(dimensions.get("height", math.nan))
    half_r1 = int(holes.get("radius1_boundary_half_holes", -1))
    full_r1 = int(holes.get("radius1_full_holes_after_overlap_and_large_hole_replacement", -1))
    full_r2 = int(holes.get("radius2_full_holes", -1))
    expected_removed_area = math.pi * (0.5 * half_r1 + full_r1 + 4.0 * full_r2)
    expected_volume = (length * width - expected_removed_area) * height
    recorded_volume = float(summary.get("analytic_volume", math.nan))

    topology_tuples = {
        (
            int(row.get("vertex_count", -1)),
            int(row.get("edge_count", -1)),
            int(row.get("face_count", -1)),
        )
        for row in api_by_name.values()
    }
    source_digests = {
        str(row.get("source_sha256") or "") for row in api_by_name.values()
    }
    step_digests = {str(row.get("step_sha256") or "") for row in api_by_name.values()}
    api_volumes = [float(row.get("volume", math.nan)) for row in api_by_name.values()]
    api_areas = [float(row.get("surface_area", math.nan)) for row in api_by_name.values()]
    api_bboxes = [row.get("bbox_size") or [] for row in api_by_name.values()]

    external_biases = [
        float(row.get("signed_volume_relative_bias", math.nan))
        for row in external_by_key.values()
    ]
    bias_spread = (
        max(external_biases) - min(external_biases)
        if external_biases and all(math.isfinite(value) for value in external_biases)
        else math.inf
    )
    diagnostics = [str(value) for value in summary.get("startup_diagnostics", [])]
    allowed_suffixes = ("/plugins", "-commandplugindir")
    startup_only_allowlisted = bool(diagnostics) and len(diagnostics) <= 2 and all(
        "Could not open file:" in row and row.rstrip().endswith(allowed_suffixes)
        for row in diagnostics
    )
    process_exit_code = int(summary.get("external_process_exit_code", -1))
    script_errors = [str(value) for value in summary.get("script_error_lines", [])]
    process_ok = process_exit_code == 0 or (
        process_exit_code in {2, 3}
        and startup_only_allowlisted
        and not script_errors
        and summary.get("external_result_artifact_fresh") is True
    )

    checks = {
        "upstream_commit_recorded": len(str(summary.get("upstream_commit") or "")) >= 7,
        "installed_version_recorded": bool(str(summary.get("version") or "").strip()),
        "builder_and_algebra_rows_recorded": set(api_by_name) == expected_apis,
        "source_examples_are_canonical_pair": {
            str(row.get("source_example") or "") for row in api_by_name.values()
        }
        == {"circuit_board.py", "circuit_board_algebra.py"},
        "source_digests_recorded_and_distinct": (
            len(source_digests) == 2 and all(len(value) == 64 for value in source_digests)
        ),
        "step_digests_recorded_and_distinct": (
            len(step_digests) == 2 and all(len(value) == 64 for value in step_digests)
        ),
        "exact_source_display_stub_mode": all(
            row.get("source_execution_mode") == "exact_source_with_display_stub"
            for row in api_by_name.values()
        ),
        "shape_valid_is_property": all(
            row.get("shape_valid_access") == "property" for row in api_by_name.values()
        ),
        "one_valid_solid_per_api": all(
            row.get("valid") is True and int(row.get("solid_count", 0)) == 1
            for row in api_by_name.values()
        ),
        "clipped_and_overlapping_hole_contract_recorded": (
            half_r1 == 31 and full_r1 == 25 and full_r2 == 4
        ),
        "analytic_volume_recomputed": (
            math.isfinite(expected_volume)
            and abs(recorded_volume - expected_volume) / max(abs(expected_volume), 1.0)
            <= float(max_analytic_relative_error)
        ),
        "api_analytic_volumes_match": all(
            float(row.get("analytic_relative_error", math.inf))
            <= float(max_analytic_relative_error)
            for row in api_by_name.values()
        ),
        "api_topology_matches": len(topology_tuples) == 1 and topology_tuples == {(190, 285, 97)},
        "api_volume_area_bbox_match": (
            len(api_volumes) == 2
            and max(api_volumes) - min(api_volumes) <= 1.0e-14 * max(abs(api_volumes[0]), 1.0)
            and max(api_areas) - min(api_areas) <= 1.0e-14 * max(abs(api_areas[0]), 1.0)
            and len(api_bboxes) == 2
            and len(api_bboxes[0]) == len(api_bboxes[1]) == 3
            and max(abs(float(a) - float(b)) for a, b in zip(api_bboxes[0], api_bboxes[1])) <= 1.0e-12
        ),
        "same_kernel_step_roundtrips_match": all(
            float(row.get("self_roundtrip_relative_error", math.inf))
            <= float(max_self_roundtrip_relative_error)
            for row in api_by_name.values()
        ),
        "two_external_import_modes_per_api": set(external_by_key) == {
            (mode, api) for mode in expected_modes for api in expected_apis
        },
        "external_single_volume_and_topology_match": all(
            int(row.get("volume_count", 0)) == 1 and int(row.get("surface_count", 0)) == 97
            for row in external_by_key.values()
        ),
        "external_volume_within_classified_tolerance": all(
            float(row.get("volume_relative_error", math.inf))
            <= float(max_external_volume_relative_error)
            for row in external_by_key.values()
        ),
        "external_area_and_bbox_match": all(
            float(row.get("surface_area_relative_error", math.inf))
            <= float(max_external_area_relative_error)
            and float(row.get("bbox_max_absolute_error", math.inf))
            <= float(max_external_bbox_absolute_error)
            for row in external_by_key.values()
        ),
        "external_bias_is_reproducible_across_api_and_import_mode": (
            len(external_biases) == 4
            and all(0.0 < value <= float(max_external_volume_relative_error) for value in external_biases)
            and bias_spread <= 1.0e-15
        ),
        "external_bias_classified_not_hidden": (
            summary.get("external_volume_bias_classification")
            == "systematic_kernel_mass_property_bias"
            and summary.get("external_volume_bias_tolerance_basis")
            == "dual_api_dual_import_mode_consistency"
        ),
        "external_cad_headless": summary.get("external_execution_mode")
        == "python_api_headless"
        and summary.get("persistent_gui_started") is False,
        "external_process_exit_semantics_acceptable": process_ok,
    }
    return {
        "policy": "build123d_dual_api_perforated_board_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "analytic_volume": expected_volume,
            "topology": list(next(iter(topology_tuples))) if len(topology_tuples) == 1 else None,
            "maximum_external_volume_relative_error": max(
                (float(row.get("volume_relative_error", math.inf)) for row in external_by_key.values()),
                default=math.inf,
            ),
            "external_volume_bias_spread": bias_spread,
            "external_process_exit_code": process_exit_code,
        },
        "external_launcher_classification": (
            "clean_exit"
            if process_exit_code == 0
            else "allowlisted_startup_diagnostic_with_clean_script"
            if process_ok
            else "execution_error"
        ),
        "notes": [
            "Builder and Algebra APIs should be equivalent in topology and mass properties, not merely similar on screen.",
            "Boundary-centered and overlapping holes require a clipped-area analytic contract rather than naive hole counting.",
            "A small external volume bias is accepted only when two source APIs and two import modes reproduce it while area, bounding box, and topology remain closed.",
            "Record the bias classification and tolerance basis; do not silently widen the CAD-volume tolerance.",
        ],
    }
