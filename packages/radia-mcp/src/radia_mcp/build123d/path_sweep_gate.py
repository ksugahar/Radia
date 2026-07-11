"""Validation gates for source-native build123d path sweeps."""

from __future__ import annotations

import math
from typing import Any, Mapping


def _positive(row: Mapping[str, Any], name: str) -> float:
    try:
        value = float(row[name])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"missing or invalid {name}") from exc
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _vector3(row: Mapping[str, Any], name: str) -> list[float]:
    try:
        values = [float(value) for value in row[name]]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"missing or invalid {name}") from exc
    if len(values) != 3 or not all(math.isfinite(value) and value > 0.0 for value in values):
        raise ValueError(f"{name} must contain three finite positive values")
    return values


def build123d_path_sweep_handoff_gate(
    result: Mapping[str, Any],
    *,
    path_length_rtol: float = 1.0e-12,
    external_volume_rtol: float = 5.0e-5,
    external_area_rtol: float = 1.0e-10,
    bbox_atol: float = 1.0e-10,
) -> dict[str, Any]:
    """Gate a curved sweep by path identity and independent STEP mass properties.

    ``section_area * path_length`` is intentionally diagnostic only. Composite
    curved paths can introduce transition overlap or trimming, so the finished
    solid volume must be checked against an independent CAD kernel.
    """

    if not isinstance(result, Mapping):
        raise ValueError("result must be a mapping")
    analytic = result.get("analytic")
    native = result.get("native")
    external = result.get("external")
    if not all(isinstance(row, Mapping) for row in (analytic, native, external)):
        raise ValueError("analytic, native, and external mappings are required")

    expected_length = _positive(analytic, "path_length_mm")
    section_area = _positive(analytic, "section_area_mm2")
    native_length = _positive(native, "path_length_mm")
    native_volume = _positive(native, "volume_mm3")
    native_area = _positive(native, "area_mm2")
    external_volume = _positive(external, "volume_mm3")
    external_area = _positive(external, "area_mm2")
    native_bbox = _vector3(native, "bbox_size_mm")
    external_bbox = _vector3(external, "bbox_size_mm")

    path_error = abs(native_length - expected_length) / expected_length
    volume_error = abs(external_volume - native_volume) / native_volume
    area_error = abs(external_area - native_area) / native_area
    bbox_error = max(abs(a - b) for a, b in zip(external_bbox, native_bbox))
    naive_volume = section_area * expected_length
    naive_gap = abs(native_volume - naive_volume) / native_volume
    oracle_policy = str(result.get("volume_oracle_policy") or "").strip().lower()
    digest = str(result.get("step_sha256") or "")

    checks = {
        "explicit_mm_units_recorded": str(result.get("length_unit") or "").lower() == "mm",
        "analytic_path_length_matches": path_error <= float(path_length_rtol),
        "native_valid_single_solid": bool(native.get("is_valid")) and int(native.get("solid_count", 0)) == 1,
        "same_step_digest_recorded": len(digest) == 64 and external.get("step_sha256") == digest,
        "external_single_volume": int(external.get("volume_count", 0)) == 1,
        "external_closed_boundary_topology": int(external.get("euler_characteristic", -999)) == 2,
        "external_volume_matches": volume_error <= float(external_volume_rtol),
        "external_area_matches": area_error <= float(external_area_rtol),
        "external_bbox_matches": bbox_error <= float(bbox_atol),
        "finished_volume_uses_cross_kernel_oracle": oracle_policy == "cross_kernel_mass_properties",
    }
    return {
        "policy": "build123d_path_sweep_handoff_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "path_length_relative_error": path_error,
            "external_volume_relative_error": volume_error,
            "external_area_relative_error": area_error,
            "external_bbox_absolute_error_mm": bbox_error,
            "naive_section_area_times_path_length_mm3": naive_volume,
            "naive_tube_volume_relative_gap": naive_gap,
        },
        "naive_tube_volume_is_oracle": False,
        "lesson": (
            "Use analytic arc-plus-line length to verify path construction, but validate the trimmed sweep "
            "solid with independent-kernel STEP mass properties; area times path length is diagnostic only."
        ),
    }


def build123d_path_sweep_source_contract_gate(result: Mapping[str, Any]) -> dict[str, Any]:
    """Gate the official algebra-mode sweep idiom and shape-validity API form."""

    if not isinstance(result, Mapping):
        raise ValueError("result must be a mapping")
    segments_raw = result.get("path_segments")
    if not isinstance(segments_raw, list):
        raise ValueError("path_segments must be a list")
    segments = [str(value).strip() for value in segments_raw]
    attribute_type = str(result.get("validity_attribute_type") or "").strip().lower()
    access = str(result.get("validity_access") or "").strip().lower()
    run_error = str(result.get("run_error") or "").strip()
    bool_called = "bool' object is not callable" in run_error.lower()
    method_called_on_bool = attribute_type == "bool" and access == "method"

    checks = {
        "build123d_source_recorded": str(result.get("ecosystem") or "").strip().lower() == "build123d",
        "official_intro_example_recorded": str(result.get("example_id") or "").strip().lower() == "introductory-ex14",
        "algebra_mode_recorded": str(result.get("api_mode") or "").strip().lower() == "algebra",
        "composite_path_segments_recorded": segments == ["JernArc", "JernArc", "Line"],
        "explicit_path_keyword_used": result.get("explicit_path_keyword") is True,
        "profile_plane_is_xz": str(result.get("profile_plane") or "").strip().upper() == "XZ",
        "rectangle_profile_recorded": str(result.get("profile_type") or "").strip().lower() == "rectangle",
        "build123d_version_recorded": bool(str(result.get("build123d_version") or "").strip()),
        "validity_attribute_form_recorded": attribute_type in {"bool", "callable"},
        "validity_access_matches_attribute": not method_called_on_bool and not bool_called,
        "source_replay_has_no_error": not run_error,
    }
    diagnosis = "source_contract_ok"
    if bool_called or method_called_on_bool:
        diagnosis = "is_valid_property_called_as_method"
    elif run_error:
        diagnosis = "source_replay_error"
    elif not all(checks.values()):
        diagnosis = "source_contract_incomplete"

    return {
        "policy": "build123d_path_sweep_source_contract_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "diagnosis": diagnosis,
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "lesson": (
            "Use an explicit path= sweep and align the profile plane with the path start. Inspect is_valid "
            "before access because build123d 0.10 exposes it as a bool property, not a method."
        ),
    }
