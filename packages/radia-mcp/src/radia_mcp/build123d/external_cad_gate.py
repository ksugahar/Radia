"""Validation gates for upstream build123d examples and external CAD kernels."""

from __future__ import annotations

import math
from typing import Any, Mapping


def _positive_metric(row: Mapping[str, Any], name: str) -> float:
    try:
        value = float(row[name])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"missing or invalid {name}") from exc
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return value


def _vector(row: Mapping[str, Any], name: str) -> list[float]:
    try:
        values = [float(value) for value in row[name]]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"missing or invalid {name}") from exc
    if len(values) != 3 or not all(math.isfinite(value) for value in values):
        raise ValueError(f"{name} must contain three finite values")
    return values


def build123d_upstream_example_roundtrip_gate(
    result: Mapping[str, Any],
    *,
    mass_property_rtol: float = 1.0e-12,
    centroid_atol: float = 1.0e-12,
) -> dict[str, Any]:
    """Gate official source identity and a build123d STEP self-roundtrip."""

    if not isinstance(result, Mapping):
        raise ValueError("result must be a mapping")
    native = result.get("native")
    imported = result.get("roundtrip")
    if not isinstance(native, Mapping) or not isinstance(imported, Mapping):
        raise ValueError("native and roundtrip mappings are required")
    reference_volume = _positive_metric(native, "volume")
    reference_area = _positive_metric(native, "area")
    imported_volume = _positive_metric(imported, "volume")
    imported_area = _positive_metric(imported, "area")
    reference_center = _vector(native, "centroid")
    imported_center = _vector(imported, "centroid")
    reference_bbox = _vector(native, "bbox_size")
    imported_bbox = _vector(imported, "bbox_size")
    topology_keys = ("vertices", "edges", "faces", "solids", "euler_characteristic")
    volume_error = abs(imported_volume - reference_volume) / reference_volume
    area_error = abs(imported_area - reference_area) / reference_area
    centroid_error = max(abs(a - b) for a, b in zip(imported_center, reference_center))
    bbox_error = max(abs(a - b) for a, b in zip(imported_bbox, reference_bbox))
    checks = {
        "upstream_native_source_recorded": result.get("source_kind") == "upstream_native_example",
        "upstream_commit_recorded": len(str(result.get("upstream_commit") or "")) == 40,
        "source_digest_recorded": len(str(result.get("source_sha256") or "")) == 64,
        "build123d_version_recorded": bool(str(result.get("build123d_version") or "").strip()),
        "step_digest_recorded": len(str(result.get("step_sha256") or "")) == 64,
        "single_solid_reference_and_roundtrip": int(native.get("solids", 0)) == int(imported.get("solids", 0)) == 1,
        "volume_roundtrip_matches": volume_error <= float(mass_property_rtol),
        "area_roundtrip_matches": area_error <= float(mass_property_rtol),
        "centroid_roundtrip_matches": centroid_error <= float(centroid_atol),
        "bbox_roundtrip_matches": bbox_error <= float(centroid_atol),
        "brep_topology_roundtrip_matches": all(native.get(key) == imported.get(key) for key in topology_keys),
        "timings_recorded": all(
            float((result.get("timings_s") or {}).get(name, -1.0)) >= 0.0
            for name in ("source_build", "step_export", "step_reimport")
        ),
    }
    return {
        "policy": "build123d_upstream_example_roundtrip_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "volume_relative_error": volume_error,
            "area_relative_error": area_error,
            "centroid_absolute_error": centroid_error,
            "bbox_absolute_error": bbox_error,
            "euler_characteristic": native.get("euler_characteristic"),
        },
        "lesson": "Bind an upstream example to its commit and source digest before treating its STEP roundtrip as durable teaching evidence.",
    }


def external_cad_mass_topology_crosscheck_gate(
    reference: Mapping[str, Any],
    external: Mapping[str, Any],
    *,
    volume_rtol: float = 2.0e-6,
    area_rtol: float = 1.0e-10,
    bbox_atol: float = 1.0e-10,
    centroid_atol: float = 1.0e-6,
) -> dict[str, Any]:
    """Compare CAD kernels while keeping entity-center and mass-centroid semantics distinct."""

    if not isinstance(reference, Mapping) or not isinstance(external, Mapping):
        raise ValueError("reference and external must be mappings")
    reference_volume = _positive_metric(reference, "volume")
    reference_area = _positive_metric(reference, "area")
    external_volume = _positive_metric(external, "volume")
    external_area = _positive_metric(external, "area")
    reference_bbox = _vector(reference, "bbox_size")
    external_bbox = _vector(external, "bbox_size")
    semantics = str(external.get("center_semantics") or "").strip().lower()
    if semantics not in {"mass_centroid", "entity_center_excluded"}:
        raise ValueError("center_semantics must be mass_centroid or entity_center_excluded")
    centroid_error = None
    if semantics == "mass_centroid":
        centroid_error = max(
            abs(a - b)
            for a, b in zip(
                _vector(reference, "center_of_mass"),
                _vector(external, "center_of_mass"),
            )
        )
    else:
        _vector(external, "representative_center")
    volume_error = abs(external_volume - reference_volume) / reference_volume
    area_error = abs(external_area - reference_area) / reference_area
    bbox_error = max(abs(a - b) for a, b in zip(external_bbox, reference_bbox))
    topology_keys = ("vertices", "edges", "faces", "solids", "euler_characteristic")
    checks = {
        "distinct_kernel_sources_recorded": (
            bool(str(reference.get("source") or "").strip())
            and bool(str(external.get("source") or "").strip())
            and reference.get("source") != external.get("source")
        ),
        "same_step_digest_recorded": (
            len(str(reference.get("step_sha256") or "")) == 64
            and reference.get("step_sha256") == external.get("step_sha256")
        ),
        "volume_matches": volume_error <= float(volume_rtol),
        "area_matches": area_error <= float(area_rtol),
        "bbox_matches": bbox_error <= float(bbox_atol),
        "brep_topology_matches": all(reference.get(key) == external.get(key) for key in topology_keys),
        "center_semantics_explicit": semantics in {"mass_centroid", "entity_center_excluded"},
        "centroid_matches_when_comparable": centroid_error is None or centroid_error <= float(centroid_atol),
    }
    return {
        "policy": "external_cad_mass_topology_crosscheck_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "volume_relative_error": volume_error,
            "area_relative_error": area_error,
            "bbox_absolute_error": bbox_error,
            "centroid_absolute_error": centroid_error,
            "center_comparison_performed": centroid_error is not None,
            "center_semantics": semantics,
        },
        "lesson": (
            "An entity center or bounding-box center is not a mass centroid. Exclude it explicitly from "
            "centroid validation while retaining volume, area, bbox, STEP digest, and Euler topology checks."
        ),
    }
