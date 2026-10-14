"""B-rep mass-property and topology invariants for CAD roundtrips."""
from __future__ import annotations

import math
from typing import Mapping


def _vector(row: Mapping, name: str) -> list[float]:
    value = row.get(name)
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError(f"{name} must contain three components")
    result = [float(component) for component in value]
    if not all(math.isfinite(component) for component in result):
        raise ValueError(f"{name} must be finite")
    return result


def _max_abs(left, right) -> float:
    return max(abs(a - b) for a, b in zip(left, right))


def _relative_error(reference: float, measured: float) -> float:
    return abs(measured - reference) / max(abs(reference), 1.0e-300)


def brep_mass_topology_roundtrip_gate(
    reference,
    measured_rows,
    *,
    expected_volume: float | None = None,
    expected_volume_rtol: float = 1.0e-6,
    volume_rtol: float = 1.0e-4,
    area_rtol: float = 2.0e-3,
    bbox_atol: float = 0.11,
    centroid_atol: float = 2.0e-4,
    required_import_modes=("heal", "noheal"),
) -> dict:
    """Gate mass properties, B-rep topology and import-mode invariance."""

    if not isinstance(reference, Mapping):
        raise ValueError("reference must be a mapping")
    rows = list(measured_rows)
    if not rows or any(not isinstance(row, Mapping) for row in rows):
        raise ValueError("measured_rows must contain mappings")
    tolerances = [expected_volume_rtol, volume_rtol, area_rtol, bbox_atol, centroid_atol]
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")

    ref_volume = float(reference.get("volume"))
    ref_area = float(reference.get("surface_area"))
    ref_bbox_min = _vector(reference, "bbox_min")
    ref_bbox_max = _vector(reference, "bbox_max")
    ref_centroid = _vector(reference, "center_of_mass")
    if not math.isfinite(ref_volume) or ref_volume <= 0.0 or not math.isfinite(ref_area) or ref_area <= 0.0:
        raise ValueError("reference volume and surface_area must be finite and positive")
    count_names = ("vertex_count", "edge_count", "face_count", "solid_count")
    ref_counts = {name: int(reference.get(name, -1)) for name in count_names}
    ref_euler = int(reference.get("boundary_euler_characteristic", -10**9))

    comparisons = []
    for row in rows:
        measured_volume = float(row.get("volume"))
        measured_area = float(row.get("surface_area"))
        bbox_min = _vector(row, "bbox_min")
        bbox_max = _vector(row, "bbox_max")
        centroid = _vector(row, "center_of_mass")
        counts = {name: int(row.get(name, -1)) for name in count_names}
        euler = int(row.get("boundary_euler_characteristic", -10**9))
        comparisons.append({
            "import_mode": str(row.get("import_mode") or "").strip().lower(),
            "center_semantics": str(row.get("center_semantics") or "").strip().lower(),
            "volume_relative_error": _relative_error(ref_volume, measured_volume),
            "area_relative_error": _relative_error(ref_area, measured_area),
            "bbox_absolute_error": max(_max_abs(ref_bbox_min, bbox_min), _max_abs(ref_bbox_max, bbox_max)),
            "centroid_absolute_error": _max_abs(ref_centroid, centroid),
            "entity_counts_match": counts == ref_counts,
            "boundary_euler_characteristic": euler,
            "euler_matches": euler == ref_euler,
            "solid_count_one": counts["solid_count"] == 1,
        })

    modes = [row["import_mode"] for row in comparisons]
    required_modes = {str(mode).strip().lower() for mode in required_import_modes}
    expected_error = None
    if expected_volume is not None:
        expected = float(expected_volume)
        if not math.isfinite(expected) or expected <= 0.0:
            raise ValueError("expected_volume must be finite and positive")
        expected_error = _relative_error(expected, ref_volume)
    checks = {
        "reference_center_is_mass_centroid": str(reference.get("center_semantics") or "").strip().lower() == "mass_centroid",
        "measured_centers_are_mass_centroids": all(row["center_semantics"] == "mass_centroid" for row in comparisons),
        "expected_reference_volume_agrees": expected_error is None or expected_error <= float(expected_volume_rtol),
        "required_import_modes_present": required_modes.issubset(set(modes)),
        "import_modes_unique": len(modes) == len(set(modes)),
        "all_volume_errors_ok": all(row["volume_relative_error"] <= float(volume_rtol) for row in comparisons),
        "all_area_errors_ok": all(row["area_relative_error"] <= float(area_rtol) for row in comparisons),
        "all_bbox_errors_ok": all(row["bbox_absolute_error"] <= float(bbox_atol) for row in comparisons),
        "all_centroid_errors_ok": all(row["centroid_absolute_error"] <= float(centroid_atol) for row in comparisons),
        "all_entity_counts_match": all(row["entity_counts_match"] for row in comparisons),
        "all_euler_characteristics_match": all(row["euler_matches"] for row in comparisons),
        "single_solid_preserved": all(row["solid_count_one"] for row in comparisons),
        "closed_genus_zero_boundary": ref_euler == 2,
    }
    return {
        "policy": "build123d_brep_mass_topology_roundtrip_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "expected_volume_relative_error": expected_error,
        "reference_entity_counts": ref_counts,
        "reference_boundary_euler_characteristic": ref_euler,
        "comparisons": comparisons,
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "lesson": (
            "Pair volume, area and bbox with mass-centroid semantics and the B-rep Euler invariant. "
            "A representative or bounding-box center is not a mass centroid, and equal volume does not prove equal topology."
        ),
    }
