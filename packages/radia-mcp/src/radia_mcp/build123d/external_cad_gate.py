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


def step_portability_diagnosis_gate(
    result: Mapping[str, Any],
    *,
    self_roundtrip_rtol: float = 1.0e-12,
    external_volume_rtol: float = 2.0e-6,
    import_mode_rtol: float = 1.0e-12,
    required_import_modes: tuple[str, ...] = ("heal", "noheal"),
) -> dict[str, Any]:
    """Locate STEP volume loss in export, an external translator, or import mode."""

    if not isinstance(result, Mapping):
        raise ValueError("result must be a mapping")
    tolerances = (self_roundtrip_rtol, external_volume_rtol, import_mode_rtol)
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and non-negative")

    native_volume = _positive_metric(result, "native_volume")
    self_volume = _positive_metric(result, "self_roundtrip_volume")
    imports_raw = result.get("external_imports")
    if not isinstance(imports_raw, list) or not imports_raw:
        raise ValueError("external_imports must be a non-empty list")

    imports: list[dict[str, Any]] = []
    for index, row in enumerate(imports_raw):
        if not isinstance(row, Mapping):
            raise ValueError(f"external_imports[{index}] must be a mapping")
        mode = str(row.get("mode") or "").strip().lower()
        volume = _positive_metric(row, "volume")
        imports.append({
            "mode": mode,
            "volume": volume,
            "volume_count": int(row.get("volume_count", -1)),
            "relative_error": abs(volume - native_volume) / native_volume,
        })

    self_error = abs(self_volume - native_volume) / native_volume
    modes = [row["mode"] for row in imports]
    required_modes = {str(mode).strip().lower() for mode in required_import_modes}
    import_spread = (
        max(row["volume"] for row in imports) - min(row["volume"] for row in imports)
    ) / native_volume
    self_ok = self_error <= float(self_roundtrip_rtol)
    external_ok = all(row["relative_error"] <= float(external_volume_rtol) for row in imports)
    modes_invariant = import_spread <= float(import_mode_rtol)

    if not self_ok:
        diagnosis = "step_self_roundtrip_loss"
    elif not modes_invariant:
        diagnosis = "external_import_mode_inconsistency"
    elif not external_ok:
        diagnosis = "external_kernel_translation_loss"
    else:
        diagnosis = "portable"

    checks = {
        "upstream_native_source_recorded": result.get("source_kind") == "upstream_native_example",
        "upstream_commit_recorded": len(str(result.get("upstream_commit") or "")) == 40,
        "source_digest_recorded": len(str(result.get("source_sha256") or "")) == 64,
        "step_digest_recorded": len(str(result.get("step_sha256") or "")) == 64,
        "self_roundtrip_matches": self_ok,
        "required_import_modes_present": required_modes.issubset(set(modes)),
        "import_modes_unique": len(modes) == len(set(modes)),
        "single_volume_preserved": all(row["volume_count"] == 1 for row in imports),
        "external_volume_matches": external_ok,
        "external_import_modes_invariant": modes_invariant,
    }
    return {
        "policy": "build123d_step_portability_diagnosis_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "diagnosis": diagnosis,
        "healing_not_root_cause": diagnosis == "external_kernel_translation_loss" and modes_invariant,
        "native_volume": native_volume,
        "self_roundtrip_relative_error": self_error,
        "external_imports": imports,
        "external_import_mode_spread_relative": import_spread,
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "lesson": (
            "A STEP file is solver-ready only after native-to-STEP self-roundtrip and independent-kernel "
            "mass properties agree. Equal heal/noheal failures implicate translation compatibility, not healing."
        ),
    }


def curved_step_topology_crosscheck_gate(
    result: Mapping[str, Any],
    *,
    self_mass_rtol: float = 1.0e-7,
    external_volume_rtol: float = 1.0e-4,
    import_mode_rtol: float = 1.0e-10,
    bbox_atol: float = 1.0e-9,
) -> dict[str, Any]:
    """Gate a curved single-solid STEP across same and independent kernels.

    A wider curved-surface volume tolerance is accepted only when source
    identity is bound, same-kernel STEP mass/topology closes, two independent
    import modes agree, and external body/face/edge counts remain exact.
    """
    if not isinstance(result, Mapping):
        raise ValueError("result must be a mapping")
    tolerances = (self_mass_rtol, external_volume_rtol, import_mode_rtol, bbox_atol)
    if any(not math.isfinite(float(value)) or float(value) < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")

    native = result.get("native")
    roundtrip = result.get("roundtrip")
    imports_raw = result.get("external_imports")
    if not isinstance(native, Mapping) or not isinstance(roundtrip, Mapping):
        raise ValueError("native and roundtrip mappings are required")
    if not isinstance(imports_raw, list) or len(imports_raw) < 2:
        raise ValueError("at least two external import modes are required")

    native_volume = _positive_metric(native, "volume")
    native_area = _positive_metric(native, "area")
    roundtrip_volume = _positive_metric(roundtrip, "volume")
    roundtrip_area = _positive_metric(roundtrip, "area")
    native_bbox = _vector(native, "bbox_size")
    roundtrip_bbox = _vector(roundtrip, "bbox_size")
    native_faces = int(native.get("faces", 0))
    native_edges = int(native.get("edges", 0))
    native_solids = int(native.get("solids", 0))
    roundtrip_faces = int(roundtrip.get("faces", 0))
    roundtrip_edges = int(roundtrip.get("edges", 0))
    roundtrip_solids = int(roundtrip.get("solids", 0))

    self_volume_error = abs(roundtrip_volume - native_volume) / native_volume
    self_area_error = abs(roundtrip_area - native_area) / native_area
    bbox_error = max(abs(left - right) for left, right in zip(native_bbox, roundtrip_bbox))

    imports = []
    for index, row in enumerate(imports_raw):
        if not isinstance(row, Mapping):
            raise ValueError(f"external_imports[{index}] must be a mapping")
        volume = _positive_metric(row, "volume")
        imports.append({
            "mode": str(row.get("mode") or "").strip().lower(),
            "volume": volume,
            "volume_relative_error": abs(volume - native_volume) / native_volume,
            "body_count": int(row.get("body_count", -1)),
            "volume_count": int(row.get("volume_count", -1)),
            "surface_count": int(row.get("surface_count", -1)),
            "curve_count": int(row.get("curve_count", -1)),
            "step_sha256": str(row.get("step_sha256") or ""),
        })
    modes = [row["mode"] for row in imports]
    external_volumes = [row["volume"] for row in imports]
    import_spread = (max(external_volumes) - min(external_volumes)) / native_volume
    max_external_error = max(row["volume_relative_error"] for row in imports)
    classified_bias = max_external_error > self_mass_rtol
    classification_ok = (
        not classified_bias
        or (
            result.get("external_volume_bias_classification")
            == "cross_kernel_curved_surface_translation"
            and bool(str(result.get("external_volume_bias_tolerance_basis") or "").strip())
        )
    )
    step_sha = str(result.get("step_sha256") or "")
    checks = {
        "upstream_native_source_recorded": result.get("source_kind")
        == "upstream_native_example",
        "upstream_commit_recorded": len(str(result.get("upstream_commit") or "")) == 40,
        "source_digest_recorded": len(str(result.get("source_sha256") or "")) == 64,
        "build123d_version_recorded": bool(str(result.get("build123d_version") or "").strip()),
        "step_digest_recorded": len(step_sha) == 64,
        "native_and_roundtrip_valid_single_solid": (
            native.get("is_valid") is True
            and roundtrip.get("is_valid") is True
            and native_solids == roundtrip_solids == 1
        ),
        "same_kernel_volume_matches": self_volume_error <= self_mass_rtol,
        "same_kernel_area_matches": self_area_error <= self_mass_rtol,
        "same_kernel_bbox_matches": bbox_error <= bbox_atol,
        "same_kernel_topology_matches": (
            native_faces == roundtrip_faces and native_edges == roundtrip_edges
        ),
        "external_import_modes_recorded_and_unique": (
            {"heal", "noheal"}.issubset(set(modes)) and len(modes) == len(set(modes))
        ),
        "same_step_digest_used_by_external_imports": all(
            row["step_sha256"] == step_sha for row in imports
        ),
        "external_single_body_single_volume": all(
            row["body_count"] == row["volume_count"] == 1 for row in imports
        ),
        "external_face_edge_topology_matches": all(
            row["surface_count"] == native_faces and row["curve_count"] == native_edges
            for row in imports
        ),
        "external_volume_within_curved_step_tolerance": (
            max_external_error <= external_volume_rtol
        ),
        "external_import_modes_volume_invariant": import_spread <= import_mode_rtol,
        "cross_kernel_volume_bias_explicitly_classified": classification_ok,
        "timings_recorded": all(
            float((result.get("timings_s") or {}).get(name, -1.0)) >= 0.0
            for name in ("source_build", "step_export", "step_reimport", "external_import")
        ),
    }
    return {
        "policy": "build123d_curved_step_topology_crosscheck_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "diagnosis": (
            "portable_with_classified_cross_kernel_bias"
            if all(checks.values()) and classified_bias
            else "portable"
            if all(checks.values())
            else "needs_attention"
        ),
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "self_volume_relative_error": self_volume_error,
            "self_area_relative_error": self_area_error,
            "bbox_absolute_error": bbox_error,
            "maximum_external_volume_relative_error": max_external_error,
            "external_import_mode_volume_spread": import_spread,
            "native_face_count": native_faces,
            "native_edge_count": native_edges,
        },
        "external_imports": imports,
        "tolerances": {
            "self_mass_rtol": self_mass_rtol,
            "external_volume_rtol": external_volume_rtol,
            "import_mode_rtol": import_mode_rtol,
            "bbox_atol": bbox_atol,
        },
        "lesson": (
            "A curved STEP tolerance may be wider only when body, face, edge, "
            "source, digest, and independent import-mode invariants stay exact."
        ),
    }
