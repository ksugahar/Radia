"""Semantic inventory for Netgen ``.vol`` exports from Cubit/Coreform.

This helper is intentionally broader than
``radia_mcp.radia_ngsolve.netgen_vol``.  The radia-ngsolve parser is the
first-order education path and rejects anything except boundary triangles and
volume tetrahedra.  Cubit, however, is the lab's hex-led and mixed-mesh lane,
so the MCP server also needs a light preflight that can *recognize* hex,
pyramid, wedge, tet, quad, and tri records before routing the file.
"""

from __future__ import annotations

from collections import Counter
from math import isfinite
from pathlib import Path
from typing import Iterable, Mapping


VOLUME_KIND_BY_NP = {
    4: "tet",
    5: "pyramid",
    6: "wedge",
    8: "hex",
}

SURFACE_KIND_BY_NP = {
    3: "triangle",
    4: "quad",
}


def read_netgen_vol_inventory(path: str | Path) -> dict[str, object]:
    """Read a Netgen ``.vol`` file and return a semantic element inventory."""

    p = Path(path)
    return summarize_netgen_vol_inventory(p.read_text(encoding="utf-8"), source=str(p))


def cubit_hex_quality_gate(
    scaled_jacobians: Iterable[float],
    *,
    expected_hex_count: int | None = None,
    min_scaled_jacobian: float = 0.2,
) -> dict[str, object]:
    """Summarize Cubit hex scaled-Jacobian values as a routing/quality gate.

    Cubit reports scaled Jacobian with 1 as ideal and nonpositive values as
    inverted or degenerate.  This helper is deliberately independent of a live
    Cubit session so archived raw JSON can be replayed by the MCP server.
    """

    values = [float(value) for value in scaled_jacobians]
    threshold = float(min_scaled_jacobian)
    if not values:
        raise ValueError("scaled_jacobians must not be empty")
    if threshold <= 0.0:
        raise ValueError("min_scaled_jacobian must be > 0")
    if expected_hex_count is not None and int(expected_hex_count) < 0:
        raise ValueError("expected_hex_count must be non-negative")

    count = len(values)
    min_value = min(values)
    max_value = max(values)
    mean_value = sum(values) / count
    bad_count = sum(1 for value in values if value < threshold)
    checks = {
        "all_positive": min_value > 0.0,
        "min_scaled_jacobian_ok": min_value >= threshold,
        "expected_hex_count_ok": True
        if expected_hex_count is None else count == int(expected_hex_count),
    }
    return {
        "policy": "cubit_hex_scaled_jacobian_quality_gate",
        "metric": "scaled jacobian",
        "count": count,
        "expected_hex_count": expected_hex_count,
        "min": min_value,
        "max": max_value,
        "mean": mean_value,
        "min_scaled_jacobian": threshold,
        "bad_count": bad_count,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
    }


def cubit_element_quality_gate(
    values: Iterable[float],
    *,
    element_type: str,
    metric: str = "scaled_jacobian",
    expected_count: int | None = None,
    min_value: float = 0.2,
) -> dict[str, object]:
    """Replay Cubit element-quality metrics without requiring a live session.

    Cubit 2026.6 added more explicit higher-order Tetra10/Tri6 Jacobian
    metrics.  The MCP-facing contract should stay solver-independent: archive
    the raw metric list, then run this lower-bound gate before routing a mesh to
    NGSolve or to a heavier Coreform validation.
    """

    element = element_type.strip().lower().replace(" ", "")
    if not element:
        raise ValueError("element_type must not be empty")
    metric_name = metric.strip().lower().replace(" ", "_")
    if not metric_name:
        raise ValueError("metric must not be empty")

    samples = [float(value) for value in values]
    threshold = float(min_value)
    if not samples:
        raise ValueError("values must not be empty")
    if threshold <= 0.0:
        raise ValueError("min_value must be > 0")
    if expected_count is not None and int(expected_count) < 0:
        raise ValueError("expected_count must be non-negative")

    count = len(samples)
    finite_count = sum(1 for value in samples if isfinite(value))
    min_sample = min(samples)
    max_sample = max(samples)
    mean_sample = sum(samples) / count
    bad_count = sum(1 for value in samples if value < threshold or not isfinite(value))
    checks = {
        "all_finite": finite_count == count,
        "min_value_ok": min_sample >= threshold,
        "expected_count_ok": True
        if expected_count is None else count == int(expected_count),
    }
    return {
        "policy": "cubit_element_quality_metric_gate",
        "element_type": element,
        "metric": metric_name,
        "count": count,
        "expected_count": expected_count,
        "min": min_sample,
        "max": max_sample,
        "mean": mean_sample,
        "min_value": threshold,
        "bad_count": bad_count,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use this archived-metric gate for Coreform Cubit 2026.6-style "
            "higher-order Tetra10/Tri6 Jacobian checks and for older Cubit "
            "quality lists when the metric is a positive lower-bound quality."
        ),
    }


def cubit_quality_distribution_gate(
    values: Iterable[float],
    *,
    element_type: str = "hex",
    metric: str = "scaled_jacobian",
    expected_count: int | None = None,
    min_value: float = 0.2,
    quantile_levels: Iterable[float] = (0.0, 0.01, 0.05, 0.5, 0.95, 0.99, 1.0),
    histogram_edges: Iterable[float] = (0.0, 0.2, 0.5, 0.8, 0.95, 1.0),
) -> dict[str, object]:
    """Replay a Cubit quality list with quantiles and a compact histogram.

    A single minimum value is useful for routing, but hex-led Cubit workflows
    also need to know whether quality loss is isolated or broad.  This helper
    keeps that knowledge in archived JSON so the MCP server can learn from
    headless Cubit runs without reopening Cubit.
    """

    element = element_type.strip().lower().replace(" ", "")
    metric_name = metric.strip().lower().replace(" ", "_")
    samples = [float(value) for value in values]
    threshold = float(min_value)
    q_levels = [float(level) for level in quantile_levels]
    edges = [float(edge) for edge in histogram_edges]

    if not element:
        raise ValueError("element_type must not be empty")
    if not metric_name:
        raise ValueError("metric must not be empty")
    if not samples:
        raise ValueError("values must not be empty")
    if threshold <= 0.0:
        raise ValueError("min_value must be > 0")
    if expected_count is not None and int(expected_count) < 0:
        raise ValueError("expected_count must be non-negative")
    if not q_levels:
        raise ValueError("quantile_levels must not be empty")
    if any(level < 0.0 or level > 1.0 or not isfinite(level) for level in q_levels):
        raise ValueError("quantile_levels must be finite values in [0, 1]")
    if len(edges) < 2:
        raise ValueError("histogram_edges must contain at least two values")
    if any(not isfinite(edge) for edge in edges):
        raise ValueError("histogram_edges must be finite")
    if any(right <= left for left, right in zip(edges, edges[1:])):
        raise ValueError("histogram_edges must be strictly increasing")

    count = len(samples)
    finite_values = sorted(value for value in samples if isfinite(value))
    finite_count = len(finite_values)
    nonfinite_count = count - finite_count
    low_count = sum(1 for value in finite_values if value < threshold)
    min_sample = finite_values[0] if finite_values else None
    max_sample = finite_values[-1] if finite_values else None
    mean_sample = sum(finite_values) / finite_count if finite_values else None
    quantiles = {
        _quantile_key(level): _linear_quantile(finite_values, level)
        for level in q_levels
    }
    histogram = _quality_histogram(finite_values, edges)
    histogram_count = (
        sum(int(bin_record["count"]) for bin_record in histogram["bins"])
        + int(histogram["underflow"])
        + int(histogram["overflow"])
    )
    checks = {
        "all_finite": nonfinite_count == 0,
        "min_value_ok": finite_count > 0 and low_count == 0,
        "expected_count_ok": True
        if expected_count is None else count == int(expected_count),
        "histogram_count_ok": histogram_count == finite_count,
    }
    return {
        "policy": "cubit_quality_distribution_gate",
        "element_type": element,
        "metric": metric_name,
        "count": count,
        "expected_count": expected_count,
        "finite_count": finite_count,
        "nonfinite_count": nonfinite_count,
        "min": min_sample,
        "max": max_sample,
        "mean": mean_sample,
        "min_value": threshold,
        "low_count": low_count,
        "bad_count": low_count + nonfinite_count,
        "quantiles": quantiles,
        "histogram": histogram,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use this replay gate when Cubit/Coreform mesh quality must be "
            "learned as a distribution, not only as a single minimum."
        ),
    }


def cubit_vol_label_metadata_gate(
    inventory: dict[str, object],
    *,
    required_materials: Iterable[str] = (),
    required_boundaries: Iterable[str] = (),
) -> dict[str, object]:
    """Check that a Cubit/Coreform ``.vol`` export preserved labels.

    Geometry and quality gates are not enough for downstream solver setup.  The
    exported mesh also has to carry material/block labels and boundary names so
    notebooks or MCP tools can apply sources and boundary conditions without
    relying on row position or private Cubit ids.
    """

    materials_raw = inventory.get("materials", {})
    boundaries_raw = inventory.get("boundary_names", {})
    if not isinstance(materials_raw, dict):
        raise ValueError("inventory['materials'] must be a mapping")
    if not isinstance(boundaries_raw, dict):
        raise ValueError("inventory['boundary_names'] must be a mapping")

    material_names = [str(value) for value in materials_raw.values()]
    boundary_names = [str(value) for value in boundaries_raw.values()]
    required_material_names = [str(value) for value in required_materials if str(value)]
    required_boundary_names = [str(value) for value in required_boundaries if str(value)]
    checks = {
        "material_labels_recorded": bool(material_names),
        "material_labels_unique": len(material_names) == len(set(material_names)),
        "required_materials_present": all(name in material_names for name in required_material_names),
        "boundary_names_recorded": bool(boundary_names),
        "boundary_names_may_be_grouped": True,
        "required_boundaries_present": all(name in boundary_names for name in required_boundary_names),
        "volume_elements_present": int(inventory.get("volume_elements", 0) or 0) > 0,
        "surface_elements_present": int(inventory.get("surface_elements", 0) or 0) > 0,
    }
    return {
        "policy": "cubit_vol_label_metadata_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "materials": material_names,
        "boundary_names": boundary_names,
        "required_materials": required_material_names,
        "required_boundaries": required_boundary_names,
        "checks": checks,
        "notes": [
            "Run this before solver setup so material and boundary conditions are not inferred from numeric ids alone.",
            "For Cubit hex-led workflows, label metadata is part of solver readiness alongside volume/area and quality gates.",
        ],
    }


def cubit_mixed_transition_metadata_gate(
    inventory: dict[str, object],
    *,
    required_volume_kinds: Iterable[str] = ("hex", "pyramid", "tet"),
    transition_kinds: Iterable[str] = ("pyramid",),
    transition_material_names: Iterable[str] = ("pyramid_transition", "pyram", "transition"),
) -> dict[str, object]:
    """Check that a mixed Cubit ``.vol`` records transition elements explicitly.

    Cubit is the lab's hex-led route.  When a hex region is connected to a tet
    region, pyramid records are meaningful topology even if a companion sidecar
    reports zero material volume for the transition block.  This replay helper
    keeps the routing contract executable from archived ``.vol`` inventory.
    """

    volume_counts_raw = inventory.get("volume_kind_counts", {})
    surface_counts_raw = inventory.get("surface_kind_counts", {})
    materials_raw = inventory.get("materials", {})
    if not isinstance(volume_counts_raw, dict):
        raise ValueError("inventory['volume_kind_counts'] must be a mapping")
    if not isinstance(surface_counts_raw, dict):
        raise ValueError("inventory['surface_kind_counts'] must be a mapping")
    if not isinstance(materials_raw, dict):
        raise ValueError("inventory['materials'] must be a mapping")

    volume_counts = {str(key): int(value) for key, value in volume_counts_raw.items()}
    surface_counts = {str(key): int(value) for key, value in surface_counts_raw.items()}
    materials = [str(value) for value in materials_raw.values()]
    required = [str(kind).strip().lower() for kind in required_volume_kinds if str(kind).strip()]
    transitions = [str(kind).strip().lower() for kind in transition_kinds if str(kind).strip()]
    transition_names = [str(name) for name in transition_material_names if str(name)]
    routing_hint = str(inventory.get("routing_hint", ""))
    checks = {
        "volume_elements_present": int(inventory.get("volume_elements", 0) or 0) > 0,
        "surface_elements_present": int(inventory.get("surface_elements", 0) or 0) > 0,
        "required_volume_kinds_present": all(volume_counts.get(kind, 0) > 0 for kind in required),
        "transition_kinds_present": all(volume_counts.get(kind, 0) > 0 for kind in transitions),
        "contains_hex_region": volume_counts.get("hex", 0) > 0,
        "contains_tet_region": volume_counts.get("tet", 0) > 0,
        "routing_hint_is_cubit_mixed": routing_hint == "cubit_hex_or_mixed_path",
        "not_tri_tet_only": not bool(inventory.get("is_tri_tet_only", False)),
        "material_labels_recorded": bool(materials),
        "transition_material_label_present": (
            True if not transition_names else any(name in materials for name in transition_names)
        ),
    }
    return {
        "policy": "cubit_mixed_transition_metadata_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "volume_kind_counts": volume_counts,
        "surface_kind_counts": surface_counts,
        "required_volume_kinds": required,
        "transition_kinds": transitions,
        "materials": materials,
        "transition_material_names": transition_names,
        "routing_hint": routing_hint,
        "checks": checks,
        "notes": [
            "Run this before importing hex+pyramid+tet meshes into radia-ngsolve or Radia-style mixed-element lanes.",
            "A pyramid transition can have zero sidecar material volume and still be required topology in the .vol inventory.",
        ],
    }


def cubit_bnd_area_interface_gate(
    *,
    external_area: float,
    material_interface_area: float,
    ngsolve_bnd_area: float,
    abs_tol: float = 1e-9,
    rel_tol: float = 1e-9,
) -> dict[str, object]:
    """Check the NGSolve ``BND`` area convention for Cubit multi-material meshes.

    Netgen ``.vol`` files exported from Cubit include material-interface faces
    in NGSolve boundary integration.  For imprinted/merged multi-volume meshes,
    ``Integrate(1, mesh, BND)`` should therefore be compared with external area
    plus each shared material-interface area once, not with the external CAD
    area alone.
    """

    external = float(external_area)
    interface = float(material_interface_area)
    observed = float(ngsolve_bnd_area)
    abs_threshold = float(abs_tol)
    rel_threshold = float(rel_tol)
    if abs_threshold < 0.0:
        raise ValueError("abs_tol must be non-negative")
    if rel_threshold < 0.0:
        raise ValueError("rel_tol must be non-negative")

    expected = external + interface
    abs_error = abs(observed - expected)
    rel_error = abs_error / max(abs(expected), 1.0)
    external_only_abs_error = abs(observed - external)
    external_only_rel_error = external_only_abs_error / max(abs(external), 1.0)
    tolerance = max(abs_threshold, rel_threshold * max(abs(expected), 1.0))
    external_only_tolerance = max(abs_threshold, rel_threshold * max(abs(external), 1.0))
    checks = {
        "all_finite": all(isfinite(value) for value in (external, interface, observed)),
        "external_area_nonnegative": external >= 0.0,
        "material_interface_area_nonnegative": interface >= 0.0,
        "ngsolve_bnd_area_nonnegative": observed >= 0.0,
        "matches_external_plus_interface": abs_error <= tolerance,
    }
    return {
        "policy": "cubit_ngsolve_bnd_area_includes_material_interfaces_once",
        "external_area": external,
        "material_interface_area": interface,
        "expected_bnd_area": expected,
        "ngsolve_bnd_area": observed,
        "abs_error": abs_error,
        "rel_error": rel_error,
        "abs_tol": abs_threshold,
        "rel_tol": rel_threshold,
        "external_only_abs_error": external_only_abs_error,
        "external_only_rel_error": external_only_rel_error,
        "matches_external_only": external_only_abs_error <= external_only_tolerance,
        "checks": checks,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "version_note": (
            "Use this gate for Cubit/Coreform imprinted multi-volume .vol exports: "
            "NGSolve BND integrates external boundary faces plus material-interface "
            "faces once."
        ),
    }


def cubit_mass_property_sidecar_gate(
    rows: Iterable[Mapping[str, object]],
    *,
    expected_total_volume: float | None = None,
    expected_total_area: float | None = None,
    expected_bbox_size: Iterable[float] | None = None,
    rel_tol: float = 1e-9,
    abs_tol: float = 1e-12,
) -> dict[str, object]:
    """Validate Cubit CAD mass-property sidecar rows before mesh routing.

    Cubit/Coreform can own the hex-led mesh route, but the CAD handoff should
    keep volume, summed surface area, and bounding-box dimensions as a small
    sidecar.  This replay gate makes those checks executable without reopening
    Cubit, and gives build123d/CST/CAD lanes a common volume/area/bbox contract.
    """

    records = [dict(row) for row in rows]
    relative_tolerance = float(rel_tol)
    absolute_tolerance = float(abs_tol)
    if relative_tolerance < 0.0:
        raise ValueError("rel_tol must be non-negative")
    if absolute_tolerance < 0.0:
        raise ValueError("abs_tol must be non-negative")
    if not records:
        raise ValueError("rows must not be empty")

    total_volume = 0.0
    total_area = 0.0
    volume_values: list[float] = []
    area_values: list[float] = []
    row_names: list[str] = []
    row_issues: list[dict[str, object]] = []
    bbox_size: list[float] | None = None

    for index, row in enumerate(records):
        name = str(row.get("name", "")).strip()
        row_names.append(name)
        try:
            volume = float(row.get("volume"))
            area = float(row.get("area"))
        except (TypeError, ValueError):
            row_issues.append({"index": index, "name": name, "reason": "missing volume or area"})
            continue
        volume_values.append(volume)
        area_values.append(area)
        total_volume += volume
        total_area += area
        bbox = row.get("bounding_box") or row.get("bbox")
        if isinstance(bbox, Mapping) and bbox.get("size") is not None:
            try:
                bbox_size = [float(value) for value in bbox["size"]]
            except (TypeError, ValueError):
                row_issues.append({"index": index, "name": name, "reason": "invalid bbox size"})

    expected_volume_ok = True
    expected_volume_rel_error = None
    if expected_total_volume is not None:
        expected_volume = float(expected_total_volume)
        expected_volume_rel_error = _relative_error(total_volume, expected_volume)
        expected_volume_ok = _within_tolerance(
            total_volume,
            expected_volume,
            rel_tol=relative_tolerance,
            abs_tol=absolute_tolerance,
        )

    expected_area_ok = True
    expected_area_rel_error = None
    if expected_total_area is not None:
        expected_area = float(expected_total_area)
        expected_area_rel_error = _relative_error(total_area, expected_area)
        expected_area_ok = _within_tolerance(
            total_area,
            expected_area,
            rel_tol=relative_tolerance,
            abs_tol=absolute_tolerance,
        )

    expected_bbox_ok = True
    bbox_size_abs_error = None
    expected_bbox_list = None
    if expected_bbox_size is not None:
        expected_bbox_list = [float(value) for value in expected_bbox_size]
        if len(expected_bbox_list) != 3:
            raise ValueError("expected_bbox_size must have three values")
        if bbox_size is None or len(bbox_size) != 3:
            expected_bbox_ok = False
        else:
            bbox_size_abs_error = max(abs(a - b) for a, b in zip(bbox_size, expected_bbox_list))
            expected_bbox_ok = bbox_size_abs_error <= absolute_tolerance

    checks = {
        "rows_present": bool(records),
        "row_names_recorded": all(bool(name) for name in row_names),
        "all_rows_have_volume_area": len(volume_values) == len(records) and len(area_values) == len(records),
        "all_values_finite": all(isfinite(value) for value in volume_values + area_values),
        "all_volumes_positive": all(value > 0.0 for value in volume_values),
        "all_areas_positive": all(value > 0.0 for value in area_values),
        "total_volume_expected_ok": expected_volume_ok,
        "total_area_expected_ok": expected_area_ok,
        "bbox_size_expected_ok": expected_bbox_ok,
    }
    return {
        "policy": "cubit_mass_property_sidecar_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "row_count": len(records),
        "row_names": row_names,
        "total_volume": total_volume,
        "total_area": total_area,
        "expected_total_volume": expected_total_volume,
        "expected_total_area": expected_total_area,
        "expected_bbox_size": expected_bbox_list,
        "bbox_size": bbox_size,
        "max_bbox_size_abs_error": bbox_size_abs_error,
        "volume_rel_error": expected_volume_rel_error,
        "area_rel_error": expected_area_rel_error,
        "rel_tol": relative_tolerance,
        "abs_tol": absolute_tolerance,
        "checks": checks,
        "issues": row_issues,
        "notes": [
            "Run this before mesh export or solver-ready promotion when Cubit owns the hex-led route.",
            "Volume alone is not enough for CAD handoff; carry summed surface area and bounding-box size when available.",
        ],
    }


def cubit_export_package_identity_gate(
    artifacts: Iterable[Mapping[str, object]],
    *,
    expected_export_id: str | None = None,
    expected_geometry_id: str | None = None,
    required_kinds: Iterable[str] = ("vol", "vol_sidecar", "raw_result"),
    expected_order: int | None = None,
    expected_routing_hint: str | None = None,
    inventory: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Check that Cubit export files travel as one identified package.

    Cubit-led validation often creates several files for the same mesh:
    ``.vol``, ``.vol.json``, raw batch JSON, and optional mass-property
    sidecars.  This gate prevents notebook or solver-ready promotion when those
    files have lost a stable export id, geometry id, order, or routing hint.
    """

    records = [dict(row) for row in artifacts]
    if not records:
        raise ValueError("artifacts must not be empty")

    required = [str(kind).strip() for kind in required_kinds if str(kind).strip()]
    kinds = [str(row.get("kind", "")).strip() for row in records]
    paths = [str(row.get("path", "")).strip() for row in records]
    export_ids = [str(row.get("export_id", "")).strip() for row in records]
    geometry_ids = [str(row.get("geometry_id", "")).strip() for row in records]
    orders: list[int] = []
    order_issues: list[dict[str, object]] = []
    for index, row in enumerate(records):
        if row.get("order") is None:
            continue
        try:
            orders.append(int(row["order"]))
        except (TypeError, ValueError):
            order_issues.append({"index": index, "kind": kinds[index], "reason": "invalid order"})

    vol_paths = [path for kind, path in zip(kinds, paths) if kind == "vol" and path]
    sidecar_paths = [path for kind, path in zip(kinds, paths) if kind in {"vol_sidecar", "sidecar"} and path]

    def normalized_path(value: str) -> str:
        return value.replace("/", "\\").rstrip("\\").lower()

    expected_sidecars = {normalized_path(f"{path}.json") for path in vol_paths}
    actual_sidecars = {normalized_path(path) for path in sidecar_paths}
    inv = dict(inventory or {})
    inv_source = str(inv.get("source", "")).strip()
    inv_routing = str(inv.get("routing_hint", "")).strip()
    expected_export = None if expected_export_id is None else str(expected_export_id).strip()
    expected_geometry = None if expected_geometry_id is None else str(expected_geometry_id).strip()
    expected_hint = None if expected_routing_hint is None else str(expected_routing_hint).strip()
    expected_order_value = None if expected_order is None else int(expected_order)

    checks = {
        "artifacts_present": bool(records),
        "kinds_recorded": all(bool(kind) for kind in kinds),
        "paths_recorded": all(bool(path) for path in paths),
        "required_kinds_present": all(kind in kinds for kind in required),
        "export_ids_recorded": all(bool(value) for value in export_ids),
        "export_id_unique": len(set(export_ids)) == 1,
        "export_id_matches_expected": expected_export is None or set(export_ids) == {expected_export},
        "geometry_ids_recorded": all(bool(value) for value in geometry_ids),
        "geometry_id_unique": len(set(geometry_ids)) == 1,
        "geometry_id_matches_expected": expected_geometry is None or set(geometry_ids) == {expected_geometry},
        "vol_path_recorded": bool(vol_paths),
        "vol_paths_have_vol_suffix": bool(vol_paths)
        and all(normalized_path(path).endswith(".vol") for path in vol_paths),
        "vol_sidecar_path_recorded": bool(sidecar_paths),
        "vol_sidecar_pairs_vol": bool(expected_sidecars) and expected_sidecars.issubset(actual_sidecars),
        "orders_valid": not order_issues,
        "order_matches_expected": expected_order_value is None or (
            bool(orders) and all(order == expected_order_value for order in orders)
        ),
        "inventory_source_matches_vol": not inv_source
        or normalized_path(inv_source) in {normalized_path(path) for path in vol_paths},
        "inventory_routing_hint_matches_expected": expected_hint is None or inv_routing == expected_hint,
    }
    return {
        "policy": "cubit_export_package_identity_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "artifact_count": len(records),
        "kinds": kinds,
        "required_kinds": required,
        "paths": paths,
        "export_ids": sorted(set(export_ids)),
        "geometry_ids": sorted(set(geometry_ids)),
        "expected_export_id": expected_export,
        "expected_geometry_id": expected_geometry,
        "orders": sorted(set(orders)),
        "expected_order": expected_order_value,
        "expected_routing_hint": expected_hint,
        "inventory_source": inv_source,
        "inventory_routing_hint": inv_routing,
        "expected_vol_sidecars": sorted(expected_sidecars),
        "actual_vol_sidecars": sorted(actual_sidecars),
        "order_issues": order_issues,
        "checks": checks,
        "notes": [
            "Use this before docs/panel notebooks or solver-ready runs consume a Cubit export package.",
            "A .vol, .vol.json, raw result, and mass-property sidecar should share export_id and geometry_id.",
        ],
    }


def summarize_netgen_vol_inventory(text: str, source: str | None = None) -> dict[str, object]:
    """Return mixed-element inventory for a Netgen ``.vol`` text.

    The result is a routing preflight, not a solver parser.  It detects whether
    a file belongs to the Netgen tri/tet-only teaching path or the Cubit
    hex/mixed path, and it deliberately refuses to split or reinterpret element
    types.
    """

    lines = text.splitlines()
    surface_section = _first_existing_section(lines, ("surfaceelements", "surfaceelementsuv"))
    surface_rows = (
        _read_counted_section(lines, surface_section, required=False)
        if surface_section is not None else []
    )
    volume_rows = _read_counted_section(lines, "volumeelements", required=False)
    point_rows = _read_counted_section(lines, "points", required=False)
    material_rows = _read_counted_section(lines, "materials", required=False)
    boundary_rows = _read_counted_section(lines, "bcnames", required=False)
    curvedelements_present = _has_section(lines, "curvedelements")

    surface_kind_counts = _count_by_np(surface_rows, 4, SURFACE_KIND_BY_NP)
    volume_kind_counts = _count_by_np(volume_rows, 1, VOLUME_KIND_BY_NP)
    materials = _parse_materials(material_rows)
    boundary_names = _parse_named_rows(boundary_rows)
    is_tri_tet_only = (
        set(surface_kind_counts).issubset({"triangle"})
        and set(volume_kind_counts).issubset({"tet"})
        and bool(volume_rows)
    )
    has_mixed_hex_transition = any(
        volume_kind_counts.get(kind, 0) > 0 for kind in ("hex", "pyramid", "wedge")
    )

    routing_hint = (
        "netgen_tri_tet_path"
        if is_tri_tet_only
        else "cubit_hex_or_mixed_path"
        if has_mixed_hex_transition
        else "inspect_before_solver_import"
    )

    return {
        "source": source,
        "surface_section": surface_section,
        "surface_elements": len(surface_rows),
        "surface_kind_counts": surface_kind_counts,
        "volume_elements": len(volume_rows),
        "volume_kind_counts": volume_kind_counts,
        "points": len(point_rows),
        "materials": materials,
        "boundary_names": boundary_names,
        "curvedelements_present": curvedelements_present,
        "is_tri_tet_only": is_tri_tet_only,
        "has_mixed_hex_transition": has_mixed_hex_transition,
        "routing_hint": routing_hint,
        "policy": (
            "Cubit/Coreform owns hex-led and mixed hex+pyramid+tet inventory; "
            "Netgen/OCC owns tet-only generation for the first-order education path. "
            "Cubit stock Netgen exports may use surfaceelementsuv; count it as "
            "the surface inventory instead of treating the boundary as empty."
        ),
        "order_series_policy": (
            "Route high-order Cubit .vol files from the first-order element arity "
            "inventory. The curvedelements section can grow with order without "
            "changing whether the mesh is tri/tet-only or hex-led mixed."
        ),
    }


def _next_payload_line(lines: list[str], start: int) -> tuple[int, str]:
    for i in range(start, len(lines)):
        line = lines[i].strip()
        if line and not line.startswith("#"):
            return i, line
    raise ValueError("unexpected end of .vol file")


def _read_counted_section(lines: list[str], name: str, *, required: bool) -> list[str]:
    for i, line in enumerate(lines):
        if line.strip().lower() == name.lower():
            count_i, count_line = _next_payload_line(lines, i + 1)
            count = int(count_line.split()[0])
            rows: list[str] = []
            j = count_i + 1
            while len(rows) < count:
                if j >= len(lines):
                    raise ValueError(f"section {name!r} ended early")
                row = lines[j].strip()
                if row and not row.startswith("#"):
                    rows.append(row)
                j += 1
            return rows
    if required:
        raise ValueError(f"section {name!r} not found")
    return []


def _has_section(lines: Iterable[str], name: str) -> bool:
    needle = name.lower()
    return any(line.strip().lower() == needle for line in lines)


def _relative_error(observed: float, expected: float) -> float:
    return abs(observed - expected) / max(abs(expected), abs(observed), 1.0)


def _within_tolerance(
    observed: float,
    expected: float,
    *,
    rel_tol: float,
    abs_tol: float,
) -> bool:
    return abs(observed - expected) <= max(abs_tol, rel_tol * max(abs(expected), 1.0))


def _first_existing_section(lines: Iterable[str], names: Iterable[str]) -> str | None:
    lowered = {line.strip().lower() for line in lines}
    for name in names:
        if name.lower() in lowered:
            return name
    return None


def _count_by_np(rows: Iterable[str], np_column: int, kind_by_np: dict[int, str]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for row in rows:
        parts = row.split()
        np_value = int(parts[np_column])
        counts[kind_by_np.get(np_value, f"np{np_value}")] += 1
    return dict(sorted(counts.items()))


def _parse_materials(rows: Iterable[str]) -> dict[int, str]:
    materials: dict[int, str] = {}
    for row in rows:
        parts = row.split(maxsplit=1)
        if not parts:
            continue
        material_id = int(parts[0])
        materials[material_id] = parts[1] if len(parts) > 1 else f"material_{material_id}"
    return materials


def _parse_named_rows(rows: Iterable[str]) -> dict[int, str]:
    names: dict[int, str] = {}
    for row in rows:
        parts = row.split(maxsplit=1)
        if not parts:
            continue
        item_id = int(parts[0])
        names[item_id] = parts[1] if len(parts) > 1 else f"name_{item_id}"
    return names


def _linear_quantile(sorted_values: list[float], level: float) -> float | None:
    if not sorted_values:
        return None
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = (len(sorted_values) - 1) * level
    lower_index = int(position)
    upper_index = min(lower_index + 1, len(sorted_values) - 1)
    weight = position - lower_index
    return (
        sorted_values[lower_index] * (1.0 - weight)
        + sorted_values[upper_index] * weight
    )


def _quantile_key(level: float) -> str:
    percent = level * 100.0
    rounded = round(percent)
    if abs(percent - rounded) <= 1e-9:
        return f"p{int(rounded):02d}"
    return f"q{level:g}"


def _quality_histogram(values: list[float], edges: list[float]) -> dict[str, object]:
    bins = [
        {"lo": left, "hi": right, "count": 0}
        for left, right in zip(edges, edges[1:])
    ]
    underflow = 0
    overflow = 0
    for value in values:
        if value < edges[0]:
            underflow += 1
            continue
        if value > edges[-1]:
            overflow += 1
            continue
        for index, bin_record in enumerate(bins):
            lo = float(bin_record["lo"])
            hi = float(bin_record["hi"])
            is_last = index == len(bins) - 1
            if lo <= value <= hi if is_last else lo <= value < hi:
                bin_record["count"] = int(bin_record["count"]) + 1
                break
    return {
        "edges": edges,
        "bins": bins,
        "underflow": underflow,
        "overflow": overflow,
    }
