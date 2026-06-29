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
from typing import Iterable


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
