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
from datetime import datetime
import math
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


_HEADLESS_STARTUP_DIAGNOSTIC_SUFFIXES = (
    "/plugins",
    "-commandplugindir",
    "-nojournal",
)


def _headless_process_evidence(summary: Mapping[str, object]) -> dict[str, object]:
    """Classify Cubit's known startup diagnostics without hiding script errors."""
    diagnostics = [str(value) for value in summary.get("startup_diagnostics", [])]
    script_errors = [str(value) for value in summary.get("script_error_lines", [])]
    normalized_suffixes = [
        next(
            (
                suffix
                for suffix in _HEADLESS_STARTUP_DIAGNOSTIC_SUFFIXES
                if row.rstrip().endswith(suffix)
            ),
            "",
        )
        for row in diagnostics
    ]
    startup_only_allowlisted = (
        bool(diagnostics)
        and len(diagnostics) <= len(_HEADLESS_STARTUP_DIAGNOSTIC_SUFFIXES)
        and all("Could not open file:" in row for row in diagnostics)
        and all(normalized_suffixes)
        and len(set(normalized_suffixes)) == len(normalized_suffixes)
    )
    process_exit_code = int(summary.get("process_exit_code", -1))
    result_fresh = summary.get("result_artifact_fresh") is True
    process_exit_acceptable = process_exit_code == 0 or (
        process_exit_code in {2, 3}
        and startup_only_allowlisted
        and not script_errors
        and result_fresh
    )
    return {
        "diagnostics": diagnostics,
        "script_errors": script_errors,
        "startup_only_allowlisted": startup_only_allowlisted,
        "process_exit_code": process_exit_code,
        "result_fresh": result_fresh,
        "process_exit_acceptable": process_exit_acceptable,
        "launcher_classification": (
            "clean_exit"
            if process_exit_code == 0
            else "allowlisted_startup_diagnostic_with_clean_script"
            if process_exit_acceptable
            else "execution_error"
        ),
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


def cubit_hex_geometry_refinement_gate(
    rows: Iterable[Mapping[str, object]],
    *,
    min_quality: float = 0.2,
    max_final_volume_relative_error: float = 0.05,
    min_error_reduction_fraction: float = 0.25,
) -> dict[str, object]:
    """Separate hex-count refinement from curved-boundary convergence.

    ``volume ... size`` can add interior hexes without changing the faceted
    boundary that controls a curved solid's mesh volume.  This replay gate
    therefore requires both a real refinement attempt and an improving
    geometry observable before it reports convergence.
    """

    records = [dict(row) for row in rows]
    quality_limit = float(min_quality)
    final_error_limit = float(max_final_volume_relative_error)
    reduction_limit = float(min_error_reduction_fraction)
    if len(records) < 2:
        raise ValueError("rows must contain at least two refinement levels")
    if quality_limit <= 0.0:
        raise ValueError("min_quality must be > 0")
    if final_error_limit < 0.0:
        raise ValueError("max_final_volume_relative_error must be >= 0")
    if reduction_limit < 0.0 or reduction_limit > 1.0:
        raise ValueError("min_error_reduction_fraction must be in [0, 1]")

    required = {
        "target_size",
        "hex_count",
        "mesh_volume",
        "expected_volume",
        "volume_relative_error",
        "quality_min",
        "quality_q05",
    }
    normalized = []
    for index, row in enumerate(records):
        missing = sorted(required - row.keys())
        if missing:
            raise ValueError(f"row {index} is missing: {', '.join(missing)}")
        item = {
            "target_size": float(row["target_size"]),
            "hex_count": int(row["hex_count"]),
            "tet_count": int(row.get("tet_count", 0)),
            "pyramid_count": int(row.get("pyramid_count", 0)),
            "wedge_count": int(row.get("wedge_count", 0)),
            "mesh_volume": float(row["mesh_volume"]),
            "expected_volume": float(row["expected_volume"]),
            "volume_relative_error": float(row["volume_relative_error"]),
            "quality_min": float(row["quality_min"]),
            "quality_q05": float(row["quality_q05"]),
        }
        if not all(isfinite(value) for key, value in item.items() if key not in {
            "hex_count", "tet_count", "pyramid_count", "wedge_count"
        }):
            raise ValueError(f"row {index} contains a non-finite value")
        if item["target_size"] <= 0.0 or item["hex_count"] <= 0:
            raise ValueError(f"row {index} target_size and hex_count must be positive")
        if item["expected_volume"] <= 0.0 or item["volume_relative_error"] < 0.0:
            raise ValueError(f"row {index} has invalid volume evidence")
        normalized.append(item)

    sizes = [row["target_size"] for row in normalized]
    counts = [row["hex_count"] for row in normalized]
    errors = [row["volume_relative_error"] for row in normalized]
    first_error = errors[0]
    final_error = errors[-1]
    error_reduction_fraction = (
        (first_error - final_error) / first_error if first_error > 0.0
        else (1.0 if final_error == 0.0 else 0.0)
    )
    if abs(error_reduction_fraction) < 1.0e-12:
        error_reduction_fraction = 0.0
    count_increased = any(right > left for left, right in zip(counts, counts[1:]))
    checks = {
        "target_sizes_decrease": all(right < left for left, right in zip(sizes, sizes[1:])),
        "hex_counts_nondecreasing": all(right >= left for left, right in zip(counts, counts[1:])),
        "hex_count_increased": count_increased,
        "hex_only": all(
            row["tet_count"] == row["pyramid_count"] == row["wedge_count"] == 0
            for row in normalized
        ),
        "quality_ok": all(
            min(row["quality_min"], row["quality_q05"]) >= quality_limit
            for row in normalized
        ),
        "geometry_error_reduced": error_reduction_fraction >= reduction_limit,
        "final_geometry_error_ok": final_error <= final_error_limit,
    }
    refinement_attempted = (
        checks["target_sizes_decrease"]
        and checks["hex_counts_nondecreasing"]
        and checks["hex_count_increased"]
    )
    plateau_detected = (
        refinement_attempted
        and checks["hex_only"]
        and checks["quality_ok"]
        and not checks["geometry_error_reduced"]
    )
    if plateau_detected:
        status = "needs_geometry_refinement"
    elif all(checks.values()):
        status = "ok"
    else:
        status = "needs_attention"

    return {
        "policy": "cubit_hex_geometry_refinement_plateau_gate_v1",
        "status": status,
        "level_count": len(normalized),
        "rows": normalized,
        "thresholds": {
            "min_quality": quality_limit,
            "max_final_volume_relative_error": final_error_limit,
            "min_error_reduction_fraction": reduction_limit,
        },
        "checks": checks,
        "refinement_attempted": refinement_attempted,
        "plateau_detected": plateau_detected,
        "initial_volume_relative_error": first_error,
        "final_volume_relative_error": final_error,
        "error_reduction_fraction": error_reduction_fraction,
        "recommendation": (
            "Refine the boundary topology or export curved/high-order geometry; "
            "interior hex-count growth alone is not evidence of curved-volume "
            "convergence."
            if plateau_detected
            else "Keep the geometry observable and quality distribution in the refinement gate."
        ),
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


def cubit_release_feature_routing_gate(
    features: Iterable[Mapping[str, object]],
    *,
    release_version: str,
    source_url: str,
    lab_cubit_role: str = "hex_led_or_mixed",
    tet_only_owner: str = "netgen_tri_tet_path",
    required_features: Iterable[str] = (
        "anisotropic_tetrahedral_meshing",
        "cohesive_element_generation",
        "higher_order_quality_metrics",
        "tri_tet_meshing_robustness",
        "sculpt_refinement_memory",
        "solver_io_compatibility",
        "python_312_runtime",
    ),
    require_headless_policy: bool = True,
) -> dict[str, object]:
    """Route public Coreform/Cubit release features into lab validation lanes.

    Release notes are useful only after they become executable lab policy.  This
    gate records which features were learned, where they should be tested, and
    whether they preserve the CAE-AI Lab split: Cubit owns hex-led/mixed mesh
    work while tet-only educational meshes stay on the Netgen/OCC path.
    """

    rows = [dict(row) for row in features]
    if not rows:
        raise ValueError("features must not be empty")

    def feature_key(value: object) -> str:
        return "_".join(
            "".join(ch.lower() if ch.isalnum() else " " for ch in str(value)).split()
        )

    version = str(release_version or "").strip()
    url = str(source_url or "").strip()
    cubit_role = feature_key(lab_cubit_role)
    tet_owner = feature_key(tet_only_owner)
    required = [feature_key(item) for item in required_features if str(item).strip()]
    feature_keys = [
        feature_key(row.get("feature_key", row.get("name", row.get("feature", ""))))
        for row in rows
    ]
    categories = [str(row.get("category", "")).strip() for row in rows]
    routes = [str(row.get("lab_route", row.get("route", ""))).strip() for row in rows]
    validation_notes = [
        str(row.get("validation_note", row.get("validation", ""))).strip()
        for row in rows
    ]
    missing = [item for item in required if item not in feature_keys]
    checks = {
        "release_version_recorded": bool(version),
        "source_url_recorded": bool(url),
        "source_url_is_coreform": "coreform.com" in url.lower(),
        "features_recorded": bool(feature_keys) and all(bool(item) for item in feature_keys),
        "categories_recorded": all(bool(item) for item in categories),
        "lab_routes_recorded": all(bool(item) for item in routes),
        "validation_notes_recorded": all(bool(item) for item in validation_notes),
        "required_features_present": not missing,
        "cubit_role_is_hex_or_mixed": cubit_role in {"hex_led_or_mixed", "cubit_hex_or_mixed_path"},
        "tet_only_owner_is_netgen": tet_owner in {"netgen_tri_tet_path", "netgen_occ", "netgen"},
        "headless_policy_recorded": bool(require_headless_policy),
    }
    return {
        "policy": "cubit_release_feature_routing_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "release_version": version,
        "source_url": url,
        "feature_keys": feature_keys,
        "required_features": required,
        "missing_features": missing,
        "lab_cubit_role": lab_cubit_role,
        "tet_only_owner": tet_only_owner,
        "feature_count": len(rows),
        "features": rows,
        "checks": checks,
        "notes": [
            "Treat release-note learning as a routing contract before it becomes a solver validation claim.",
            "Cubit remains the hex-led and mixed hex+pyramid+tet lane; tet-only education remains on Netgen/OCC unless a slot explicitly says otherwise.",
            "Run Coreform/Cubit from headless batch scripts in this lab loop; GUI-facing improvements are documentation, not automation evidence.",
        ],
    }


def cubit_headless_installation_route_gate(
    installation: Mapping[str, object],
    *,
    required_headless_flags: Iterable[str] = ("-nographics", "-batch"),
) -> dict[str, object]:
    """Check that Cubit execution evidence matches the installed headless lane.

    Release-note learning and local execution evidence should not be conflated.
    This replay gate records the installed version, binary discovery, required
    headless flags, and any newer release-note watchlist version in one package
    so a validation slot cannot claim live evidence from a release that is not
    installed on the host.
    """

    record = dict(installation)
    version = str(record.get("installed_version", record.get("version", ""))).strip()
    binary_path = str(record.get("binary_path", "")).strip()
    binary_exists = record.get("binary_exists")
    flags_raw = record.get("headless_flags", record.get("required_flags", ())) or ()
    if isinstance(flags_raw, (str, bytes)):
        flags = [str(flags_raw)]
    else:
        flags = [str(flag) for flag in flags_raw]
    normalized_flags = {flag.strip().lower() for flag in flags if str(flag).strip()}
    required_flags = [
        str(flag).strip().lower() for flag in required_headless_flags if str(flag).strip()
    ]
    command_line = str(record.get("command_line", "")).strip().lower()
    gui_policy = str(record.get("gui_policy", record.get("gui_daemon_policy", ""))).strip().lower()
    allow_gui_daemon = bool(record.get("allow_gui_daemon", False))
    release_note_version = str(record.get("release_note_version", "")).strip()
    release_note_status = str(record.get("release_note_status", "")).strip().lower()
    live_claim_version = str(
        record.get("live_claimed_release_version", record.get("live_claim_version", version))
    ).strip()
    license_status = str(record.get("license_status", "")).strip()
    license_status_norm = license_status.lower().replace(" ", "")
    version_probe_command = str(record.get("version_probe_command", "")).strip()
    version_probe_command_norm = version_probe_command.lower().replace("\\", "/")
    version_probe_summary = record.get("version_probe_summary", {})
    if version_probe_summary is None:
        version_probe_summary = {}
    if not isinstance(version_probe_summary, Mapping):
        raise ValueError("version_probe_summary must be a mapping when present")
    version_probe_summary_text = " ".join(
        str(value) for value in version_probe_summary.values()
        if value is not None
    )
    binary_path_norm = binary_path.lower().replace("\\", "/")
    binary_is_console = binary_path_norm.endswith("coreform_cubit.com")

    has_flag = {
        flag: flag in normalized_flags or flag in command_line
        for flag in required_flags
    }
    release_note_matches_install = (
        bool(release_note_version) and bool(version) and release_note_version == version
    )
    release_note_is_watchlist = (
        not release_note_version
        or release_note_matches_install
        or release_note_status in {"watchlist", "not_installed", "documentation_only"}
    )
    checks = {
        "installed_version_recorded": bool(version),
        "binary_path_recorded": bool(binary_path),
        "binary_path_is_console_com": binary_is_console,
        "binary_exists_recorded": isinstance(binary_exists, bool),
        "binary_exists": binary_exists is True,
        "headless_flags_recorded": bool(flags) or bool(command_line),
        "required_headless_flags_present": all(has_flag.values()) if required_flags else True,
        "gui_daemon_disabled_by_default": (not allow_gui_daemon)
        and ("no_gui" in gui_policy or "headless" in gui_policy or "disabled" in gui_policy),
        "live_claim_matches_installed_version": bool(version) and live_claim_version == version,
        "release_note_watchlist_not_live_claim": release_note_is_watchlist,
    }
    if "license_status" in record:
        checks["license_status_recorded"] = bool(license_status)
        checks["license_status_allows_headless_probe"] = (
            license_status_norm.startswith("valid")
            and "expired" not in license_status_norm
            and "noavailableseats" not in license_status_norm
        )
    if version_probe_command:
        checks["version_probe_is_synchronous_console"] = (
            "coreform_cubit.com" in version_probe_command_norm
            and "-version" in version_probe_command_norm
        )
        checks["version_probe_uses_recorded_binary"] = (
            bool(binary_path_norm) and binary_path_norm in version_probe_command_norm
        )
    if version_probe_summary:
        checks["version_probe_summary_records_installed_version"] = (
            bool(version) and version in version_probe_summary_text
        )
        checks["version_probe_summary_records_license_status"] = (
            not license_status or license_status in version_probe_summary_text
        )
    return {
        "policy": "cubit_headless_installation_route_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "installed_version": version,
        "binary_path": binary_path,
        "binary_exists": binary_exists,
        "headless_flags": flags,
        "required_headless_flags": required_flags,
        "license_status": license_status,
        "version_probe_command": version_probe_command,
        "version_probe_summary": dict(version_probe_summary),
        "release_note_version": release_note_version,
        "release_note_status": release_note_status,
        "live_claimed_release_version": live_claim_version,
        "checks": checks,
        "notes": [
            "Use installed-version evidence for live headless Cubit claims.",
            "Treat newer release notes as a watchlist until that release is installed and replayed.",
            "Keep GUI daemon usage out of the default validation lane.",
            "Record synchronous .com -version probes and license status separately from batch-mesh evidence.",
            "Use coreform_cubit.com for console evidence; coreform_cubit.exe can be a GUI stub and is not enough for the live headless lane.",
        ],
    }


def cubit_curvilinear_handoff_manifest_gate(
    manifest: Mapping[str, object],
    *,
    expected_mesh_id: str | None = None,
    expected_export_id: str | None = None,
    min_order: int = 2,
    min_scaled_jacobian: float = 0.2,
) -> dict[str, object]:
    """Check a Cubit high-order curvilinear mesh handoff manifest.

    The lesson from third-party curvilinear mesh workflows is that the mesh
    order alone is not enough.  A reusable handoff needs the imported/source
    mesh identity, CAD or geometry association, projection-error evidence,
    curved export order, quality metric, negative-Jacobian count, and routing
    policy in one package before downstream FEM/BEM code treats the file as
    solver-ready.
    """

    if not isinstance(manifest, Mapping):
        raise ValueError("manifest must be a mapping")
    min_order_int = int(min_order)
    min_quality = float(min_scaled_jacobian)
    if min_order_int < 2:
        raise ValueError("min_order must be at least 2 for a curvilinear handoff")
    if min_quality <= 0.0:
        raise ValueError("min_scaled_jacobian must be > 0")

    mesh_id = str(manifest.get("mesh_id", "")).strip()
    export_id = str(manifest.get("export_id", "")).strip()
    source_mesh = dict(manifest.get("source_mesh", {}) or {})
    geometry = dict(manifest.get("geometry_association", {}) or {})
    curved_export = dict(manifest.get("curved_export", {}) or {})
    quality = dict(manifest.get("quality", {}) or {})
    provenance = dict(manifest.get("provenance", {}) or {})

    source_kind = str(source_mesh.get("kind", "")).strip().lower().replace("-", "_")
    volume_kinds = {
        str(kind).strip().lower()
        for kind in source_mesh.get("volume_kinds", source_mesh.get("volume_element_kinds", []))
        if str(kind).strip()
    }
    surface_kinds = {
        str(kind).strip().lower()
        for kind in source_mesh.get("surface_kinds", source_mesh.get("surface_element_kinds", []))
        if str(kind).strip()
    }
    order = int(curved_export.get("order", 0) or 0)
    export_format = str(curved_export.get("format", "")).strip().lower()
    routing_hint = str(curved_export.get("routing_hint", manifest.get("routing_hint", ""))).strip().lower()
    implicit_conversion = bool(
        curved_export.get("implicit_element_conversion", curved_export.get("implicit_tetization", False))
    )
    quality_metric = str(quality.get("metric", "")).strip().lower().replace(" ", "_")
    quality_min = float(quality.get("min", quality.get("min_scaled_jacobian", float("nan"))))
    quality_count = int(quality.get("count", 0) or 0)
    negative_jacobian_raw = quality.get(
        "negative_jacobian_count",
        quality.get("negative_jacobians", quality.get("negative_count")),
    )
    negative_jacobian_count: int | None = None
    if negative_jacobian_raw is not None:
        try:
            negative_jacobian_count = int(negative_jacobian_raw)
        except (TypeError, ValueError):
            negative_jacobian_count = None
    geometry_policy = str(
        geometry.get("projection_policy", geometry.get("association_policy", geometry.get("policy", "")))
    ).strip()
    cad_source = str(geometry.get("cad_source", geometry.get("geometry_source", ""))).strip()
    projection_raw = geometry.get(
        "projection_quality",
        geometry.get("boundary_conformity", manifest.get("projection_quality", {})),
    ) or {}
    projection_quality = dict(projection_raw) if isinstance(projection_raw, Mapping) else {}

    def finite_float(value: object) -> float:
        try:
            return float(value)
        except (TypeError, ValueError):
            return float("nan")

    projection_distance = finite_float(
        projection_quality.get(
            "max_distance",
            projection_quality.get("max_cad_distance", projection_quality.get("distance_max")),
        )
    )
    projection_tolerance = finite_float(
        projection_quality.get(
            "tolerance",
            projection_quality.get("distance_tolerance", projection_quality.get("cad_distance_tolerance")),
        )
    )
    literature_note = str(
        provenance.get("literature_note", manifest.get("literature_note", ""))
    ).strip()

    expected_mesh_ok = True if expected_mesh_id is None else mesh_id == str(expected_mesh_id)
    expected_export_ok = True if expected_export_id is None else export_id == str(expected_export_id)
    hex_or_mixed = bool(volume_kinds.intersection({"hex", "pyramid", "wedge"}))
    tet_only = volume_kinds == {"tet"}
    checks = {
        "mesh_id_recorded": bool(mesh_id),
        "export_id_recorded": bool(export_id),
        "expected_mesh_id_ok": expected_mesh_ok,
        "expected_export_id_ok": expected_export_ok,
        "source_mesh_is_imported": source_kind in {"third_party_mesh", "imported_mesh", "external_mesh"},
        "volume_kinds_recorded": bool(volume_kinds),
        "surface_kinds_recorded": bool(surface_kinds),
        "hex_or_mixed_cubit_route": hex_or_mixed and not tet_only,
        "geometry_association_recorded": bool(cad_source) and bool(geometry_policy),
        "boundary_ids_preserved": bool(geometry.get("boundary_ids_preserved", False)),
        "projection_error_recorded": isfinite(projection_distance) and isfinite(projection_tolerance),
        "projection_error_within_tolerance": (
            isfinite(projection_distance)
            and isfinite(projection_tolerance)
            and projection_distance <= projection_tolerance
        ),
        "curved_export_order_ok": order >= min_order_int,
        "curved_export_format_recorded": bool(export_format),
        "routing_hint_is_cubit_hex_or_mixed": routing_hint == "cubit_hex_or_mixed_path",
        "no_implicit_element_conversion": not implicit_conversion,
        "quality_metric_is_scaled_jacobian": quality_metric in {"scaled_jacobian", "jacobian"},
        "quality_count_positive": quality_count > 0,
        "quality_min_ok": isfinite(quality_min) and quality_min >= min_quality,
        "negative_jacobian_count_recorded": negative_jacobian_count is not None,
        "negative_jacobian_count_zero": negative_jacobian_count == 0,
        "literature_context_recorded": bool(literature_note),
    }
    return {
        "policy": "cubit_curvilinear_handoff_manifest_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "mesh_id": mesh_id,
        "export_id": export_id,
        "source_kind": source_kind,
        "volume_kinds": sorted(volume_kinds),
        "surface_kinds": sorted(surface_kinds),
        "cad_source": cad_source,
        "geometry_policy": geometry_policy,
        "order": order,
        "format": export_format,
        "routing_hint": routing_hint,
        "quality_metric": quality_metric,
        "quality_min": quality_min,
        "quality_count": quality_count,
        "negative_jacobian_count": negative_jacobian_count,
        "projection_distance": projection_distance,
        "projection_tolerance": projection_tolerance,
        "min_order": min_order_int,
        "min_scaled_jacobian": min_quality,
        "checks": checks,
        "notes": [
            "Do not relax the first-order tri/tet .vol parser for Cubit high-order hex or mixed meshes.",
            "Keep CAD/geometry association and curved export order together before solver-ready promotion.",
            "Record CAD projection error and zero negative-Jacobian count before treating a curved mesh as solver-ready.",
            "Use Cubit for hex-led or mixed curvilinear handoff; tet-only education remains on the Netgen/OCC lane.",
        ],
    }


def cubit_mixed_solver_route_manifest_gate(
    inventory: Mapping[str, object],
    manifest: Mapping[str, object],
    *,
    expected_package_id: str | None = None,
    expected_solver_contract_artifact_id: str | None = None,
    expected_solver_contract_digest: str | None = None,
    expected_solver_contract_path: str | None = None,
    expected_solver_route_convention_schema_id: str | None = None,
    require_solver_contract_artifact: bool = False,
    require_solver_route_convention_schema: bool = False,
    expected_routing_hint: str = "cubit_hex_or_mixed_path",
    required_volume_kinds: Iterable[str] = ("hex", "pyramid", "tet"),
    required_surface_kinds: Iterable[str] = ("quad", "triangle"),
) -> dict[str, object]:
    """Check the downstream solver route for a mixed Cubit mesh package.

    Cubit is the hex-led route in the CAE-AI Lab.  A mixed
    hex+pyramid+tet package is not solver-ready until the route manifest says
    which element families are primary volume cells, which are transition
    bridge cells, and which cells are compatibility/subregion cells.
    """

    inv = dict(inventory)
    spec = dict(manifest)

    def norm(value: object) -> str:
        return str(value or "").strip().lower().replace(" ", "_").replace("-", "_")

    def norm_list(values: Iterable[object]) -> list[str]:
        return [text for value in values if (text := norm(value))]

    volume_counts_raw = inv.get("volume_kind_counts", {})
    surface_counts_raw = inv.get("surface_kind_counts", {})
    if not isinstance(volume_counts_raw, Mapping):
        raise ValueError("inventory['volume_kind_counts'] must be a mapping")
    if not isinstance(surface_counts_raw, Mapping):
        raise ValueError("inventory['surface_kind_counts'] must be a mapping")
    volume_counts = {norm(key): int(value) for key, value in volume_counts_raw.items()}
    surface_counts = {norm(key): int(value) for key, value in surface_counts_raw.items()}
    active_volume_kinds = sorted(kind for kind, count in volume_counts.items() if count > 0)
    active_surface_kinds = sorted(kind for kind, count in surface_counts.items() if count > 0)

    volume_routes_raw = spec.get("volume_routes", [])
    surface_routes_raw = spec.get("surface_routes", [])
    if not isinstance(volume_routes_raw, Iterable) or isinstance(volume_routes_raw, (str, bytes)):
        raise ValueError("manifest['volume_routes'] must be an iterable of mappings")
    if not isinstance(surface_routes_raw, Iterable) or isinstance(surface_routes_raw, (str, bytes)):
        raise ValueError("manifest['surface_routes'] must be an iterable of mappings")
    volume_routes = [dict(row) for row in volume_routes_raw]
    surface_routes = [dict(row) for row in surface_routes_raw]

    route_volume_kinds = sorted({
        kind
        for row in volume_routes
        if (kind := norm(row.get("volume_kind") or row.get("kind")))
    })
    route_surface_kinds = sorted({
        kind
        for row in surface_routes
        if (kind := norm(row.get("surface_kind") or row.get("kind")))
    })
    required_volumes = norm_list(required_volume_kinds)
    required_surfaces = norm_list(required_surface_kinds)
    package_id = str(spec.get("solver_route_package_id", spec.get("package_id", ""))).strip()
    expected_id = None if expected_package_id is None else str(expected_package_id).strip()
    routing_hint = str(spec.get("routing_hint", inv.get("routing_hint", ""))).strip()
    expected_hint = str(expected_routing_hint).strip()
    route_policy = str(spec.get("route_policy", "")).strip()
    downstream_solver = str(spec.get("downstream_solver", spec.get("solver", ""))).strip()
    convention_aliases = (
        "solver_route_convention_schema_id",
        "solverRouteConventionSchemaId",
        "route_convention_schema_id",
        "routeConventionSchemaId",
        "mesh_route_convention_schema_id",
        "meshRouteConventionSchemaId",
        "mixed_solver_route_convention_schema_id",
        "mixedSolverRouteConventionSchemaId",
    )
    solver_route_convention_schema_ids = sorted({
        str(spec[key]).strip()
        for key in convention_aliases
        if key in spec and str(spec[key]).strip()
    })
    solver_route_convention_schema_id = (
        solver_route_convention_schema_ids[0]
        if solver_route_convention_schema_ids
        else None
    )
    solver_contract_artifact_id = str(
        spec.get(
            "solver_contract_artifact_id",
            spec.get(
                "downstream_solver_contract_artifact_id",
                spec.get("solver_reader_contract_artifact_id", ""),
            ),
        )
    ).strip()
    solver_contract_digest = str(
        spec.get(
            "solver_contract_digest",
            spec.get(
                "downstream_solver_contract_digest",
                spec.get("solver_reader_contract_digest", ""),
            ),
        )
    ).strip()
    solver_contract_path = str(
        spec.get(
            "solver_contract_path",
            spec.get(
                "downstream_solver_contract_path",
                spec.get("solver_reader_contract_path", ""),
            ),
        )
    ).strip()
    expected_contract_id = (
        None
        if expected_solver_contract_artifact_id is None
        else str(expected_solver_contract_artifact_id).strip()
    )
    expected_contract_digest = (
        None if expected_solver_contract_digest is None else str(expected_solver_contract_digest).strip()
    )
    expected_contract_path = (
        None if expected_solver_contract_path is None else str(expected_solver_contract_path).strip()
    )
    expected_route_convention_schema = (
        None
        if expected_solver_route_convention_schema_id is None
        else str(expected_solver_route_convention_schema_id).strip()
    )
    route_convention_schema_required = bool(
        require_solver_route_convention_schema
        or expected_route_convention_schema is not None
    )
    tet_only_owner = norm(spec.get("tet_only_owner", ""))
    no_implicit_tetization = bool(
        spec.get("no_implicit_tetization", spec.get("no_implicit_element_conversion", False))
    )

    def has_role(kind: str, role_words: Iterable[str], *, not_primary: bool | None = None) -> bool:
        role_needles = [norm(word) for word in role_words]
        for row in volume_routes:
            if norm(row.get("volume_kind") or row.get("kind")) != kind:
                continue
            role = norm(row.get("solver_role", row.get("role", "")))
            route = norm(row.get("route", ""))
            haystack = f"{role} {route}"
            if not all(word in haystack for word in role_needles):
                continue
            if not_primary is not None and bool(row.get("not_primary_region", False)) is not not_primary:
                continue
            return True
        return False

    checks = {
        "inventory_is_mixed_hex_pyramid_tet": (
            volume_counts.get("hex", 0) > 0
            and volume_counts.get("pyramid", 0) > 0
            and volume_counts.get("tet", 0) > 0
        ),
        "routing_hint_matches_expected": routing_hint == expected_hint,
        "package_id_recorded": bool(package_id),
        "expected_package_id_matches": expected_id is None or package_id == expected_id,
        "route_policy_recorded": bool(route_policy),
        "downstream_solver_recorded": bool(downstream_solver),
        "solver_route_convention_schema_id_consistent_when_present": (
            len(solver_route_convention_schema_ids) <= 1
        ),
        "solver_route_convention_schema_id_recorded_when_required": (
            not route_convention_schema_required
            or bool(solver_route_convention_schema_id)
        ),
        "solver_route_convention_schema_id_recorded_when_expected": (
            expected_route_convention_schema is None
            or bool(solver_route_convention_schema_id)
        ),
        "expected_solver_route_convention_schema_id_matches": (
            expected_route_convention_schema is None
            or solver_route_convention_schema_id == expected_route_convention_schema
        ),
        "solver_contract_artifact_id_recorded_when_required": (
            not require_solver_contract_artifact or bool(solver_contract_artifact_id)
        ),
        "solver_contract_digest_recorded_when_required": (
            not require_solver_contract_artifact or bool(solver_contract_digest)
        ),
        "solver_contract_path_recorded_when_required": (
            not require_solver_contract_artifact or bool(solver_contract_path)
        ),
        "expected_solver_contract_artifact_id_matches": (
            expected_contract_id is None or solver_contract_artifact_id == expected_contract_id
        ),
        "expected_solver_contract_digest_matches": (
            expected_contract_digest is None or solver_contract_digest == expected_contract_digest
        ),
        "expected_solver_contract_path_matches": (
            expected_contract_path is None or solver_contract_path == expected_contract_path
        ),
        "volume_route_kinds_cover_inventory": set(active_volume_kinds).issubset(set(route_volume_kinds)),
        "required_volume_routes_present": set(required_volumes).issubset(set(route_volume_kinds)),
        "surface_route_kinds_cover_inventory": set(active_surface_kinds).issubset(set(route_surface_kinds)),
        "required_surface_routes_present": set(required_surfaces).issubset(set(route_surface_kinds)),
        "hex_primary_volume_role_recorded": has_role("hex", ("primary", "volume")),
        "pyramid_transition_role_recorded": has_role("pyramid", ("transition",), not_primary=True),
        "tet_compatibility_role_recorded": (
            has_role("tet", ("compatibility",)) or has_role("tet", ("subregion",))
        ),
        "no_implicit_tetization_recorded": no_implicit_tetization,
        "tet_only_owner_is_netgen": tet_only_owner == "netgen_tri_tet_path",
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_mixed_solver_route_manifest_gate",
        "status": "ok" if not issues else "needs_attention",
        "solver_route_package_id": package_id,
        "expected_package_id": expected_id,
        "routing_hint": routing_hint,
        "expected_routing_hint": expected_hint,
        "route_policy": route_policy,
        "downstream_solver": downstream_solver,
        "solver_route_convention_schema_id": solver_route_convention_schema_id,
        "solver_route_convention_schema_ids": solver_route_convention_schema_ids,
        "expected_solver_route_convention_schema_id": expected_route_convention_schema,
        "solver_contract_artifact_id": solver_contract_artifact_id,
        "solver_contract_digest": solver_contract_digest,
        "solver_contract_path": solver_contract_path,
        "expected_solver_contract_artifact_id": expected_contract_id,
        "expected_solver_contract_digest": expected_contract_digest,
        "expected_solver_contract_path": expected_contract_path,
        "solver_contract_artifact_required": bool(require_solver_contract_artifact),
        "solver_route_convention_schema_required": route_convention_schema_required,
        "tet_only_owner": tet_only_owner,
        "active_volume_kinds": active_volume_kinds,
        "route_volume_kinds": route_volume_kinds,
        "required_volume_kinds": required_volumes,
        "active_surface_kinds": active_surface_kinds,
        "route_surface_kinds": route_surface_kinds,
        "required_surface_kinds": required_surfaces,
        "volume_routes": volume_routes,
        "surface_routes": surface_routes,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use this after Cubit inventory proves a hex+pyramid+tet route and before solver-ready promotion.",
            "Pyramid cells should be explicit transition bridge cells, not silently tetized or treated as a primary solver region.",
            "Tet-only educational .vol files remain on the Netgen/OCC tri/tet route; mixed Cubit packages need their own solver route manifest.",
            "The solver-route convention schema id distinguishes the recorded hex/pyramid/tet role mapping from a value-only or layout-only route manifest.",
            "When promoting a mixed route to solver-ready evidence, bind it to the downstream solver/reader contract artifact that accepts the recorded element families.",
        ],
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
    required_surface_kinds: Iterable[str] = ("quad", "triangle"),
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
    required_surfaces = [
        str(kind).strip().lower()
        for kind in required_surface_kinds
        if str(kind).strip()
    ]
    transitions = [str(kind).strip().lower() for kind in transition_kinds if str(kind).strip()]
    transition_names = [str(name) for name in transition_material_names if str(name)]
    routing_hint = str(inventory.get("routing_hint", ""))
    checks = {
        "volume_elements_present": int(inventory.get("volume_elements", 0) or 0) > 0,
        "surface_elements_present": int(inventory.get("surface_elements", 0) or 0) > 0,
        "required_volume_kinds_present": all(volume_counts.get(kind, 0) > 0 for kind in required),
        "required_surface_kinds_present": all(surface_counts.get(kind, 0) > 0 for kind in required_surfaces),
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
        "required_surface_kinds": required_surfaces,
        "transition_kinds": transitions,
        "materials": materials,
        "transition_material_names": transition_names,
        "routing_hint": routing_hint,
        "checks": checks,
        "notes": [
            "Run this before importing hex+pyramid+tet meshes into radia-ngsolve or Radia-style mixed-element lanes.",
            "A pyramid transition can have zero sidecar material volume and still be required topology in the .vol inventory.",
            "Mixed Cubit handoff should expose both quad and triangle surface families before solver boundary labels are trusted.",
        ],
    }


def cubit_live_mixed_mesh_python_gate(
    summary: Mapping[str, object],
    *,
    expected_total_volume: float | None = None,
    volume_relative_tolerance: float = 1.0e-9,
    min_scaled_jacobian: float = 0.0,
) -> dict[str, object]:
    """Gate source-journal mixed-mesh evidence from a headless Python run.

    This complements the exported-mesh inventory gates: it checks that a
    source-native ``.jou``/``.py`` case was actually replayed without a GUI,
    produced the hex-pyramid-tet transition, retained positive element
    quality, and preserved independently known CAD volume.
    """

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    if volume_relative_tolerance < 0.0:
        raise ValueError("volume_relative_tolerance must be non-negative")

    def finite_number(value: object, field: str) -> float:
        try:
            number = float(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{field} must be a finite number") from exc
        if not isfinite(number):
            raise ValueError(f"{field} must be a finite number")
        return number

    counts_raw = summary.get("element_counts", {})
    quality_raw = summary.get("quality", {})
    volumes_raw = summary.get("volumes", {})
    if not isinstance(counts_raw, Mapping):
        raise ValueError("summary['element_counts'] must be a mapping")
    if not isinstance(quality_raw, Mapping):
        raise ValueError("summary['quality'] must be a mapping")
    if not isinstance(volumes_raw, Mapping):
        raise ValueError("summary['volumes'] must be a mapping")

    counts = {str(key).lower(): int(value) for key, value in counts_raw.items()}
    volumes = {
        str(key): finite_number(value, f"volumes[{key!r}]")
        for key, value in volumes_raw.items()
    }
    total_volume = finite_number(summary.get("total_volume"), "total_volume")
    mesh_s = finite_number(summary.get("mesh_s"), "mesh_s")
    source_journal = Path(str(summary.get("source_journal", ""))).name
    execution_mode = str(summary.get("execution_mode", "")).strip().lower()
    version = str(summary.get("version", "")).strip()

    quality_minima: dict[str, float] = {}
    quality_metrics: dict[str, str] = {}
    for kind in ("hex", "tet"):
        row = quality_raw.get(kind, {})
        if not isinstance(row, Mapping):
            raise ValueError(f"summary['quality']['{kind}'] must be a mapping")
        quality_minima[kind] = finite_number(row.get("min"), f"quality.{kind}.min")
        quality_metrics[kind] = str(row.get("metric", "")).strip().lower()

    reference_volume = (
        total_volume
        if expected_total_volume is None
        else finite_number(expected_total_volume, "expected_total_volume")
    )
    denominator = max(abs(reference_volume), 1.0)
    relative_error = abs(total_volume - reference_volume) / denominator
    component_sum_error = abs(sum(volumes.values()) - total_volume) / max(abs(total_volume), 1.0)
    checks = {
        "source_native_journal": Path(source_journal).suffix.lower() in {".jou", ".py"},
        "headless_python_api": execution_mode == "python_api_headless",
        "version_recorded": bool(version),
        "mesh_time_nonnegative": mesh_s >= 0.0,
        "hex_present": counts.get("hex", 0) > 0,
        "pyramid_present": counts.get("pyramid", 0) > 0,
        "tet_present": counts.get("tet", 0) > 0,
        "scaled_jacobian_metrics": all(
            quality_metrics[kind] == "scaled jacobian" for kind in ("hex", "tet")
        ),
        "quality_above_threshold": all(
            quality_minima[kind] > min_scaled_jacobian for kind in ("hex", "tet")
        ),
        "cad_volumes_positive": bool(volumes) and all(value > 0.0 for value in volumes.values()),
        "cad_volume_sum_matches": component_sum_error <= volume_relative_tolerance,
        "expected_volume_matches": relative_error <= volume_relative_tolerance,
    }
    return {
        "policy": "cubit_live_mixed_mesh_python_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "source_journal": source_journal,
        "execution_mode": execution_mode,
        "version": version,
        "element_counts": counts,
        "quality_minima": quality_minima,
        "total_volume": total_volume,
        "expected_total_volume": reference_volume,
        "volume_relative_error": relative_error,
        "component_volume_sum_relative_error": component_sum_error,
        "checks": checks,
        "notes": [
            "Cubit owns the hex-led mixed route; a tet-only case belongs to the Netgen route.",
            "Use the documented Python entity families for inventory queries; unsupported aliases must not be treated as empty mesh sets.",
            "This gate validates execution evidence, while exported .vol/Gmsh gates validate downstream topology and labels.",
        ],
    }


def cubit_mapped_boundary_layer_shell_gate(
    summary: Mapping[str, object],
    *,
    min_scaled_jacobian: float = 0.2,
    min_shape: float = 0.0,
    shell_radius_tolerance: float = 1.0e-9,
    volume_relative_tolerance: float = 1.0e-12,
) -> dict[str, object]:
    """Gate a mapped all-hex boundary layer by its nodal radial shells.

    Curved boundary-layer hex centroids lie inside the arithmetic radial
    midpoint because straight chords approximate the cylinder.  Therefore the
    layer contract is checked from the min/max nodal radii of each shell, not
    from centroid-distance thresholds.
    """
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    tolerances = (
        float(min_scaled_jacobian),
        float(min_shape),
        float(shell_radius_tolerance),
        float(volume_relative_tolerance),
    )
    if not all(math.isfinite(value) for value in tolerances):
        raise ValueError("all tolerances must be finite")
    if min_scaled_jacobian <= 0.0 or min_shape < 0.0:
        raise ValueError("quality thresholds must be nonnegative and scaled Jacobian positive")
    if shell_radius_tolerance < 0.0 or volume_relative_tolerance < 0.0:
        raise ValueError("geometry tolerances must be nonnegative")

    counts_raw = summary.get("element_counts", {})
    quality_raw = summary.get("quality", {})
    boundary_raw = summary.get("boundary_layer", {})
    if not all(isinstance(row, Mapping) for row in (counts_raw, quality_raw, boundary_raw)):
        raise ValueError("element_counts, quality, and boundary_layer must be mappings")
    counts = {str(key).lower(): int(value) for key, value in counts_raw.items()}
    hex_count = counts.get("hex", 0)
    first_height = float(boundary_raw.get("first_height", math.nan))
    growth = float(boundary_raw.get("growth", math.nan))
    layer_count = int(boundary_raw.get("layers", 0))
    outer_radius = float(boundary_raw.get("outer_radius", math.nan))
    radial_levels = sorted(float(value) for value in boundary_raw.get("radial_node_levels", []))
    shell_counts_raw = boundary_raw.get("radial_shell_element_counts", {})
    if not isinstance(shell_counts_raw, Mapping):
        raise ValueError("radial_shell_element_counts must be a mapping")
    shell_counts = [
        int(shell_counts_raw.get(f"wall_layer_{index}", 0))
        for index in range(1, layer_count + 1)
    ]
    core_count = int(shell_counts_raw.get("core", 0))

    cumulative = []
    total = 0.0
    if first_height > 0.0 and growth > 0.0 and layer_count > 0:
        for index in range(layer_count):
            total += first_height * growth**index
            cumulative.append(total)
    expected_levels = sorted([outer_radius] + [outer_radius - value for value in cumulative])
    level_errors = [
        min((abs(observed - expected) for observed in radial_levels), default=math.inf)
        for expected in expected_levels
    ]

    scaled = quality_raw.get("scaled_jacobian", {})
    shape = quality_raw.get("shape", {})
    if not isinstance(scaled, Mapping) or not isinstance(shape, Mapping):
        raise ValueError("scaled_jacobian and shape quality rows must be mappings")
    scaled_count = int(scaled.get("count", 0))
    scaled_min = float(scaled.get("min", math.nan))
    shape_count = int(shape.get("count", 0))
    shape_minimum = float(shape.get("min", math.nan))

    volume_before = float(summary.get("cad_volume_before_scale", math.nan))
    analytic_volume = float(summary.get("analytic_volume_before_scale", math.nan))
    volume_after = float(summary.get("cad_volume_after_scale", math.nan))
    unit_scale = float(summary.get("unit_scale", math.nan))
    coordinate_scale_error = float(summary.get("coordinate_scale_max_abs_error", math.nan))
    volume_error = abs(volume_before - analytic_volume) / max(abs(analytic_volume), 1.0e-300)
    scale_error = abs(volume_after / volume_before - unit_scale**3) if volume_before else math.inf

    process = _headless_process_evidence(summary)
    script_errors = process["script_errors"]
    checks = {
        "source_native_journal_digest_recorded": (
            Path(str(summary.get("source_journal", ""))).suffix.lower() == ".jou"
            and len(str(summary.get("source_sha256") or "")) == 64
        ),
        "headless_python_api": str(summary.get("execution_mode", "")).lower()
        == "python_api_headless",
        "headless_flags_recorded": {"-nographics", "-batch"}.issubset(
            {str(value).lower() for value in summary.get("headless_flags", [])}
        ),
        "persistent_gui_not_started": summary.get("persistent_gui_started") is False,
        "single_line_compile_wrapper_recorded": summary.get("batch_wrapper_mode")
        == "single_line_compile_wrapper",
        "direct_multiline_batch_failure_recorded": summary.get(
            "direct_multiline_batch_rejected"
        )
        is True,
        "fresh_result_artifact": process["result_fresh"] is True,
        "script_error_lines_empty": not script_errors,
        "process_exit_semantics_acceptable": process["process_exit_acceptable"] is True,
        "hex_only_volume_mesh": hex_count > 0 and all(
            counts.get(kind, 0) == 0 for kind in ("pyramid", "wedge", "tet")
        ),
        "boundary_layer_parameters_valid": (
            first_height > 0.0
            and growth >= 1.0
            and layer_count > 0
            and outer_radius > sum(first_height * growth**index for index in range(layer_count))
        ),
        "every_requested_radial_interface_present": (
            len(level_errors) == layer_count + 1
            and all(error <= shell_radius_tolerance for error in level_errors)
        ),
        "every_boundary_layer_shell_occupied": (
            len(shell_counts) == layer_count and all(value > 0 for value in shell_counts)
        ),
        "interior_core_occupied": core_count > 0,
        "all_hexes_classified_by_shell_or_core": sum(shell_counts) + core_count == hex_count,
        "scaled_jacobian_above_threshold": (
            scaled_count == hex_count
            and math.isfinite(scaled_min)
            and scaled_min >= min_scaled_jacobian
        ),
        "shape_above_threshold": (
            shape_count == hex_count
            and math.isfinite(shape_minimum)
            and shape_minimum > min_shape
        ),
        "cad_volume_matches_analytic_geometry": volume_error <= volume_relative_tolerance,
        "uniform_scale_preserves_cubic_volume_law": (
            unit_scale > 0.0 and scale_error <= volume_relative_tolerance
        ),
        "node_coordinates_follow_uniform_scale": (
            math.isfinite(coordinate_scale_error)
            and coordinate_scale_error <= shell_radius_tolerance
        ),
    }
    return {
        "policy": "cubit_mapped_boundary_layer_shell_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "hex_count": hex_count,
            "boundary_layer_shell_counts": shell_counts,
            "core_hex_count": core_count,
            "expected_radial_levels": expected_levels,
            "radial_level_absolute_errors": level_errors,
            "scaled_jacobian_min": scaled_min,
            "shape_min": shape_minimum,
            "cad_volume_relative_error": volume_error,
            "cubic_scale_absolute_error": scale_error,
            "process_exit_code": process["process_exit_code"],
        },
        "launcher_classification": process["launcher_classification"],
        "notes": [
            "Validate curved boundary-layer thickness from nodal radial interfaces; polygon-chord centroids are not radial midpoints.",
            "Cubit is the hex-led boundary-layer route; tet-only educational meshes belong to the Netgen route.",
            "A nonzero launcher exit is acceptable only for a fresh passing artifact, no script errors, and exact allowlisted startup diagnostics.",
            "Cubit batch reads Python line by line, so multiline implementations require a one-line compile/exec wrapper.",
        ],
    }


def cubit_sweep_along_curve_gate(
    summary: Mapping[str, object],
    *,
    min_scaled_jacobian: float = 0.2,
    volume_relative_tolerance: float = 1.0e-12,
) -> dict[str, object]:
    """Gate a source-quad x sweep-interval all-hex replay and launcher status."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    if min_scaled_jacobian <= 0.0 or volume_relative_tolerance < 0.0:
        raise ValueError("quality threshold must be positive and volume tolerance nonnegative")

    counts_raw = summary.get("element_counts", {})
    quality_raw = summary.get("quality", {})
    volumes_raw = summary.get("cad_volume_by_body", {})
    if not isinstance(counts_raw, Mapping) or not isinstance(quality_raw, Mapping):
        raise ValueError("element_counts and quality must be mappings")
    if not isinstance(volumes_raw, Mapping):
        raise ValueError("cad_volume_by_body must be a mapping")

    counts = {str(key).lower(): int(value) for key, value in counts_raw.items()}
    source_quads = int(summary.get("source_quad_count", 0))
    sweep_intervals = int(summary.get("sweep_interval_count", 0))
    total_volume = float(summary.get("total_cad_volume", math.nan))
    analytic_volume = float(summary.get("analytic_volume", math.nan))
    scaled = quality_raw.get("scaled_jacobian", {})
    shape = quality_raw.get("shape", {})
    if not isinstance(scaled, Mapping) or not isinstance(shape, Mapping):
        raise ValueError("quality rows must be mappings")
    scaled_min = float(scaled.get("min", math.nan))
    shape_min = float(shape.get("min", math.nan))

    process = _headless_process_evidence(summary)
    script_errors = process["script_errors"]
    process_exit_code = int(process["process_exit_code"])
    result_fresh = process["result_fresh"] is True
    process_exit_acceptable = process["process_exit_acceptable"] is True
    volume_error = abs(total_volume - analytic_volume) / max(abs(analytic_volume), 1.0)
    volume_sum_error = abs(sum(float(value) for value in volumes_raw.values()) - total_volume) / max(abs(total_volume), 1.0)
    checks = {
        "source_native_journal": Path(str(summary.get("source_journal", ""))).suffix.lower() == ".jou",
        "headless_python_api": str(summary.get("execution_mode", "")).lower() == "python_api_headless",
        "persistent_gui_not_started": summary.get("persistent_gui_started") is False,
        "version_recorded": bool(str(summary.get("version", "")).strip()),
        "fresh_result_artifact": result_fresh,
        "script_error_lines_empty": not script_errors,
        "process_exit_semantics_acceptable": process_exit_acceptable,
        "hex_only_volume_mesh": counts.get("hex", 0) > 0 and all(
            counts.get(kind, 0) == 0 for kind in ("pyramid", "wedge", "tet")
        ),
        "source_quads_and_sweep_intervals_positive": source_quads > 0 and sweep_intervals > 0,
        "sweep_layer_count_conserved": counts.get("hex", 0) == source_quads * sweep_intervals,
        "scaled_jacobian_above_threshold": scaled_min >= float(min_scaled_jacobian),
        "shape_positive": shape_min > 0.0,
        "cad_body_volumes_positive": bool(volumes_raw) and all(float(value) > 0.0 for value in volumes_raw.values()),
        "cad_volume_sum_matches": volume_sum_error <= volume_relative_tolerance,
        "analytic_volume_matches": volume_error <= volume_relative_tolerance,
    }
    return {
        "policy": "cubit_sweep_along_curve_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "hex_count": counts.get("hex", 0),
            "expected_hex_count": source_quads * sweep_intervals,
            "scaled_jacobian_min": scaled_min,
            "shape_min": shape_min,
            "volume_relative_error": volume_error,
            "component_volume_sum_relative_error": volume_sum_error,
            "process_exit_code": process_exit_code,
        },
        "launcher_classification": process["launcher_classification"],
        "notes": [
            "For a mesh-carrying sweep, hex count must equal mapped source quads times path intervals.",
            "Do not treat an unsupported Python entity alias as an empty element family.",
            "A nonzero launcher exit is accepted only for allowlisted startup option/path diagnostics, with no script errors and a fresh passing artifact.",
        ],
    }


def cubit_partitioned_sweep_compatibility_gate(
    summary: Mapping[str, object],
    *,
    min_scaled_jacobian: float = 0.2,
    volume_relative_tolerance: float = 1.0e-12,
) -> dict[str, object]:
    """Gate a legacy partition/webcut journal promoted to an all-hex sweep."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    if min_scaled_jacobian <= 0.0 or volume_relative_tolerance < 0.0:
        raise ValueError("quality threshold must be positive and volume tolerance nonnegative")
    counts_raw = summary.get("element_counts", {})
    quality_raw = summary.get("quality", {})
    volumes_raw = summary.get("cad_volume_by_body", {})
    transforms_raw = summary.get("compatibility_transforms", {})
    if not all(isinstance(row, Mapping) for row in (counts_raw, quality_raw, volumes_raw, transforms_raw)):
        raise ValueError("counts, quality, volumes, and compatibility transforms must be mappings")
    counts = {str(key).lower(): int(value) for key, value in counts_raw.items()}
    transforms = {str(key).strip().lower(): str(value).strip().lower() for key, value in transforms_raw.items()}
    volumes = {int(key): float(value) for key, value in volumes_raw.items()}
    ids = sorted(volumes)
    base_volume_count = int(summary.get("base_volume_count", 0))
    final_volume_count = int(summary.get("volume_count", 0))
    base_hex_count = int(summary.get("base_hex_count", 0))
    copy_factor = int(summary.get("mesh_copy_factor", 0))
    command_count = int(summary.get("source_command_count", 0))
    executed_count = int(summary.get("executed_command_count", 0))
    webcut_count = int(summary.get("webcut_count", 0))
    scaled = quality_raw.get("scaled_jacobian", {})
    shape = quality_raw.get("shape", {})
    if not isinstance(scaled, Mapping) or not isinstance(shape, Mapping):
        raise ValueError("quality rows must be mappings")
    scaled_min = float(scaled.get("min", math.nan))
    shape_min = float(shape.get("min", math.nan))
    total_volume = float(summary.get("total_cad_volume", math.nan))
    volume_sum_error = abs(sum(volumes.values()) - total_volume) / max(abs(total_volume), 1.0e-300)
    base_ids = ids[:base_volume_count]
    copy_ids = ids[base_volume_count:]
    base_volume = sum(volumes[key] for key in base_ids)
    copy_volume = sum(volumes[key] for key in copy_ids)
    copy_volume_error = abs(copy_volume - base_volume * max(copy_factor - 1, 0)) / max(abs(base_volume), 1.0e-300)
    process = _headless_process_evidence(summary)
    script_errors = process["script_errors"]
    process_exit_code = int(process["process_exit_code"])
    result_fresh = process["result_fresh"] is True
    process_exit_acceptable = process["process_exit_acceptable"] is True
    checks = {
        "source_native_journal_digest_recorded": Path(str(summary.get("source_journal", ""))).suffix.lower() == ".jou"
        and bool(str(summary.get("source_sha256") or "").strip()),
        "headless_python_api": str(summary.get("execution_mode", "")).lower() == "python_api_headless",
        "persistent_gui_not_started": summary.get("persistent_gui_started") is False,
        "fresh_result_artifact": result_fresh,
        "script_error_lines_empty": not script_errors,
        "process_exit_semantics_acceptable": process_exit_acceptable,
        "all_source_commands_replayed": command_count > 0 and command_count == executed_count,
        "legacy_quad_dominant_promoted_to_pave": transforms.get(
            "surface 22 scheme quad_dominant"
        ) == "surface 22 scheme pave",
        "partition_webcuts_recorded": webcut_count == 11,
        "mesh_copy_factor_recorded": copy_factor == 2,
        "volume_copy_count_conserved": base_volume_count > 0
        and final_volume_count == base_volume_count * copy_factor
        and len(volumes) == final_volume_count,
        "hex_copy_count_conserved": base_hex_count > 0
        and counts.get("hex", 0) == base_hex_count * copy_factor,
        "hex_only_volume_mesh": counts.get("hex", 0) > 0 and all(
            counts.get(kind, 0) == 0 for kind in ("pyramid", "wedge", "tet")
        ),
        "scaled_jacobian_above_threshold": scaled_min >= float(min_scaled_jacobian),
        "shape_positive": shape_min > 0.0,
        "cad_body_volumes_positive": bool(volumes) and all(value > 0.0 for value in volumes.values()),
        "cad_volume_sum_matches": volume_sum_error <= volume_relative_tolerance,
        "copied_cad_volume_matches_base": copy_volume_error <= volume_relative_tolerance,
    }
    return {
        "policy": "cubit_partitioned_sweep_compatibility_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "hex_count": counts.get("hex", 0),
            "base_hex_count": base_hex_count,
            "volume_count": final_volume_count,
            "scaled_jacobian_min": scaled_min,
            "shape_min": shape_min,
            "cad_volume_sum_relative_error": volume_sum_error,
            "copy_volume_relative_error": copy_volume_error,
            "process_exit_code": process_exit_code,
        },
        "launcher_classification": process["launcher_classification"],
        "notes": [
            "Promote an obsolete meshing scheme through an explicit compatibility map; never silently skip it.",
            "A mesh-carrying copy must conserve both hex count and CAD volume.",
            "A nonzero launcher exit is acceptable only for allowlisted startup diagnostics with a fresh artifact and no script errors.",
        ],
    }


def cubit_pyramid_degenerate_hex_export_gate(
    summary: Mapping[str, object],
    *,
    min_hex_scaled_jacobian: float = 0.2,
    min_pyramid_geometric_volume: float = 0.0,
) -> dict[str, object]:
    """Gate CPYRAM and degenerate-CHEXA exports from one Cubit mesh.

    The source database can contain more element families than the registered
    export blocks.  The ``nopyramid`` route also has an order-specific detail:
    an order-2 PYRAMID13 is emitted as a linear degenerate CHEXA8 while native
    hexes remain CHEXA20.  This gate makes both facts explicit before a deck is
    called solver-ready.
    """
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    if min_hex_scaled_jacobian <= 0.0 or min_pyramid_geometric_volume < 0.0:
        raise ValueError("quality thresholds must be nonnegative and hex threshold positive")

    def mapping(name: str) -> Mapping[str, object]:
        row = summary.get(name, {})
        if not isinstance(row, Mapping):
            raise ValueError(f"summary[{name!r}] must be a mapping")
        return row

    def integer_map(row: Mapping[str, object]) -> dict[str, int]:
        return {str(key).strip().lower(): int(value) for key, value in row.items()}

    counts = integer_map(mapping("element_counts"))
    quality = mapping("quality")
    blocks = mapping("block_inventory")
    true_order1 = mapping("pyramid_card_deck")
    false_order1 = mapping("nopyramid_deck")
    true_order2 = mapping("pyramid_card_deck_order2")
    false_order2 = mapping("nopyramid_deck_order2")

    def deck_counts(deck: Mapping[str, object]) -> dict[str, int]:
        row = deck.get("card_counts", {})
        if not isinstance(row, Mapping):
            raise ValueError("deck card_counts must be a mapping")
        return {str(key).upper(): int(value) for key, value in row.items()}

    def int_list(deck: Mapping[str, object], name: str) -> list[int]:
        row = deck.get(name, [])
        if isinstance(row, (str, bytes)) or not isinstance(row, Iterable):
            raise ValueError(f"deck {name} must be an iterable")
        return [int(value) for value in row]

    block_rows = []
    for key, value in blocks.items():
        if not isinstance(value, Mapping):
            raise ValueError(f"block_inventory[{key!r}] must be a mapping")
        block_rows.append(integer_map({
            kind: value.get(kind, 0) for kind in ("hex", "pyramid", "tet", "wedge")
        }))
    block_totals = {
        kind: sum(row.get(kind, 0) for row in block_rows)
        for kind in ("hex", "pyramid", "tet", "wedge")
    }

    hex_quality = quality.get("hex", {})
    pyramid_quality = quality.get("pyramid", {})
    if not isinstance(hex_quality, Mapping) or not isinstance(pyramid_quality, Mapping):
        raise ValueError("hex and pyramid quality rows must be mappings")
    hex_minimum = float(hex_quality.get("minimum", math.nan))
    pyramid_volume_minimum = float(
        pyramid_quality.get("geometric_volume_minimum", math.nan)
    )

    process = _headless_process_evidence(summary)
    script_errors = process["script_errors"]
    process_exit_code = int(process["process_exit_code"])
    result_fresh = process["result_fresh"] is True
    process_exit_acceptable = process["process_exit_acceptable"] is True

    hex_count = counts.get("hex", 0)
    pyramid_count = counts.get("pyramid", 0)
    exported_count = hex_count + pyramid_count
    true1_counts = deck_counts(true_order1)
    false1_counts = deck_counts(false_order1)
    true2_counts = deck_counts(true_order2)
    false2_counts = deck_counts(false_order2)
    expected_order2_nodes = sorted([20] * hex_count + [8] * pyramid_count)
    expected_order2_unique = sorted([20] * hex_count + [5] * pyramid_count)
    digests = [
        str(deck.get("sha256") or "").strip()
        for deck in (true_order1, false_order1, true_order2, false_order2)
    ]

    checks = {
        "source_native_journal_digest_recorded": (
            Path(str(summary.get("source_journal", ""))).suffix.lower() == ".jou"
            and len(str(summary.get("source_sha256") or "")) == 64
        ),
        "headless_python_api": str(summary.get("execution_mode", "")).lower()
        == "python_api_headless",
        "headless_flags_recorded": {"-nographics", "-batch"}.issubset(
            {str(value).lower() for value in summary.get("headless_flags", [])}
        ),
        "persistent_gui_not_started": summary.get("persistent_gui_started") is False,
        "single_line_compile_wrapper_recorded": summary.get("batch_wrapper_mode")
        == "single_line_compile_wrapper",
        "direct_multiline_batch_failure_recorded": summary.get(
            "direct_multiline_batch_rejected"
        )
        is True,
        "fresh_result_artifact": result_fresh,
        "script_error_lines_empty": not script_errors,
        "process_exit_semantics_acceptable": process_exit_acceptable,
        "full_database_inventory_is_mixed": (
            hex_count > 0
            and pyramid_count > 0
            and counts.get("tet", 0) > 0
            and counts.get("wedge", 0) > 0
        ),
        "registered_blocks_export_only_hex_and_pyramid": (
            block_totals.get("hex") == hex_count
            and block_totals.get("pyramid") == pyramid_count
            and block_totals.get("tet") == 0
            and block_totals.get("wedge") == 0
        ),
        "export_scope_is_registered_blocks_only": summary.get("export_scope_claim")
        == "registered_blocks_only",
        "hex_scaled_jacobian_above_threshold": (
            int(hex_quality.get("count", 0)) == hex_count
            and isfinite(hex_minimum)
            and hex_minimum >= float(min_hex_scaled_jacobian)
        ),
        "pyramid_quality_fallback_is_positive_geometry": (
            int(pyramid_quality.get("api_value_count", -1)) == 0
            and int(pyramid_quality.get("geometric_volume_count", 0)) == pyramid_count
            and isfinite(pyramid_volume_minimum)
            and pyramid_volume_minimum > float(min_pyramid_geometric_volume)
        ),
        "order1_cpyram_deck_matches_registered_blocks": (
            true1_counts.get("CHEXA") == hex_count
            and true1_counts.get("CPYRAM") == pyramid_count
            and sorted(int_list(true_order1, "chexa_node_counts")) == [8] * hex_count
            and sorted(int_list(true_order1, "cpyram_node_counts")) == [5] * pyramid_count
        ),
        "order1_nopyramid_is_degenerate_chexa8": (
            false1_counts.get("CHEXA") == exported_count
            and false1_counts.get("CPYRAM") == 0
            and int(false_order1.get("regular_chexa_count", -1)) == hex_count
            and int(false_order1.get("degenerate_chexa_count", -1)) == pyramid_count
            and sorted(int_list(false_order1, "chexa_unique_node_counts"))
            == sorted([8] * hex_count + [5] * pyramid_count)
        ),
        "order2_cpyram_deck_is_hex20_and_pyramid13": (
            true2_counts.get("CHEXA") == hex_count
            and true2_counts.get("CPYRAM") == pyramid_count
            and sorted(int_list(true_order2, "chexa_node_counts")) == [20] * hex_count
            and sorted(int_list(true_order2, "cpyram_node_counts")) == [13] * pyramid_count
        ),
        "order2_nopyramid_linearization_recorded": (
            false2_counts.get("CHEXA") == exported_count
            and false2_counts.get("CPYRAM") == 0
            and sorted(int_list(false_order2, "chexa_node_counts"))
            == expected_order2_nodes
            and sorted(int_list(false_order2, "chexa_unique_node_counts"))
            == expected_order2_unique
            and summary.get("pyramid_conversion_order_policy")
            == "linearized_degenerate_chexa8"
            and summary.get("order2_nopyramid_uniform_order_claimed") is False
        ),
        "four_export_digests_recorded_and_distinct": (
            all(len(value) == 64 for value in digests) and len(set(digests)) == 4
        ),
        "cad_volume_positive": float(summary.get("total_cad_volume", math.nan)) > 0.0,
    }
    return {
        "policy": "cubit_pyramid_degenerate_hex_export_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "database_element_counts": counts,
            "registered_block_totals": block_totals,
            "exported_volume_element_count": exported_count,
            "order2_nopyramid_chexa_node_counts": int_list(
                false_order2, "chexa_node_counts"
            ),
            "order2_nopyramid_chexa_unique_node_counts": int_list(
                false_order2, "chexa_unique_node_counts"
            ),
            "hex_scaled_jacobian_minimum": hex_minimum,
            "pyramid_geometric_volume_minimum": pyramid_volume_minimum,
            "process_exit_code": process_exit_code,
        },
        "launcher_classification": process["launcher_classification"],
        "notes": [
            "Database inventory and registered export-block inventory are different contracts; do not claim unregistered tet or wedge cells are present in the deck.",
            "nopyramid preserves first-order pyramid count by repeated-node CHEXA8 cards rather than deleting the transition cells.",
            "At order 2, native hexes remain CHEXA20 but converted pyramids are linearized to degenerate CHEXA8; reject a uniform-order claim.",
            "Cubit batch reads Python line by line, so multiline implementations require a one-line compile/exec wrapper.",
        ],
    }


def cubit_mixed_order_series_inventory_gate(
    rows: Iterable[Mapping[str, object]],
    *,
    required_volume_kinds: Iterable[str] = ("hex", "pyramid", "tet"),
    required_surface_kinds: Iterable[str] = ("quad", "triangle"),
    expected_routing_hint: str = "cubit_hex_or_mixed_path",
    require_first_order_inventory: bool = True,
) -> dict[str, object]:
    """Check that Cubit mixed-mesh routing is invariant across export orders.

    High-order Cubit ``.vol`` files may add or enlarge ``curvedelements``
    records, but the first-order arity inventory still decides whether the mesh
    belongs to the Cubit hex/mixed lane or the Netgen tri/tet education lane.
    This gate replays an order series and rejects topology/routing drift.
    """

    def normalize_counts(raw: object, field: str) -> dict[str, int]:
        if not isinstance(raw, dict):
            raise ValueError(f"{field} must be a mapping")
        return {str(key): int(value) for key, value in raw.items()}

    series: list[dict[str, object]] = []
    for index, row in enumerate(rows):
        if not isinstance(row, Mapping):
            raise ValueError("each order-series row must be a mapping")
        inventory = row.get("inventory", row)
        if not isinstance(inventory, Mapping):
            raise ValueError("row inventory must be a mapping")
        order_value = row.get("order", inventory.get("order", index + 1))
        order = int(order_value) if order_value is not None else index + 1
        volume_counts = normalize_counts(inventory.get("volume_kind_counts", {}), "volume_kind_counts")
        surface_counts = normalize_counts(inventory.get("surface_kind_counts", {}), "surface_kind_counts")
        series.append({
            "order": order,
            "source": str(inventory.get("source", row.get("source", ""))),
            "volume_kind_counts": volume_counts,
            "surface_kind_counts": surface_counts,
            "routing_hint": str(inventory.get("routing_hint", "")),
            "curvedelements_present": bool(inventory.get("curvedelements_present", False)),
            "is_tri_tet_only": bool(inventory.get("is_tri_tet_only", False)),
            "has_mixed_hex_transition": bool(inventory.get("has_mixed_hex_transition", False)),
        })

    required_volume = [str(kind).strip().lower() for kind in required_volume_kinds if str(kind).strip()]
    required_surface = [str(kind).strip().lower() for kind in required_surface_kinds if str(kind).strip()]
    orders = [int(row["order"]) for row in series]
    baseline = series[0] if series else None
    baseline_volume = baseline["volume_kind_counts"] if baseline else {}
    baseline_surface = baseline["surface_kind_counts"] if baseline else {}
    first_order_rows = [row for row in series if row["order"] == 1]
    checks = {
        "series_rows_present": bool(series),
        "orders_recorded": bool(orders),
        "orders_unique": len(orders) == len(set(orders)),
        "first_order_inventory_present": (
            True if not require_first_order_inventory else bool(first_order_rows)
        ),
        "first_order_inventory_not_curved": (
            True if not require_first_order_inventory or not first_order_rows
            else not any(row["curvedelements_present"] for row in first_order_rows)
        ),
        "volume_kind_counts_invariant": all(row["volume_kind_counts"] == baseline_volume for row in series),
        "surface_kind_counts_invariant": all(row["surface_kind_counts"] == baseline_surface for row in series),
        "routing_hint_invariant": all(row["routing_hint"] == expected_routing_hint for row in series),
        "routing_hint_is_cubit_mixed": all(row["routing_hint"] == expected_routing_hint for row in series),
        "required_volume_kinds_present": (
            bool(baseline_volume)
            and all(int(baseline_volume.get(kind, 0)) > 0 for kind in required_volume)
        ),
        "required_surface_kinds_present": (
            bool(baseline_surface)
            and all(int(baseline_surface.get(kind, 0)) > 0 for kind in required_surface)
        ),
        "has_mixed_hex_transition": all(row["has_mixed_hex_transition"] for row in series),
        "not_tri_tet_only": all(not row["is_tri_tet_only"] for row in series),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_mixed_order_series_inventory_gate",
        "status": "ok" if not issues else "needs_attention",
        "orders": orders,
        "expected_routing_hint": expected_routing_hint,
        "baseline_volume_kind_counts": baseline_volume,
        "baseline_surface_kind_counts": baseline_surface,
        "required_volume_kinds": required_volume,
        "required_surface_kinds": required_surface,
        "series": series,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Route high-order Cubit .vol files from the first-order volume/surface arity inventory.",
            "The curvedelements section may grow with order, but hex/pyramid/tet and quad/triangle routing must not drift.",
            "Mixed Cubit order series are not Netgen tri/tet education inputs even when the file extension is .vol.",
        ],
    }


def cubit_mixed_interface_adjacency_gate(
    rows: Iterable[Mapping[str, object]],
    *,
    required_roles: Iterable[str] = ("hex_to_transition", "transition_to_tet"),
    transition_material_names: Iterable[str] = ("pyramid_transition", "pyram", "transition"),
    expected_role_surface_kinds: Mapping[str, str] | None = None,
    expected_role_volume_kind_pairs: Mapping[str, Iterable[str]] | None = None,
) -> dict[str, object]:
    """Check that a mixed Cubit interface ledger names both sides of transitions.

    Hex/pyramid/tet inventories prove that all topology families exist, but they
    do not prove which surface is the hex-to-transition interface and which one
    is the transition-to-tet interface.  This gate keeps a small replayable
    ledger of interface surface id/name, surface kind, adjacent materials, and
    adjacent volume kinds before BND rows or submodel boundaries are trusted.
    """

    def text(value) -> str:
        return str(value or "").strip()

    def normalize_kind(value) -> str:
        key = text(value).lower()
        aliases = {
            "tri": "triangle",
            "tri3": "triangle",
            "face3": "triangle",
            "quad4": "quad",
            "quadrilateral": "quad",
            "tetra": "tet",
            "tetrahedron": "tet",
            "tet4": "tet",
            "hexahedron": "hex",
            "hex8": "hex",
            "pyram": "pyramid",
            "pyramid5": "pyramid",
        }
        return aliases.get(key, key)

    def normalize_list(value) -> list[str]:
        if value is None:
            return []
        if isinstance(value, str):
            return [item.strip() for item in value.replace(";", ",").split(",") if item.strip()]
        try:
            iterator = iter(value)
        except TypeError:
            iterator = iter((value,))
        out = []
        for item in iterator:
            if isinstance(item, Mapping):
                raw = item.get("name", item.get("label", item.get("id", "")))
            else:
                raw = item
            item_text = text(raw)
            if item_text:
                out.append(item_text)
        return out

    role_surface_defaults = {
        "hex_to_transition": "quad",
        "transition_to_tet": "triangle",
    }
    role_volume_defaults = {
        "hex_to_transition": ("hex", "pyramid"),
        "transition_to_tet": ("pyramid", "tet"),
    }
    role_surface = {
        text(role).lower(): normalize_kind(kind)
        for role, kind in (expected_role_surface_kinds or role_surface_defaults).items()
        if text(role)
    }
    role_volume = {
        text(role).lower(): sorted({normalize_kind(kind) for kind in normalize_list(kinds) if text(kind)})
        for role, kinds in (expected_role_volume_kind_pairs or role_volume_defaults).items()
        if text(role)
    }
    required = [text(role).lower() for role in required_roles if text(role)]
    transition_names = {text(name) for name in transition_material_names if text(name)}

    normalized_rows = []
    missing_identity_rows = []
    missing_adjacency_rows = []
    transition_missing_rows = []
    non_transition_neighbor_missing_rows = []
    role_surface_mismatch_rows = []
    role_volume_mismatch_rows = []
    for index, row in enumerate(rows, start=1):
        record = dict(row)
        role = text(
            record.get("role")
            or record.get("interface_role")
            or record.get("surface_role")
        ).lower()
        surface_kind = normalize_kind(record.get("surface_kind") or record.get("kind"))
        surface_id = text(record.get("surface_id") or record.get("id"))
        surface_name = text(record.get("surface_name") or record.get("name") or record.get("boundary_name"))
        adjacent_materials = normalize_list(
            record.get("adjacent_material_names")
            or record.get("adjacent_materials")
            or record.get("materials")
        )
        adjacent_kinds = [
            normalize_kind(kind)
            for kind in normalize_list(
                record.get("adjacent_volume_kinds")
                or record.get("adjacent_volume_types")
                or record.get("volume_kinds")
            )
        ]
        transition_material_present = (
            True if not transition_names else any(name in transition_names for name in adjacent_materials)
        )
        non_transition_neighbor_present = bool(
            [name for name in adjacent_materials if name not in transition_names]
            or [kind for kind in adjacent_kinds if kind != "pyramid"]
        )
        expected_surface_kind = role_surface.get(role)
        expected_volume_kinds = role_volume.get(role, [])
        surface_kind_matches = (
            True if not expected_surface_kind else surface_kind == expected_surface_kind
        )
        volume_kinds_match = (
            True
            if not expected_volume_kinds
            else set(expected_volume_kinds).issubset(set(adjacent_kinds))
        )
        if not (surface_id or surface_name):
            missing_identity_rows.append(index)
        if len(adjacent_materials) < 2 and len(adjacent_kinds) < 2:
            missing_adjacency_rows.append(index)
        if not transition_material_present:
            transition_missing_rows.append(index)
        if not non_transition_neighbor_present:
            non_transition_neighbor_missing_rows.append(index)
        if not surface_kind_matches:
            role_surface_mismatch_rows.append(index)
        if not volume_kinds_match:
            role_volume_mismatch_rows.append(index)
        normalized_rows.append(
            {
                "row": index,
                "role": role,
                "surface_id": surface_id,
                "surface_name": surface_name,
                "surface_kind": surface_kind,
                "adjacent_material_names": adjacent_materials,
                "adjacent_volume_kinds": adjacent_kinds,
                "expected_surface_kind": expected_surface_kind,
                "expected_volume_kinds": expected_volume_kinds,
            }
        )

    roles_present = sorted({row["role"] for row in normalized_rows if row["role"]})
    surface_ids = [row["surface_id"] for row in normalized_rows if row["surface_id"]]
    checks = {
        "rows_present": bool(normalized_rows),
        "interface_identity_recorded": not missing_identity_rows,
        "interface_roles_recorded": all(bool(row["role"]) for row in normalized_rows),
        "required_roles_present": all(role in roles_present for role in required),
        "surface_ids_unique_when_present": len(surface_ids) == len(set(surface_ids)),
        "adjacency_recorded": not missing_adjacency_rows,
        "transition_material_touches_every_interface": not transition_missing_rows,
        "non_transition_neighbor_recorded": not non_transition_neighbor_missing_rows,
        "role_surface_kinds_match": not role_surface_mismatch_rows,
        "role_volume_kind_pairs_match": not role_volume_mismatch_rows,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_mixed_interface_adjacency_gate",
        "status": "ok" if not issues else "needs_attention",
        "required_roles": required,
        "roles_present": roles_present,
        "transition_material_names": sorted(transition_names),
        "expected_role_surface_kinds": role_surface,
        "expected_role_volume_kind_pairs": role_volume,
        "rows": normalized_rows,
        "missing_identity_rows": missing_identity_rows,
        "missing_adjacency_rows": missing_adjacency_rows,
        "transition_missing_rows": transition_missing_rows,
        "non_transition_neighbor_missing_rows": non_transition_neighbor_missing_rows,
        "role_surface_mismatch_rows": role_surface_mismatch_rows,
        "role_volume_mismatch_rows": role_volume_mismatch_rows,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Run this after the mixed-transition inventory gate and before BND or submodel boundary rows are promoted.",
            "The interface ledger should say which surface connects hex to pyramid and which connects pyramid to tet.",
            "A mixed mesh can have the right element counts while still losing the adjacency role needed for solver-ready boundary setup.",
        ],
    }


def cubit_meshing_scheme_trace_gate(
    trace: Mapping[str, object],
    *,
    expected_trace_id: str | None = None,
    expected_command_digest: str | None = None,
    expected_volume_schemes: Mapping[str, str] | None = None,
    required_command_fragments: Iterable[str] = ("imprint all", "merge all", "export netgen"),
    expected_export_order: int | None = None,
    expected_export_output_artifact_id: str | None = None,
    expected_export_output_digest: str | None = None,
    expected_export_output_path: str | None = None,
    require_export_output_artifact: bool = False,
) -> dict[str, object]:
    """Check Cubit meshing scheme and export-command provenance.

    A mixed ``.vol`` inventory can prove that hex, pyramid, and tet records
    exist, but reproducibility also depends on the Cubit commands that selected
    the volume schemes and exported the file.  Keep this trace next to the mesh
    package so stale ``scheme`` or ``export netgen`` commands fail before solver
    rows are trusted.  The concrete exported ``.vol`` artifact should also be
    named when a trace is promoted to solver-ready evidence; otherwise a fresh
    command trace can accidentally point at an old mesh file.
    """

    if not isinstance(trace, Mapping):
        raise ValueError("trace must be a mapping")

    def text(value) -> str:
        return str(value or "").strip()

    def normalize_commands(value) -> list[str]:
        if value is None:
            return []
        if isinstance(value, str):
            return [line.strip() for line in value.splitlines() if line.strip()]
        try:
            iterator = iter(value)
        except TypeError:
            iterator = iter((value,))
        commands = []
        for item in iterator:
            item_text = text(item)
            if item_text:
                commands.append(item_text)
        return commands

    def normalize_schemes(value) -> dict[str, str]:
        if value is None:
            return {}
        if isinstance(value, Mapping):
            return {text(k): text(v).lower() for k, v in value.items() if text(k) and text(v)}
        try:
            iterator = iter(value)
        except TypeError:
            iterator = iter((value,))
        schemes = {}
        for item in iterator:
            if not isinstance(item, Mapping):
                continue
            volume_id = text(
                item.get("volume_id")
                or item.get("volume")
                or item.get("id")
                or item.get("block_id")
            )
            scheme = text(item.get("scheme") or item.get("mesh_scheme")).lower()
            if volume_id and scheme:
                schemes[volume_id] = scheme
        return schemes

    trace_id = text(
        trace.get("trace_id")
        or trace.get("meshing_trace_id")
        or trace.get("command_trace_id")
    )
    digest = text(
        trace.get("command_digest")
        or trace.get("commands_digest")
        or trace.get("meshing_command_digest")
    )
    commands = normalize_commands(
        trace.get("commands")
        or trace.get("command_trace")
        or trace.get("journal_commands")
    )
    schemes = normalize_schemes(
        trace.get("volume_schemes")
        or trace.get("schemes")
        or trace.get("mesh_schemes")
    )
    required_fragments = [text(fragment).lower() for fragment in required_command_fragments if text(fragment)]
    command_text = "\n".join(command.lower() for command in commands)
    expected_schemes = {
        text(volume): text(scheme).lower()
        for volume, scheme in (expected_volume_schemes or {}).items()
        if text(volume) and text(scheme)
    }
    export_order = trace.get("export_order", trace.get("order"))
    export_order_int = None if export_order is None or text(export_order) == "" else int(export_order)
    export_output_artifact_id = text(
        trace.get("export_output_artifact_id")
        or trace.get("output_artifact_id")
        or trace.get("vol_artifact_id")
        or trace.get("mesh_artifact_id")
    )
    export_output_digest = text(
        trace.get("export_output_digest")
        or trace.get("output_digest")
        or trace.get("vol_digest")
        or trace.get("mesh_digest")
        or trace.get("export_output_sha256")
        or trace.get("output_sha256")
    )
    export_output_path = text(
        trace.get("export_output_path")
        or trace.get("output_path")
        or trace.get("vol_path")
        or trace.get("mesh_path")
        or trace.get("export_path")
    )
    expected_order = None if expected_export_order is None else int(expected_export_order)
    expected_trace = None if expected_trace_id is None else text(expected_trace_id)
    expected_digest = None if expected_command_digest is None else text(expected_command_digest)
    expected_output_artifact = (
        None if expected_export_output_artifact_id is None else text(expected_export_output_artifact_id)
    )
    expected_output_digest = (
        None if expected_export_output_digest is None else text(expected_export_output_digest)
    )
    expected_output_path = (
        None if expected_export_output_path is None else text(expected_export_output_path)
    )
    output_required = bool(
        require_export_output_artifact
        or expected_output_artifact is not None
        or expected_output_digest is not None
        or expected_output_path is not None
    )

    missing_expected_scheme_volumes = [
        volume for volume, expected in expected_schemes.items() if schemes.get(volume) != expected
    ]
    checks = {
        "trace_id_recorded": bool(trace_id),
        "expected_trace_id_matches": expected_trace is None or trace_id == expected_trace,
        "command_digest_recorded": bool(digest),
        "expected_command_digest_matches": expected_digest is None or digest == expected_digest,
        "commands_recorded": bool(commands),
        "required_command_fragments_present": all(fragment in command_text for fragment in required_fragments),
        "volume_schemes_recorded": bool(schemes),
        "expected_volume_schemes_match": not missing_expected_scheme_volumes,
        "export_order_recorded": expected_order is None or export_order_int is not None,
        "expected_export_order_matches": expected_order is None or export_order_int == expected_order,
        "export_output_artifact_id_recorded_when_required": (
            not output_required or bool(export_output_artifact_id)
        ),
        "export_output_digest_recorded_when_required": (
            not output_required or bool(export_output_digest)
        ),
        "export_output_path_recorded_when_required": (
            not output_required or bool(export_output_path)
        ),
        "export_output_artifact_id_recorded_when_expected": (
            expected_output_artifact is None or bool(export_output_artifact_id)
        ),
        "expected_export_output_artifact_id_matches": (
            expected_output_artifact is None or export_output_artifact_id == expected_output_artifact
        ),
        "export_output_digest_recorded_when_expected": (
            expected_output_digest is None or bool(export_output_digest)
        ),
        "expected_export_output_digest_matches": (
            expected_output_digest is None or export_output_digest == expected_output_digest
        ),
        "export_output_path_recorded_when_expected": (
            expected_output_path is None or bool(export_output_path)
        ),
        "expected_export_output_path_matches": (
            expected_output_path is None or export_output_path == expected_output_path
        ),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_meshing_scheme_trace_gate",
        "status": "ok" if not issues else "needs_attention",
        "trace_id": trace_id or None,
        "expected_trace_id": expected_trace,
        "command_digest": digest or None,
        "expected_command_digest": expected_digest,
        "commands": commands,
        "required_command_fragments": required_fragments,
        "volume_schemes": schemes,
        "expected_volume_schemes": expected_schemes,
        "missing_expected_scheme_volumes": missing_expected_scheme_volumes,
        "export_order": export_order_int,
        "expected_export_order": expected_order,
        "export_output_artifact_id": export_output_artifact_id or None,
        "export_output_digest": export_output_digest or None,
        "export_output_path": export_output_path or None,
        "expected_export_output_artifact_id": expected_output_artifact,
        "expected_export_output_digest": expected_output_digest,
        "expected_export_output_path": expected_output_path,
        "require_export_output_artifact": output_required,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Bind Cubit scheme commands and export netgen order to the same package as the .vol inventory.",
            "Record the concrete export output artifact id, digest, and path with the scheme trace before reusing a mesh.",
            "Hex-led and mixed meshes should preserve the scheme trace instead of relying on element counts alone.",
            "Tet-only validation remains the Netgen/OCC education route unless the slot explicitly overrides it.",
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
    expected_length_unit: str | None = None,
    expected_volume_unit: str | None = None,
    expected_area_unit: str | None = None,
    expected_material_names: Iterable[str] = (),
    allow_zero_measurement_names: Iterable[str] = (),
    rel_tol: float = 1e-9,
    abs_tol: float = 1e-12,
) -> dict[str, object]:
    """Validate Cubit CAD mass-property sidecar rows before mesh routing.

    Cubit/Coreform can own the hex-led mesh route, but the CAD handoff should
    keep volume, summed surface area, and bounding-box dimensions as a small
    sidecar.  This replay gate makes those checks executable without reopening
    Cubit, and gives build123d/CST/CAD lanes a common volume/area/bbox contract.
    For mixed hex/pyramid/tet routes, a transition material row may legitimately
    report zero CAD volume in a material sidecar, but its label must still match
    the `.vol` inventory.
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
    length_units: list[str] = []
    volume_units: list[str] = []
    area_units: list[str] = []
    value_rows: list[tuple[str, float, float]] = []

    def norm_unit(value: object) -> str:
        return str(value or "").strip().lower().replace(" ", "")

    def norm_name(value: object) -> str:
        return str(value or "").strip().lower()

    def row_unit(row: Mapping[str, object], kind: str) -> str:
        key_groups = {
            "length": ("length_unit", "unit_length", "lengthUnit"),
            "volume": ("volume_unit", "unit_volume", "volumeUnit"),
            "area": ("area_unit", "unit_area", "areaUnit"),
        }
        units = row.get("units")
        for key in key_groups[kind]:
            value = row.get(key)
            if value:
                return norm_unit(value)
        if isinstance(units, Mapping):
            for key in key_groups[kind]:
                value = units.get(key) or units.get(kind)
                if value:
                    return norm_unit(value)
        return ""

    for index, row in enumerate(records):
        name = str(row.get("name", "")).strip()
        row_names.append(name)
        length_units.append(row_unit(row, "length"))
        volume_units.append(row_unit(row, "volume"))
        area_units.append(row_unit(row, "area"))
        try:
            volume = float(row.get("volume"))
            area = float(row.get("area"))
        except (TypeError, ValueError):
            row_issues.append({"index": index, "name": name, "reason": "missing volume or area"})
            continue
        volume_values.append(volume)
        area_values.append(area)
        value_rows.append((name, volume, area))
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

    expected_length_unit_norm = norm_unit(expected_length_unit)
    expected_volume_unit_norm = norm_unit(expected_volume_unit)
    expected_area_unit_norm = norm_unit(expected_area_unit)
    unit_sets = {
        "length": sorted({unit for unit in length_units if unit}),
        "volume": sorted({unit for unit in volume_units if unit}),
        "area": sorted({unit for unit in area_units if unit}),
    }
    expected_material_list = [str(name).strip() for name in expected_material_names if str(name).strip()]
    expected_material_norm = [norm_name(name) for name in expected_material_list]
    row_name_norm = [norm_name(name) for name in row_names]
    allowed_zero_norm = {
        norm_name(name)
        for name in allow_zero_measurement_names
        if str(name).strip()
    }
    missing_material_names = [
        name
        for name, normalized in zip(expected_material_list, expected_material_norm)
        if normalized not in set(row_name_norm)
    ]
    zero_volume_row_names = [name for name, volume, _area in value_rows if volume == 0.0]
    zero_area_row_names = [name for name, _volume, area in value_rows if area == 0.0]
    disallowed_zero_volume_row_names = [
        name
        for name in zero_volume_row_names
        if norm_name(name) not in allowed_zero_norm
    ]
    disallowed_zero_area_row_names = [
        name
        for name in zero_area_row_names
        if norm_name(name) not in allowed_zero_norm
    ]

    checks = {
        "rows_present": bool(records),
        "row_names_recorded": all(bool(name) for name in row_names),
        "row_names_unique": len(row_name_norm) == len(set(row_name_norm)),
        "all_rows_have_volume_area": len(volume_values) == len(records) and len(area_values) == len(records),
        "all_values_finite": all(isfinite(value) for value in volume_values + area_values),
        "all_volumes_positive": all(value > 0.0 for value in volume_values)
        or (
            bool(allowed_zero_norm)
            and all(value >= 0.0 for value in volume_values)
            and not disallowed_zero_volume_row_names
        ),
        "all_areas_positive": all(value > 0.0 for value in area_values)
        or (
            bool(allowed_zero_norm)
            and all(value >= 0.0 for value in area_values)
            and not disallowed_zero_area_row_names
        ),
        "expected_material_names_present": not expected_material_list or not missing_material_names,
        "total_volume_expected_ok": expected_volume_ok,
        "total_area_expected_ok": expected_area_ok,
        "bbox_size_expected_ok": expected_bbox_ok,
    }
    if any(unit_sets.values()):
        checks["unit_metadata_unique_when_present"] = all(len(units) <= 1 for units in unit_sets.values())
    if expected_length_unit_norm:
        checks["length_unit_expected_ok"] = unit_sets["length"] == [expected_length_unit_norm]
    if expected_volume_unit_norm:
        checks["volume_unit_expected_ok"] = unit_sets["volume"] == [expected_volume_unit_norm]
    if expected_area_unit_norm:
        checks["area_unit_expected_ok"] = unit_sets["area"] == [expected_area_unit_norm]
    return {
        "policy": "cubit_mass_property_sidecar_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "row_count": len(records),
        "row_names": row_names,
        "expected_material_names": expected_material_list,
        "missing_material_names": missing_material_names,
        "allowed_zero_measurement_names": sorted(allowed_zero_norm),
        "zero_volume_row_names": zero_volume_row_names,
        "zero_area_row_names": zero_area_row_names,
        "disallowed_zero_volume_row_names": disallowed_zero_volume_row_names,
        "disallowed_zero_area_row_names": disallowed_zero_area_row_names,
        "total_volume": total_volume,
        "total_area": total_area,
        "units": unit_sets,
        "expected_length_unit": expected_length_unit_norm or None,
        "expected_volume_unit": expected_volume_unit_norm or None,
        "expected_area_unit": expected_area_unit_norm or None,
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
            "When comparing Cubit, build123d, CST, or other CAD lanes, record units so mm^3 and m^3 rows cannot be mixed.",
            "For mixed transition meshes, keep material/block row names even when an allowed transition row has zero sidecar volume.",
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
    expected_export_output_artifact_id: str | None = None,
    expected_export_output_digest: str | None = None,
    require_export_output_artifact: bool = False,
    expected_export_observable_id: str | None = None,
    expected_export_observable_family: str | None = None,
    require_export_observable: bool = False,
    expected_vol_sidecar_schema_id: str | None = None,
    require_vol_sidecar_schema: bool = False,
    require_vol_sidecar_inventory_counts: bool = False,
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
    sidecar_rows = [row for kind, row in zip(kinds, records) if kind in {"vol_sidecar", "sidecar"}]

    def normalized_path(value: str) -> str:
        return value.replace("/", "\\").rstrip("\\").lower()

    def collect_values(rows: list[dict[str, object]], *names: str) -> list[str]:
        values: list[str] = []
        for row in rows:
            for name in names:
                if name in row and row[name] is not None:
                    text = str(row[name]).strip()
                    if text:
                        values.append(text)
        return values

    def collect_int_values(rows: list[dict[str, object]], *names: str) -> tuple[list[int], list[dict[str, object]]]:
        values: list[int] = []
        issues: list[dict[str, object]] = []
        for index, row in enumerate(rows):
            for name in names:
                if name not in row or row[name] is None:
                    continue
                try:
                    values.append(int(row[name]))
                except (TypeError, ValueError):
                    issues.append({"index": index, "field": name, "reason": "invalid integer"})
                break
        return values, issues

    expected_sidecars = {normalized_path(f"{path}.json") for path in vol_paths}
    actual_sidecars = {normalized_path(path) for path in sidecar_paths}
    inv = dict(inventory or {})
    inv_source = str(inv.get("source", "")).strip()
    inv_routing = str(inv.get("routing_hint", "")).strip()
    inventory_volume_elements = inv.get("volume_elements")
    inventory_points = inv.get("points")
    expected_export = None if expected_export_id is None else str(expected_export_id).strip()
    expected_geometry = None if expected_geometry_id is None else str(expected_geometry_id).strip()
    expected_hint = None if expected_routing_hint is None else str(expected_routing_hint).strip()
    expected_order_value = None if expected_order is None else int(expected_order)
    expected_output_artifact = (
        None if expected_export_output_artifact_id is None else str(expected_export_output_artifact_id).strip()
    )
    expected_output_digest = (
        None if expected_export_output_digest is None else str(expected_export_output_digest).strip()
    )
    expected_observable = (
        None if expected_export_observable_id is None else str(expected_export_observable_id).strip()
    )
    expected_observable_family = (
        None if expected_export_observable_family is None else str(expected_export_observable_family).strip()
    )
    expected_sidecar_schema = (
        None if expected_vol_sidecar_schema_id is None else str(expected_vol_sidecar_schema_id).strip()
    )
    output_artifact_ids = sorted(set(collect_values(
        records,
        "export_output_artifact_id",
        "package_output_artifact_id",
        "vol_output_artifact_id",
        "output_artifact_id",
    )))
    output_digests = sorted(set(collect_values(
        records,
        "export_output_digest",
        "package_output_digest",
        "vol_output_digest",
        "output_digest",
        "export_output_sha256",
        "vol_sha256",
    )))
    output_paths = sorted(set(collect_values(
        records,
        "export_output_path",
        "package_output_path",
        "vol_output_path",
        "output_path",
    )))
    observable_rows = records + ([inv] if inv else [])
    observable_ids = sorted(set(collect_values(
        observable_rows,
        "export_observable_id",
        "inventory_observable_id",
        "output_observable_id",
        "observable_id",
    )))
    observable_families = sorted(set(collect_values(
        observable_rows,
        "export_observable_family",
        "inventory_observable_family",
        "output_observable_family",
        "observable_family",
    )))
    output_artifact_id = output_artifact_ids[0] if output_artifact_ids else ""
    output_digest = output_digests[0] if output_digests else ""
    output_path = output_paths[0] if output_paths else ""
    observable_id = observable_ids[0] if observable_ids else ""
    observable_family = observable_families[0] if observable_families else ""
    sidecar_schema_ids = sorted(set(collect_values(
        sidecar_rows,
        "vol_sidecar_schema_id",
        "sidecar_schema_id",
        "inventory_schema_id",
        "schema_id",
        "schema",
    )))
    sidecar_schema_id = sidecar_schema_ids[0] if sidecar_schema_ids else ""
    sidecar_element_counts, sidecar_element_count_issues = collect_int_values(
        sidecar_rows,
        "n_elements",
        "element_count",
        "volume_elements",
    )
    sidecar_point_counts, sidecar_point_count_issues = collect_int_values(
        sidecar_rows,
        "n_points",
        "point_count",
        "points",
    )
    sidecar_orders, sidecar_order_issues = collect_int_values(sidecar_rows, "order")
    try:
        inventory_volume_count_value = (
            None if inventory_volume_elements is None else int(inventory_volume_elements)
        )
    except (TypeError, ValueError):
        inventory_volume_count_value = None
    try:
        inventory_point_count_value = None if inventory_points is None else int(inventory_points)
    except (TypeError, ValueError):
        inventory_point_count_value = None
    sidecar_count_required = bool(require_vol_sidecar_inventory_counts)
    sidecar_schema_required = bool(require_vol_sidecar_schema or expected_sidecar_schema is not None)

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
        "vol_sidecar_schema_id_consistent_when_present": len(sidecar_schema_ids) <= 1,
        "vol_sidecar_schema_id_recorded_when_required": (
            not sidecar_schema_required or bool(sidecar_schema_id)
        ),
        "vol_sidecar_schema_id_recorded_when_expected": (
            expected_sidecar_schema is None or bool(sidecar_schema_id)
        ),
        "expected_vol_sidecar_schema_id_matches": (
            expected_sidecar_schema is None or sidecar_schema_id == expected_sidecar_schema
        ),
        "orders_valid": not order_issues,
        "order_matches_expected": expected_order_value is None or (
            bool(orders) and all(order == expected_order_value for order in orders)
        ),
        "inventory_source_matches_vol": not inv_source
        or normalized_path(inv_source) in {normalized_path(path) for path in vol_paths},
        "inventory_routing_hint_matches_expected": expected_hint is None or inv_routing == expected_hint,
        "vol_sidecar_inventory_counts_recorded_when_required": (
            not sidecar_count_required
            or (bool(sidecar_element_counts) and bool(sidecar_point_counts))
        ),
        "vol_sidecar_inventory_count_fields_valid": (
            not sidecar_element_count_issues
            and not sidecar_point_count_issues
            and not sidecar_order_issues
        ),
        "vol_sidecar_element_count_consistent_when_present": len(set(sidecar_element_counts)) <= 1,
        "vol_sidecar_point_count_consistent_when_present": len(set(sidecar_point_counts)) <= 1,
        "vol_sidecar_order_consistent_when_present": len(set(sidecar_orders)) <= 1,
        "vol_sidecar_element_count_matches_inventory": (
            not sidecar_element_counts
            or inventory_volume_count_value is None
            or set(sidecar_element_counts) == {inventory_volume_count_value}
        ),
        "vol_sidecar_point_count_matches_inventory": (
            not sidecar_point_counts
            or inventory_point_count_value is None
            or set(sidecar_point_counts) == {inventory_point_count_value}
        ),
        "vol_sidecar_order_matches_expected": (
            expected_order_value is None
            or not sidecar_orders
            or set(sidecar_orders) == {expected_order_value}
        ),
        "export_output_artifact_id_consistent_when_present": len(output_artifact_ids) <= 1,
        "export_output_digest_consistent_when_present": len(output_digests) <= 1,
        "export_output_path_consistent_when_present": len(output_paths) <= 1,
        "export_output_artifact_id_recorded_when_required": (
            not require_export_output_artifact or bool(output_artifact_id)
        ),
        "export_output_digest_recorded_when_required": (
            not require_export_output_artifact or bool(output_digest)
        ),
        "export_output_path_recorded_when_required": (
            not require_export_output_artifact or bool(output_path)
        ),
        "export_output_artifact_id_recorded_when_expected": (
            expected_output_artifact is None or bool(output_artifact_id)
        ),
        "expected_export_output_artifact_id_matches": (
            expected_output_artifact is None or output_artifact_id == expected_output_artifact
        ),
        "export_output_digest_recorded_when_expected": (
            expected_output_digest is None or bool(output_digest)
        ),
        "expected_export_output_digest_matches": (
            expected_output_digest is None or output_digest == expected_output_digest
        ),
        "export_observable_id_consistent_when_present": len(observable_ids) <= 1,
        "export_observable_family_consistent_when_present": len(observable_families) <= 1,
        "export_observable_id_recorded_when_required": (
            not require_export_observable or bool(observable_id)
        ),
        "export_observable_family_recorded_when_required": (
            not require_export_observable or bool(observable_family)
        ),
        "export_observable_id_recorded_when_expected": (
            expected_observable is None or bool(observable_id)
        ),
        "expected_export_observable_id_matches": (
            expected_observable is None or observable_id == expected_observable
        ),
        "export_observable_family_recorded_when_expected": (
            expected_observable_family is None or bool(observable_family)
        ),
        "expected_export_observable_family_matches": (
            expected_observable_family is None or observable_family == expected_observable_family
        ),
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
        "export_output_artifact_id": output_artifact_id or None,
        "export_output_digest": output_digest or None,
        "export_output_path": output_path or None,
        "export_output_artifact_ids": output_artifact_ids,
        "export_output_digests": output_digests,
        "export_output_paths": output_paths,
        "expected_export_output_artifact_id": expected_output_artifact,
        "expected_export_output_digest": expected_output_digest,
        "require_export_output_artifact": bool(require_export_output_artifact),
        "export_observable_id": observable_id or None,
        "export_observable_family": observable_family or None,
        "export_observable_ids": observable_ids,
        "export_observable_families": observable_families,
        "expected_export_observable_id": expected_observable,
        "expected_export_observable_family": expected_observable_family,
        "require_export_observable": bool(require_export_observable),
        "inventory_source": inv_source,
        "inventory_routing_hint": inv_routing,
        "inventory_volume_elements": inventory_volume_count_value,
        "inventory_points": inventory_point_count_value,
        "vol_sidecar_element_counts": sorted(set(sidecar_element_counts)),
        "vol_sidecar_point_counts": sorted(set(sidecar_point_counts)),
        "vol_sidecar_orders": sorted(set(sidecar_orders)),
        "vol_sidecar_schema_id": sidecar_schema_id or None,
        "vol_sidecar_schema_ids": sidecar_schema_ids,
        "expected_vol_sidecar_schema_id": expected_sidecar_schema,
        "require_vol_sidecar_schema": sidecar_schema_required,
        "require_vol_sidecar_inventory_counts": sidecar_count_required,
        "expected_vol_sidecars": sorted(expected_sidecars),
        "actual_vol_sidecars": sorted(actual_sidecars),
        "order_issues": order_issues + sidecar_order_issues,
        "vol_sidecar_count_issues": sidecar_element_count_issues + sidecar_point_count_issues,
        "checks": checks,
        "notes": [
            "Use this before docs/panel notebooks or solver-ready runs consume a Cubit export package.",
            "A .vol, .vol.json, raw result, and mass-property sidecar should share export_id and geometry_id.",
            "The .vol.json sidecar schema id should be recorded so old inventory JSON layouts cannot masquerade as current sidecars.",
            "The .vol.json n_elements/n_points/order fields should match the parsed .vol inventory when those fields are recorded.",
            "The emitted export output artifact id, digest, and path identify the concrete package consumed by notebooks or solver-ready steps.",
            "The export observable id/family says what the package measures, such as a .vol inventory, quality distribution, or solver-ready mesh contract.",
        ],
    }


def cubit_headless_batch_quality_package_gate(
    batch_result: Mapping[str, object],
    quality_result: Mapping[str, object],
    *,
    expected_export_id: str | None = None,
    expected_geometry_id: str | None = None,
    expected_element_type: str | None = "hex",
    export_inventory: Mapping[str, object] | None = None,
    expected_routing_hint: str | None = "cubit_hex_or_mixed_path",
    expected_process_exit_policy: str | None = None,
) -> dict[str, object]:
    """Check that a headless Cubit run and quality gate describe one mesh.

    A live Coreform/Cubit batch can end with noisy process status even when the
    archived JSON says the mesh was produced and quality checks passed.  This
    replay gate keeps the session/run identity next to the quality summary so a
    notebook or validation case does not pair a mesh with stale quality rows.
    """

    batch = dict(batch_result)
    quality = dict(quality_result)
    export_id = str(batch.get("export_id", "")).strip()
    quality_export_id = str(quality.get("export_id", "")).strip()
    geometry_id = str(batch.get("geometry_id", "")).strip()
    quality_geometry_id = str(quality.get("geometry_id", "")).strip()
    command_line = str(batch.get("command_line", "")).strip()
    version = str(batch.get("version", batch.get("cubit_version", ""))).strip()
    process_mode = str(batch.get("process_mode", batch.get("execution_mode", ""))).strip().lower().replace("-", "_")
    batch_script = str(batch.get("batch_script", batch.get("script_path", ""))).strip()
    journal_policy = str(batch.get("journal_policy", "")).strip().lower().replace("-", "_")
    exit_code_note = str(batch.get("process_exit_note", batch.get("exit_code_note", ""))).strip().lower()
    process_exit_policy = str(
        batch.get("process_exit_policy", batch.get("exit_code_policy", ""))
    ).strip().lower().replace("-", "_").replace(" ", "_")
    gui_daemon_raw = batch.get("gui_daemon", batch.get("allow_gui_daemon", None))
    exit_code_raw = batch.get("exit_code", batch.get("process_exit_code", None))
    solver_ready_claimed = bool(
        batch.get("solver_ready_claimed", batch.get("solver_ready", False))
    )
    element = str(quality.get("element_type", "")).strip().lower().replace(" ", "")
    expected_export = None if expected_export_id is None else str(expected_export_id).strip()
    expected_geometry = None if expected_geometry_id is None else str(expected_geometry_id).strip()
    expected_element = None if expected_element_type is None else str(expected_element_type).strip().lower().replace(" ", "")
    expected_hint = None if expected_routing_hint is None else str(expected_routing_hint).strip()
    expected_exit_policy = (
        None
        if expected_process_exit_policy is None
        else str(expected_process_exit_policy).strip().lower().replace("-", "_").replace(" ", "_")
    )
    output_paths = batch.get("output_paths", batch.get("paths", ())) or ()
    if isinstance(output_paths, (str, bytes)):
        paths = [str(output_paths)]
    else:
        paths = [str(path) for path in output_paths]
    inventory = None if export_inventory is None else dict(export_inventory)
    inventory_volume_counts = {}
    inventory_surface_counts = {}
    inventory_source = ""
    inventory_routing_hint = ""
    inventory_is_tri_tet_only = False
    if inventory is not None:
        counts = inventory.get("volume_kind_counts", {})
        if not isinstance(counts, Mapping):
            raise ValueError("export_inventory['volume_kind_counts'] must be a mapping when provided")
        inventory_volume_counts = {str(key).strip().lower(): int(value) for key, value in counts.items()}
        surface_counts = inventory.get("surface_kind_counts", {})
        if isinstance(surface_counts, Mapping):
            inventory_surface_counts = {
                str(key).strip().lower(): int(value) for key, value in surface_counts.items()
            }
        inventory_source = str(inventory.get("source", "")).strip()
        inventory_routing_hint = str(inventory.get("routing_hint", "")).strip()
        inventory_is_tri_tet_only = bool(inventory.get("is_tri_tet_only", False)) or (
            bool(inventory_volume_counts)
            and set(inventory_volume_counts).issubset({"tet"})
            and set(inventory_surface_counts).issubset({"triangle"})
        )
    inventory_volume_total = sum(inventory_volume_counts.values())
    expected_inventory_count = (
        int(quality.get("count", 0) or 0)
        if expected_element and expected_element in inventory_volume_counts
        else None
    )
    try:
        exit_code = None if exit_code_raw is None else int(exit_code_raw)
    except (TypeError, ValueError):
        exit_code = None
    if isinstance(gui_daemon_raw, str):
        gui_daemon_enabled = gui_daemon_raw.strip().lower() in {"1", "true", "yes", "enabled", "on"}
    else:
        gui_daemon_enabled = bool(gui_daemon_raw) if gui_daemon_raw is not None else None
    process_evidence_requested = any(
        value is not None and value != ""
        for value in (
            process_mode,
            batch_script,
            journal_policy,
            exit_code_note,
            gui_daemon_raw,
            exit_code_raw,
            process_exit_policy,
        )
    )
    documented_headless_warning_exit = (
        exit_code in {3}
        and "headless" in exit_code_note
        and ("startup" in exit_code_note or "warning" in exit_code_note)
    )
    documented_nonzero_exit = (
        exit_code == 0
        or documented_headless_warning_exit
        or (
            exit_code in {1, 3}
            and process_exit_policy == "artifact_evidence_over_process_exit"
            and bool(exit_code_note)
            and not solver_ready_claimed
        )
    )

    checks = {
        "batch_passed": bool(batch.get("pass", False)) or str(batch.get("status", "")).strip().lower() in {"ok", "pass", "passed"},
        "quality_passed": str(quality.get("status", "")).strip().lower() == "ok",
        "export_id_recorded": bool(export_id and quality_export_id),
        "export_id_matches_quality": bool(export_id) and export_id == quality_export_id,
        "geometry_id_recorded": bool(geometry_id and quality_geometry_id),
        "geometry_id_matches_quality": bool(geometry_id) and geometry_id == quality_geometry_id,
        "export_id_matches_expected": expected_export is None or export_id == expected_export,
        "geometry_id_matches_expected": expected_geometry is None or geometry_id == expected_geometry,
        "headless_command_recorded": "-nographics" in command_line.lower() or bool(batch.get("headless", False)),
        "version_recorded": bool(version),
        "output_paths_recorded": bool(paths),
        "quality_element_type_matches_expected": expected_element is None or element == expected_element,
        "quality_count_positive": int(quality.get("count", 0) or 0) > 0,
    }
    if process_evidence_requested:
        command_lower = command_line.lower()
        checks.update({
            "process_mode_is_headless_batch": process_mode in {"headless_batch", "batch_headless"},
            "nographics_flag_present": "-nographics" in command_lower,
            "batch_flag_present": "-batch" in command_lower,
            "gui_daemon_disabled": gui_daemon_enabled is False,
            "batch_script_recorded": bool(batch_script),
            "batch_script_appears_in_command": bool(batch_script) and batch_script.lower() in command_lower,
            "process_exit_code_recorded": exit_code is not None,
            "process_exit_policy_recorded_when_nonzero": (
                exit_code is None
                or exit_code == 0
                or documented_headless_warning_exit
                or bool(process_exit_policy)
            ),
            "expected_process_exit_policy_matches": (
                expected_exit_policy is None or process_exit_policy == expected_exit_policy
            ),
            "nonzero_exit_does_not_claim_solver_ready": (
                exit_code is None or exit_code == 0 or not solver_ready_claimed
            ),
            "process_exit_code_success_or_documented": documented_nonzero_exit,
        })
        if journal_policy:
            checks["journal_policy_records_batch_not_gui_daemon"] = (
                "batch" in journal_policy or "nojournal" in journal_policy
            )
    if inventory is not None:
        checks.update({
            "export_inventory_recorded": True,
            "export_inventory_source_in_output_paths": bool(inventory_source) and inventory_source in paths,
            "export_inventory_volume_elements_positive": inventory_volume_total > 0,
            "export_inventory_routing_hint_matches_expected": expected_hint is None or inventory_routing_hint == expected_hint,
            "export_inventory_count_matches_quality": (
                expected_inventory_count is None
                or inventory_volume_counts.get(expected_element, -1) == expected_inventory_count
            ),
            "export_inventory_contains_quality_element": (
                expected_element is None
                or inventory_volume_counts.get(expected_element, 0) > 0
            ),
            "export_inventory_not_tri_tet_only_for_cubit_hex_route": (
                expected_hint != "cubit_hex_or_mixed_path"
                or not inventory_is_tri_tet_only
            ),
        })
    return {
        "policy": "cubit_headless_batch_quality_package_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "export_id": export_id,
        "geometry_id": geometry_id,
        "expected_export_id": expected_export,
        "expected_geometry_id": expected_geometry,
        "command_line": command_line,
        "version": version,
        "process_mode": process_mode,
        "batch_script": batch_script,
        "journal_policy": journal_policy,
        "exit_code_note": exit_code_note,
        "process_exit_policy": process_exit_policy,
        "expected_process_exit_policy": expected_exit_policy,
        "solver_ready_claimed": solver_ready_claimed,
        "gui_daemon": gui_daemon_enabled,
        "exit_code": exit_code,
        "output_paths": paths,
        "quality_policy": quality.get("policy"),
        "quality_element_type": element,
        "quality_count": int(quality.get("count", 0) or 0),
        "export_inventory_source": inventory_source,
        "export_inventory_volume_kind_counts": inventory_volume_counts,
        "export_inventory_routing_hint": inventory_routing_hint,
        "export_inventory_is_tri_tet_only": inventory_is_tri_tet_only,
        "expected_routing_hint": expected_hint,
        "checks": checks,
        "notes": [
            "Use this after a headless Cubit run so quality rows cannot be reused with a stale mesh package.",
            "The gate treats headless execution and no-GUI policy as part of mesh evidence provenance.",
            "When a .vol inventory is available, bind its volume-element count and routing hint to the same headless package.",
            "Cubit hex-led evidence must not be silently routed into the Netgen/MATLAB tri-tet-only education path.",
            "When process evidence is recorded, require headless_batch mode, -nographics -batch, disabled GUI daemon, batch script identity, and zero exit code or a documented process_exit_policy that keeps solver-ready claims separate.",
        ],
    }


def cubit_mesh_quality_ledger_identity_gate(
    quality_ledger: Mapping[str, object],
    *,
    quality_gate: Mapping[str, object] | None = None,
    export_package_gate: Mapping[str, object] | None = None,
    headless_batch_quality_gate: Mapping[str, object] | None = None,
    inventory: Mapping[str, object] | None = None,
    expected_quality_artifact_id: str | None = None,
    expected_quality_digest: str | None = None,
    expected_metric_set_id: str | None = None,
    expected_export_id: str | None = None,
    expected_geometry_id: str | None = None,
    expected_mesh_artifact_id: str | None = None,
    expected_mesh_digest: str | None = None,
    expected_version: str | None = None,
    expected_parameter_set_artifact_id: str | None = None,
    expected_parameter_set_digest: str | None = None,
    expected_parameter_set_path: str | None = None,
    expected_objective_observable_id: str | None = None,
    expected_objective_observable_family: str | None = None,
    expected_quality_postprocess_row_convention_schema_id: str | None = None,
    expected_quality_component_basis_schema_id: str | None = None,
    expected_routing_hint: str | None = "cubit_hex_or_mixed_path",
    expected_element_type_counts: Mapping[str, int] | None = None,
    min_scaled_jacobian: float = 0.2,
    require_hex_or_mixed_route: bool = True,
    require_execution_metadata: bool = False,
    require_parameter_set_artifact: bool = False,
    require_quality_postprocess_row_convention_schema: bool = False,
    require_quality_component_basis_schema: bool = False,
    require_timing_breakdown: bool = False,
    min_timing_stages: int = 4,
) -> dict[str, object]:
    """Bind Cubit mesh-quality evidence to the concrete exported mesh.

    Cubit can produce several valid-looking rows for one slot: the exported
    ``.vol``, its package identity, a headless batch result, and one or more
    quality distributions.  This gate makes the reusable quality ledger carry
    a stable artifact id, digest, metric-set id, mesh id/digest, and element
    inventory before notebooks or solver-ready steps reuse the quality row.
    For parametric mesh studies, the meshing/design parameter set and objective
    observable identity should travel with the quality ledger too.
    """

    ledger = dict(quality_ledger)
    quality = dict(quality_gate or {})
    export_package = dict(export_package_gate or {})
    headless = dict(headless_batch_quality_gate or {})
    inv = dict(inventory or {})

    def first(row: Mapping[str, object], *names: str) -> str:
        for name in names:
            value = row.get(name)
            if value is not None:
                text = str(value).strip()
                if text:
                    return text
        return ""

    def normalized(value: object) -> str:
        return str(value or "").strip().lower().replace(" ", "_").replace("-", "_")

    def counts_from(row: Mapping[str, object], *names: str) -> dict[str, int]:
        for name in names:
            value = row.get(name)
            if isinstance(value, Mapping):
                return {normalized(key): int(count) for key, count in value.items()}
        return {}

    def list_or_scalar_first(row: Mapping[str, object], *names: str) -> str:
        for name in names:
            value = row.get(name)
            if isinstance(value, (list, tuple)) and value:
                text = str(value[0]).strip()
                if text:
                    return text
            if value is not None and not isinstance(value, (list, tuple)):
                text = str(value).strip()
                if text:
                    return text
        return ""

    def nested_mapping(row: Mapping[str, object], *names: str) -> dict[str, object]:
        for name in names:
            value = row.get(name)
            if isinstance(value, Mapping):
                return dict(value)
        return {}

    def string_values(rows: Iterable[Mapping[str, object]], *names: str) -> list[str]:
        values = []
        for row in rows:
            for name in names:
                value = row.get(name)
                if value is None:
                    continue
                if isinstance(value, (list, tuple)):
                    raw_values = value
                else:
                    raw_values = (value,)
                for raw in raw_values:
                    text = str(raw).strip()
                    if text:
                        values.append(text)
        return values

    def unique_strings(values: Iterable[str]) -> list[str]:
        return sorted({str(value).strip() for value in values if str(value).strip()})

    def parseable_time(value: object) -> bool:
        text = str(value or "").strip()
        if not text:
            return False
        try:
            datetime.fromisoformat(text.replace("Z", "+00:00"))
            return True
        except ValueError:
            return False

    def nonnegative_float(value: object) -> float | None:
        try:
            number = float(value)
        except (TypeError, ValueError):
            return None
        return number if isfinite(number) and number >= 0.0 else None

    quality_artifact_id = first(
        ledger,
        "mesh_quality_artifact_id",
        "quality_artifact_id",
        "quality_output_artifact_id",
        "artifact_id",
    )
    quality_digest = first(
        ledger,
        "mesh_quality_digest",
        "quality_digest",
        "quality_output_digest",
        "digest",
    )
    quality_path = first(
        ledger,
        "mesh_quality_path",
        "quality_path",
        "quality_output_path",
        "path",
    )
    metric_set_id = first(ledger, "quality_metric_set_id", "metric_set_id", "quality_metric_id")
    export_id = first(ledger, "export_id")
    geometry_id = first(ledger, "geometry_id")
    mesh_artifact_id = first(
        ledger,
        "mesh_artifact_id",
        "export_output_artifact_id",
        "vol_artifact_id",
    )
    mesh_digest = first(
        ledger,
        "mesh_digest",
        "export_output_digest",
        "vol_digest",
        "vol_sha256",
    )
    routing_hint = first(ledger, "routing_hint", "mesh_route")
    version = first(ledger, "version", "cubit_version", "coreform_cubit_version")
    created_at_utc = first(ledger, "created_at_utc", "created_at", "run_date_utc")
    elapsed_s = nonnegative_float(ledger.get("elapsed_s", ledger.get("run_duration_s")))
    timing_raw = ledger.get("timing_breakdown_s", ledger.get("timing_breakdown", {}))
    timing_breakdown = {}
    if isinstance(timing_raw, Mapping):
        for key, value in timing_raw.items():
            number = nonnegative_float(value)
            if number is not None:
                timing_breakdown[str(key)] = number
    metric = normalized(first(ledger, "metric", "quality_metric"))
    parameter_set = nested_mapping(ledger, "parameter_set", "parameterSet", "mesh_parameter_set")
    optimization = nested_mapping(ledger, "optimization", "optimizer", "mesh_optimization")
    objective = nested_mapping(ledger, "objective", "objective_metadata", "objectiveMetadata")
    parameter_rows = (ledger, parameter_set, optimization)
    objective_rows = (ledger, optimization, objective)
    parameter_set_artifact_ids = unique_strings(string_values(
        parameter_rows,
        "parameter_set_artifact_id",
        "parameterSetArtifactId",
        "mesh_parameter_set_artifact_id",
        "meshParameterSetArtifactId",
        "sizing_parameter_set_artifact_id",
        "sizingParameterSetArtifactId",
        "optimization_parameter_set_artifact_id",
        "optimizationParameterSetArtifactId",
    ))
    parameter_set_digests = unique_strings(string_values(
        parameter_rows,
        "parameter_set_digest",
        "parameterSetDigest",
        "parameter_set_sha256",
        "parameterSetSha256",
        "mesh_parameter_set_digest",
        "meshParameterSetDigest",
        "sizing_parameter_set_digest",
        "sizingParameterSetDigest",
        "optimization_parameter_set_digest",
        "optimizationParameterSetDigest",
    ))
    parameter_set_paths = unique_strings(string_values(
        parameter_rows,
        "parameter_set_path",
        "parameterSetPath",
        "mesh_parameter_set_path",
        "meshParameterSetPath",
        "sizing_parameter_set_path",
        "sizingParameterSetPath",
        "optimization_parameter_set_path",
        "optimizationParameterSetPath",
    ))
    objective_observable_ids = unique_strings(string_values(
        objective_rows,
        "objective_observable_id",
        "objectiveObservableId",
        "objective_id",
        "objectiveId",
        "mesh_objective_observable_id",
        "meshObjectiveObservableId",
        "quality_objective_id",
        "qualityObjectiveId",
        "target_observable_id",
        "targetObservableId",
    ))
    objective_observable_families = sorted({
        normalized(value)
        for value in string_values(
            objective_rows,
            "objective_observable_family",
            "objectiveObservableFamily",
            "objective_family",
            "objectiveFamily",
            "mesh_objective_observable_family",
            "meshObjectiveObservableFamily",
            "quality_objective_family",
            "qualityObjectiveFamily",
            "target_observable_family",
            "targetObservableFamily",
        )
    })
    parameter_set_artifact_id = (
        parameter_set_artifact_ids[0] if parameter_set_artifact_ids else ""
    )
    parameter_set_digest = parameter_set_digests[0] if parameter_set_digests else ""
    parameter_set_path = parameter_set_paths[0] if parameter_set_paths else ""
    objective_observable_id = (
        objective_observable_ids[0] if objective_observable_ids else ""
    )
    objective_observable_family = (
        objective_observable_families[0] if objective_observable_families else ""
    )
    row_convention_rows = (ledger, quality, headless)
    quality_postprocess_row_convention_schema_ids = unique_strings(string_values(
        row_convention_rows,
        "mesh_quality_postprocess_row_convention_schema_id",
        "meshQualityPostprocessRowConventionSchemaId",
        "quality_postprocess_row_convention_schema_id",
        "qualityPostprocessRowConventionSchemaId",
        "postprocess_row_convention_schema_id",
        "postprocessRowConventionSchemaId",
        "quality_row_convention_schema_id",
        "qualityRowConventionSchemaId",
    ))
    quality_postprocess_row_convention_schema_id = (
        quality_postprocess_row_convention_schema_ids[0]
        if quality_postprocess_row_convention_schema_ids
        else ""
    )
    expected_quality_postprocess_row_convention_schema = (
        ""
        if expected_quality_postprocess_row_convention_schema_id is None
        else str(expected_quality_postprocess_row_convention_schema_id).strip()
    )
    quality_postprocess_row_convention_schema_required = bool(
        require_quality_postprocess_row_convention_schema
        or expected_quality_postprocess_row_convention_schema
    )
    component_basis_rows = (ledger, quality, headless)
    quality_component_basis_schema_ids = unique_strings(string_values(
        component_basis_rows,
        "mesh_quality_component_basis_schema_id",
        "meshQualityComponentBasisSchemaId",
        "quality_component_basis_schema_id",
        "qualityComponentBasisSchemaId",
        "element_quality_component_basis_schema_id",
        "elementQualityComponentBasisSchemaId",
        "component_basis_schema_id",
        "componentBasisSchemaId",
    ))
    quality_component_basis_schema_id = (
        quality_component_basis_schema_ids[0]
        if quality_component_basis_schema_ids
        else ""
    )
    expected_quality_component_basis_schema = (
        ""
        if expected_quality_component_basis_schema_id is None
        else str(expected_quality_component_basis_schema_id).strip()
    )
    quality_component_basis_schema_required = bool(
        require_quality_component_basis_schema
        or expected_quality_component_basis_schema
    )
    min_j_raw = ledger.get("min_scaled_jacobian", ledger.get("min_jacobian", ledger.get("min")))
    negative_j_raw = ledger.get("negative_jacobian_count", ledger.get("negative_jacobian_count_after", 0))
    element_counts = counts_from(ledger, "element_type_counts", "volume_kind_counts")
    expected_counts = {
        normalized(key): int(value)
        for key, value in dict(expected_element_type_counts or {}).items()
    }

    try:
        min_j = float(min_j_raw) if min_j_raw is not None else None
    except (TypeError, ValueError):
        min_j = None
    try:
        negative_j = int(negative_j_raw)
    except (TypeError, ValueError):
        negative_j = -1

    inv_counts = counts_from(inv, "volume_kind_counts", "element_type_counts")
    inv_routing = first(inv, "routing_hint", "mesh_route")
    inv_is_tri_tet_only = bool(inv.get("is_tri_tet_only", False)) or (
        bool(inv_counts) and set(inv_counts).issubset({"tet"})
    )
    quality_count = int(quality.get("count", 0) or 0) if quality else 0
    quality_metric = normalized(first(quality, "metric", "quality_metric"))
    quality_element = normalized(first(quality, "element_type"))
    export_package_export_id = list_or_scalar_first(export_package, "export_id", "export_ids")
    export_package_geometry_id = list_or_scalar_first(export_package, "geometry_id", "geometry_ids")
    export_package_mesh_artifact_id = first(export_package, "export_output_artifact_id", "mesh_artifact_id")
    export_package_mesh_digest = first(export_package, "export_output_digest", "mesh_digest")
    headless_export_id = first(headless, "export_id")
    headless_geometry_id = first(headless, "geometry_id")

    threshold = float(min_scaled_jacobian)
    if threshold <= 0.0:
        raise ValueError("min_scaled_jacobian must be > 0")

    checks = {
        "quality_artifact_id_recorded": bool(quality_artifact_id),
        "quality_digest_recorded": bool(quality_digest),
        "quality_path_recorded": bool(quality_path),
        "metric_set_id_recorded": bool(metric_set_id),
        "export_id_recorded": bool(export_id),
        "geometry_id_recorded": bool(geometry_id),
        "mesh_artifact_id_recorded": bool(mesh_artifact_id),
        "mesh_digest_recorded": bool(mesh_digest),
        "metric_is_quality_lower_bound": metric in {"scaled_jacobian", "jacobian"},
        "min_scaled_jacobian_recorded": min_j is not None and isfinite(min_j),
        "min_scaled_jacobian_above_threshold": min_j is not None and min_j >= threshold,
        "negative_jacobian_count_recorded": negative_j >= 0,
        "negative_jacobian_count_zero": negative_j == 0,
        "element_type_counts_recorded": bool(element_counts),
        "element_type_counts_positive": bool(element_counts) and all(value > 0 for value in element_counts.values()),
        "routing_hint_matches_expected": expected_routing_hint is None or routing_hint == expected_routing_hint,
        "expected_quality_artifact_id_matches": (
            expected_quality_artifact_id is None or quality_artifact_id == expected_quality_artifact_id
        ),
        "expected_quality_digest_matches": (
            expected_quality_digest is None or quality_digest == expected_quality_digest
        ),
        "expected_metric_set_id_matches": (
            expected_metric_set_id is None or metric_set_id == expected_metric_set_id
        ),
        "expected_export_id_matches": expected_export_id is None or export_id == expected_export_id,
        "expected_geometry_id_matches": expected_geometry_id is None or geometry_id == expected_geometry_id,
        "expected_mesh_artifact_id_matches": (
            expected_mesh_artifact_id is None or mesh_artifact_id == expected_mesh_artifact_id
        ),
        "expected_mesh_digest_matches": expected_mesh_digest is None or mesh_digest == expected_mesh_digest,
        "expected_version_matches": expected_version is None or version == expected_version,
        "parameter_set_artifact_id_consistent_when_present": (
            len(parameter_set_artifact_ids) <= 1
        ),
        "parameter_set_digest_consistent_when_present": len(parameter_set_digests) <= 1,
        "parameter_set_path_consistent_when_present": len(parameter_set_paths) <= 1,
        "objective_observable_id_consistent_when_present": (
            len(objective_observable_ids) <= 1
        ),
        "objective_observable_family_consistent_when_present": (
            len(objective_observable_families) <= 1
        ),
        "parameter_set_artifact_id_recorded_when_required": (
            not require_parameter_set_artifact or bool(parameter_set_artifact_id)
        ),
        "parameter_set_digest_recorded_when_required": (
            not require_parameter_set_artifact or bool(parameter_set_digest)
        ),
        "parameter_set_path_recorded_when_required": (
            not require_parameter_set_artifact or bool(parameter_set_path)
        ),
        "parameter_set_artifact_id_recorded_when_expected": (
            expected_parameter_set_artifact_id is None or bool(parameter_set_artifact_id)
        ),
        "expected_parameter_set_artifact_id_matches": (
            expected_parameter_set_artifact_id is None
            or parameter_set_artifact_id == str(expected_parameter_set_artifact_id).strip()
        ),
        "parameter_set_digest_recorded_when_expected": (
            expected_parameter_set_digest is None or bool(parameter_set_digest)
        ),
        "expected_parameter_set_digest_matches": (
            expected_parameter_set_digest is None
            or parameter_set_digest == str(expected_parameter_set_digest).strip()
        ),
        "parameter_set_path_recorded_when_expected": (
            expected_parameter_set_path is None or bool(parameter_set_path)
        ),
        "expected_parameter_set_path_matches": (
            expected_parameter_set_path is None
            or parameter_set_path == str(expected_parameter_set_path).strip()
        ),
        "objective_observable_id_recorded_when_expected": (
            expected_objective_observable_id is None or bool(objective_observable_id)
        ),
        "expected_objective_observable_id_matches": (
            expected_objective_observable_id is None
            or objective_observable_id == str(expected_objective_observable_id).strip()
        ),
        "objective_observable_family_recorded_when_expected": (
            expected_objective_observable_family is None or bool(objective_observable_family)
        ),
        "expected_objective_observable_family_matches": (
            expected_objective_observable_family is None
            or objective_observable_family == normalized(expected_objective_observable_family)
        ),
        "quality_postprocess_row_convention_schema_id_consistent_when_present": (
            len(quality_postprocess_row_convention_schema_ids) <= 1
        ),
        "quality_postprocess_row_convention_schema_id_recorded_when_required": (
            not quality_postprocess_row_convention_schema_required
            or bool(quality_postprocess_row_convention_schema_id)
        ),
        "quality_postprocess_row_convention_schema_id_recorded_when_expected": (
            not expected_quality_postprocess_row_convention_schema
            or bool(quality_postprocess_row_convention_schema_id)
        ),
        "expected_quality_postprocess_row_convention_schema_id_matches": (
            not expected_quality_postprocess_row_convention_schema
            or quality_postprocess_row_convention_schema_id
            == expected_quality_postprocess_row_convention_schema
        ),
        "quality_component_basis_schema_id_consistent_when_present": (
            len(quality_component_basis_schema_ids) <= 1
        ),
        "quality_component_basis_schema_id_recorded_when_required": (
            not quality_component_basis_schema_required
            or bool(quality_component_basis_schema_id)
        ),
        "quality_component_basis_schema_id_recorded_when_expected": (
            not expected_quality_component_basis_schema
            or bool(quality_component_basis_schema_id)
        ),
        "expected_quality_component_basis_schema_id_matches": (
            not expected_quality_component_basis_schema
            or quality_component_basis_schema_id
            == expected_quality_component_basis_schema
        ),
        "expected_element_type_counts_match": not expected_counts or element_counts == expected_counts,
    }
    if require_execution_metadata or version or created_at_utc or elapsed_s is not None:
        checks.update({
            "created_at_utc_recorded_when_required": (
                not require_execution_metadata or bool(created_at_utc)
            ),
            "created_at_utc_parseable_when_present": (
                not created_at_utc or parseable_time(created_at_utc)
            ),
            "version_recorded_when_required": (
                not require_execution_metadata or bool(version)
            ),
            "elapsed_s_recorded_when_required": (
                not require_execution_metadata or elapsed_s is not None
            ),
            "elapsed_s_finite_nonnegative_when_present": elapsed_s is None or elapsed_s >= 0.0,
        })
    if require_timing_breakdown or timing_breakdown:
        checks.update({
            "timing_breakdown_recorded_when_required": (
                not require_timing_breakdown or bool(timing_breakdown)
            ),
            "timing_breakdown_has_required_stage_count": (
                len(timing_breakdown) >= int(min_timing_stages)
            ),
            "timing_breakdown_values_finite_nonnegative": (
                bool(timing_breakdown)
                and all(isfinite(value) and value >= 0.0 for value in timing_breakdown.values())
            ),
            "timing_breakdown_total_within_elapsed_when_present": (
                elapsed_s is None or sum(timing_breakdown.values()) <= elapsed_s * 1.5 + 1e-12
            ),
        })
    if require_hex_or_mixed_route:
        checks["hex_or_mixed_volume_family_present"] = any(
            element_counts.get(kind, 0) > 0 for kind in ("hex", "pyramid", "wedge")
        )
        checks["not_tri_tet_only_for_cubit_quality_ledger"] = not inv_is_tri_tet_only

    if quality:
        checks.update({
            "quality_distribution_gate_ok": (
                quality.get("policy") == "cubit_quality_distribution_gate"
                and quality.get("status") == "ok"
            ),
            "quality_metric_matches_ledger": not quality_metric or quality_metric == metric,
            "quality_count_matches_ledger_counts": quality_count == sum(element_counts.values()),
            "quality_bad_count_zero": int(quality.get("bad_count", 0) or 0) == 0,
            "quality_element_present_in_ledger": not quality_element or element_counts.get(quality_element, 0) > 0,
        })
    if export_package:
        checks.update({
            "export_package_gate_ok": (
                export_package.get("policy") == "cubit_export_package_identity_gate"
                and export_package.get("status") == "ok"
            ),
            "export_package_export_id_matches_ledger": (
                not export_package_export_id or export_package_export_id == export_id
            ),
            "export_package_geometry_id_matches_ledger": (
                not export_package_geometry_id or export_package_geometry_id == geometry_id
            ),
            "mesh_artifact_matches_export_package_when_present": (
                not export_package_mesh_artifact_id or export_package_mesh_artifact_id == mesh_artifact_id
            ),
            "mesh_digest_matches_export_package_when_present": (
                not export_package_mesh_digest or export_package_mesh_digest == mesh_digest
            ),
        })
    if headless:
        checks.update({
            "headless_batch_quality_gate_ok": (
                headless.get("policy") == "cubit_headless_batch_quality_package_gate"
                and headless.get("status") == "ok"
            ),
            "headless_export_id_matches_ledger": not headless_export_id or headless_export_id == export_id,
            "headless_geometry_id_matches_ledger": not headless_geometry_id or headless_geometry_id == geometry_id,
        })
    if inv:
        checks.update({
            "inventory_volume_kind_counts_recorded": bool(inv_counts),
            "inventory_routing_hint_matches_ledger": not inv_routing or inv_routing == routing_hint,
            "inventory_counts_match_ledger_counts": (
                bool(inv_counts)
                and all(inv_counts.get(kind) == count for kind, count in element_counts.items())
            ),
        })

    return {
        "policy": "cubit_mesh_quality_ledger_identity_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "quality_artifact_id": quality_artifact_id,
        "quality_digest": quality_digest,
        "quality_path": quality_path,
        "metric_set_id": metric_set_id,
        "export_id": export_id,
        "geometry_id": geometry_id,
        "mesh_artifact_id": mesh_artifact_id,
        "mesh_digest": mesh_digest,
        "routing_hint": routing_hint,
        "version": version,
        "created_at_utc": created_at_utc,
        "elapsed_s": elapsed_s,
        "timing_breakdown_s": timing_breakdown,
        "parameter_set_artifact_id": parameter_set_artifact_id or None,
        "parameter_set_digest": parameter_set_digest or None,
        "parameter_set_path": parameter_set_path or None,
        "parameter_set_artifact_ids": parameter_set_artifact_ids,
        "parameter_set_digests": parameter_set_digests,
        "parameter_set_paths": parameter_set_paths,
        "objective_observable_id": objective_observable_id or None,
        "objective_observable_family": objective_observable_family or None,
        "objective_observable_ids": objective_observable_ids,
        "objective_observable_families": objective_observable_families,
        "expected_parameter_set_artifact_id": (
            None if expected_parameter_set_artifact_id is None else str(expected_parameter_set_artifact_id).strip()
        ),
        "expected_parameter_set_digest": (
            None if expected_parameter_set_digest is None else str(expected_parameter_set_digest).strip()
        ),
        "expected_parameter_set_path": (
            None if expected_parameter_set_path is None else str(expected_parameter_set_path).strip()
        ),
        "expected_objective_observable_id": (
            None if expected_objective_observable_id is None else str(expected_objective_observable_id).strip()
        ),
        "expected_objective_observable_family": (
            None if expected_objective_observable_family is None else normalized(expected_objective_observable_family)
        ),
        "quality_postprocess_row_convention_schema_id": (
            quality_postprocess_row_convention_schema_id or None
        ),
        "quality_postprocess_row_convention_schema_ids": (
            quality_postprocess_row_convention_schema_ids
        ),
        "expected_quality_postprocess_row_convention_schema_id": (
            expected_quality_postprocess_row_convention_schema or None
        ),
        "quality_postprocess_row_convention_schema_required": (
            quality_postprocess_row_convention_schema_required
        ),
        "quality_component_basis_schema_id": (
            quality_component_basis_schema_id or None
        ),
        "quality_component_basis_schema_ids": quality_component_basis_schema_ids,
        "expected_quality_component_basis_schema_id": (
            expected_quality_component_basis_schema or None
        ),
        "quality_component_basis_schema_required": quality_component_basis_schema_required,
        "parameter_set_artifact_required": bool(require_parameter_set_artifact),
        "metric": metric,
        "min_scaled_jacobian": min_j,
        "min_scaled_jacobian_threshold": threshold,
        "negative_jacobian_count": negative_j,
        "element_type_counts": element_counts,
        "expected_element_type_counts": expected_counts,
        "inventory_volume_kind_counts": inv_counts,
        "inventory_routing_hint": inv_routing,
        "inventory_is_tri_tet_only": inv_is_tri_tet_only,
        "quality_gate_policy": quality.get("policy"),
        "export_package_policy": export_package.get("policy"),
        "headless_batch_quality_policy": headless.get("policy"),
        "checks": checks,
        "notes": [
            "Use this after Cubit quality distribution replay but before notebook or solver-ready reuse.",
            "A mesh-quality row is reusable only with its own artifact id, digest, metric-set id, mesh id/digest, route, and element counts.",
            "For Cubit hex-led lanes, reject tri/tet-only inventories even when a stale quality row looks healthy.",
            "When execution metadata is required, keep version, run timestamp, elapsed time, and the dominant timing stages with the quality ledger.",
            "For parametric meshing or optimization, keep the mesh parameter-set artifact and objective observable identity with the quality ledger.",
            "The mesh-quality postprocess-row convention schema distinguishes how quality rows were selected, aggregated, and reduced from the metric-set id itself.",
            "The mesh-quality component-basis schema distinguishes which element family, quality component, coordinate basis, and normalization a row represents.",
        ],
    }


def cubit_mixed_solver_ready_package_gate(
    inventory: Mapping[str, object],
    transition_gate: Mapping[str, object],
    export_package_gate: Mapping[str, object],
    bnd_area_gate: Mapping[str, object],
    quality_gate: Mapping[str, object],
    *,
    expected_routing_hint: str = "cubit_hex_or_mixed_path",
    routing_policy_gate: Mapping[str, object] | None = None,
    interface_adjacency_gate: Mapping[str, object] | None = None,
    scheme_trace_gate: Mapping[str, object] | None = None,
    headless_batch_quality_gate: Mapping[str, object] | None = None,
    curvilinear_handoff_gate: Mapping[str, object] | None = None,
    solver_route_gate: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Check a Cubit mixed hex+pyramid+tet package before solver promotion.

    The individual gates catch important local errors.  This package-level
    replay keeps the mixed topology inventory, pyramid-transition metadata,
    export package identity, NGSolve BND area convention, and quality
    distribution together as one solver-ready handoff.
    """

    inv = dict(inventory)
    transition = dict(transition_gate)
    export_package = dict(export_package_gate)
    bnd = dict(bnd_area_gate)
    quality = dict(quality_gate)
    routing_policy = dict(routing_policy_gate or {})
    interface = dict(interface_adjacency_gate or {})
    scheme_trace = dict(scheme_trace_gate or {})
    headless_batch = dict(headless_batch_quality_gate or {})
    curvilinear_handoff = dict(curvilinear_handoff_gate or {})
    solver_route = dict(solver_route_gate or {})
    volume_counts_raw = inv.get("volume_kind_counts", {})
    if not isinstance(volume_counts_raw, Mapping):
        raise ValueError("inventory['volume_kind_counts'] must be a mapping")
    volume_counts = {str(key): int(value) for key, value in volume_counts_raw.items()}
    routing_hint = str(inv.get("routing_hint", "")).strip()
    expected_hint = str(expected_routing_hint).strip()
    transition_checks = transition.get("checks", {})
    export_checks = export_package.get("checks", {})
    bnd_checks = bnd.get("checks", {})
    quality_checks = quality.get("checks", {})
    if not isinstance(transition_checks, Mapping):
        transition_checks = {}
    if not isinstance(export_checks, Mapping):
        export_checks = {}
    if not isinstance(bnd_checks, Mapping):
        bnd_checks = {}
    if not isinstance(quality_checks, Mapping):
        quality_checks = {}
    routing_checks = routing_policy.get("checks", {})
    if not isinstance(routing_checks, Mapping):
        routing_checks = {}
    interface_checks = interface.get("checks", {})
    if not isinstance(interface_checks, Mapping):
        interface_checks = {}
    scheme_checks = scheme_trace.get("checks", {})
    if not isinstance(scheme_checks, Mapping):
        scheme_checks = {}
    headless_batch_checks = headless_batch.get("checks", {})
    if not isinstance(headless_batch_checks, Mapping):
        headless_batch_checks = {}
    curvilinear_checks = curvilinear_handoff.get("checks", {})
    if not isinstance(curvilinear_checks, Mapping):
        curvilinear_checks = {}
    solver_route_checks = solver_route.get("checks", {})
    if not isinstance(solver_route_checks, Mapping):
        solver_route_checks = {}

    checks = {
        "inventory_is_mixed_hex_pyramid_tet": (
            volume_counts.get("hex", 0) > 0
            and volume_counts.get("pyramid", 0) > 0
            and volume_counts.get("tet", 0) > 0
        ),
        "routing_hint_is_cubit_mixed": routing_hint == expected_hint,
        "not_tri_tet_only": inv.get("is_tri_tet_only") is not True,
        "transition_gate_ok": (
            transition.get("policy") == "cubit_mixed_transition_metadata_gate"
            and transition.get("status") == "ok"
            and transition_checks.get("transition_kinds_present") is True
        ),
        "export_package_gate_ok": (
            export_package.get("policy") == "cubit_export_package_identity_gate"
            and export_package.get("status") == "ok"
            and export_checks.get("inventory_routing_hint_matches_expected") is True
        ),
        "bnd_area_gate_ok": (
            bnd.get("policy") == "cubit_ngsolve_bnd_area_includes_material_interfaces_once"
            and bnd.get("status") == "ok"
            and bnd_checks.get("matches_external_plus_interface") is True
        ),
        "quality_distribution_gate_ok": (
            quality.get("policy") == "cubit_quality_distribution_gate"
            and quality.get("status") == "ok"
            and quality_checks.get("min_value_ok") is True
            and int(quality.get("count", 0) or 0) > 0
        ),
    }
    if routing_policy:
        checks["routing_policy_gate_ok"] = (
            routing_policy.get("policy") == "cubit_release_feature_routing_gate"
            and routing_policy.get("status") == "ok"
            and routing_checks.get("cubit_role_is_hex_or_mixed") is True
            and routing_checks.get("tet_only_owner_is_netgen") is True
            and routing_checks.get("headless_policy_recorded") is True
        )
    if interface:
        checks["interface_adjacency_gate_ok"] = (
            interface.get("policy") == "cubit_mixed_interface_adjacency_gate"
            and interface.get("status") == "ok"
            and interface_checks.get("required_roles_present") is True
            and interface_checks.get("role_surface_kinds_match") is True
            and interface_checks.get("role_volume_kind_pairs_match") is True
        )
    if scheme_trace:
        checks["scheme_trace_gate_ok"] = (
            scheme_trace.get("policy") == "cubit_meshing_scheme_trace_gate"
            and scheme_trace.get("status") == "ok"
            and scheme_checks.get("expected_command_digest_matches") is True
            and scheme_checks.get("expected_volume_schemes_match") is True
            and scheme_checks.get("required_command_fragments_present") is True
        )
    if headless_batch:
        export_package_id = str(
            export_package.get("expected_export_id")
            or (export_package.get("export_ids") or [""])[0]
        ).strip()
        export_package_geometry = str(
            export_package.get("expected_geometry_id")
            or (export_package.get("geometry_ids") or [""])[0]
        ).strip()
        checks["headless_batch_quality_gate_ok"] = (
            headless_batch.get("policy") == "cubit_headless_batch_quality_package_gate"
            and headless_batch.get("status") == "ok"
            and headless_batch_checks.get("process_mode_is_headless_batch") is True
            and headless_batch_checks.get("nographics_flag_present") is True
            and headless_batch_checks.get("batch_flag_present") is True
            and headless_batch_checks.get("gui_daemon_disabled") is True
            and headless_batch_checks.get("batch_script_recorded") is True
            and headless_batch_checks.get("process_exit_code_success_or_documented") is True
            and headless_batch.get("export_id") == export_package_id
            and headless_batch.get("geometry_id") == export_package_geometry
        )
    if curvilinear_handoff:
        checks["curvilinear_handoff_gate_ok"] = (
            curvilinear_handoff.get("policy") == "cubit_curvilinear_handoff_manifest_gate"
            and curvilinear_handoff.get("status") == "ok"
            and curvilinear_checks.get("source_mesh_is_imported") is True
            and curvilinear_checks.get("geometry_association_recorded") is True
            and curvilinear_checks.get("boundary_ids_preserved") is True
            and curvilinear_checks.get("projection_error_within_tolerance") is True
            and curvilinear_checks.get("curved_export_order_ok") is True
            and curvilinear_checks.get("routing_hint_is_cubit_hex_or_mixed") is True
            and curvilinear_checks.get("no_implicit_element_conversion") is True
            and curvilinear_checks.get("negative_jacobian_count_zero") is True
        )
    if solver_route:
        checks["solver_route_manifest_gate_ok"] = (
            solver_route.get("policy") == "cubit_mixed_solver_route_manifest_gate"
            and solver_route.get("status") == "ok"
            and solver_route_checks.get("volume_route_kinds_cover_inventory") is True
            and solver_route_checks.get("surface_route_kinds_cover_inventory") is True
            and solver_route_checks.get("pyramid_transition_role_recorded") is True
            and solver_route_checks.get("no_implicit_tetization_recorded") is True
            and solver_route_checks.get("tet_only_owner_is_netgen") is True
        )
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_mixed_solver_ready_package_gate",
        "status": "ok" if not issues else "needs_attention",
        "routing_hint": routing_hint,
        "expected_routing_hint": expected_hint,
        "volume_kind_counts": volume_counts,
        "transition_policy": transition.get("policy"),
        "export_package_policy": export_package.get("policy"),
        "bnd_area_policy": bnd.get("policy"),
        "quality_policy": quality.get("policy"),
        "routing_policy": routing_policy.get("policy"),
        "interface_adjacency_policy": interface.get("policy"),
        "scheme_trace_policy": scheme_trace.get("policy"),
        "headless_batch_quality_policy": headless_batch.get("policy"),
        "curvilinear_handoff_policy": curvilinear_handoff.get("policy"),
        "solver_route_policy": solver_route.get("policy"),
        "quality_min": quality.get("min"),
        "quality_count": quality.get("count"),
        "expected_bnd_area": bnd.get("expected_bnd_area"),
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use this after Cubit creates a mixed hex+pyramid+tet .vol package and before solver-ready promotion.",
            "Tet-only .vol files should stay on the Netgen/OCC education route, not this mixed Cubit package gate.",
            "The NGSolve BND check must include external plus material-interface area once.",
            "When release/route learning is available, bind it to the package so Cubit remains the hex-led mixed path.",
            "When interface-adjacency learning is available, bind it so hex-pyramid and pyramid-tet faces cannot be swapped silently.",
            "When scheme-trace learning is available, bind the Cubit volume schemes and export netgen command to the package.",
            "When headless-batch learning is available, bind process mode, no-GUI flags, batch script identity, and exit policy to the mixed solver-ready package.",
            "When curvilinear-handoff learning is available, bind CAD association, boundary preservation, projection tolerance, curved order, and zero negative-Jacobian evidence to the same mixed package.",
            "When solver-route learning is available, bind hex primary volume, pyramid transition, tet compatibility, surface trace routes, and no implicit tetization to the same mixed package.",
        ],
    }


def cubit_submodel_boundary_handoff_mesh_package_gate(
    inventory: Mapping[str, object],
    handoff: Mapping[str, object],
    *,
    expected_routing_hint: str = "cubit_hex_or_mixed_path",
    expected_boundary_name: str | None = None,
    max_boundary_transfer_error: float | None = None,
    expected_volume_kinds: Iterable[str] = ("hex",),
    transition_kinds_requiring_policy: Iterable[str] = ("pyramid",),
) -> dict[str, object]:
    """Bind Cubit ``.vol`` inventory to a submodel boundary handoff contract.

    Cubit is the lab's hex-led and mixed-mesh route.  When a mesh package is
    used for a local/zoomed submodel, the archived ``.vol`` inventory and the
    parent-to-local boundary transfer metadata must travel together.
    """

    inv = dict(inventory)
    h = dict(handoff)
    routing_hint = str(inv.get("routing_hint", "")).strip()
    expected_hint = str(expected_routing_hint or "").strip()
    source = str(inv.get("source", "")).strip()
    volume_elements = int(inv.get("volume_elements", 0) or 0)
    surface_elements = int(inv.get("surface_elements", 0) or 0)
    is_tri_tet_only = bool(inv.get("is_tri_tet_only", False))
    volume_counts_raw = inv.get("volume_kind_counts", {})
    if not isinstance(volume_counts_raw, Mapping):
        raise ValueError("inventory['volume_kind_counts'] must be a mapping when provided")
    volume_counts = {
        str(key).strip().lower(): int(value)
        for key, value in volume_counts_raw.items()
        if str(key).strip()
    }
    expected_kinds = [
        str(kind).strip().lower()
        for kind in expected_volume_kinds
        if str(kind).strip()
    ]
    transition_kinds = [
        str(kind).strip().lower()
        for kind in transition_kinds_requiring_policy
        if str(kind).strip()
    ]
    present_volume_kinds = sorted(
        kind for kind, count in volume_counts.items() if int(count) > 0
    )
    present_transition_kinds = sorted(
        kind for kind in transition_kinds if volume_counts.get(kind, 0) > 0
    )
    boundary_names_raw = inv.get("boundary_names", {})
    if not isinstance(boundary_names_raw, Mapping):
        boundary_names_raw = {}
    boundary_names = {
        str(key): str(value).strip()
        for key, value in boundary_names_raw.items()
        if str(value).strip()
    }
    boundary_name_values = set(boundary_names.values())

    def text_field(name: str) -> str:
        return str(h.get(name, "") or "").strip()

    parent_model_id = text_field("parent_model_id")
    submodel_region_id = text_field("submodel_region_id")
    zoom_boundary_id = text_field("zoom_boundary_id")
    boundary_transfer_quantity = text_field("boundary_transfer_quantity")
    boundary_transfer_error_unit = text_field("boundary_transfer_error_unit")
    local_refinement_rule = text_field("local_refinement_rule")
    target_observable_id = text_field("target_observable_id")
    parent_mesh_id = text_field("parent_mesh_id")
    local_mesh_id = text_field("local_mesh_id")
    boundary_trace_id = text_field("boundary_trace_id")
    transition_policy = (
        text_field("transition_policy")
        or text_field("pyramid_transition_policy")
        or str(inv.get("transition_policy", "") or "").strip()
    )
    expected_boundary = (
        str(expected_boundary_name).strip()
        if expected_boundary_name is not None
        else zoom_boundary_id
    )

    def finite_float_or_none(value):
        if value is None:
            return None
        try:
            converted = float(value)
        except (TypeError, ValueError):
            return None
        return converted if isfinite(converted) else None

    transfer_error = finite_float_or_none(h.get("boundary_transfer_error_estimate"))
    max_error = finite_float_or_none(max_boundary_transfer_error)

    checks = {
        "inventory_source_recorded": bool(source),
        "inventory_source_is_vol": source.lower().endswith(".vol"),
        "inventory_volume_elements_positive": volume_elements > 0,
        "inventory_surface_elements_positive": surface_elements > 0,
        "volume_kind_counts_recorded": bool(volume_counts),
        "expected_volume_kinds_present": all(
            volume_counts.get(kind, 0) > 0 for kind in expected_kinds
        ),
        "hex_family_present_for_cubit_submodel": (
            expected_hint != "cubit_hex_or_mixed_path"
            or volume_counts.get("hex", 0) > 0
        ),
        "transition_policy_recorded_when_present": (
            not present_transition_kinds or bool(transition_policy)
        ),
        "routing_hint_matches_expected": routing_hint == expected_hint,
        "not_tri_tet_only_for_cubit_submodel": not is_tri_tet_only,
        "boundary_names_recorded": bool(boundary_names),
        "parent_model_id_recorded": bool(parent_model_id),
        "submodel_region_id_recorded": bool(submodel_region_id),
        "zoom_boundary_id_recorded": bool(zoom_boundary_id),
        "boundary_transfer_quantity_recorded": bool(boundary_transfer_quantity),
        "boundary_transfer_error_estimate_recorded": transfer_error is not None,
        "boundary_transfer_error_nonnegative": (
            transfer_error is not None and transfer_error >= 0.0
        ),
        "boundary_transfer_error_unit_recorded": bool(boundary_transfer_error_unit),
        "local_refinement_rule_recorded": bool(local_refinement_rule),
        "target_observable_id_recorded": bool(target_observable_id),
        "boundary_handoff_not_value_only": all(
            bool(item)
            for item in (
                parent_model_id,
                submodel_region_id,
                zoom_boundary_id,
                boundary_transfer_quantity,
                local_refinement_rule,
                target_observable_id,
            )
        ),
    }
    if expected_boundary:
        checks["expected_boundary_name_recorded"] = bool(expected_boundary)
        checks["zoom_boundary_present_in_vol_inventory"] = expected_boundary in boundary_name_values
    if max_boundary_transfer_error is not None:
        checks["max_boundary_transfer_error_recorded"] = max_error is not None and max_error >= 0.0
        checks["boundary_transfer_error_within_limit"] = (
            transfer_error is not None and max_error is not None and transfer_error <= max_error
        )
    if parent_mesh_id or local_mesh_id:
        checks["parent_mesh_id_recorded"] = bool(parent_mesh_id)
        checks["local_mesh_id_recorded"] = bool(local_mesh_id)
        checks["parent_local_mesh_identity_separated"] = (
            bool(parent_mesh_id) and bool(local_mesh_id) and parent_mesh_id != local_mesh_id
        )
    if boundary_trace_id:
        checks["boundary_trace_id_recorded"] = True

    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_submodel_boundary_handoff_mesh_package_gate",
        "status": "ok" if not issues else "needs_attention",
        "source": source,
        "routing_hint": routing_hint,
        "expected_routing_hint": expected_hint,
        "volume_elements": volume_elements,
        "surface_elements": surface_elements,
        "volume_kind_counts": volume_counts,
        "present_volume_kinds": present_volume_kinds,
        "expected_volume_kinds": expected_kinds,
        "transition_kinds_requiring_policy": transition_kinds,
        "transition_policy": transition_policy,
        "is_tri_tet_only": is_tri_tet_only,
        "boundary_names": boundary_names,
        "parent_model_id": parent_model_id,
        "submodel_region_id": submodel_region_id,
        "zoom_boundary_id": zoom_boundary_id,
        "expected_boundary_name": expected_boundary,
        "boundary_transfer_quantity": boundary_transfer_quantity,
        "boundary_transfer_error_estimate": transfer_error,
        "boundary_transfer_error_unit": boundary_transfer_error_unit,
        "max_boundary_transfer_error": max_error,
        "local_refinement_rule": local_refinement_rule,
        "target_observable_id": target_observable_id,
        "parent_mesh_id": parent_mesh_id,
        "local_mesh_id": local_mesh_id,
        "boundary_trace_id": boundary_trace_id,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Use this for Cubit hex-led or mixed .vol packages that feed a local/zoomed submodel.",
            "The mesh inventory and parent-to-local boundary handoff must be verified together.",
            "Tet-only .vol files remain on the Netgen/OCC education route, not the Cubit submodel route.",
            "For Cubit submodels, keep the hex/mixed volume family inventory and pyramid-transition policy with the boundary handoff.",
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


def cubit_webcut_conformal_hex_gate(summary: Mapping[str, object], *, max_webcut_volume_relative_drift: float = 1.0e-5, max_partition_relative_spread: float = 1.0e-12, min_scaled_jacobian: float = 0.2) -> dict[str, object]:
    """Gate a webcut decomposition by conservation and conformal interfaces."""
    counts = summary.get("element_counts") or {}
    interfaces = summary.get("interfaces") or []
    drift = float(summary.get("webcut_volume_relative_drift", math.inf))
    spread = float(summary.get("quarter_volume_relative_spread", math.inf))
    quality = summary.get("quality") or {}
    sj = float((quality.get("scaled_jacobian") or {}).get("min", -math.inf))
    shape = float((quality.get("shape") or {}).get("min", -math.inf))
    checks = {
        "four_partition_volumes": len(summary.get("volume_ids") or []) == 4,
        "all_hex_mesh": int(counts.get("hex", 0)) > 0 and all(int(counts.get(name, 0)) == 0 for name in ("tet", "wedge", "pyramid")),
        "partition_volume_spread_bounded": spread <= float(max_partition_relative_spread),
        "webcut_volume_drift_bounded": drift <= float(max_webcut_volume_relative_drift),
        "four_internal_interfaces": len(interfaces) == 4,
        "interfaces_shared_and_meshed": all(len(row.get("adjacent_volumes") or []) == 2 and int(row.get("face_count", 0)) > 0 for row in interfaces),
        "interface_area_positive": all(float(row.get("area", 0.0)) > 0.0 for row in interfaces),
        "boundary_face_block_occupied": int(summary.get("boundary_block_face_count", 0)) > 0,
        "scaled_jacobian_ok": sj >= float(min_scaled_jacobian),
        "shape_positive": shape > 0.0,
    }
    return {"policy":"cubit_webcut_conformal_hex_gate_v1","status":"ok" if all(checks.values()) else "needs_attention","checks":checks,"issues":[name for name,ok in checks.items() if not ok],"metrics":{"hex_count":int(counts.get("hex",0)),"webcut_volume_relative_drift":drift,"partition_relative_spread":spread,"interface_count":len(interfaces),"interface_face_count":sum(int(row.get("face_count",0)) for row in interfaces),"minimum_scaled_jacobian":sj,"minimum_shape":shape}}


def cubit_webcut_journal_execution_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Gate source-journal ordering and the shared-owner headless process lane."""
    commands = [str(value).lower() for value in summary.get("commands", [])]
    def first(token):
        return next((i for i,row in enumerate(commands) if token in row), -1)
    process = _headless_process_evidence(summary)
    checks = {
        "source_journal_digest_recorded": len(str(summary.get("source_sha256", ""))) == 64,
        "source_native_journal": summary.get("source_kind") == "source_native_local_journal",
        "headless_python_api": summary.get("execution_mode") == "python_api_headless" and {"-nographics","-batch"}.issubset(set(summary.get("headless_flags") or [])),
        "persistent_gui_disabled": summary.get("gui_daemon_enabled") is False,
        "two_orthogonal_webcuts": first("xplane") >= 0 and first("yplane") >= 0,
        "imprint_merge_before_mesh": 0 <= first("imprint all") < first("merge all") < first("mesh volume"),
        "boundary_selection_is_geometric": any("surface with area" in row and "z_min" in row and "z_max" in row for row in commands),
        "process_exit_acceptable": process["process_exit_acceptable"],
        "owned_processes_closed": int(summary.get("owned_processes_remaining", 0)) == 0,
        "public_geometry_gate_passed": summary.get("public_gate_status") == "ok",
    }
    return {"policy":"cubit_webcut_journal_execution_gate_v1","status":"ok" if all(checks.values()) else "needs_attention","checks":checks,"issues":[name for name,ok in checks.items() if not ok],"process":process}


def cubit_helical_partition_mesh_gate(
    summary: Mapping[str, object],
    *,
    max_webcut_volume_relative_drift: float = 1.0e-5,
    max_analytic_volume_relative_error: float = 5.0e-4,
    min_scaled_jacobian: float = 0.2,
    min_partitioned_volume_count: int = 1000,
) -> dict[str, object]:
    """Gate a many-body helical webcut mesh and its exported volume inventory."""

    counts = summary.get("element_counts") or {}
    quality = summary.get("quality") or {}
    hex_quality = (quality.get("hex") or {}).get("scaled_jacobian") or {}
    inventory = summary.get("export_inventory") or {}
    inventory_counts = inventory.get("volume_kind_counts") or {}
    hex_count = int(counts.get("hex", 0))
    volume_element_count = int(inventory.get("volume_elements", 0))
    drift = float(summary.get("webcut_volume_relative_drift", math.inf))
    analytic_error = float(summary.get("analytic_volume_relative_error", math.inf))
    minimum_scaled_jacobian = float(hex_quality.get("min", -math.inf))
    checks = {
        "many_body_partition_present": int(summary.get("volume_count", 0)) >= int(min_partitioned_volume_count),
        "hex_led_mesh_present": hex_count > 0,
        "webcut_volume_drift_bounded": drift <= float(max_webcut_volume_relative_drift),
        "analytic_enclosure_volume_bounded": analytic_error <= float(max_analytic_volume_relative_error),
        "scaled_jacobian_acceptable": minimum_scaled_jacobian >= float(min_scaled_jacobian),
        "shared_interfaces_present": int(summary.get("shared_surface_count", 0)) > 0,
        "shared_interfaces_meshed": int(summary.get("shared_meshed_surface_count", 0)) > 0,
        "export_volume_elements_present": volume_element_count > 0,
        "export_volume_count_matches_cubit": volume_element_count == sum(int(value) for value in counts.values()),
        "export_hex_count_matches_cubit": int(inventory_counts.get("hex", 0)) == hex_count,
        "export_routes_to_cubit_hex_or_mixed": inventory.get("routing_hint") == "cubit_hex_or_mixed_path",
        "export_points_present": int(inventory.get("points", 0)) > 0,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_helical_partition_mesh_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "volume_count": int(summary.get("volume_count", 0)),
            "hex_count": hex_count,
            "webcut_volume_relative_drift": drift,
            "analytic_volume_relative_error": analytic_error,
            "minimum_scaled_jacobian": minimum_scaled_jacobian,
            "shared_surface_count": int(summary.get("shared_surface_count", 0)),
            "shared_meshed_surface_count": int(summary.get("shared_meshed_surface_count", 0)),
            "export_volume_elements": volume_element_count,
            "export_hex_count": int(inventory_counts.get("hex", 0)),
        },
        "notes": [
            "A positive Cubit hex count is not solver-ready evidence by itself.",
            "Bind Boolean/webcut conservation, non-inverted quality, shared interfaces, and parsed .vol inventory in one gate.",
            "A surface-only .vol export must be rejected even when Cubit still holds volume hexes in memory.",
        ],
    }


def cubit_source_journal_replay_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Gate safe synchronous replay of a heavy source-native Cubit journal."""

    process = _headless_process_evidence(summary)
    operations = summary.get("operations") or {}
    selection_counts = [int(value) for value in operations.get("webcut_intersection_selection_counts", [])]
    timing = summary.get("timing_breakdown_s") or {}
    async_probe = summary.get("playback_async_probe") or {}
    expected_status = str(summary.get("expected_public_gate_status", ""))
    actual_status = str(summary.get("public_gate_status", ""))
    checks = {
        "source_digest_recorded": len(str(summary.get("source_sha256", ""))) == 64,
        "source_native_seed": str(summary.get("source_kind", "")).startswith("source_native_local_journal"),
        "synchronous_python_replay": summary.get("execution_mode") == "python_api_headless_synchronous_commands",
        "headless_flags_present": {"-nographics", "-batch"}.issubset(set(summary.get("headless_flags") or [])),
        "persistent_gui_disabled": summary.get("gui_daemon_enabled") is False,
        "early_playback_artifact_rejected": async_probe.get("unsafe_zero_entity_artifact_rejected") is True,
        "playback_replaced_by_synchronous_commands": "synchronous" in str(async_probe.get("replacement", "")).lower(),
        "all_cut_planes_use_nonempty_selection": len(selection_counts) == 15 and min(selection_counts, default=0) > 0,
        "four_stage_timing_recorded": len(timing) == 4 and all(float(value) >= 0.0 for value in timing.values()),
        "process_exit_acceptable": process["process_exit_acceptable"],
        "owned_processes_closed": int(summary.get("owned_processes_remaining", 0)) == 0,
        "expected_mesh_disposition_recorded": expected_status in {"ok", "needs_attention"},
        "expected_mesh_disposition_matched": expected_status == actual_status,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_source_journal_replay_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "process": process,
        "selection_count_range": [min(selection_counts, default=0), max(selection_counts, default=0)],
        "public_gate_status": actual_status,
        "notes": [
            "Do not query or write result JSON immediately after queuing a playback command; use synchronous command execution or a completion sentinel.",
            "Treat an expected mesh rejection as a valid source-workflow lesson only when the public quality/export gate independently explains it.",
            "Select only volumes intersecting each cut plane and record geometry, cut/merge, mesh, and export timing separately.",
        ],
    }
