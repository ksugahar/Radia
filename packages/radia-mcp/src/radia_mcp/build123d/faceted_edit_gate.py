"""Gates for an upstream faceted-solid edit and external CAD replay."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def _rows(value: object, name: str) -> list[Mapping[str, object]]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ValueError(f"{name} must be a sequence")
    result = []
    for index, row in enumerate(value):
        result.append(_mapping(row, f"{name}[{index}]"))
    return result


def _relative(left: object, right: object) -> float:
    a = float(left)
    b = float(right)
    return abs(a - b) / max(abs(a), abs(b), 1.0)


def build123d_faceted_edit_portability_gate(
    summary: Mapping[str, object],
    *,
    minimum_solver_ready_scaled_jacobian: float = 0.05,
) -> dict[str, object]:
    """Separate CAD portability from downstream tetrahedral mesh readiness."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    if minimum_solver_ready_scaled_jacobian <= 0.0:
        raise ValueError("minimum_solver_ready_scaled_jacobian must be positive")
    native_runs = _rows(summary.get("native_runs"), "native_runs")
    external_rows = _rows(summary.get("external_imports"), "external_imports")
    mesh = _mapping(summary.get("mesh_evidence"), "mesh_evidence")

    native_metrics = [_mapping(row.get("native"), "native") for row in native_runs]
    source_metrics = [
        _mapping(row.get("imported_stl_before_edit"), "imported_stl_before_edit")
        for row in native_runs
    ]
    roundtrips = [
        _mapping(row.get("self_roundtrip"), "self_roundtrip") for row in native_runs
    ]
    native_topology = [
        tuple(int(row.get(name, -1)) for name in ("solid_count", "face_count", "edge_count", "vertex_count"))
        for row in native_metrics
    ]
    external_topology = [
        tuple(int(row.get(name, -1)) for name in ("surface_count", "curve_count", "vertex_count"))
        for row in external_rows
    ]
    external_volume_errors = [
        float(row.get("native_volume_relative_error", math.inf)) for row in external_rows
    ]
    external_area_errors = [
        float(row.get("native_area_relative_error", math.inf)) for row in external_rows
    ]
    self_volume_errors = [
        _relative(native["volume_mm3"], replay["volume_mm3"])
        for native, replay in zip(native_metrics, roundtrips)
    ]
    self_area_errors = [
        _relative(native["area_mm2"], replay["area_mm2"])
        for native, replay in zip(native_metrics, roundtrips)
    ]
    min_scaled = float(mesh.get("min_scaled_jacobian", math.nan))
    observed_solver_ready = math.isfinite(min_scaled) and min_scaled >= float(
        minimum_solver_ready_scaled_jacobian
    )
    claimed_solver_ready = summary.get("solver_ready_claimed") is True
    expected_classification = "solver_ready" if observed_solver_ready else "diagnostic_only"

    checks = {
        "two_native_runs_recorded": len(native_runs) == 2,
        "native_runs_are_valid_single_solids": all(
            row.get("is_valid") is True and int(row.get("solid_count", 0)) == 1
            for row in native_metrics
        ),
        "faceted_source_was_materially_edited": all(
            _relative(native["volume_mm3"], source["volume_mm3"]) > 1.0e-4
            and int(native.get("face_count", 0)) != int(source.get("face_count", 0))
            for native, source in zip(native_metrics, source_metrics)
        ),
        "native_volume_area_and_topology_replay": len(set(native_topology)) == 1
        and _relative(native_metrics[0]["volume_mm3"], native_metrics[1]["volume_mm3"])
        <= 1.0e-12
        and _relative(native_metrics[0]["area_mm2"], native_metrics[1]["area_mm2"])
        <= 1.0e-12,
        "self_step_roundtrip_closes": max(self_volume_errors, default=math.inf) <= 1.0e-10
        and max(self_area_errors, default=math.inf) <= 1.0e-10
        and all(row.get("is_valid") is True for row in roundtrips),
        "four_external_heal_noheal_rows": len(external_rows) == 4
        and {(int(row.get("run_index", 0)), str(row.get("import_mode", ""))) for row in external_rows}
        == {(1, "heal"), (1, "noheal"), (2, "heal"), (2, "noheal")},
        "external_imports_are_single_positive_volumes": all(
            int(row.get("volume_count", 0)) == 1
            and int(row.get("positive_volume_count", 0)) == 1
            for row in external_rows
        ),
        "heal_noheal_topology_is_invariant": len(set(external_topology)) == 1,
        "external_mass_properties_close": max(external_volume_errors, default=math.inf)
        <= 1.0e-4
        and max(external_area_errors, default=math.inf) <= 1.0e-4,
        "tet_mesh_volume_is_diagnostic_consistent": list(
            mesh.get("connectivity_orders", [])
        )
        == [4]
        and float(mesh.get("cad_volume_relative_error", math.inf)) <= 5.0e-3,
        "gmsh_launch_companions_recorded": all(
            _mapping(mesh.get("gmsh_companions"), "gmsh_companions").get(name) is True
            for name in ("geo", "geo_opt", "msh_opt")
        ),
        "solver_ready_claim_matches_quality": claimed_solver_ready
        == observed_solver_ready,
        "mesh_quality_classification_matches": summary.get("mesh_quality_classification")
        == expected_classification,
    }
    return {
        "policy": "build123d_faceted_edit_portability_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "solver_ready": observed_solver_ready,
        "mesh_quality_classification": expected_classification,
        "metrics": {
            "native_volume_mm3": float(native_metrics[0]["volume_mm3"]),
            "native_area_mm2": float(native_metrics[0]["area_mm2"]),
            "native_face_count": int(native_metrics[0]["face_count"]),
            "source_face_count": int(source_metrics[0]["face_count"]),
            "maximum_self_roundtrip_volume_relative_error": max(self_volume_errors),
            "maximum_external_volume_relative_error": max(external_volume_errors),
            "maximum_external_area_relative_error": max(external_area_errors),
            "tet_count": int(mesh.get("tet_count", 0)),
            "tet_integrated_volume_relative_error": float(
                mesh.get("cad_volume_relative_error", math.inf)
            ),
            "minimum_tet_scaled_jacobian": min_scaled,
            "minimum_solver_ready_scaled_jacobian": float(
                minimum_solver_ready_scaled_jacobian
            ),
        },
        "lesson": (
            "A valid STEP with excellent cross-kernel volume closure is not automatically "
            "solver-ready. Keep CAD portability and mesh readiness as separate gates: "
            "faceted source edits can preserve mass properties while producing sliver "
            "tetrahedra that require defeaturing, local sizing, or remeshing."
        ),
    }


def build123d_faceted_source_replay_gate(
    summary: Mapping[str, object],
) -> dict[str, object]:
    """Gate upstream source, dependent STL, exact execution, and external replay."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    process = _mapping(summary.get("external_process"), "external_process")
    public = _mapping(summary.get("public_gate"), "public_gate")
    timing = _mapping(summary.get("timing_breakdown_s"), "timing_breakdown_s")
    flags = {str(value).lower() for value in summary.get("headless_flags", [])}
    checks = {
        "tagged_upstream_source_identity_recorded": summary.get("source_kind")
        == "upstream_native_example_tag_v0.10.0"
        and summary.get("upstream_tag") == "v0.10.0"
        and len(str(summary.get("upstream_commit", ""))) == 40
        and len(str(summary.get("source_sha256", ""))) == 64,
        "dependent_stl_identity_recorded": summary.get("dependent_asset_bound") is True
        and len(str(summary.get("stl_sha256", ""))) == 64,
        "source_and_dependency_preserved": summary.get("source_preserved") is True
        and summary.get("stl_preserved") is True,
        "exact_source_viewer_stub_only": summary.get("execution_mode")
        == "exact_source_with_viewer_stub_only"
        and summary.get("viewer_stub_only") is True,
        "installed_version_and_two_runs_recorded": summary.get("build123d_version")
        == "0.10.0"
        and int(summary.get("source_run_count", 0)) == 2,
        "faceted_import_edit_semantics_recorded": summary.get(
            "mesh_imported_solid_edited"
        )
        is True
        and int(summary.get("face_count_delta", 0)) != 0,
        "external_cad_headless_runtime_preserved": summary.get("external_cad_version")
        == "2025.12"
        and {"-nographics", "-batch"}.issubset(flags)
        and summary.get("gui_daemon_enabled") is False
        and process.get("acceptable") is True
        and process.get("result_artifact_fresh") is True
        and int(process.get("owned_processes_remaining", -1)) == 0
        and not process.get("unexpected_error_lines"),
        "public_portability_quality_gate_passed": public.get("policy")
        == "build123d_faceted_edit_portability_gate_v1"
        and public.get("status") == "ok"
        and public.get("solver_ready") is False,
        "public_negative_control_rejected": summary.get("public_negative_status")
        == "needs_attention",
        "four_stage_timing_recorded": len(timing) == 4
        and all(float(value) >= 0.0 for value in timing.values()),
    }
    return {
        "policy": "build123d_faceted_source_replay_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "lesson": (
            "Bind every upstream dependent asset, not only the Python source. Execute the "
            "tagged example with a viewer-only stub, preserve both digests, and require a "
            "headless external-CAD replay before promoting a mesh-import/edit example."
        ),
    }
