"""Jointed assembly STEP closure and source-replay validation gates."""

from __future__ import annotations

import math
from typing import Mapping


def _relative_error(measured: float, reference: float) -> float:
    if not math.isfinite(reference) or reference <= 0.0 or not math.isfinite(measured):
        return math.inf
    return abs(measured - reference) / reference


def jointed_assembly_step_closure_gate(
    summary: Mapping[str, object],
    *,
    tessellated_volume_rtol: float = 2.0e-4,
    self_roundtrip_volume_rtol: float = 1.0e-6,
    external_volume_rtol: float = 2.0e-5,
) -> dict[str, object]:
    """Diagnose component-level solid loss in a jointed assembly STEP handoff.

    ``status=ok`` means the diagnosis is supported by a portable component,
    a rejected component, and independent volume evidence.  It does not mean
    the assembly is solver-ready; that is reported separately.
    """

    tolerances = (
        float(tessellated_volume_rtol),
        float(self_roundtrip_volume_rtol),
        float(external_volume_rtol),
    )
    if any(value < 0.0 or not math.isfinite(value) for value in tolerances):
        raise ValueError("volume tolerances must be finite and nonnegative")

    raw_components = summary.get("components") or []
    if not isinstance(raw_components, list) or len(raw_components) < 2:
        raise ValueError("components must contain at least two rows")
    external_rows = summary.get("external_rows") or []
    if not isinstance(external_rows, list):
        raise ValueError("external_rows must be a list")
    external_by_name = {str(row.get("name", "")): row for row in external_rows}

    components = []
    portable = []
    rejected = []
    for raw in raw_components:
        name = str(raw.get("name", "")).strip()
        native_volume = float(raw.get("native_volume_mm3", math.nan))
        tessellated_volume = float(raw.get("tessellated_volume_mm3", math.nan))
        self_roundtrip = raw.get("self_roundtrip") or {}
        self_volume = float(self_roundtrip.get("total_volume_mm3", math.nan))
        self_solid_count = int(self_roundtrip.get("solid_count", -1))
        external = external_by_name.get(f"{name}.step", {})
        external_volume = float(external.get("total_volume_mm3", math.nan))
        external_volume_count = int(external.get("volume_count", -1))
        tess_error = _relative_error(tessellated_volume, native_volume)
        self_error = _relative_error(self_volume, native_volume)
        external_error = _relative_error(external_volume, native_volume)
        native_supported = (
            bool(name)
            and raw.get("native_valid") is True
            and int(raw.get("native_solid_count", 0)) == 1
            and native_volume > 0.0
            and tess_error <= tolerances[0]
        )
        is_portable = (
            native_supported
            and self_solid_count == 1
            and self_error <= tolerances[1]
            and external_volume_count >= 1
            and external_volume > 0.0
            and external_error <= tolerances[2]
        )
        is_closure_loss = (
            native_supported
            and self_solid_count == 0
            and self_volume == 0.0
            and external_volume_count >= 1
            and external_volume == 0.0
        )
        expected = str(raw.get("expected_disposition", ""))
        disposition = (
            "portable_control"
            if is_portable
            else "reject_solid_closure_loss"
            if is_closure_loss
            else "unresolved"
        )
        row = {
            "name": name,
            "expected_disposition": expected,
            "observed_disposition": disposition,
            "native_volume_mm3": native_volume,
            "tessellated_volume_relative_error": tess_error,
            "self_roundtrip_volume_relative_error": self_error,
            "external_volume_relative_error": external_error,
            "self_roundtrip_solid_count": self_solid_count,
            "external_volume_count": external_volume_count,
            "external_volume_mm3": external_volume,
            "disposition_matches": disposition == expected,
        }
        components.append(row)
        if is_portable:
            portable.append(row)
        if is_closure_loss:
            rejected.append(row)

    assembly = summary.get("assembly") or {}
    assembly_external = external_by_name.get(str(assembly.get("step_name", "")), {})
    native_total = float(assembly.get("native_total_volume_mm3", math.nan))
    self_total = float((assembly.get("self_roundtrip") or {}).get("total_volume_mm3", math.nan))
    external_total = float(assembly_external.get("total_volume_mm3", math.nan))
    rejected_native_total = sum(float(row["native_volume_mm3"]) for row in rejected)
    self_lost = native_total - self_total
    external_lost = native_total - external_total
    names = [str(row["name"]) for row in components]
    checks = {
        "component_names_are_nonempty_and_unique": all(names) and len(set(names)) == len(names),
        "external_rows_cover_every_component_and_assembly": all(
            f"{name}.step" in external_by_name for name in names
        )
        and str(assembly.get("step_name", "")) in external_by_name,
        "native_brep_supported_by_tessellated_volume": all(
            float(row["tessellated_volume_relative_error"]) <= tolerances[0]
            for row in components
        ),
        "portable_component_control_present": bool(portable),
        "solid_closure_loss_component_present": bool(rejected),
        "component_dispositions_match_expectations": all(
            bool(row["disposition_matches"]) for row in components
        ),
        "assembly_self_loss_matches_rejected_components": rejected_native_total > 0.0
        and _relative_error(self_lost, rejected_native_total) <= tolerances[1],
        "assembly_external_loss_matches_rejected_components": rejected_native_total > 0.0
        and _relative_error(external_lost, rejected_native_total) <= tolerances[2],
    }
    issues = [name for name, ok in checks.items() if not ok]
    diagnosis = "component_solid_closure_loss" if not issues else "incomplete_component_evidence"
    return {
        "policy": "build123d_jointed_assembly_step_closure_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "diagnosis": diagnosis,
        "solver_ready": False,
        "checks": checks,
        "issues": issues,
        "components": components,
        "assembly": {
            "native_total_volume_mm3": native_total,
            "self_roundtrip_total_volume_mm3": self_total,
            "external_total_volume_mm3": external_total,
            "rejected_component_native_volume_mm3": rejected_native_total,
            "self_roundtrip_lost_volume_mm3": self_lost,
            "external_lost_volume_mm3": external_lost,
        },
        "notes": [
            "Assembly total volume alone cannot identify which source component lost solid closure.",
            "Require a portable component control and an independently tessellated source volume before blaming STEP translation.",
            "Do not promote an assembly when any positive-volume source component returns as a zero-volume external entity.",
        ],
    }


def jointed_assembly_source_replay_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Gate immutable source replay, joint identity, and headless CAD diagnosis."""

    components = summary.get("components") or []
    connections = summary.get("joint_connections") or []
    execution = summary.get("external_execution") or {}
    timing = summary.get("timing_breakdown_s") or {}
    joint_names = {
        str(name)
        for row in components
        for name in (row.get("joint_names") or [])
        if str(name)
    }
    connection_endpoints = {
        str(connection.get(side, ""))
        for connection in connections
        for side in ("from", "to")
    }
    checks = {
        "upstream_source_identity_recorded": str(summary.get("source_kind", "")).startswith(
            "upstream_source_native_example"
        )
        and len(str(summary.get("source_sha256", ""))) == 64
        and "/v0.10.0/" in str(summary.get("source_url", "")),
        "source_preserved_except_display_stub": summary.get("source_preserved") is True
        and summary.get("display_stubbed_only") is True,
        "component_inventory_and_joint_names_recorded": len(components) >= 2
        and all(str(row.get("name", "")).strip() and row.get("joint_names") for row in components),
        "joint_connection_endpoints_resolve": bool(connections)
        and bool(connection_endpoints)
        and connection_endpoints.issubset(joint_names),
        "headless_external_cad_replay": execution.get("mode")
        == "python_api_headless_synchronous_commands"
        and {"-nographics", "-batch"}.issubset(set(execution.get("headless_flags") or []))
        and execution.get("gui_daemon_enabled") is False,
        "fresh_result_and_owned_process_cleanup": execution.get("result_artifact_fresh") is True
        and int(execution.get("owned_processes_remaining", -1)) == 0,
        "component_closure_diagnosis_verified": summary.get("diagnosis_gate_status") == "ok"
        and summary.get("diagnosis") == "component_solid_closure_loss"
        and summary.get("solver_ready") is False,
        "four_dominant_timing_stages_recorded": len(timing) == 4
        and all(float(value) >= 0.0 for value in timing.values()),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_jointed_assembly_source_replay_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "joint_names": sorted(joint_names),
        "connection_endpoints": sorted(connection_endpoints),
        "notes": [
            "Build123d joint semantics are source metadata; persist their names and connection graph beside a neutral CAD export.",
            "A display stub is acceptable for headless replay, but geometry and joint-building statements must remain unchanged.",
            "Record component-level external entities so a zero-volume body cannot hide behind a surviving assembly member.",
        ],
    }
