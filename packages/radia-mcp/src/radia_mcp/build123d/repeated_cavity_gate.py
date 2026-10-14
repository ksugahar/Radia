"""Dual-API and source-replay gates for repeated-feature cavity solids."""

from __future__ import annotations

import math
from collections.abc import Mapping


def _relative(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def _metrics(row: Mapping[str, object]) -> dict[str, object]:
    return {
        "volume": float(row.get("volume", math.nan)),
        "area": float(row.get("area", math.nan)),
        "bbox_extent": [float(value) for value in row.get("bbox_extent") or []],
        "solid_count": int(row.get("solid_count", 0)),
        "shell_count": int(row.get("shell_count", 0)),
        "face_count": int(row.get("face_count", 0)),
        "edge_count": int(row.get("edge_count", 0)),
        "vertex_count": int(row.get("vertex_count", 0)),
    }


def build123d_repeated_cavity_dual_api_gate(
    summary: Mapping[str, object],
    *,
    native_relative_tolerance: float = 1.0e-12,
    external_relative_tolerance: float = 2.0e-6,
) -> dict[str, object]:
    """Gate dual authoring APIs and four external STEP import paths."""

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    native_tol = float(native_relative_tolerance)
    external_tol = float(external_relative_tolerance)
    if not all(math.isfinite(value) and value >= 0.0 for value in (native_tol, external_tol)):
        raise ValueError("tolerances must be finite and nonnegative")

    native = summary.get("native")
    topology = summary.get("topology_contract")
    features = summary.get("feature_contract")
    external = summary.get("external_imports")
    if not isinstance(native, Mapping) or set(native) != {"builder", "algebra"}:
        raise ValueError("native must contain builder and algebra records")
    if not isinstance(topology, Mapping) or not isinstance(features, Mapping):
        raise ValueError("topology_contract and feature_contract must be mappings")
    if not isinstance(external, list) or len(external) != 4:
        raise ValueError("external_imports must contain four records")

    expected_topology = {
        key: int(topology.get(key, 0))
        for key in ("solid_count", "shell_count", "face_count", "edge_count", "vertex_count")
    }
    expected_bbox = [float(value) for value in topology.get("bbox_extent") or []]
    if len(expected_bbox) != 3 or not all(math.isfinite(value) and value > 0.0 for value in expected_bbox):
        raise ValueError("topology_contract.bbox_extent must contain three positive values")
    parsed_native = {name: _metrics(row) for name, row in native.items() if isinstance(row, Mapping)}
    if set(parsed_native) != {"builder", "algebra"}:
        raise ValueError("native records must be mappings")
    builder = parsed_native["builder"]
    algebra = parsed_native["algebra"]
    official = float(summary.get("official_expected_volume", math.nan))

    native_topology_ok = all(
        all(int(row[key]) == expected_topology[key] for key in expected_topology)
        for row in parsed_native.values()
    )
    native_bbox_ok = all(row["bbox_extent"] == expected_bbox for row in parsed_native.values())
    self_roundtrip_ok = True
    for name, row in native.items():
        replay = row.get("self_roundtrip") if isinstance(row, Mapping) else None
        if not isinstance(replay, Mapping):
            self_roundtrip_ok = False
            continue
        parsed_replay = _metrics(replay)
        self_roundtrip_ok = self_roundtrip_ok and (
            _relative(float(parsed_native[name]["volume"]), float(parsed_replay["volume"])) <= native_tol
            and _relative(float(parsed_native[name]["area"]), float(parsed_replay["area"])) <= native_tol
            and parsed_replay["bbox_extent"] == expected_bbox
            and all(int(parsed_replay[key]) == expected_topology[key] for key in expected_topology)
        )

    expected_paths = {
        ("builder", "noheal"),
        ("builder", "heal"),
        ("algebra", "noheal"),
        ("algebra", "heal"),
    }
    observed_paths = {
        (str(row.get("source_mode")), str(row.get("import_mode")))
        for row in external if isinstance(row, Mapping)
    }
    external_topology_ok = all(
        int(row.get("volume_count", 0)) == expected_topology["solid_count"]
        and int(row.get("surface_count", 0)) == expected_topology["face_count"]
        and int(row.get("curve_count", 0)) == expected_topology["edge_count"]
        and int(row.get("vertex_count", 0)) == expected_topology["vertex_count"]
        for row in external if isinstance(row, Mapping)
    ) and len(external) == 4
    external_mass_ok = all(
        float(row.get("volume_relative_error", math.inf)) <= external_tol
        and float(row.get("area_relative_error", math.inf)) <= external_tol
        and float(row.get("bbox_extent_relative_error", math.inf)) <= native_tol
        for row in external if isinstance(row, Mapping)
    ) and len(external) == 4

    checks = {
        "feature_roles_recorded": int(features.get("repeated_top_feature_count", 0)) > 1
        and int(features.get("internal_support_count", 0)) > 0
        and features.get("internal_cavity_present") is True,
        "topology_contract_is_nontrivial_single_solid": expected_topology["solid_count"] == 1
        and expected_topology["shell_count"] == 1
        and expected_topology["face_count"] > 20
        and expected_topology["edge_count"] > expected_topology["face_count"]
        and expected_topology["vertex_count"] > expected_topology["face_count"],
        "official_volume_reproduced_by_both_apis": math.isfinite(official)
        and _relative(float(builder["volume"]), official) <= native_tol
        and _relative(float(algebra["volume"]), official) <= native_tol,
        "native_dual_api_mass_matches": _relative(float(builder["volume"]), float(algebra["volume"])) <= native_tol
        and _relative(float(builder["area"]), float(algebra["area"])) <= native_tol,
        "native_dual_api_bbox_and_topology_match": native_bbox_ok and native_topology_ok,
        "both_step_roundtrips_preserve_mass_and_topology": self_roundtrip_ok,
        "four_external_import_paths_present": observed_paths == expected_paths,
        "external_imports_preserve_topology": external_topology_ok,
        "external_imports_preserve_mass_and_bbox": external_mass_ok,
        "heal_noheal_is_invariant_per_source": all(
            len({
                (float(row.get("volume", math.nan)), int(row.get("surface_count", 0)), int(row.get("curve_count", 0)), int(row.get("vertex_count", 0)))
                for row in external
                if row.get("source_mode") == source_mode
            }) == 1
            for source_mode in ("builder", "algebra")
        ),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_repeated_cavity_dual_api_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "native_dual_api_volume_relative_error": _relative(float(builder["volume"]), float(algebra["volume"])),
            "native_dual_api_area_relative_error": _relative(float(builder["area"]), float(algebra["area"])),
            "maximum_external_volume_relative_error": max(float(row.get("volume_relative_error", math.inf)) for row in external),
            "maximum_external_area_relative_error": max(float(row.get("area_relative_error", math.inf)) for row in external),
            "expected_topology": expected_topology,
        },
        "lesson": (
            "Repeated external features and an internal cavity need an explicit topology contract. "
            "Native API parity, self-roundtrip, and heal/noheal external imports must preserve "
            "mass, extents, and entity counts together; equal volume alone can hide feature loss."
        ),
    }


def build123d_repeated_cavity_source_replay_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Gate immutable dual sources, viewer suppression, STEP artifacts, and external replay."""

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    sources = summary.get("sources")
    artifacts = summary.get("step_artifacts")
    process = summary.get("external_process")
    timing = summary.get("timing_breakdown_s")
    public_gate = summary.get("public_gate")
    if not isinstance(sources, list) or len(sources) != 2:
        raise ValueError("sources must contain two records")
    if not isinstance(artifacts, list) or len(artifacts) != 2:
        raise ValueError("step_artifacts must contain two records")
    if not isinstance(process, Mapping) or not isinstance(timing, Mapping) or not isinstance(public_gate, Mapping):
        raise ValueError("external_process, timing_breakdown_s, and public_gate must be mappings")
    source_modes = {str(row.get("mode")) for row in sources if isinstance(row, Mapping)}
    artifact_modes = {str(row.get("mode")) for row in artifacts if isinstance(row, Mapping)}
    checks = {
        "upstream_tagged_sources_preserved": summary.get("source_kind") == "upstream_native_examples"
        and len(str(summary.get("source_commit", ""))) == 40
        and summary.get("source_files_preserved") is True,
        "builder_and_algebra_source_digests_recorded": source_modes == {"builder", "algebra"}
        and all(len(str(row.get("sha256", ""))) == 64 for row in sources),
        "viewer_suppression_is_stub_only": summary.get("viewer_stub_only") is True,
        "cubit_batch_python_uses_exec_compile_entry": summary.get("cubit_batch_entry_mode")
        == "exec_compile_wrapper",
        "source_parameter_and_feature_contract_bound": int(summary.get("pip_count_x", 0)) > 1
        and int(summary.get("pip_count_y", 0)) > 1
        and int(summary.get("repeated_top_feature_count", 0))
        == int(summary.get("pip_count_x", 0)) * int(summary.get("pip_count_y", 0))
        and summary.get("internal_cavity_present") is True,
        "official_source_assertion_reproduced": summary.get("official_volume_assertion_reproduced") is True,
        "two_step_artifacts_bound": artifact_modes == {"builder", "algebra"}
        and all(len(str(row.get("sha256", ""))) == 64 and int(row.get("bytes", 0)) > 0 for row in artifacts),
        "headless_external_replay_completed": process.get("execution_mode") == "python_api_headless"
        and {"-nographics", "-batch"}.issubset(set(process.get("headless_flags") or []))
        and process.get("gui_daemon_enabled") is False
        and process.get("acceptable") is True,
        "fresh_result_and_no_owned_process_leak": process.get("result_artifact_fresh") is True
        and int(process.get("owned_processes_remaining", -1)) == 0,
        "exactly_four_timing_stages_recorded": len(timing) == 4
        and all(math.isfinite(float(value)) and float(value) >= 0.0 for value in timing.values()),
        "repeated_cavity_public_gate_passed": public_gate.get("policy")
        == "build123d_repeated_cavity_dual_api_gate_v1"
        and public_gate.get("status") == "ok",
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_repeated_cavity_source_replay_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "notes": [
            "Execute upstream Builder and Algebra sources without editing geometry; suppress only optional visualization.",
            "Pass multiline Cubit Python through a one-line exec/compile entry because direct batch playback is line-oriented.",
            "Bind both STEP digests and require the external replay to preserve the repeated-feature topology and internal cavity contract.",
        ],
    }
