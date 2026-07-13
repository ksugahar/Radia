"""Gates for partial-volume hex diagnosis in Cubit Power Tools replays."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence


_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_VOLUME_KINDS = ("hex", "tet", "pyramid", "wedge")


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def _sequence(value: object, name: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be a sequence")
    return value


def _counts(value: object, name: str = "element_counts") -> dict[str, int]:
    row = _mapping(value, name)
    return {str(key).strip().lower(): int(count) for key, count in row.items()}


def _ids(value: object, name: str) -> list[int]:
    return sorted({int(item) for item in _sequence(value, name)})


def _normalized(value: object) -> str:
    return " ".join(str(value).strip().lower().split())


def _subsequence(commands: Sequence[str], required: Sequence[str]) -> bool:
    cursor = 0
    for command in commands:
        if cursor < len(required) and required[cursor] in command:
            cursor += 1
    return cursor == len(required)


def cubit_partial_volume_hex_diagnosis_gate(
    summary: Mapping[str, object],
    *,
    min_scaled_jacobian: float = 0.2,
) -> dict[str, object]:
    """Accept a truthful rejection of a partial or low-quality all-hex mesh.

    This gate deliberately separates *diagnosis success* from solver readiness.
    A replay can be accepted as useful evidence while ``solver_ready`` remains
    false.  Cubit block membership and command return values are not treated as
    proof that every CAD volume owns three-dimensional elements.
    """

    summary = _mapping(summary, "summary")
    if min_scaled_jacobian <= 0.0:
        raise ValueError("min_scaled_jacobian must be positive")

    source_volume_count = int(summary.get("source_volume_count", 0))
    prepared_ids = _ids(summary.get("prepared_volume_ids", []), "prepared_volume_ids")
    expected_unmeshed = _ids(
        summary.get("expected_unmeshed_volume_ids", []),
        "expected_unmeshed_volume_ids",
    )
    ownership_raw = _mapping(
        summary.get("volume_element_counts", {}), "volume_element_counts"
    )
    ownership = {
        int(volume_id): _counts(row, f"volume_element_counts[{volume_id}]")
        for volume_id, row in ownership_raw.items()
    }
    totals = _counts(summary.get("element_counts", {}))
    actual_unmeshed = sorted(
        volume_id
        for volume_id in prepared_ids
        if sum(ownership.get(volume_id, {}).get(kind, 0) for kind in _VOLUME_KINDS)
        == 0
    )
    owned_totals = {
        kind: sum(ownership.get(volume_id, {}).get(kind, 0) for volume_id in prepared_ids)
        for kind in _VOLUME_KINDS
    }

    quality = _mapping(summary.get("quality", {}), "quality")
    scaled = _mapping(quality.get("scaled_jacobian", {}), "scaled_jacobian")
    scaled_min = float(scaled.get("minimum", scaled.get("min", math.nan)))
    gmsh = _mapping(summary.get("gmsh_export", {}), "gmsh_export")
    block_ids = _ids(
        summary.get("exported_block_volume_ids", []), "exported_block_volume_ids"
    )
    reported_solver_ready = summary.get("solver_ready") is True
    disposition = str(summary.get("disposition", "")).strip().lower()

    checks = {
        "decomposition_created_multiple_volumes": source_volume_count == 1
        and len(prepared_ids) >= 2,
        "volume_ownership_rows_complete": set(ownership) == set(prepared_ids),
        "partial_mesh_is_explicit": bool(actual_unmeshed)
        and actual_unmeshed == expected_unmeshed
        and len(actual_unmeshed) < len(prepared_ids),
        "produced_volume_elements_are_hex_only": owned_totals["hex"] > 0
        and all(owned_totals[kind] == 0 for kind in ("tet", "pyramid", "wedge")),
        "owned_counts_match_global_counts": all(
            owned_totals[kind] == totals.get(kind, 0) for kind in _VOLUME_KINDS
        ),
        "quality_rejection_is_measured": math.isfinite(scaled_min)
        and scaled_min < min_scaled_jacobian,
        "block_membership_does_not_hide_missing_mesh": set(block_ids)
        == set(prepared_ids)
        and bool(set(block_ids).intersection(actual_unmeshed)),
        "gmsh_ascii_v41_matches_produced_mesh": str(gmsh.get("mesh_format", ""))
        == "4.1"
        and gmsh.get("binary") is False
        and int(gmsh.get("hex_count", -1)) == totals.get("hex", 0)
        and int(gmsh.get("other_volume_count", -1)) == 0,
        "solver_handoff_is_suppressed": not reported_solver_ready
        and disposition == "reject_partial_or_low_quality",
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_partial_volume_hex_diagnosis_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "solver_ready": False,
        "checks": checks,
        "issues": issues,
        "metrics": {
            "prepared_volume_count": len(prepared_ids),
            "meshed_volume_count": len(prepared_ids) - len(actual_unmeshed),
            "unmeshed_volume_ids": actual_unmeshed,
            "hex_count": totals.get("hex", 0),
            "scaled_jacobian_minimum": scaled_min,
            "scaled_jacobian_threshold": min_scaled_jacobian,
        },
        "notes": [
            "A successful command return or block assignment does not prove that every CAD volume owns 3-D elements.",
            "Accepting the diagnostic does not promote the mesh: partial ownership or sub-threshold quality keeps solver_ready false.",
            "Compare per-volume ownership with the independent Gmsh inventory before solver handoff.",
        ],
    }


def cubit_power_tools_cleanup_source_replay_gate(
    summary: Mapping[str, object],
) -> dict[str, object]:
    """Gate an official Power Tools cleanup trace and its headless diagnosis."""

    summary = _mapping(summary, "summary")
    process = _mapping(summary.get("process", {}), "process")
    raw = _mapping(summary.get("raw_control", {}), "raw_control")
    cleanup = _mapping(summary.get("cleanup_replay", {}), "cleanup_replay")
    explicit = _mapping(summary.get("explicit_sweep_attempt", {}), "explicit_sweep_attempt")
    public_gate = _mapping(summary.get("public_gate", {}), "public_gate")
    timing = _mapping(summary.get("timing_breakdown_s", {}), "timing_breakdown_s")
    docs = _mapping(summary.get("source_doc_sha256", {}), "source_doc_sha256")

    commands = [
        _normalized(value)
        for value in _sequence(summary.get("command_trace", []), "command_trace")
    ]
    required = [
        "import acis",
        "healer autoheal body all",
        "split surface 22",
        "split surface 10",
        "webcut volume 1 with plane normal to curve 35 close_to vertex 51",
        "webcut volume 1 with plane normal to curve 35 close_to vertex 49",
        "remove surface 17 extend",
        "remove surface 15 extend",
        "tweak surface 16 offset -0.9",
        "imprint volume all",
        "merge volume all",
        "composite create surface 9 27",
        "composite create surface 6 24",
        "composite create surface 11 25 35",
        "volume all scheme auto",
        "mesh volume all",
    ]
    raw_counts = _counts(raw.get("element_counts", {}), "raw_control.element_counts")
    cleanup_counts = _counts(
        cleanup.get("element_counts", {}), "cleanup_replay.element_counts"
    )
    flags = {
        _normalized(value)
        for value in _sequence(summary.get("headless_flags", []), "headless_flags")
    }
    diagnostics = "\n".join(
        str(value).lower()
        for value in _sequence(
            explicit.get("console_diagnostics", []),
            "explicit_sweep_attempt.console_diagnostics",
        )
    )
    source_sha256 = str(summary.get("source_sha256", "")).strip().lower()
    replay = _mapping(summary.get("deterministic_replay", {}), "deterministic_replay")

    checks = {
        "official_power_tools_asset_recorded": str(summary.get("source_kind", ""))
        == "installed-official-help-power-tools-cad"
        and str(summary.get("source_name", "")).lower() == "knuckle.sat",
        "source_and_tutorial_digests_recorded": bool(
            _SHA256_RE.fullmatch(source_sha256)
        )
        and len(docs) >= 9
        and all(_SHA256_RE.fullmatch(str(value).lower()) for value in docs.values()),
        "headless_console_execution": str(summary.get("binary_name", "")).lower()
        == "coreform_cubit.com"
        and {"-nographics", "-batch"}.issubset(flags)
        and summary.get("gui_daemon_enabled") is False,
        "official_cleanup_order_preserved": _subsequence(commands, required),
        "raw_autoscheme_rejects_unsplit_volume": int(raw.get("volume_count", 0))
        == 1
        and not str(raw.get("volume_scheme", "")).strip()
        and all(raw_counts.get(kind, 0) == 0 for kind in _VOLUME_KINDS),
        "cleanup_reaches_partial_all_hex_state": int(
            cleanup.get("volume_count", 0)
        )
        == 3
        and cleanup.get("volume_schemes") == {"1": "sweep", "2": "sweep", "3": ""}
        and cleanup_counts.get("hex", 0) > 0
        and all(cleanup_counts.get(kind, 0) == 0 for kind in ("tet", "pyramid", "wedge")),
        "explicit_sweep_failure_is_observed_not_inferred": explicit.get(
            "command_returned"
        )
        is True
        and int(explicit.get("volume_id", 0)) == 3
        and int(explicit.get("volume_hex_count", -1)) == 0
        and "internal loop(s) do not have a corner or end node" in diagnostics
        and "meshing unsuccessful" in diagnostics,
        "process_errors_are_classified": process.get("acceptable") is True
        and process.get("result_artifact_fresh") is True
        and not list(process.get("unexpected_error_lines", [])),
        "public_diagnosis_passed_without_promotion": public_gate.get("status")
        == "ok"
        and public_gate.get("solver_ready") is False,
        "deterministic_replay_matches": int(replay.get("repeat_count", 0)) >= 2
        and replay.get("stable_fields_match") is True,
        "four_stage_timing_recorded": len(timing) == 4
        and all(float(value) >= 0.0 for value in timing.values()),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_power_tools_cleanup_source_replay_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "solver_ready": False,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Capture the console stream: cubit.cmd returning True is not sufficient evidence that a volume meshed.",
            "Power Tools cleanup improves meshability, but every decomposed volume still needs explicit 3-D element ownership and quality checks.",
            "Treat partial export as diagnostic evidence, never as a solver-ready package.",
        ],
    }
