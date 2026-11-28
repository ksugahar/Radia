"""Validation gates for level-set extraction and Sculpt all-hex workflows."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence


_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def _sequence(value: object, name: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be a sequence")
    return value


def _normalized(value: object) -> str:
    return " ".join(str(value).strip().lower().split())


def _subsequence(commands: Sequence[str], required: Sequence[str]) -> bool:
    cursor = 0
    for command in commands:
        if cursor < len(required) and required[cursor] in command:
            cursor += 1
    return cursor == len(required)


def cubit_levelset_sculpt_hex_validation_gate(
    summary: Mapping[str, object],
    *,
    quality_floor: float = 0.2,
) -> dict[str, object]:
    """Gate a coarse-to-fine Sculpt all-hex package.

    The source Cubit quality metric owns the configured quality floor.
    """

    summary = _mapping(summary, "summary")
    if quality_floor <= 0.0:
        raise ValueError("quality_floor must be positive")

    mbg_volume = float(summary.get("mbg_volume", math.nan))
    rows = []
    for index, raw in enumerate(
        _sequence(summary.get("mesh_series", []), "mesh_series")
    ):
        row = _mapping(raw, f"mesh_series[{index}]")
        rows.append(
            {
                "label": str(row.get("label", "")).strip().lower(),
                "cell_size": float(row.get("cell_size", math.nan)),
                "hex_count": int(row.get("hex_count", 0)),
                "other_volume_count": int(row.get("other_volume_count", 0)),
                "source_minimum_quality": float(
                    row.get("source_minimum_quality", math.nan)
                ),
            }
        )

    coarse = rows[0] if len(rows) == 2 else {}
    fine = rows[1] if len(rows) == 2 else {}
    checks = {
        "levelset_tet_and_iso_surface_recorded": int(
            summary.get("source_tet_count", 0)
        )
        > 0
        and int(summary.get("iso_triangle_count", 0)) > 0,
        "closed_mbg_volume_recorded": math.isfinite(mbg_volume)
        and mbg_volume > 0.0,
        "coarse_and_fine_rows_present": len(rows) == 2
        and coarse.get("label") == "coarse"
        and fine.get("label") == "fine",
        "refinement_reduces_cell_size_and_adds_hexes": len(rows) == 2
        and coarse["cell_size"] > fine["cell_size"] > 0.0
        and fine["hex_count"] > coarse["hex_count"] > 0,
        "both_meshes_are_all_hex": len(rows) == 2
        and all(
            row["hex_count"] > 0 and row["other_volume_count"] == 0
            for row in rows
        ),
        "source_quality_crosses_configured_floor": len(rows) == 2
        and math.isfinite(coarse["source_minimum_quality"])
        and math.isfinite(fine["source_minimum_quality"])
        and coarse["source_minimum_quality"] < quality_floor
        and fine["source_minimum_quality"] >= quality_floor
        and fine["source_minimum_quality"] > coarse["source_minimum_quality"],
        "mesh_ready_without_solver_ready_overclaim": summary.get("mesh_ready")
        is True
        and summary.get("solver_ready") is False
        and str(summary.get("disposition", "")).strip().lower()
        == "accept_mesh_require_physics_label_handoff",
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_levelset_sculpt_hex_validation_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "mesh_ready": not issues,
        "solver_ready": False,
        "checks": checks,
        "issues": issues,
        "metrics": {
            "source_tet_count": int(summary.get("source_tet_count", 0)),
            "iso_triangle_count": int(summary.get("iso_triangle_count", 0)),
            "mbg_volume": mbg_volume,
            "quality_floor": quality_floor,
            "mesh_series": rows,
        },
        "notes": [
            "A level-set tet mesh is source data; the promoted solver mesh is the Sculpt all-hex result.",
            "A coarse mesh below the quality floor is diagnostic evidence, not a solver handoff.",
            "Mesh readiness does not imply solver readiness until material and boundary labels are attached.",
        ],
    }


def cubit_ato_levelset_sculpt_source_replay_gate(
    summary: Mapping[str, object],
) -> dict[str, object]:
    """Gate official ATO provenance, MBG migration, and headless replay."""

    summary = _mapping(summary, "summary")
    source = _mapping(summary.get("source_mesh", {}), "source_mesh")
    iso = _mapping(summary.get("iso_surface", {}), "iso_surface")
    obsolete = _mapping(
        summary.get("obsolete_stl_negative_control", {}),
        "obsolete_stl_negative_control",
    )
    process = _mapping(summary.get("process", {}), "process")
    public_gate = _mapping(summary.get("public_gate", {}), "public_gate")
    replay = _mapping(summary.get("deterministic_replay", {}), "deterministic_replay")
    timing = _mapping(summary.get("timing_breakdown_s", {}), "timing_breakdown_s")
    flags = {
        _normalized(value)
        for value in _sequence(summary.get("headless_flags", []), "headless_flags")
    }
    commands = [
        _normalized(value)
        for value in _sequence(summary.get("command_trace", []), "command_trace")
    ]
    required_commands = [
        "set dev on",
        'nodal_var "lsd" no_geom',
        'create tri iso tet all nodal_var "lsd"',
        "block 3 4 overwrite",
        "reset",
        'small_bracket_iso.e" no_geom',
        "node in tri in block 3 position fixed",
        "smooth tri in block 4 target free mesh iteration 5",
        "create mesh geom tri all",
        "create vol surf all",
        "save cub5",
        "open",
        "size 0.5",
        "open",
        "size 0.125",
    ]
    source_digest = str(summary.get("source_sha256", "")).strip().lower()
    doc_digest = str(summary.get("source_doc_sha256", "")).strip().lower()
    obsolete_diagnostics = "\n".join(
        str(value).lower()
        for value in _sequence(
            obsolete.get("console_diagnostics", []),
            "obsolete_stl_negative_control.console_diagnostics",
        )
    )
    checks = {
        "official_ato_asset_and_help_recorded": str(summary.get("source_kind", ""))
        == "installed-official-help-ato-levelset-exodus"
        and str(summary.get("source_name", "")).lower() == "small_bracket.exo"
        and str(summary.get("source_doc_name", "")).lower() == "ato_to_mesh.htm",
        "source_and_help_digests_recorded": bool(
            _SHA256_RE.fullmatch(source_digest)
        )
        and bool(_SHA256_RE.fullmatch(doc_digest)),
        "headless_execution_recorded": str(summary.get("binary_name", "")).lower()
        == "coreform_cubit.com"
        and {"-nographics", "-batch"}.issubset(flags)
        and summary.get("gui_daemon_enabled") is False,
        "source_levelset_counts_recorded": int(source.get("node_count", 0))
        == 10021
        and int(source.get("tet_count", 0)) == 52472,
        "documented_iso_blocks_recorded": int(iso.get("block_3_tri_count", 0))
        == 2704
        and int(iso.get("block_4_tri_count", 0)) == 8148,
        "mbg_sculpt_migration_order_preserved": _subsequence(
            commands, required_commands
        ),
        "obsolete_stl_failure_is_observed": obsolete.get("command_returned")
        is True
        and obsolete.get("artifact_exists") is False
        and "could not export the data" in obsolete_diagnostics
        and "stl is not a valid type" in obsolete_diagnostics
        and str(obsolete.get("migration", "")).strip()
        == "export iso blocks to Exodus, reconstruct MBG, then sculpt volume",
        "headless_process_errors_are_classified": process.get("acceptable") is True
        and process.get("result_artifact_fresh") is True
        and not list(process.get("unexpected_error_lines", []))
        and int(process.get("owned_processes_remaining", -1)) == 0,
        "public_mesh_gate_passed_without_solver_promotion": public_gate.get(
            "status"
        )
        == "ok"
        and public_gate.get("mesh_ready") is True
        and public_gate.get("solver_ready") is False,
        "two_replays_match": int(replay.get("repeat_count", 0)) >= 2
        and replay.get("stable_fields_match") is True,
        "four_stage_timing_recorded": len(timing) == 4
        and all(float(value) >= 0.0 for value in timing.values()),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "cubit_ato_levelset_sculpt_source_replay_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "mesh_ready": not issues,
        "solver_ready": False,
        "checks": checks,
        "issues": issues,
        "notes": [
            "Preserve nodal_var during the source import; without it, the zero level-set cannot be reconstructed.",
            "On current Cubit, a documented free-triangle STL command can return true while producing no file; console evidence is mandatory.",
            "The durable migration is iso-block Exodus, mesh-based geometry reconstruction, then Sculpt volume.",
            "Run Cubit headless and retain two deterministic replays before accepting the mesh lesson.",
        ],
    }
