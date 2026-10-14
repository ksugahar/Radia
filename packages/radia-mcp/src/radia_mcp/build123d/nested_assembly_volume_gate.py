"""Nested assembly mass-property and source-replay validation gates."""

from __future__ import annotations

import math
from typing import Mapping


_STUD_WALL_SOURCE_SHA256 = (
    "f1f707eaa299d9d858c30ccdf0f01f6e9374ac47e5195b2b935d3d1901fcdc00"
)
_STUD_WALL_SOURCE_COMMIT = "fa8e93687c2e6069d0eae0e4b0b8ae128e33de1f"


def _finite(row: Mapping[str, object], key: str) -> float:
    value = float(row.get(key, math.nan))
    if not math.isfinite(value):
        raise ValueError(f"{key} must be finite")
    return value


def _integer(row: Mapping[str, object], key: str) -> int:
    value = int(row.get(key, -1))
    return value


def _relative_error(measured: float, reference: float) -> float:
    if not math.isfinite(measured) or not math.isfinite(reference) or reference <= 0.0:
        return math.inf
    return abs(measured - reference) / reference


def nested_assembly_volume_gate(
    summary: Mapping[str, object],
    *,
    same_kernel_rtol: float = 1.0e-9,
    external_rtol: float = 2.0e-6,
    union_rtol: float = 1.0e-12,
) -> dict[str, object]:
    """Distinguish a zero-valued parent Compound from an empty CAD model.

    A nested build123d ``Compound`` can report ``volume == 0`` while its leaf
    solids all have positive volume.  Acceptance therefore comes from the leaf
    inventory, same-kernel roundtrips, and an independent CAD import.  A CAD
    handoff may pass without claiming that a volume mesh or solver input exists.
    """

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    tolerances = (float(same_kernel_rtol), float(external_rtol), float(union_rtol))
    if any(value < 0.0 or not math.isfinite(value) for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")

    native = summary.get("native") or {}
    roundtrips = summary.get("same_kernel_roundtrips") or {}
    external = summary.get("external") or {}
    if not isinstance(native, Mapping):
        raise ValueError("native must be a mapping")
    if not isinstance(roundtrips, Mapping):
        raise ValueError("same_kernel_roundtrips must be a mapping")
    if not isinstance(external, Mapping):
        raise ValueError("external must be a mapping")

    native_leaf_sum = _finite(native, "leaf_volume_sum_mm3")
    native_parent_volume = _finite(native, "top_level_volume_mm3")
    native_solid_count = _integer(native, "solid_count")
    native_topology = {
        "surface_count": _integer(native, "face_count"),
        "curve_count": _integer(native, "edge_count"),
        "vertex_count": _integer(native, "vertex_count"),
    }

    same_kernel_errors: dict[str, float] = {}
    same_kernel_rows_ok: dict[str, bool] = {}
    for name in ("step", "brep"):
        row = roundtrips.get(name) or {}
        if not isinstance(row, Mapping):
            raise ValueError(f"same_kernel_roundtrips.{name} must be a mapping")
        row_sum = _finite(row, "leaf_volume_sum_mm3")
        error = _relative_error(row_sum, native_leaf_sum)
        same_kernel_errors[name] = error
        same_kernel_rows_ok[name] = (
            row.get("is_valid") is True
            and _integer(row, "solid_count") == native_solid_count
            and _integer(row, "leaf_positive_volume_count") == native_solid_count
            and _integer(row, "face_count") == native_topology["surface_count"]
            and _integer(row, "edge_count") == native_topology["curve_count"]
            and _integer(row, "vertex_count") == native_topology["vertex_count"]
            and error <= tolerances[0]
        )

    external_replays = external.get("replays") or []
    if not isinstance(external_replays, list):
        raise ValueError("external.replays must be a list")
    observed = []
    for replay in external_replays:
        if not isinstance(replay, Mapping):
            raise ValueError("each external replay must be a mapping")
        imported = replay.get("imported") or {}
        united = replay.get("united") or {}
        import_command = replay.get("import_command") or {}
        unite_command = replay.get("unite_command") or {}
        imported_sum = _finite(imported, "volume_sum_mm3")
        united_sum = _finite(united, "volume_sum_mm3")
        observed.append(
            {
                "index": int(replay.get("index", -1)),
                "import_command_ok": import_command.get("returned") is True
                and import_command.get("exception") is None,
                "unite_command_ok": unite_command.get("returned") is True
                and unite_command.get("exception") is None,
                "imported_volume_count": _integer(imported, "volume_count"),
                "imported_positive_volume_count": _integer(
                    imported, "positive_volume_count"
                ),
                "imported_surface_count": _integer(imported, "surface_count"),
                "imported_curve_count": _integer(imported, "curve_count"),
                "imported_vertex_count": _integer(imported, "vertex_count"),
                "imported_volume_sum_mm3": imported_sum,
                "external_volume_relative_error": _relative_error(
                    imported_sum, native_leaf_sum
                ),
                "united_volume_count": _integer(united, "volume_count"),
                "united_positive_volume_count": _integer(
                    united, "positive_volume_count"
                ),
                "united_volume_sum_mm3": united_sum,
                "union_volume_relative_error": _relative_error(
                    united_sum, imported_sum
                ),
            }
        )

    replay_projections = [
        {
            key: value
            for key, value in row.items()
            if key not in {"index"}
        }
        for row in observed
    ]
    checks = {
        "native_is_valid_nested_compound": str(native.get("type", "")).lower()
        == "compound"
        and native.get("is_valid") is True
        and _integer(native, "direct_child_count") >= 2,
        "parent_scalar_volume_is_zero": native_parent_volume == 0.0,
        "leaf_inventory_is_positive_and_complete": native_leaf_sum > 0.0
        and native_solid_count >= 2
        and _integer(native, "leaf_positive_volume_count") == native_solid_count,
        "same_kernel_step_and_brep_preserve_leaf_measure_and_topology": set(
            roundtrips
        )
        >= {"step", "brep"}
        and all(same_kernel_rows_ok.values()),
        "two_external_replays_recorded": len(observed) >= 2,
        "external_import_commands_succeeded": bool(observed)
        and all(row["import_command_ok"] for row in observed),
        "external_import_preserves_positive_leaf_ownership": bool(observed)
        and all(
            row["imported_volume_count"] == native_solid_count
            and row["imported_positive_volume_count"] == native_solid_count
            for row in observed
        ),
        "external_import_preserves_topology": bool(observed)
        and all(
            row["imported_surface_count"] == native_topology["surface_count"]
            and row["imported_curve_count"] == native_topology["curve_count"]
            and row["imported_vertex_count"] == native_topology["vertex_count"]
            for row in observed
        ),
        "external_leaf_volume_sum_matches": bool(observed)
        and all(
            row["external_volume_relative_error"] <= tolerances[1]
            for row in observed
        ),
        "external_unite_closes_without_volume_loss": bool(observed)
        and all(
            row["unite_command_ok"]
            and row["united_volume_count"] == 1
            and row["united_positive_volume_count"] == 1
            and row["union_volume_relative_error"] <= tolerances[2]
            for row in observed
        ),
        "external_replays_are_deterministic": len(replay_projections) >= 2
        and all(row == replay_projections[0] for row in replay_projections[1:]),
        "cad_handoff_ready_is_explicit": external.get("cad_handoff_ready") is True,
        "mesh_and_solver_readiness_are_not_overclaimed": external.get(
            "mesh_attempted"
        )
        is False
        and external.get("solver_ready") is False,
    }
    issues = [name for name, ok in checks.items() if not ok]
    diagnosis = (
        "valid_nested_compound_zero_parent_scalar"
        if not issues
        else "empty_or_incompletely_verified_nested_assembly"
    )
    return {
        "policy": "build123d_nested_assembly_volume_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "diagnosis": diagnosis,
        "parent_volume_zero_is_not_empty": all(
            checks[name]
            for name in (
                "native_is_valid_nested_compound",
                "parent_scalar_volume_is_zero",
                "leaf_inventory_is_positive_and_complete",
            )
        ),
        "cad_handoff_ready": not issues,
        "solver_ready": False,
        "checks": checks,
        "issues": issues,
        "metrics": {
            "native_parent_volume_mm3": native_parent_volume,
            "native_leaf_volume_sum_mm3": native_leaf_sum,
            "native_leaf_solid_count": native_solid_count,
            "same_kernel_relative_errors": same_kernel_errors,
            "maximum_external_volume_relative_error": max(
                (row["external_volume_relative_error"] for row in observed),
                default=math.inf,
            ),
            "maximum_union_volume_relative_error": max(
                (row["union_volume_relative_error"] for row in observed),
                default=math.inf,
            ),
        },
        "external_replays": observed,
        "notes": [
            "For nested compounds, sum positive leaf solids; do not treat the parent scalar volume as the physical assembly volume.",
            "Require STEP/BREP roundtrips and an independent CAD volume sum before accepting the handoff.",
            "A successful unite proves CAD closure for this replay, not volume-mesh or solver readiness.",
        ],
    }


def stud_wall_source_replay_gate(summary: Mapping[str, object]) -> dict[str, object]:
    """Gate the tagged stud-wall source, joints, and headless CAD replay."""

    if not isinstance(summary, Mapping):
        raise ValueError("summary must be a mapping")
    build = summary.get("build") or {}
    external = summary.get("external") or {}
    timing = summary.get("timing_breakdown_s") or {}
    if not isinstance(build, Mapping) or not isinstance(external, Mapping):
        raise ValueError("build and external must be mappings")
    source = build.get("source") or {}
    replays = build.get("replays") or []
    files = build.get("files") or {}
    if not isinstance(source, Mapping):
        raise ValueError("build.source must be a mapping")
    if not isinstance(replays, list):
        raise ValueError("build.replays must be a list")
    if not isinstance(timing, Mapping):
        raise ValueError("timing_breakdown_s must be a mapping")

    replay_equal = len(replays) == 2 and replays[0] == replays[1]
    first = replays[0] if replays else {}
    x_wall = first.get("x_wall") or {}
    y_wall = first.get("y_wall") or {}
    assembly = first.get("assembly") or {}
    child_types = first.get("wall_child_types") or {}
    external_process = external.get("process") or {}
    external_replays = external.get("replays") or []
    checks = {
        "tagged_upstream_source_identity_is_exact": source.get("repository")
        == "gumyr/build123d"
        and source.get("tag") == "v0.10.0"
        and source.get("commit") == _STUD_WALL_SOURCE_COMMIT
        and source.get("path") == "examples/stud_wall.py"
        and source.get("sha256") == _STUD_WALL_SOURCE_SHA256,
        "source_copy_is_immutable_and_viewer_only_stubbed": source.get(
            "copy_sha256"
        )
        == _STUD_WALL_SOURCE_SHA256
        and source.get("source_preserved") is True
        and source.get("display_stubbed_only") is True,
        "two_exact_source_replays_match": replay_equal,
        "wall_inventory_is_thirteen_plus_ten_studs": int(
            x_wall.get("solid_count", -1)
        )
        == 13
        and int(y_wall.get("solid_count", -1)) == 10
        and int(assembly.get("solid_count", -1)) == 23
        and set(child_types.get("x_wall") or []) == {"Stud"}
        and set(child_types.get("y_wall") or []) == {"Stud"},
        "wall_rigid_joint_names_are_preserved": first.get(
            "x_wall_joint_names"
        )
        == ["end0", "inside0"]
        and first.get("y_wall_joint_names") == ["end0", "inside0"],
        "neutral_cad_artifact_digests_are_bound": all(
            len(str((files.get(name) or {}).get("sha256", ""))) == 64
            for name in ("step", "brep")
        ),
        "nested_volume_diagnosis_is_verified": summary.get(
            "nested_gate_status"
        )
        == "ok"
        and summary.get("nested_gate_diagnosis")
        == "valid_nested_compound_zero_parent_scalar",
        "headless_external_replay_is_fresh_and_clean": external.get(
            "execution_mode"
        )
        == "headless_python_api_synchronous_commands"
        and {"-nographics", "-batch"}.issubset(
            set(external.get("headless_flags") or [])
        )
        and external.get("gui_daemon_enabled") is False
        and len(external_replays) == 2
        and external_process.get("acceptable") is True
        and external_process.get("result_artifact_fresh") is True
        and int(external_process.get("owned_processes_remaining", -1)) == 0,
        "cad_handoff_not_solver_readiness": external.get("cad_handoff_ready")
        is True
        and external.get("mesh_attempted") is False
        and external.get("solver_ready") is False,
        "four_dominant_timing_stages_recorded": len(timing) == 4
        and all(float(value) >= 0.0 and math.isfinite(float(value)) for value in timing.values()),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_stud_wall_source_replay_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "cad_handoff_ready": not issues,
        "solver_ready": False,
        "notes": [
            "RigidJoint names and assembly hierarchy are source metadata and must accompany neutral CAD exports.",
            "The parent Compound scalar may be zero even when all 23 Stud leaves are positive and portable.",
            "Keep CAD handoff acceptance separate from mesh and solver acceptance.",
        ],
    }
