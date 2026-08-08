"""Evidence gates for dense face-first perforation and external CAD replay."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import math


_CONSTRUCTION_MODE = "face_from_outer_wire_and_hole_wires_then_single_extrude"
_SOURCE_EXECUTION_MODES = {"exact_source", "exact_source_with_display_stub"}
_HEADLESS_FLAGS = {"-nographics", "-batch"}
_IMPORT_MODES = {"heal", "noheal"}
_NONZERO_EXIT_POLICY = "fresh_pass_artifact_plus_allowlisted_startup_diagnostics"


def _mapping(value: object, name: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def _sequence(value: object, name: str) -> Sequence[object]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be a sequence")
    return value


def _count(value: object, name: str, *, minimum: int = 0) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _positive_number(value: object, name: str) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ValueError(f"{name} must be numeric")
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _digest(value: object, length: int) -> bool:
    if not isinstance(value, str) or len(value) != length:
        return False
    return all(character in "0123456789abcdef" for character in value.lower())


def _relative_error(measured: float, reference: float) -> float:
    return abs(measured - reference) / reference


def face_first_perforation_handoff_gate(
    summary: Mapping[str, object], *, volume_rtol: float = 1.0e-9
) -> dict[str, object]:
    """Bind generated hole counts to native, STEP/BREP, and external topology."""

    summary = _mapping(summary, "summary")
    tolerance = float(volume_rtol)
    if not math.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("volume_rtol must be finite and >= 0")

    construction = _mapping(summary.get("construction"), "construction")
    requested = _count(
        construction.get("requested_hole_count"),
        "construction.requested_hole_count",
    )
    generated_locations = _count(
        construction.get("generated_location_count"),
        "construction.generated_location_count",
    )
    hole_wires = _count(
        construction.get("hole_wire_count"), "construction.hole_wire_count"
    )
    face_wires = _count(
        construction.get("face_wire_count"), "construction.face_wire_count", minimum=1
    )
    hole_sides = _count(
        construction.get("hole_side_count"),
        "construction.hole_side_count",
        minimum=3,
    )
    outer_sides = _count(
        construction.get("outer_side_count"),
        "construction.outer_side_count",
        minimum=3,
    )
    expected_surfaces = 2 + outer_sides + requested * hole_sides

    native = _mapping(summary.get("native"), "native")
    native_body_count = _count(native.get("body_count"), "native.body_count")
    native_surface_count = _count(
        native.get("surface_count"), "native.surface_count"
    )
    native_volume = _positive_number(native.get("volume"), "native.volume")

    self_roundtrips = _sequence(
        summary.get("self_roundtrips"), "self_roundtrips"
    )
    external_imports = _sequence(
        summary.get("external_imports"), "external_imports"
    )

    replay_rows: list[dict[str, object]] = []
    for index, value in enumerate(self_roundtrips):
        row = _mapping(value, f"self_roundtrips[{index}]")
        replay_rows.append(
            {
                "format": str(row.get("format") or "").lower(),
                "body_count": _count(
                    row.get("body_count"), f"self_roundtrips[{index}].body_count"
                ),
                "surface_count": _count(
                    row.get("surface_count"),
                    f"self_roundtrips[{index}].surface_count",
                ),
                "volume": _positive_number(
                    row.get("volume"), f"self_roundtrips[{index}].volume"
                ),
            }
        )

    external_rows: list[dict[str, object]] = []
    for index, value in enumerate(external_imports):
        row = _mapping(value, f"external_imports[{index}]")
        external_rows.append(
            {
                "mode": str(row.get("mode") or "").lower(),
                "body_count": _count(
                    row.get("body_count"), f"external_imports[{index}].body_count"
                ),
                "surface_count": _count(
                    row.get("surface_count"),
                    f"external_imports[{index}].surface_count",
                ),
                "volume": _positive_number(
                    row.get("volume"), f"external_imports[{index}].volume"
                ),
            }
        )

    replay_formats = {str(row["format"]) for row in replay_rows}
    external_modes = {str(row["mode"]) for row in external_rows}
    all_exchange_rows = replay_rows + external_rows
    volume_errors = [
        _relative_error(float(row["volume"]), native_volume)
        for row in all_exchange_rows
    ]
    checks = {
        "face_first_single_extrude": construction.get("mode") == _CONSTRUCTION_MODE,
        "requested_locations_match": generated_locations == requested,
        "generated_hole_wires_match": hole_wires == requested,
        "face_contains_outer_and_hole_wires": face_wires == requested + 1,
        "native_single_body": native_body_count == 1,
        "native_surface_topology_matches": native_surface_count == expected_surfaces,
        "step_and_brep_replayed": replay_formats == {"step", "brep"},
        "self_roundtrip_topology_matches": bool(replay_rows)
        and all(
            row["body_count"] == 1 and row["surface_count"] == expected_surfaces
            for row in replay_rows
        ),
        "heal_and_noheal_imported": external_modes == _IMPORT_MODES,
        "external_topology_matches": bool(external_rows)
        and all(
            row["body_count"] == 1 and row["surface_count"] == expected_surfaces
            for row in external_rows
        ),
        "all_volumes_match": bool(all_exchange_rows)
        and max(volume_errors, default=math.inf) <= tolerance,
    }
    return {
        "policy": "build123d_face_first_perforation_handoff_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
        "expected_surface_count": expected_surfaces,
        "metrics": {
            "requested_hole_count": requested,
            "generated_location_count": generated_locations,
            "hole_wire_count": hole_wires,
            "face_wire_count": face_wires,
            "native_surface_count": native_surface_count,
            "self_roundtrip_surface_counts": [
                row["surface_count"] for row in replay_rows
            ],
            "external_surface_counts": [
                row["surface_count"] for row in external_rows
            ],
            "max_volume_relative_error": max(volume_errors, default=math.inf),
            "volume_rtol": tolerance,
        },
        "lesson": (
            "A dense perforation handoff must bind requested locations, generated "
            "hole wires, the face wire inventory, and downstream wall topology; "
            "caller-declared hole count and volume alone are insufficient."
        ),
    }


def face_first_perforation_source_replay_gate(
    summary: Mapping[str, object],
) -> dict[str, object]:
    """Gate immutable source replay and its owned headless external-CAD process."""

    summary = _mapping(summary, "summary")
    source = _mapping(summary.get("source"), "source")
    execution = _mapping(summary.get("execution"), "execution")
    exports = _mapping(summary.get("exports"), "exports")
    external = _mapping(summary.get("external"), "external")

    expected_digest = source.get("expected_sha256")
    observed_digest = source.get("observed_sha256")
    flags = {
        str(value).lower()
        for value in _sequence(external.get("headless_flags"), "external.headless_flags")
    }
    import_modes = {
        str(value).lower()
        for value in _sequence(external.get("import_modes"), "external.import_modes")
    }
    run_count = _count(
        execution.get("native_run_count"), "execution.native_run_count"
    )
    owned_processes = _count(
        external.get("owned_processes_remaining"),
        "external.owned_processes_remaining",
    )
    exit_code = external.get("process_exit_code")
    if not isinstance(exit_code, int) or isinstance(exit_code, bool):
        raise ValueError("external.process_exit_code must be an integer")
    nonzero_exit_classified = (
        exit_code != 0
        and external.get("process_exit_policy") == _NONZERO_EXIT_POLICY
        and external.get("allowlisted_startup_diagnostics_only") is True
        and external.get("artifact_fresh") is True
    )
    checks = {
        "upstream_source_bound": source.get("kind") == "upstream_native_example",
        "source_tag_recorded": bool(str(source.get("tag") or "")),
        "source_commit_bound": _digest(source.get("commit"), 40),
        "source_digest_bound": _digest(expected_digest, 64)
        and expected_digest == observed_digest,
        "source_preserved": source.get("preserved") is True,
        "exact_source_executed": execution.get("mode") in _SOURCE_EXECUTION_MODES,
        "source_replayed_at_least_twice": run_count >= 2,
        "native_runs_deterministic": execution.get("native_runs_deterministic")
        is True,
        "source_counts_observed": execution.get("source_counts_observed") is True,
        "step_digest_bound": _digest(exports.get("step_sha256"), 64),
        "brep_digest_bound": _digest(exports.get("brep_sha256"), 64),
        "headless_external_execution": external.get("execution_mode")
        == "python_api_headless"
        and _HEADLESS_FLAGS.issubset(flags),
        "heal_and_noheal_recorded": import_modes == _IMPORT_MODES,
        "no_gui_daemon": external.get("gui_daemon_enabled") is False,
        "owned_processes_cleaned": owned_processes == 0,
        "fresh_external_artifact": external.get("artifact_fresh") is True,
        "process_outcome_classified": exit_code == 0 or nonzero_exit_classified,
        "public_handoff_gate_passed": summary.get("public_gate_status") == "ok",
    }
    return {
        "policy": "build123d_face_first_perforation_source_replay_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, passed in checks.items() if not passed],
        "process_exit": {
            "code": exit_code,
            "classified_nonzero": nonzero_exit_classified,
        },
        "lesson": (
            "Source-native learning requires an immutable source digest, repeated "
            "native execution, fresh headless heal/noheal replay, and owned-process "
            "cleanup in addition to the public geometry gate."
        ),
    }
