"""Validation gates for a lofted shell crossing CAD kernels."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


def _mapping(value: object, field: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{field} must be a mapping")
    return value


def _rows(value: object, field: str) -> list[Mapping[str, object]]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ValueError(f"{field} must be a sequence of mappings")
    rows = list(value)
    if not all(isinstance(row, Mapping) for row in rows):
        raise ValueError(f"{field} must contain mappings")
    return rows


def _finite(value: object, field: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} must be a finite number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{field} must be a finite number")
    return result


def _positive(value: object, field: str) -> float:
    result = _finite(value, field)
    if result <= 0.0:
        raise ValueError(f"{field} must be positive")
    return result


def _relative_error(value: float, reference: float) -> float:
    return abs(value - reference) / abs(reference)


def _native_topology(row: Mapping[str, object]) -> tuple[int, int, int, int]:
    return (
        int(row.get("solid_count", -1)),
        int(row.get("face_count", -1)),
        int(row.get("edge_count", -1)),
        int(row.get("vertex_count", -1)),
    )


def _external_topology(row: Mapping[str, object]) -> tuple[int, int, int, int]:
    return (
        int(row.get("volume_count", -1)),
        int(row.get("surface_count", -1)),
        int(row.get("curve_count", -1)),
        int(row.get("vertex_count", -1)),
    )


def build123d_lofted_shell_handoff_gate(
    summary: Mapping[str, object],
    *,
    max_step_volume_relative_error: float = 1.0e-7,
    max_step_area_relative_error: float = 1.0e-8,
    max_brep_relative_error: float = 1.0e-12,
    min_external_volume_relative_error: float = 1.0e-4,
    max_external_volume_relative_error: float = 5.0e-3,
    max_external_area_relative_error: float = 1.0e-4,
    max_mode_relative_spread: float = 1.0e-12,
) -> dict[str, object]:
    """Accept a bounded curved-shell CAD handoff without solver overclaim.

    The non-zero lower bound is deliberate: this policy describes a measured,
    repeatable cross-kernel approximation, not exact mass equality.  A larger
    error remains a portability diagnosis and a smaller error should use a
    stricter exact-handoff policy.
    """

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    thresholds = {
        "max_step_volume_relative_error": _positive(
            max_step_volume_relative_error, "max_step_volume_relative_error"
        ),
        "max_step_area_relative_error": _positive(
            max_step_area_relative_error, "max_step_area_relative_error"
        ),
        "max_brep_relative_error": _positive(
            max_brep_relative_error, "max_brep_relative_error"
        ),
        "min_external_volume_relative_error": _positive(
            min_external_volume_relative_error,
            "min_external_volume_relative_error",
        ),
        "max_external_volume_relative_error": _positive(
            max_external_volume_relative_error,
            "max_external_volume_relative_error",
        ),
        "max_external_area_relative_error": _positive(
            max_external_area_relative_error, "max_external_area_relative_error"
        ),
        "max_mode_relative_spread": _positive(
            max_mode_relative_spread, "max_mode_relative_spread"
        ),
    }
    if (
        thresholds["min_external_volume_relative_error"]
        >= thresholds["max_external_volume_relative_error"]
    ):
        raise ValueError(
            "min_external_volume_relative_error must be smaller than the maximum"
        )

    build = _mapping(summary.get("build"), "build")
    external = _mapping(summary.get("external"), "external")
    native = _mapping(build.get("native"), "build.native")
    same = _mapping(build.get("same_kernel_roundtrips"), "build.same_kernel_roundtrips")
    step = _mapping(same.get("step"), "build.same_kernel_roundtrips.step")
    brep = _mapping(same.get("brep"), "build.same_kernel_roundtrips.brep")
    replays = _rows(external.get("replays"), "external.replays")
    snapshots = [
        _mapping(row.get("snapshot"), f"external.replays[{index}].snapshot")
        for index, row in enumerate(replays)
    ]

    native_volume = _positive(native.get("volume_mm3"), "build.native.volume_mm3")
    native_area = _positive(
        native.get("surface_area_mm2"), "build.native.surface_area_mm2"
    )
    step_volume_error = _relative_error(
        _positive(step.get("volume_mm3"), "step.volume_mm3"), native_volume
    )
    step_area_error = _relative_error(
        _positive(step.get("surface_area_mm2"), "step.surface_area_mm2"),
        native_area,
    )
    brep_volume_error = _relative_error(
        _positive(brep.get("volume_mm3"), "brep.volume_mm3"), native_volume
    )
    brep_area_error = _relative_error(
        _positive(brep.get("surface_area_mm2"), "brep.surface_area_mm2"),
        native_area,
    )
    external_volume_errors = [
        _relative_error(
            _positive(row.get("volume_sum_mm3"), f"external.snapshots[{index}].volume_sum_mm3"),
            native_volume,
        )
        for index, row in enumerate(snapshots)
    ]
    external_area_errors = [
        _relative_error(
            _positive(
                row.get("surface_area_sum_mm2"),
                f"external.snapshots[{index}].surface_area_sum_mm2",
            ),
            native_area,
        )
        for index, row in enumerate(snapshots)
    ]
    modes = [str(row.get("mode", "")) for row in replays]
    native_topology = _native_topology(native)
    step_topology = _native_topology(step)
    brep_topology = _native_topology(brep)
    external_topologies = [_external_topology(row) for row in snapshots]
    external_positive_entity_counts = [
        (
            int(row.get("positive_volume_count", -1)),
            int(row.get("positive_surface_count", -1)),
        )
        for row in snapshots
    ]

    by_mode: dict[str, list[Mapping[str, object]]] = {
        mode: [
            snapshot
            for replay, snapshot in zip(replays, snapshots, strict=True)
            if str(replay.get("mode", "")) == mode
        ]
        for mode in ("heal", "noheal")
    }

    def mass_pair(row: Mapping[str, object]) -> tuple[float, float]:
        return (
            _finite(row.get("volume_sum_mm3"), "snapshot.volume_sum_mm3"),
            _finite(row.get("surface_area_sum_mm2"), "snapshot.surface_area_sum_mm2"),
        )

    mode_repeat_match = all(
        len(rows) == 2 and mass_pair(rows[0]) == mass_pair(rows[1])
        and _external_topology(rows[0]) == _external_topology(rows[1])
        for rows in by_mode.values()
    )
    cross_mode_volume_spread = _relative_error(
        mass_pair(by_mode["heal"][0])[0], mass_pair(by_mode["noheal"][0])[0]
    ) if all(by_mode.values()) else math.inf
    cross_mode_area_spread = _relative_error(
        mass_pair(by_mode["heal"][0])[1], mass_pair(by_mode["noheal"][0])[1]
    ) if all(by_mode.values()) else math.inf

    checks = {
        "native_valid_positive_single_solid": native.get("is_valid") is True
        and native_topology[0] == 1
        and min(native_topology) > 0,
        "step_and_brep_topology_match_native": step_topology == native_topology
        and brep_topology == native_topology,
        "same_kernel_mass_is_within_strict_budgets": step_volume_error
        <= thresholds["max_step_volume_relative_error"]
        and step_area_error <= thresholds["max_step_area_relative_error"]
        and brep_volume_error <= thresholds["max_brep_relative_error"]
        and brep_area_error <= thresholds["max_brep_relative_error"],
        "four_external_replays_cover_heal_and_noheal_twice": len(replays) == 4
        and modes.count("heal") == 2
        and modes.count("noheal") == 2,
        "external_commands_succeeded": all(
            _mapping(row.get("command"), "external.replay.command").get("returned")
            is True
            and _mapping(row.get("command"), "external.replay.command").get("exception")
            is None
            for row in replays
        ),
        "external_topology_matches_native": all(
            topology == native_topology for topology in external_topologies
        ),
        "external_positive_entity_counts_match_topology": all(
            positive_volume_count == topology[0]
            and positive_surface_count == topology[1]
            for (positive_volume_count, positive_surface_count), topology in zip(
                external_positive_entity_counts,
                external_topologies,
                strict=True,
            )
        ),
        "external_replays_are_repeatable_and_heal_invariant": mode_repeat_match
        and cross_mode_volume_spread <= thresholds["max_mode_relative_spread"]
        and cross_mode_area_spread <= thresholds["max_mode_relative_spread"],
        "external_volume_drift_is_nonzero_but_bounded": all(
            thresholds["min_external_volume_relative_error"]
            <= error
            <= thresholds["max_external_volume_relative_error"]
            for error in external_volume_errors
        ),
        "external_area_drift_is_bounded": all(
            error <= thresholds["max_external_area_relative_error"]
            for error in external_area_errors
        ),
        "source_and_external_replays_passed": build.get("pass") is True
        and external.get("pass") is True,
        "cad_handoff_without_mesh_or_solver_overclaim": external.get(
            "cad_handoff_ready"
        )
        is True
        and external.get("mesh_attempted") is False
        and external.get("solver_ready") is False,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_lofted_shell_handoff_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "cad_handoff_ready": not issues,
        "mesh_ready": False,
        "solver_ready": False,
        "checks": checks,
        "issues": issues,
        "metrics": {
            "native_topology": list(native_topology),
            "native_volume_mm3": native_volume,
            "native_surface_area_mm2": native_area,
            "step_volume_relative_error": step_volume_error,
            "step_area_relative_error": step_area_error,
            "brep_volume_relative_error": brep_volume_error,
            "brep_area_relative_error": brep_area_error,
            "external_volume_relative_errors": external_volume_errors,
            "external_area_relative_errors": external_area_errors,
            "external_positive_entity_counts": [
                list(counts) for counts in external_positive_entity_counts
            ],
            "cross_mode_volume_relative_spread": cross_mode_volume_spread,
            "cross_mode_area_relative_spread": cross_mode_area_spread,
            "thresholds": thresholds,
        },
        "notes": [
            "Topology equality is necessary but mass properties still need independent budgets.",
            "A bounded deterministic CAD-kernel drift can permit CAD handoff without proving mesh or solver readiness.",
            "Heal/noheal equality is evidence of stable interpretation, not exact geometric equality.",
        ],
    }


def build123d_loft_example_source_replay_gate(
    summary: Mapping[str, object],
) -> dict[str, object]:
    """Bind the bounded handoff to the immutable upstream loft example."""

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be a mapping")
    build = _mapping(summary.get("build"), "build")
    external = _mapping(summary.get("external"), "external")
    source = _mapping(build.get("source"), "build.source")
    replays = _rows(build.get("replays"), "build.replays")
    source_checks = _mapping(build.get("checks"), "build.checks")
    process = _mapping(external.get("process"), "external.process")
    timing = _mapping(summary.get("timing_breakdown_s"), "timing_breakdown_s")
    source_sha = str(source.get("sha256", "")).lower()
    commit = str(source.get("commit", "")).lower()
    expected_volume = 1306.3405290344635
    first = replays[0] if replays else {}
    first_metrics = _mapping(first.get("metrics"), "build.replays[0].metrics")
    observed_volume = _positive(
        first_metrics.get("volume_mm3"), "build.replays[0].metrics.volume_mm3"
    )
    official_tolerance = _positive(
        first.get("official_assertion_tolerance_mm3"),
        "build.replays[0].official_assertion_tolerance_mm3",
    )
    public_gate = build123d_lofted_shell_handoff_gate(summary)

    checks = {
        "upstream_source_identity_bound": summary.get("source_kind")
        == "upstream-tagged-example-exact-replay-plus-headless-external-cad"
        and source.get("repository") == "gumyr/build123d"
        and source.get("tag") == "v0.10.0"
        and source.get("path") == "examples/loft.py"
        and source.get("license") == "Apache-2.0"
        and len(source_sha) == 64
        and all(character in "0123456789abcdef" for character in source_sha)
        and len(commit) == 40
        and all(character in "0123456789abcdef" for character in commit),
        "viewer_stub_only_and_source_preserved": source.get("source_preserved")
        is True
        and source.get("display_stubbed_only") is True,
        "two_exact_source_replays_match": len(replays) == 2
        and replays[0] == replays[1]
        and source_checks.get("source_replays_are_deterministic") is True,
        "eleven_profiles_two_openings_and_half_mm_shell_recorded": int(
            first.get("slice_count", -1)
        )
        == 10
        and int(first.get("top_bottom_face_count", -1)) == 2
        and math.isclose(
            _finite(first.get("shell_offset_mm"), "shell_offset_mm"),
            0.5,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
        and math.isclose(
            _finite(first.get("official_expected_volume_mm3"), "official_expected_volume_mm3"),
            expected_volume,
            rel_tol=1.0e-12,
            abs_tol=1.0e-9,
        ),
        "official_volume_assertion_reproduced": abs(observed_volume - expected_volume)
        < official_tolerance
        and source_checks.get("official_volume_assertion_reproduced") is True,
        "all_source_roundtrip_checks_passed": build.get("pass") is True
        and all(value is True for value in source_checks.values()),
        "public_handoff_gate_passed": public_gate["status"] == "ok"
        and public_gate["cad_handoff_ready"] is True
        and public_gate["solver_ready"] is False,
        "external_replay_is_headless": external.get("execution_mode")
        == "headless_python_api_synchronous_commands"
        and {"-nographics", "-batch"}.issubset(
            set(external.get("headless_flags") or [])
        )
        and external.get("gui_daemon_enabled") is False,
        "external_process_is_fresh_classified_and_clean": process.get("acceptable")
        is True
        and process.get("result_artifact_fresh") is True
        and not list(process.get("unexpected_error_lines") or [])
        and int(process.get("owned_processes_remaining", -1)) == 0,
        "exactly_four_timing_stages_recorded": len(timing) == 4
        and all(
            _finite(value, f"timing_breakdown_s.{name}") >= 0.0
            for name, value in timing.items()
        ),
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "build123d_loft_example_source_replay_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "cad_handoff_ready": not issues,
        "mesh_ready": False,
        "solver_ready": False,
        "checks": checks,
        "issues": issues,
        "source": "examples/loft.py",
        "notes": [
            "The source example contains eleven loft profiles because slice_count=10 includes both endpoints.",
            "Stub only display integration; preserve the geometry source and official volume assertion.",
            "CAD handoff remains distinct from mesh and physics handoff.",
        ],
    }
