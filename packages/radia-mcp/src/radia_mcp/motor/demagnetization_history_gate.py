"""Solver-independent history gate for irreversible permanent-magnet demagnetization."""
from __future__ import annotations

import json
import math


def _demagnetization_family_gate(
    summary: dict,
    state_tolerance: float,
    minimum_damage_fraction: float,
) -> dict:
    """Validate elementwise state memory across cases and load blocks."""

    cases = summary.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("cases must be a nonempty list")

    case_checks: dict[str, dict[str, bool]] = {}
    case_metrics: dict[str, dict[str, float | int]] = {}
    case_ids: list[str] = []
    total_blocks = 0
    damage_blocks = 0
    no_additional_damage_blocks = 0

    for case_position, case in enumerate(cases):
        if not isinstance(case, dict):
            raise ValueError(f"cases[{case_position}] must be an object")
        case_id = str(case.get("case_id", "")).strip()
        if not case_id:
            raise ValueError(f"cases[{case_position}].case_id must be nonempty")
        case_ids.append(case_id)

        raw_states = case.get("states")
        if not isinstance(raw_states, list) or len(raw_states) < 3:
            raise ValueError(f"cases[{case_position}].states must contain at least three records")

        state_ids: list[str] = []
        step_ids: list[int] = []
        state_vectors: dict[str, list[float]] = {}
        field_vectors: dict[str, list[float]] = {}
        expected_elements: tuple[int, ...] | None = None
        expected_materials: tuple[int, ...] | None = None
        element_identity_stable = True
        state_values_finite = True
        field_values_finite = True

        for state_position, state in enumerate(raw_states):
            if not isinstance(state, dict):
                raise ValueError(
                    f"cases[{case_position}].states[{state_position}] must be an object"
                )
            state_id = str(state.get("state_id", "")).strip()
            if not state_id:
                raise ValueError(
                    f"cases[{case_position}].states[{state_position}].state_id must be nonempty"
                )
            state_ids.append(state_id)
            step_ids.append(int(state.get("step", state_position)))

            raw_elements = state.get("element_ids")
            raw_materials = state.get("material_ids")
            raw_remanence = state.get("remanence_ratio")
            raw_field = state.get("flux_density_T")
            if not isinstance(raw_elements, list) or not raw_elements:
                raise ValueError("element_ids must be a nonempty list")
            if not isinstance(raw_materials, list) or len(raw_materials) != len(raw_elements):
                raise ValueError("material_ids must match element_ids")
            if not isinstance(raw_remanence, list) or len(raw_remanence) != len(raw_elements):
                raise ValueError("remanence_ratio must match element_ids")
            if not isinstance(raw_field, list) or len(raw_field) != len(raw_elements):
                raise ValueError("flux_density_T must match element_ids")

            elements = tuple(int(value) for value in raw_elements)
            materials = tuple(int(value) for value in raw_materials)
            remanence = [float(value) for value in raw_remanence]
            field = [float(value) for value in raw_field]
            if expected_elements is None:
                expected_elements = elements
                expected_materials = materials
            else:
                element_identity_stable &= elements == expected_elements
                element_identity_stable &= materials == expected_materials
            state_values_finite &= all(math.isfinite(value) for value in remanence)
            field_values_finite &= all(math.isfinite(value) and value >= 0.0 for value in field)
            state_vectors[state_id] = remanence
            field_vectors[state_id] = field

        state_never_recovers = all(
            later[element] <= earlier[element] + state_tolerance
            for earlier, later in zip(
                (state_vectors[state_id] for state_id in state_ids),
                (state_vectors[state_id] for state_id in state_ids[1:]),
            )
            for element in range(len(earlier))
        )

        raw_blocks = case.get("history_blocks")
        if not isinstance(raw_blocks, list) or not raw_blocks:
            raise ValueError(f"cases[{case_position}].history_blocks must be a nonempty list")
        block_references_valid = True
        stress_never_increases_state = True
        unload_preserves_stress_state = True
        block_expectations_explicit = True
        expected_damage_behavior = True
        field_changes_on_unload = True
        maximum_damage = 0.0
        maximum_unload_drift = 0.0
        state_order = {state_id: position for position, state_id in enumerate(state_ids)}

        for block in raw_blocks:
            if not isinstance(block, dict):
                raise ValueError("history_blocks entries must be objects")
            pre_id = str(block.get("pre_state", ""))
            stress_id = str(block.get("stress_state", ""))
            unload_id = str(block.get("unloaded_state", ""))
            valid = (
                pre_id in state_vectors
                and stress_id in state_vectors
                and unload_id in state_vectors
                and state_order[pre_id] < state_order[stress_id] < state_order[unload_id]
            )
            block_references_valid &= valid
            if not valid:
                continue

            pre = state_vectors[pre_id]
            stressed = state_vectors[stress_id]
            unloaded = state_vectors[unload_id]
            damage = [before - after for before, after in zip(pre, stressed)]
            unload_drift = [after - during for during, after in zip(stressed, unloaded)]
            block_damage = max(damage, default=0.0)
            block_unload_drift = max((abs(value) for value in unload_drift), default=0.0)
            maximum_damage = max(maximum_damage, block_damage)
            maximum_unload_drift = max(maximum_unload_drift, block_unload_drift)
            stress_never_increases_state &= min(damage, default=0.0) >= -state_tolerance
            unload_preserves_stress_state &= block_unload_drift <= state_tolerance

            expectation = block.get("expect_additional_demagnetization")
            block_expectations_explicit &= isinstance(expectation, bool)
            expects_damage = expectation is True
            if expects_damage:
                damage_blocks += 1
                expected_damage_behavior &= block_damage > minimum_damage_fraction
            else:
                no_additional_damage_blocks += 1
                expected_damage_behavior &= max((abs(value) for value in damage), default=0.0) <= state_tolerance
            total_blocks += 1

            stress_field = field_vectors[stress_id]
            unloaded_field = field_vectors[unload_id]
            field_changes_on_unload &= max(
                (abs(after - during) for during, after in zip(stress_field, unloaded_field)),
                default=0.0,
            ) > state_tolerance

        replay_count = int(case.get("replay_count", 0))
        replay_max_abs = float(case.get("replay_max_abs", math.inf))
        bounded = all(
            -state_tolerance <= value <= 1.0 + state_tolerance
            for vector in state_vectors.values()
            for value in vector
        )
        checks = {
            "state_ids_unique": len(state_ids) == len(set(state_ids)),
            "step_ids_strictly_increase": all(a < b for a, b in zip(step_ids, step_ids[1:])),
            "element_and_material_identity_stable": element_identity_stable,
            "remanence_state_finite_and_bounded": state_values_finite and bounded,
            "instantaneous_field_finite_and_nonnegative": field_values_finite,
            "state_never_spontaneously_recovers": state_never_recovers,
            "history_block_references_ordered": block_references_valid,
            "stress_never_increases_remanence_state": stress_never_increases_state,
            "unload_preserves_stressed_remanence_state": unload_preserves_stress_state,
            "damage_expectation_is_explicit_for_each_block": block_expectations_explicit,
            "damage_expectation_matches_each_block": expected_damage_behavior,
            "instantaneous_field_changes_without_state_healing": field_changes_on_unload,
            "two_or_more_replays_recorded": replay_count >= 2,
            "replays_match_within_state_tolerance": (
                math.isfinite(replay_max_abs) and 0.0 <= replay_max_abs <= state_tolerance
            ),
        }
        case_checks[case_id] = checks
        case_metrics[case_id] = {
            "state_count": len(state_ids),
            "partition_count": len(expected_elements or ()),
            "history_block_count": len(raw_blocks),
            "maximum_damage_fraction": maximum_damage,
            "maximum_unload_state_drift": maximum_unload_drift,
            "replay_max_abs": replay_max_abs,
        }

    checks = {
        "schema_identifies_history_family": summary.get("schema")
        == "permanent-magnet-demagnetization-history/v1",
        "result_authority_recorded": summary.get("result_authority") == ".mao",
        "state_variable_is_elementwise_remanence": summary.get("state_variable")
        == "elementwise_remanence_ratio",
        "field_observable_is_instantaneous_flux_density": summary.get(
            "instantaneous_field_observable"
        )
        == "elementwise_flux_density_magnitude_T",
        "case_ids_unique": len(case_ids) == len(set(case_ids)),
        "all_case_histories_pass": all(
            ok for per_case in case_checks.values() for ok in per_case.values()
        ),
    }
    return {
        "policy": "permanent_magnet_demagnetization_history_gate_v2",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "case_checks": case_checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "case_count": len(cases),
            "history_block_count": total_blocks,
            "damage_block_count": damage_blocks,
            "no_additional_damage_block_count": no_additional_damage_blocks,
            "cases": case_metrics,
        },
        "notes": [
            "instantaneous flux density may recover after unloading while the remanence state must not heal",
            "a preconditioned history block may validly show no additional damage when its state remains locked",
            "replay evidence is part of the state-history claim",
        ],
    }


def permanent_magnet_demagnetization_history_gate(
    summary_json: str,
    state_tolerance: float = 1.0e-9,
    minimum_damage_fraction: float = 1.0e-3,
) -> dict:
    """Check that a load excursion leaves a bounded, persistent magnet state."""

    if state_tolerance < 0.0 or minimum_damage_fraction <= 0.0:
        raise ValueError("state_tolerance must be nonnegative and minimum_damage_fraction positive")
    try:
        summary = json.loads(summary_json)
    except json.JSONDecodeError as exc:
        raise ValueError(f"summary_json must be valid JSON: {exc.msg}") from exc
    if not isinstance(summary, dict):
        raise ValueError("summary_json must decode to an object")
    if "cases" in summary and "steps" not in summary:
        return _demagnetization_family_gate(
            summary,
            state_tolerance,
            minimum_damage_fraction,
        )

    steps = summary.get("steps")
    if not isinstance(steps, list) or len(steps) < 3:
        raise ValueError("steps must contain at least initial, stress, and recovery records")
    stress_index = int(summary.get("stress_step_index", -1))
    recovery_index = int(summary.get("recovery_step_index", -1))
    if not (0 < stress_index < recovery_index < len(steps)):
        raise ValueError("stress_step_index and recovery_step_index must identify an ordered load history")

    vectors: list[list[float]] = []
    peak_fields: list[float] = []
    step_ids: list[int] = []
    expected_size: int | None = None
    for position, step in enumerate(steps):
        if not isinstance(step, dict):
            raise ValueError(f"steps[{position}] must be an object")
        step_ids.append(int(step.get("step", position)))
        raw_state = step.get("magnet_state_fraction")
        if not isinstance(raw_state, list) or not raw_state:
            raise ValueError(f"steps[{position}].magnet_state_fraction must be a nonempty list")
        state = [float(value) for value in raw_state]
        if not all(math.isfinite(value) for value in state):
            raise ValueError("magnet state fractions must be finite")
        if expected_size is None:
            expected_size = len(state)
        elif len(state) != expected_size:
            raise ValueError("all state vectors must use the same magnet partition")
        vectors.append(state)
        peak = float(step.get("peak_flux_density_T", 0.0))
        if not math.isfinite(peak) or peak < 0.0:
            raise ValueError("peak_flux_density_T must be finite and nonnegative")
        peak_fields.append(peak)

    initial = vectors[0]
    stressed = vectors[stress_index]
    recovered = vectors[recovery_index]
    state_increases = [
        later[element] - earlier[element]
        for earlier, later in zip(vectors, vectors[1:])
        for element in range(len(initial))
    ]
    damage = [before - after for before, after in zip(initial, stressed)]
    recovery_drift = [after - during for during, after in zip(stressed, recovered)]
    irreversible_elements = sum(value > minimum_damage_fraction for value in damage)
    stress_to_recovery_drift = max(
        abs(vectors[position][element] - stressed[element])
        for position in range(stress_index + 1, recovery_index + 1)
        for element in range(len(initial))
    )

    checks = {
        "state_unit_is_fraction": summary.get("state_unit") == "fraction_of_reference_remanence",
        "field_unit_is_T": summary.get("field_unit") == "T",
        "step_ids_strictly_increase": all(a < b for a, b in zip(step_ids, step_ids[1:])),
        "state_fractions_bounded": all(
            -state_tolerance <= value <= 1.0 + state_tolerance
            for vector in vectors
            for value in vector
        ),
        "state_never_spontaneously_recovers": max(state_increases, default=0.0) <= state_tolerance,
        "stress_causes_resolved_damage": irreversible_elements > 0,
        "recovery_retains_stressed_state": max((abs(value) for value in recovery_drift), default=0.0) <= state_tolerance,
        "stress_to_recovery_state_is_history_locked": stress_to_recovery_drift <= state_tolerance,
        "field_changes_between_stress_and_recovery": (
            max(peak_fields[stress_index : recovery_index + 1])
            - min(peak_fields[stress_index : recovery_index + 1])
            > 0.0
        ),
    }
    return {
        "policy": "permanent_magnet_demagnetization_history_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "partition_count": len(initial),
            "irreversibly_damaged_partition_count": irreversible_elements,
            "maximum_damage_fraction": max(damage),
            "minimum_stressed_state_fraction": min(stressed),
            "maximum_state_increase": max(state_increases, default=0.0),
            "maximum_recovery_state_drift": max((abs(value) for value in recovery_drift), default=0.0),
            "maximum_stress_to_recovery_state_drift": stress_to_recovery_drift,
            "stress_to_recovery_peak_flux_excursion_T": (
                max(peak_fields[stress_index : recovery_index + 1])
                - min(peak_fields[stress_index : recovery_index + 1])
            ),
        },
        "notes": [
            "demagnetization is a path-dependent state check, not a final-field-only threshold",
            "the state vector must retain the worst excursion even when the magnetic field changes after unloading",
        ],
    }
