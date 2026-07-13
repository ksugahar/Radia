"""Solver-neutral transient conductor replay and constitutive identities."""

from __future__ import annotations

import math
from typing import Any


_SERIES_KEYS = (
    "times_s",
    "current_a",
    "joule_loss_w",
    "circuit_power_w",
    "flux_linkage_wb",
    "resistance_ohm",
    "inductance_h",
)
_OBSERVABLE_KEYS = _SERIES_KEYS[1:]


def _series(row: dict[str, Any], key: str) -> list[float]:
    value = row.get(key)
    if not isinstance(value, list):
        raise ValueError(f"{key} must be a list")
    parsed = [float(item) for item in value]
    if not parsed or not all(math.isfinite(item) for item in parsed):
        raise ValueError(f"{key} must contain finite values")
    return parsed


def _relative_error(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-30)


def transient_conductor_replay_identity_gate(
    summary: dict[str, Any],
    *,
    identity_rtol: float = 1.0e-10,
    equivalence_rtol: float = 1.0e-12,
) -> dict[str, Any]:
    """Gate transient conductor waveforms by replay and constitutive closure.

    Each replay must contain two independently defined but physically equivalent
    variants. The gate checks the full sampled histories, not only extrema.
    """

    if not isinstance(summary, dict):
        raise ValueError("summary must be an object")
    original_steps = int(summary.get("original_step_count", 0))
    bounded_steps = int(summary.get("bounded_step_count", 0))
    time_step = float(summary.get("time_step_s", math.nan))
    replays = summary.get("replays")
    timing = summary.get("timing_breakdown_s")
    if not isinstance(replays, list) or len(replays) < 2:
        raise ValueError("at least two replays are required")

    parsed_replays: list[list[dict[str, Any]]] = []
    waveform_lengths: list[int] = []
    time_monotonic = True
    passive = True
    joule_errors: list[float] = []
    flux_errors: list[float] = []
    mesh_pairs: list[tuple[int, int]] = []
    for replay in replays:
        variants = replay.get("variants") if isinstance(replay, dict) else None
        if not isinstance(variants, list) or len(variants) != 2:
            raise ValueError("each replay must contain exactly two variants")
        parsed_variants = []
        for variant in variants:
            if not isinstance(variant, dict):
                raise ValueError("variant entries must be objects")
            series = {key: _series(variant, key) for key in _SERIES_KEYS}
            lengths = {len(values) for values in series.values()}
            if len(lengths) != 1:
                raise ValueError("all waveform series must have equal length")
            count = lengths.pop()
            waveform_lengths.append(count)
            time_monotonic &= all(
                right > left
                for left, right in zip(series["times_s"], series["times_s"][1:])
            )
            for current, joule, resistance, flux, inductance in zip(
                series["current_a"],
                series["joule_loss_w"],
                series["resistance_ohm"],
                series["flux_linkage_wb"],
                series["inductance_h"],
                strict=True,
            ):
                if abs(current) > 1.0e-12:
                    passive &= resistance > 0.0 and inductance > 0.0 and joule >= 0.0
                    joule_errors.append(_relative_error(joule, current * current * resistance))
                    flux_errors.append(_relative_error(flux, current * inductance))
            mesh_pairs.append(
                (int(variant.get("mesh_elements", 0)), int(variant.get("mesh_vertices", 0)))
            )
            parsed_variants.append({"label": str(variant.get("label", "")), **series})
        parsed_replays.append(parsed_variants)

    variant_errors: list[float] = []
    for variants in parsed_replays:
        left, right = variants
        for key in _OBSERVABLE_KEYS:
            variant_errors.extend(
                _relative_error(a, b)
                for a, b in zip(left[key], right[key], strict=True)
            )

    replay_errors: list[float] = []
    reference = parsed_replays[0]
    for replay in parsed_replays[1:]:
        for left, right in zip(reference, replay, strict=True):
            for key in _OBSERVABLE_KEYS:
                replay_errors.extend(
                    _relative_error(a, b)
                    for a, b in zip(left[key], right[key], strict=True)
                )

    timing_ok = False
    if isinstance(timing, dict) and len(timing) == 4:
        try:
            timing_ok = all(
                math.isfinite(float(value)) and float(value) >= 0.0
                for value in timing.values()
            )
        except (TypeError, ValueError):
            timing_ok = False

    maximum_joule_error = max(joule_errors, default=math.inf)
    maximum_flux_error = max(flux_errors, default=math.inf)
    maximum_variant_error = max(variant_errors, default=math.inf)
    maximum_replay_error = max(replay_errors, default=math.inf)
    checks = {
        "bounded_replay_preserves_positive_time_step": original_steps > bounded_steps >= 3
        and math.isfinite(time_step)
        and time_step > 0.0,
        "complete_waveform_lengths_match_bounded_steps": bool(waveform_lengths)
        and all(length == bounded_steps for length in waveform_lengths),
        "time_axes_are_strictly_increasing": time_monotonic,
        "passive_rl_observables": passive and bool(joule_errors),
        "joule_loss_closes_i_squared_r": maximum_joule_error <= float(identity_rtol),
        "flux_linkage_closes_l_times_i": maximum_flux_error <= float(identity_rtol),
        "equivalent_definitions_match_full_waveforms": maximum_variant_error
        <= float(equivalence_rtol),
        "independent_replays_match_full_waveforms": maximum_replay_error
        <= float(equivalence_rtol),
        "mesh_inventory_is_replay_invariant": bool(mesh_pairs)
        and mesh_pairs[0][0] > 0
        and mesh_pairs[0][1] > 0
        and len(set(mesh_pairs)) == 1,
        "exactly_four_timing_stages": timing_ok,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "transient_conductor_replay_identity_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            "replay_count": len(parsed_replays),
            "waveform_points": bounded_steps,
            "maximum_joule_i2r_relative_error": maximum_joule_error,
            "maximum_flux_li_relative_error": maximum_flux_error,
            "maximum_equivalent_definition_relative_error": maximum_variant_error,
            "maximum_independent_replay_relative_error": maximum_replay_error,
        },
        "lesson": (
            "For a transient conductor, validate complete sampled histories with "
            "Joule=I^2 R and flux-linkage=L I, then require independent replay and "
            "equivalent material-definition closure. Record staging/load separately "
            "from native solve time because model preparation can dominate runtime."
        ),
    }
