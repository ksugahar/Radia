"""Solver-neutral identities for an ideal two-winding transformer."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any


def _finite(value: object, name: str, *, positive: bool = False) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not math.isfinite(parsed) or (positive and parsed <= 0.0):
        requirement = "positive and finite" if positive else "finite"
        raise ValueError(f"{name} must be {requirement}")
    return parsed


def _complex_pair(value: object, name: str) -> complex:
    if (
        not isinstance(value, Sequence)
        or isinstance(value, (str, bytes))
        or len(value) != 2
    ):
        raise ValueError(f"{name} must contain [real, imag]")
    return complex(_finite(value[0], name), _finite(value[1], name))


def _relative_error(actual: complex | float, expected: complex | float) -> float:
    return float(abs(actual - expected) / max(abs(expected), 1.0e-300))


def _phasor_replay_generations_ok(positive: Mapping[str, object]) -> bool:
    evidence = positive.get("phasor_replay_evidence")
    if evidence is None:
        return True
    if (
        not isinstance(evidence, Sequence)
        or isinstance(evidence, (str, bytes))
        or len(evidence) < 2
    ):
        return False
    rows = [row for row in evidence if isinstance(row, Mapping)]
    return len(rows) == len(evidence) and all(
        bool(row.get("fit_id"))
        and bool(row.get("raw_generation_id"))
        and row.get("trace_group_generation_id") == row.get("raw_generation_id")
        and len(str(row.get("scalar_phasor_digest") or "")) == 64
        for row in rows
    )


def _fit_window_stays_in_one_segment(
    positive: Mapping[str, object], fit_start: float, fit_stop: float
) -> bool:
    contract = positive.get("fit_window_segment_contract")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    generation = str(contract.get("fit_run_generation_id") or "")
    sample_generations = contract.get("sample_run_generation_ids")
    restarts = contract.get("restart_discontinuities_s")
    if (
        not generation
        or not isinstance(sample_generations, Sequence)
        or isinstance(sample_generations, (str, bytes))
        or not sample_generations
        or not isinstance(restarts, Sequence)
        or isinstance(restarts, (str, bytes))
        or any(str(item) != generation for item in sample_generations)
    ):
        return False
    try:
        restart_times = [float(value) for value in restarts]
    except (TypeError, ValueError):
        return False
    return all(
        math.isfinite(value) and not (fit_start < value < fit_stop)
        for value in restart_times
    )


def _phasor_basis_contract_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get("phasor_basis_contract")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    bases = contract.get("quantity_basis")
    factors = contract.get("normalization_factor_to_rms")
    names = {
        "source_voltage_phasor_rms_v",
        "primary_voltage_phasor_rms_v",
        "secondary_voltage_phasor_rms_v",
        "source_delivery_current_phasor_rms_a",
        "primary_current_phasor_rms_a",
        "secondary_current_phasor_rms_a",
    }
    if not isinstance(bases, Mapping) or not isinstance(factors, Mapping):
        return False
    try:
        factors_are_rms = all(
            math.isclose(float(factors.get(name)), 1.0, rel_tol=0.0, abs_tol=1.0e-15)
            for name in names
        )
    except (TypeError, ValueError):
        return False
    return (
        set(bases) == names
        and set(factors) == names
        and all(bases.get(name) == "rms" for name in names)
        and factors_are_rms
        and contract.get("complex_power_formula")
        == "rms_voltage_times_conjugate_rms_current"
    )


def _phase_unwrap_contract_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get("phase_unwrap_contract")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    try:
        period = float(contract.get("branch_period_rad"))
    except (TypeError, ValueError):
        return False
    fitted = str(contract.get("fit_branch_sign_digest") or "")
    replayed = str(contract.get("replay_branch_sign_digest") or "")
    return (
        contract.get("phase_unit") == "radian"
        and math.isclose(period, 2.0 * math.pi, rel_tol=0.0, abs_tol=1.0e-12)
        and contract.get("unwrap_convention") == "continuous_signed_phase"
        and contract.get("reference_trace") == "source_voltage_phasor_rms_v"
        and len(fitted) == 64
        and replayed == fitted
    )


def _power_sign_convention_ok(positive: Mapping[str, object]) -> bool:
    contract = positive.get("power_sign_convention")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    return (
        contract.get("source_power_role") == "delivered_positive"
        and contract.get("passive_power_role") == "absorbed_positive"
        and contract.get("source_current_reference")
        == "leaving_positive_terminal"
        and contract.get("passive_current_reference")
        == "entering_positive_terminal"
        and contract.get("balance_equation")
        == "source_delivered_equals_passive_absorbed"
        and contract.get("recorded_sign_transform") == "none"
    )


def _ac_frequency_interpolation_contract_ok(
    positive: Mapping[str, object],
) -> bool:
    contract = positive.get("ac_frequency_interpolation_contract")
    if contract is None:
        return True
    if not isinstance(contract, Mapping):
        return False
    source_coordinate = str(contract.get("source_frequency_coordinate") or "")
    target_coordinate = str(contract.get("target_frequency_coordinate") or "")
    interpolation_coordinate = str(contract.get("interpolation_coordinate") or "")
    source_generation = str(contract.get("source_grid_generation") or "")
    trace_generation = str(
        contract.get("interpolated_trace_source_grid_generation") or ""
    )
    return (
        source_coordinate in {"linear_hz", "log10_hz", "natural_log_hz"}
        and source_coordinate == target_coordinate == interpolation_coordinate
        and contract.get("frequency_unit") == "Hz"
        and bool(source_generation)
        and trace_generation == source_generation
    )


def ideal_transformer_identity_gate(summary: Mapping[str, object]) -> dict[str, Any]:
    """Gate turns ratio, reflected impedance, network closure, power, and replay."""
    if not isinstance(summary, Mapping):
        raise ValueError("summary must be an object")
    model = summary.get("model_contract")
    metrics = summary.get("metrics")
    timing = summary.get("timing_breakdown_s")
    if not isinstance(model, Mapping) or not isinstance(metrics, Mapping):
        raise ValueError("model_contract and metrics must be objects")
    positive = metrics.get("positive")
    if not isinstance(positive, Mapping):
        raise ValueError("metrics.positive must be an object")

    turns_ratio = _finite(
        model.get("turns_ratio_primary_to_secondary"),
        "turns_ratio_primary_to_secondary",
        positive=True,
    )
    source_offset = _finite(model.get("source_offset_v"), "source_offset_v")
    source_peak = _finite(model.get("source_peak_v"), "source_peak_v", positive=True)
    frequency = _finite(model.get("frequency_hz"), "frequency_hz", positive=True)
    series_resistance = _finite(
        model.get("series_resistance_ohm"), "series_resistance_ohm", positive=True
    )
    load_resistance = _finite(
        model.get("load_resistance_ohm"), "load_resistance_ohm", positive=True
    )

    source_rms = source_peak / math.sqrt(2.0)
    reflected_load = turns_ratio**2 * load_resistance
    primary_current_expected = source_rms / (series_resistance + reflected_load)
    primary_voltage_expected = primary_current_expected * reflected_load
    secondary_voltage_expected = primary_voltage_expected / turns_ratio
    secondary_current_expected = secondary_voltage_expected / load_resistance
    load_power_expected = secondary_voltage_expected**2 / load_resistance
    series_loss_expected = primary_current_expected**2 * series_resistance
    derived = {
        "expected_source_voltage_rms_v": source_rms,
        "expected_reflected_load_ohm": reflected_load,
        "expected_primary_voltage_rms_v": primary_voltage_expected,
        "expected_secondary_voltage_rms_v": secondary_voltage_expected,
        "expected_primary_current_rms_a": primary_current_expected,
        "expected_secondary_current_rms_a": secondary_current_expected,
        "expected_load_power_w": load_power_expected,
        "expected_series_loss_w": series_loss_expected,
    }
    declared_expectation_error = max(
        _relative_error(_finite(model.get(name), name, positive=True), value)
        for name, value in derived.items()
    )

    v_source = _complex_pair(
        positive.get("source_voltage_phasor_rms_v"),
        "source_voltage_phasor_rms_v",
    )
    v_primary = _complex_pair(
        positive.get("primary_voltage_phasor_rms_v"),
        "primary_voltage_phasor_rms_v",
    )
    v_secondary = _complex_pair(
        positive.get("secondary_voltage_phasor_rms_v"),
        "secondary_voltage_phasor_rms_v",
    )
    i_source = _complex_pair(
        positive.get("source_delivery_current_phasor_rms_a"),
        "source_delivery_current_phasor_rms_a",
    )
    i_primary = _complex_pair(
        positive.get("primary_current_phasor_rms_a"),
        "primary_current_phasor_rms_a",
    )
    i_secondary = _complex_pair(
        positive.get("secondary_current_phasor_rms_a"),
        "secondary_current_phasor_rms_a",
    )
    if min(abs(v_source), abs(v_primary), abs(v_secondary), abs(i_primary), abs(i_secondary)) <= 1.0e-300:
        raise ValueError("transformer phasors must be nonzero")

    transformer_input_power = v_primary * i_primary.conjugate()
    transformer_output_power = v_secondary * i_secondary.conjugate()
    source_delivery_power = v_source * i_source.conjugate()
    series_loss = series_resistance * abs(i_primary) ** 2
    recomputed = {
        "source_voltage_relative_error": _relative_error(abs(v_source), source_rms),
        "primary_voltage_relative_error": _relative_error(
            abs(v_primary), primary_voltage_expected
        ),
        "secondary_voltage_relative_error": _relative_error(
            abs(v_secondary), secondary_voltage_expected
        ),
        "primary_current_relative_error": _relative_error(
            abs(i_primary), primary_current_expected
        ),
        "secondary_current_relative_error": _relative_error(
            abs(i_secondary), secondary_current_expected
        ),
        "voltage_turns_identity_relative_error": _relative_error(
            v_secondary, -v_primary / turns_ratio
        ),
        "current_turns_identity_relative_error": _relative_error(
            i_secondary, -turns_ratio * i_primary
        ),
        "secondary_load_ohm_law_relative_error": _relative_error(
            v_secondary, load_resistance * i_secondary
        ),
        "reflected_load_relative_error": _relative_error(
            v_primary / i_primary, reflected_load
        ),
        "source_series_kvl_relative_error": _relative_error(
            v_source, v_primary + series_resistance * i_primary
        ),
        "source_primary_current_relative_error": _relative_error(i_source, i_primary),
        "transformer_complex_power_relative_error": _relative_error(
            transformer_input_power, transformer_output_power
        ),
        "source_power_balance_relative_error": _relative_error(
            source_delivery_power, transformer_output_power + series_loss
        ),
    }
    reported_error_drift = max(
        abs(_finite(positive.get(name), name) - value)
        for name, value in recomputed.items()
    )

    reported_powers = {
        "transformer_input_complex_power_va": _complex_pair(
            positive.get("transformer_input_complex_power_va"),
            "transformer_input_complex_power_va",
        ),
        "transformer_output_complex_power_va": _complex_pair(
            positive.get("transformer_output_complex_power_va"),
            "transformer_output_complex_power_va",
        ),
        "source_delivery_complex_power_va": _complex_pair(
            positive.get("source_delivery_complex_power_va"),
            "source_delivery_complex_power_va",
        ),
    }
    power_report_drift = max(
        _relative_error(
            reported_powers["transformer_input_complex_power_va"],
            transformer_input_power,
        ),
        _relative_error(
            reported_powers["transformer_output_complex_power_va"],
            transformer_output_power,
        ),
        _relative_error(
            reported_powers["source_delivery_complex_power_va"],
            source_delivery_power,
        ),
        _relative_error(_finite(positive.get("series_loss_w"), "series_loss_w"), series_loss),
    )

    instantaneous_transformer_error = _finite(
        positive.get("instantaneous_transformer_power_relative_error"),
        "instantaneous_transformer_power_relative_error",
    )
    instantaneous_source_error = _finite(
        positive.get("instantaneous_source_balance_relative_error"),
        "instantaneous_source_balance_relative_error",
    )
    fit_residual = _finite(
        positive.get("maximum_phasor_fit_relative_residual"),
        "maximum_phasor_fit_relative_residual",
    )
    fit_frequency_matches_source = True
    if "phasor_fit_frequency_hz" in positive:
        fit_frequency = _finite(
            positive.get("phasor_fit_frequency_hz"),
            "phasor_fit_frequency_hz",
            positive=True,
        )
        fit_frequency_matches_source = (
            abs(fit_frequency - frequency) <= 1.0e-12 * max(frequency, 1.0)
        )

    expected_current_roles = {
        "source_delivery_current_phasor_rms_a": "source_delivery_into_network",
        "primary_current_phasor_rms_a": "transformer_primary_absorption",
        "secondary_current_phasor_rms_a": "transformer_secondary_delivery_to_load",
    }
    role_contract = positive.get("current_role_contract")
    observed_roles = positive.get("phasor_current_roles")
    role_evidence_present = role_contract is not None or observed_roles is not None
    current_roles_match_contract = not role_evidence_present or (
        isinstance(role_contract, Mapping)
        and isinstance(observed_roles, Mapping)
        and dict(role_contract) == expected_current_roles
        and dict(observed_roles) == expected_current_roles
    )
    point_count = int(_finite(positive.get("point_count"), "point_count", positive=True))
    fit_point_count = int(
        _finite(positive.get("fit_point_count"), "fit_point_count", positive=True)
    )
    fit_start = _finite(positive.get("fit_window_start_s"), "fit_window_start_s")
    fit_stop = _finite(
        positive.get("fit_window_stop_s"), "fit_window_stop_s", positive=True
    )
    replay_error = _finite(
        metrics.get("maximum_phasor_replay_relative_error"),
        "maximum_phasor_replay_relative_error",
    )

    timing_ok = False
    if isinstance(timing, Mapping) and len(timing) == 4:
        try:
            timing_ok = all(_finite(value, "timing") >= 0.0 for value in timing.values())
        except ValueError:
            timing_ok = False

    checks = {
        "ideal_two_winding_transformer_contract": model.get("topology")
        == "ideal_two_winding_transformer_with_series_source_resistance"
        and abs(source_offset) <= 1.0e-15,
        "declared_analytic_expectations_are_recomputed": declared_expectation_error
        <= 1.0e-12,
        "source_and_analytic_rms_values_close": max(
            recomputed["source_voltage_relative_error"],
            recomputed["primary_voltage_relative_error"],
            recomputed["secondary_voltage_relative_error"],
            recomputed["primary_current_relative_error"],
            recomputed["secondary_current_relative_error"],
        )
        <= 2.0e-5,
        "voltage_and_current_turns_identities_close": max(
            recomputed["voltage_turns_identity_relative_error"],
            recomputed["current_turns_identity_relative_error"],
        )
        <= 2.0e-6,
        "load_reflection_and_series_source_network_close": max(
            recomputed["secondary_load_ohm_law_relative_error"],
            recomputed["reflected_load_relative_error"],
            recomputed["source_series_kvl_relative_error"],
            recomputed["source_primary_current_relative_error"],
        )
        <= 2.0e-6,
        "complex_and_instantaneous_power_are_conserved": max(
            recomputed["transformer_complex_power_relative_error"],
            recomputed["source_power_balance_relative_error"],
            instantaneous_transformer_error,
            instantaneous_source_error,
        )
        <= 2.0e-6,
        "reported_errors_and_complex_powers_match_recomputation": max(
            reported_error_drift, power_report_drift
        )
        <= 1.0e-12,
        "two_period_phasor_fit_is_dense_and_accurate": point_count >= 1000
        and fit_point_count >= 500
        and fit_stop - fit_start >= 1.9 / frequency
        and fit_residual <= 2.0e-5,
        "phasor_fit_frequency_matches_source_contract": fit_frequency_matches_source,
        "current_phasor_roles_match_terminal_contract": current_roles_match_contract,
        "positive_phasor_replay_is_deterministic": replay_error <= 1.0e-12,
        "phasor_replays_bind_each_fit_to_one_raw_generation": (
            _phasor_replay_generations_ok(positive)
        ),
        "phasor_fit_window_does_not_cross_a_restart_segment": (
            _fit_window_stays_in_one_segment(positive, fit_start, fit_stop)
        ),
        "all_phasors_share_one_rms_normalization_basis": (
            _phasor_basis_contract_ok(positive)
        ),
        "phase_replay_preserves_unwrap_branch_orientation": (
            _phase_unwrap_contract_ok(positive)
        ),
        "source_and_passive_power_use_recorded_sign_conventions": (
            _power_sign_convention_ok(positive)
        ),
        "ac_traces_share_frequency_interpolation_coordinate": (
            _ac_frequency_interpolation_contract_ok(positive)
        ),
        "exactly_four_timing_stages": timing_ok,
    }
    return {
        "schema": "radia-spice-lab.ideal-transformer-identity.v1",
        "policy": "ideal_transformer_identity_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            **derived,
            "maximum_declared_expectation_relative_error": declared_expectation_error,
            "maximum_recomputed_identity_relative_error": max(recomputed.values()),
            "maximum_reported_error_absolute_drift": reported_error_drift,
            "maximum_reported_power_relative_drift": power_report_drift,
            "transformer_input_active_power_w": transformer_input_power.real,
            "transformer_output_active_power_w": transformer_output_power.real,
            "source_delivery_active_power_w": source_delivery_power.real,
        },
        "notes": [
            "For the chosen dot/sign convention, V_secondary=-V_primary/N and I_secondary=-N*I_primary.",
            "The load reflected to the primary is N^2*R_load; include source resistance in the analytic RMS reference.",
            "Turns-ratio agreement alone is insufficient: require instantaneous and complex-power conservation plus deterministic replay.",
            "When phasor frequency or terminal-current role metadata is supplied, bind it to the source and sign contract before accepting scalar power closure.",
            "A matching scalar phasor is not replay evidence when its fitted traces mix RAW generations or cross a transient restart.",
            "Peak, RMS, and phasor quantities must declare one normalization basis before turns-ratio or complex-power closure is compared.",
            "Phase replay must retain the fitted unwrap branch orientation; equal magnitudes do not resolve a sign-changing branch alias.",
        ],
    }
