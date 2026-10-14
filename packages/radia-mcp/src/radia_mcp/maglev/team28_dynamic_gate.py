"""Scope and evidence gate for the TEAM 28 cycle-averaged motion model."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence


_CYCLE_AVERAGED_SCOPE = "cycle_averaged_mechanical_motion"
_FULL_TRANSIENT_SCOPE = "full_electromagnetic_transient"
_ABSOLUTE_WINDOWS_PATH = re.compile(r"^[A-Za-z]:[\\/]")


def _mapping(value: object) -> Mapping[str, object]:
    return value if isinstance(value, Mapping) else {}


def _finite_number(value: object) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _portable_relative_path(value: object) -> bool:
    if not isinstance(value, str) or not value.strip():
        return False
    path = value.strip()
    return not (
        _ABSOLUTE_WINDOWS_PATH.match(path)
        or path.startswith("/")
        or path.startswith("\\\\")
        or ".." in path.replace("\\", "/").split("/")
    )


def _all_portable_paths(value: object) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and bool(value)
        and all(_portable_relative_path(item) for item in value)
    )


def team28_cycle_averaged_motion_gate(
    summary: Mapping[str, object],
    *,
    claim_scope: str = _CYCLE_AVERAGED_SCOPE,
    expected_frequency_hz: float = 50.0,
) -> dict[str, object]:
    """Validate a slow mechanical replay without calling it a full EM transient.

    The accepted model interpolates a validated fixed-frequency,
    cycle-averaged force-height family while Simulink advances the mechanical
    state.  It intentionally excludes carrier-resolved electromagnetic states,
    motion-induced EMF, and experimentally identified damping.
    """

    if not isinstance(summary, Mapping):
        raise TypeError("summary must be an object")
    if claim_scope not in {_CYCLE_AVERAGED_SCOPE, _FULL_TRANSIENT_SCOPE}:
        raise ValueError(
            "claim_scope must be 'cycle_averaged_mechanical_motion' or "
            "'full_electromagnetic_transient'"
        )
    expected_frequency = _finite_number(expected_frequency_hz)
    if expected_frequency is None or expected_frequency <= 0.0:
        raise ValueError("expected_frequency_hz must be finite and positive")

    model = _mapping(summary.get("model_contract"))
    checks_section = _mapping(summary.get("checks"))
    details = _mapping(checks_section.get("details"))
    errors = _mapping(summary.get("errors"))
    tolerances = _mapping(summary.get("tolerances"))
    observables = _mapping(summary.get("observables"))

    frequency = _finite_number(model.get("excitation_frequency_hz"))
    terminal_height_error = _finite_number(errors.get("terminal_height_abs_m"))
    terminal_force_error = _finite_number(errors.get("terminal_force_balance_abs_N"))
    terminal_speed = _finite_number(errors.get("terminal_speed_abs_m_per_s"))
    height_tolerance = _finite_number(tolerances.get("terminal_height_abs_m"))
    force_tolerance = _finite_number(tolerances.get("terminal_force_balance_abs_N"))
    speed_tolerance = _finite_number(tolerances.get("terminal_speed_abs_m_per_s"))
    terminal_lift = _finite_number(observables.get("terminal_upward_lift_N"))
    disk_weight = _finite_number(model.get("disk_weight_N"))
    sample_count = _finite_number(observables.get("sample_count"))
    snapshot_count = _finite_number(model.get("force_family_snapshot_count"))
    state_order = _finite_number(model.get("eddy_state_order"))

    required_details = (
        "all_outputs_are_finite",
        "initial_height_is_preserved",
        "trajectory_stays_inside_validated_height_family",
        "terminal_height_within_0p1_mm_of_equilibrium",
        "terminal_force_balance_below_0p02_N",
        "terminal_speed_below_1_mm_per_s",
        "coilbuilder_source_validation_is_retained",
        "common_rank_three_eddy_basis_is_retained",
    )
    terminal_errors_close = (
        terminal_height_error is not None
        and height_tolerance is not None
        and terminal_height_error <= height_tolerance
        and terminal_force_error is not None
        and force_tolerance is not None
        and terminal_force_error <= force_tolerance
        and terminal_speed is not None
        and speed_tolerance is not None
        and terminal_speed <= speed_tolerance
    )
    lift_balance_closes = (
        terminal_lift is not None
        and disk_weight is not None
        and force_tolerance is not None
        and abs(terminal_lift - disk_weight) <= force_tolerance
    )
    checks = {
        "passing_solver_run_artifact": (
            summary.get("schema") == "cae-ai-lab.solver-run.v1"
            and summary.get("pass") is True
            and checks_section.get("validation_passed") is True
        ),
        "fixed_frequency_cycle_averaged_backend": (
            summary.get("solver_backend")
            == "fixed-50Hz-cycle-average-lut-plus-mechanical-plant"
            and model.get("electromagnetic_model_class")
            == "fixed_frequency_cycle_averaged_force_height_lut"
            and model.get("height_coupling") == "quasi_steady_interpolation"
        ),
        "excitation_frequency_matches": (
            frequency is not None
            and math.isclose(frequency, expected_frequency, rel_tol=0.0, abs_tol=1.0e-12)
        ),
        "electromagnetic_transient_is_explicitly_excluded": (
            model.get("electromagnetic_state_transient_included") is False
        ),
        "motional_emf_is_explicitly_excluded": (
            model.get("motional_emf_included") is False
        ),
        "damping_is_not_claimed_as_identified_data": (
            model.get("damping_identified_from_measurement") is False
            and "not identified" in str(model.get("damping_provenance", "")).lower()
        ),
        "validated_force_family_is_retained": (
            details.get("coilbuilder_source_validation_is_retained") is True
            and details.get("common_rank_three_eddy_basis_is_retained") is True
            and snapshot_count == 25.0
            and state_order == 3.0
            and _portable_relative_path(model.get("force_family_source"))
        ),
        "all_dynamic_validation_details_hold": all(
            details.get(name) is True for name in required_details
        ),
        "terminal_mechanical_equilibrium_closes": terminal_errors_close
        and lift_balance_closes,
        "trajectory_has_multiple_finite_samples": (
            sample_count is not None and sample_count >= 2.0
        ),
        "artifact_paths_are_portable_and_relative": (
            _portable_relative_path(summary.get("source_artifact"))
            and _all_portable_paths(summary.get("result_files"))
        ),
        "claimed_scope_matches_artifact": claim_scope == _CYCLE_AVERAGED_SCOPE,
    }
    issues = [name for name, accepted in checks.items() if not accepted]
    return {
        "policy": "team28_cycle_averaged_motion_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "validated_scope": _CYCLE_AVERAGED_SCOPE,
        "requested_claim_scope": claim_scope,
        "checks": checks,
        "issues": issues,
        "metrics": {
            "excitation_frequency_hz": frequency,
            "force_family_snapshot_count": snapshot_count,
            "eddy_state_order": state_order,
            "sample_count": sample_count,
            "terminal_height_abs_m": terminal_height_error,
            "terminal_force_balance_abs_n": terminal_force_error,
            "terminal_speed_abs_m_per_s": terminal_speed,
        },
        "unsupported_claims": [
            _FULL_TRANSIENT_SCOPE,
            "carrier_resolved_electromagnetic_waveform",
            "motion_induced_emf",
            "experimentally_identified_damping",
        ],
        "next_action_for_full_transient": (
            "Advance position-dependent electromagnetic states together with the mechanical state, "
            "including the motion derivative terms of the R/L/P family; then validate current, "
            "force, energy, and motion against an independent transient reference."
        ),
        "lesson": (
            "A force-height LUT can validate slow mechanical settling at its fixed carrier frequency. "
            "It is not evidence of a carrier-resolved or motion-EMF electromagnetic transient."
        ),
    }
