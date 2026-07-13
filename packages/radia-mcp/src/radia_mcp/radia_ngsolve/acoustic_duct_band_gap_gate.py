"""Solver-neutral acoustic duct band-gap and confinement gate."""

from __future__ import annotations

import math
from typing import Any


def _finite(summary: dict[str, Any], name: str) -> float:
    value = float(summary.get(name, math.nan))
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def acoustic_duct_band_gap_gate(
    summary: dict[str, Any],
    *,
    max_empty_lattice_relative_error: float = 1.0e-3,
    max_empty_duct_transparency_error: float = 1.0e-2,
    min_gap_width: float = 1.0,
    min_passband_transmission: float = 0.1,
    max_gap_transmission: float = 1.0e-2,
    min_pass_to_gap_contrast: float = 50.0,
    max_free_space_insertion_loss_spread_db: float = 1.5,
    max_replay_relative_error: float = 1.0e-7,
) -> dict[str, Any]:
    """Separate a confined duct stop band from free-space attenuation."""

    if not isinstance(summary, dict):
        raise ValueError("summary must be an object")
    model = summary.get("model_contract")
    duct = summary.get("duct_result")
    free = summary.get("free_space_control")
    if not all(isinstance(value, dict) for value in (model, duct, free)):
        raise ValueError("model_contract, duct_result, and free_space_control must be objects")

    pitch = _finite(model, "cell_pitch")
    radius = _finite(model, "inclusion_radius")
    side = _finite(model, "duct_side")
    count = int(_finite(model, "finite_crystal_cell_count"))
    sweep_max = _finite(model, "maximum_wavenumber")
    cutoff = math.pi / side
    bragg = math.pi / pitch

    values = {
        "empty_lattice_relative_error": _finite(duct, "empty_lattice_relative_error"),
        "empty_duct_transparency_error": _finite(duct, "empty_duct_transparency_error"),
        "gap_low": _finite(duct, "gap_low"),
        "gap_high": _finite(duct, "gap_high"),
        "first_band_minimum": _finite(duct, "first_band_minimum"),
        "maximum_passband_transmission": _finite(duct, "maximum_passband_transmission"),
        "maximum_gap_transmission": _finite(duct, "maximum_gap_transmission"),
        "maximum_below_band_transmission": _finite(duct, "maximum_below_band_transmission"),
        "pass_to_gap_contrast": _finite(duct, "pass_to_gap_contrast"),
        "replay_relative_error": _finite(duct, "replay_relative_error"),
        "free_space_minimum_insertion_loss_db": _finite(free, "minimum_insertion_loss_db"),
        "free_space_maximum_insertion_loss_db": _finite(free, "maximum_insertion_loss_db"),
        "free_space_insertion_loss_spread_db": _finite(free, "insertion_loss_spread_db"),
    }
    wavenumbers = free.get("wavenumbers")
    if not isinstance(wavenumbers, list) or len(wavenumbers) < 3:
        raise ValueError("free_space_control.wavenumbers must contain at least three values")
    wavenumbers = [float(value) for value in wavenumbers]
    if not all(math.isfinite(value) and value > 0.0 for value in wavenumbers):
        raise ValueError("free-space wavenumbers must be finite and positive")

    timing = summary.get("timing_breakdown_s")
    timing_ok = False
    if isinstance(timing, dict) and len(timing) == 4:
        try:
            timing_ok = all(math.isfinite(float(value)) and float(value) >= 0.0 for value in timing.values())
        except (TypeError, ValueError):
            timing_ok = False

    checks = {
        "same_spherical_inclusion_family_is_compared": (
            model.get("inclusion") == "sound_soft_sphere"
            and model.get("confined_geometry") == "rigid_duct_periodic_cell"
            and model.get("free_space_geometry") == "finite_linear_chain"
            and model.get("same_inclusion_family") is True
            and radius > 0.0
            and count >= 3
        ),
        "duct_sweep_stays_below_first_transverse_cutoff": sweep_max < cutoff,
        "empty_lattice_matches_analytic_dispersion": (
            values["empty_lattice_relative_error"] <= max_empty_lattice_relative_error
        ),
        "empty_duct_is_transparent": (
            values["empty_duct_transparency_error"] <= max_empty_duct_transparency_error
        ),
        "bloch_gap_has_positive_ordered_width": (
            values["gap_low"] > values["first_band_minimum"]
            and values["gap_high"] - values["gap_low"] >= min_gap_width
        ),
        "finite_crystal_has_passband": (
            values["maximum_passband_transmission"] >= min_passband_transmission
        ),
        "finite_crystal_attenuates_in_gap": (
            0.0 <= values["maximum_gap_transmission"] <= max_gap_transmission
            and 0.0 <= values["maximum_below_band_transmission"] <= max_gap_transmission
        ),
        "pass_to_gap_contrast_is_resolved": (
            values["pass_to_gap_contrast"] >= min_pass_to_gap_contrast
        ),
        "free_space_control_brackets_bragg_wavenumber": (
            min(wavenumbers) < bragg < max(wavenumbers)
        ),
        "free_space_chain_has_attenuation_but_no_deep_stop_band": (
            values["free_space_minimum_insertion_loss_db"] > 0.0
            and values["free_space_maximum_insertion_loss_db"] < 5.0
            and values["free_space_insertion_loss_spread_db"]
            <= max_free_space_insertion_loss_spread_db
        ),
        "fresh_reference_replays_saved_observables": (
            values["replay_relative_error"] <= max_replay_relative_error
        ),
        "exactly_four_timing_stages": timing_ok,
    }
    issues = [name for name, ok in checks.items() if not ok]
    return {
        "policy": "acoustic_duct_band_gap_gate_v1",
        "status": "ok" if not issues else "needs_attention",
        "checks": checks,
        "issues": issues,
        "metrics": {
            **values,
            "gap_width": values["gap_high"] - values["gap_low"],
            "first_transverse_cutoff": cutoff,
            "bragg_wavenumber": bragg,
        },
        "lesson": (
            "Calibrate the empty lattice and empty duct before claiming a sonic-crystal gap. "
            "A confined Bloch gap must align with finite-crystal attenuation, while the same "
            "sparse inclusions in free space provide a negative control rather than a stop band."
        ),
    }
