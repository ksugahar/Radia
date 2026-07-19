"""Frequency-sweep and contact-result identity checks for v50 artifacts."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


_FREQUENCY = "frequency_sweep_complex_branch_phase_unit_dataset_interpolation_owner_identity"
_CONTACT = "contact_pair_augmented_lagrange_penalty_gap_pressure_frame_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation_closed(row: Mapping[str, object], *names: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _result_identity_ok(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_numbers(value: object) -> list[float] | None:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        return None
    try:
        values = [float(item) for item in value]
    except (TypeError, ValueError):
        return None
    return values if all(math.isfinite(item) for item in values) else None


def _strictly_increasing(values: Sequence[float]) -> bool:
    return all(left < right for left, right in zip(values, values[1:]))


def _frequency_ok(row: Mapping[str, object]) -> bool:
    frequencies = _finite_numbers(row.get("frequency_hz"))
    units = row.get("field_units")
    return (
        _generation_closed(
            row,
            "frequency_generation",
            "branch_generation",
            "phase_generation",
            "unit_generation",
            "dataset_generation",
            "interpolation_generation",
            "solution_generation",
            "result_generation",
        )
        and frequencies is not None
        and all(value > 0.0 for value in frequencies)
        and _strictly_increasing(frequencies)
        and row.get("result_frequency_hz") == row.get("frequency_hz")
        and row.get("complex_branch") == "positive_frequency"
        and row.get("result_complex_branch") == row.get("complex_branch")
        and row.get("phase_convention") == "exp(+jomega_t)"
        and row.get("result_phase_convention") == row.get("phase_convention")
        and units == {"electric_field": "V/m", "magnetic_field": "A/m"}
        and row.get("result_field_units") == units
        and str(row.get("dataset_tag") or "").startswith("dataset:")
        and row.get("result_dataset_tag") == row.get("dataset_tag")
        and row.get("dataset_interpolation") == "linear_complex"
        and row.get("result_dataset_interpolation") == row.get("dataset_interpolation")
        and str(row.get("solution_owner") or "").startswith("solution:")
        and row.get("result_solution_owner") == row.get("solution_owner")
        and _result_identity_ok(row)
    )


def _contact_ok(row: Mapping[str, object]) -> bool:
    gaps = _finite_numbers(row.get("gap_m"))
    pressures = _finite_numbers(row.get("pressure_pa"))
    penalty = row.get("penalty_factor")
    return (
        _generation_closed(
            row,
            "pair_generation",
            "method_generation",
            "penalty_generation",
            "gap_generation",
            "pressure_generation",
            "frame_generation",
            "owner_generation",
            "result_generation",
        )
        and str(row.get("contact_pair_id") or "").startswith("pair:")
        and row.get("result_contact_pair_id") == row.get("contact_pair_id")
        and row.get("contact_method") == "augmented_lagrange"
        and row.get("result_contact_method") == row.get("contact_method")
        and isinstance(penalty, (int, float))
        and math.isfinite(float(penalty))
        and float(penalty) > 0.0
        and row.get("result_penalty_factor") == penalty
        and gaps is not None
        and pressures is not None
        and len(gaps) == len(pressures)
        and all(value >= 0.0 for value in gaps)
        and all(value >= 0.0 for value in pressures)
        and all(left <= right for left, right in zip(gaps, gaps[1:]))
        and all(left >= right for left, right in zip(pressures, pressures[1:]))
        and row.get("result_gap_m") == row.get("gap_m")
        and row.get("result_pressure_pa") == row.get("pressure_pa")
        and row.get("coordinate_frame") == "spatial"
        and row.get("result_coordinate_frame") == row.get("coordinate_frame")
        and str(row.get("contact_owner") or "").startswith("contact:")
        and row.get("result_contact_owner") == row.get("contact_owner")
        and _result_identity_ok(row)
    )


def validate_public_v50_identity(payload: object) -> dict[str, object]:
    """Validate optional frequency-sweep and augmented-contact records."""
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    frequency = payload.get(_FREQUENCY)
    contact = payload.get(_CONTACT)
    if frequency is not None:
        checks["v50_frequency_branch_phase_unit_dataset_interpolation_owner"] = (
            isinstance(frequency, Mapping) and _frequency_ok(frequency)
        )
    if contact is not None:
        checks["v50_contact_pair_method_penalty_gap_pressure_frame_owner"] = (
            isinstance(contact, Mapping) and _contact_ok(contact)
        )
    if not checks:
        return {}
    return {
        "policy": "frequency_contact_identity_v50",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
    }
