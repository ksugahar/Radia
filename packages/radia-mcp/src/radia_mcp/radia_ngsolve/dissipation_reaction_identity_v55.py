"""Thermoelastic-damping and electrochemical-reaction identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


THERMO = (
    "thermoelastic_damping_complexeigenfrequency_energy_dissipation_"
    "normalization_owner_identity"
)
ELECTROCHEM = (
    "electrochem_current_species_stoichiometry_boundaryflux_time_owner_identity"
)
_FARADAY_C_PER_MOL = 96485.33212


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], *names: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _number(value: object, *, positive: bool = False, nonnegative: bool = False) -> bool:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return False
    number = float(value)
    return (
        math.isfinite(number)
        and (not positive or number > 0.0)
        and (not nonnegative or number >= 0.0)
    )


def _close(left: object, right: object) -> bool:
    return _number(left) and _number(right) and math.isclose(
        float(left), float(right), rel_tol=1.0e-10, abs_tol=1.0e-12
    )


def _result_identity(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get(
        "accepted_result_sha256"
    ) == row.get("result_sha256")


def _thermoelastic_ok(row: Mapping[str, object]) -> bool:
    eigenfrequency = row.get("complex_eigenfrequency_hz")
    if (
        not isinstance(eigenfrequency, Sequence)
        or isinstance(eigenfrequency, (str, bytes))
        or len(eigenfrequency) != 2
        or not all(_number(value) for value in eigenfrequency)
    ):
        return False
    frequency_hz, decay_hz = (float(value) for value in eigenfrequency)
    stored_energy = row.get("stored_energy_j")
    cycle_dissipation = row.get("cycle_dissipation_j")
    quality_factor = row.get("quality_factor")
    if frequency_hz <= 0.0 or decay_hz >= 0.0:
        return False
    expected_quality_factor = frequency_hz / (-2.0 * decay_hz)
    expected_dissipation = (
        2.0 * math.pi * float(stored_energy) / expected_quality_factor
        if _number(stored_energy, positive=True)
        else math.nan
    )
    names = (
        "complex_eigenfrequency_hz",
        "stored_energy_j",
        "cycle_dissipation_j",
        "quality_factor",
        "modal_normalization",
    )
    return (
        _generation(
            row,
            "eigenfrequency_generation",
            "energy_generation",
            "dissipation_generation",
            "normalization_generation",
            "owner_generation",
            "result_generation",
        )
        and _number(stored_energy, positive=True)
        and _number(cycle_dissipation, positive=True)
        and _number(quality_factor, positive=True)
        and _close(quality_factor, expected_quality_factor)
        and _close(cycle_dissipation, expected_dissipation)
        and row.get("modal_normalization") == "unit_total_stored_energy"
        and all(row.get("result_" + name) == row.get(name) for name in names)
        and str(row.get("solution_owner") or "").startswith("solution:")
        and row.get("result_solution_owner") == row.get("solution_owner")
        and _result_identity(row)
    )


def _electrochem_ok(row: Mapping[str, object]) -> bool:
    rates = row.get("species_rate_mol_s")
    charges = row.get("species_charge_number")
    boundary_flux = row.get("boundary_species_flux_mol_s")
    if not all(isinstance(value, Mapping) for value in (rates, charges, boundary_flux)):
        return False
    if (
        len(rates) < 2
        or set(rates) != set(charges)
        or not all(isinstance(name, str) and name and _number(rate) for name, rate in rates.items())
        or not all(
            isinstance(charge, int) and not isinstance(charge, bool)
            for charge in charges.values()
        )
        or not set(boundary_flux).issubset(rates)
        or not boundary_flux
        or not all(_number(value) for value in boundary_flux.values())
    ):
        return False
    ionic_current = _FARADAY_C_PER_MOL * sum(
        int(charges[name]) * float(rates[name])
        for name in rates
        if int(charges[name]) > 0
    )
    electronic_current = _FARADAY_C_PER_MOL * sum(
        int(charges[name]) * float(rates[name])
        for name in rates
        if int(charges[name]) < 0
    )
    terminal_current = row.get("terminal_current_a")
    charge_conserved = _close(ionic_current + electronic_current, 0.0)
    flux_closed = all(
        _close(float(boundary_flux[name]) + float(rates[name]), 0.0)
        for name in boundary_flux
    )
    names = (
        "terminal_current_a",
        "species_rate_mol_s",
        "species_charge_number",
        "boundary_species_flux_mol_s",
        "time_s",
    )
    return (
        _generation(
            row,
            "current_generation",
            "species_generation",
            "stoichiometry_generation",
            "flux_generation",
            "time_generation",
            "owner_generation",
            "result_generation",
        )
        and _number(terminal_current)
        and _close(terminal_current, ionic_current)
        and charge_conserved
        and flux_closed
        and _number(row.get("time_s"), nonnegative=True)
        and all(row.get("result_" + name) == row.get(name) for name in names)
        and str(row.get("solution_owner") or "").startswith("solution:")
        and row.get("result_solution_owner") == row.get("solution_owner")
        and _result_identity(row)
    )


def validate_public_v55_identity(payload: object) -> dict[str, object]:
    """Validate optional v55 damping and reaction identities."""
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    thermo = payload.get(THERMO)
    electrochem = payload.get(ELECTROCHEM)
    if thermo is not None:
        checks["v55_thermoelastic_damping_energy_normalization_owner"] = (
            isinstance(thermo, Mapping) and _thermoelastic_ok(thermo)
        )
    if electrochem is not None:
        checks["v55_electrochem_current_stoichiometry_flux_time_owner"] = (
            isinstance(electrochem, Mapping) and _electrochem_ok(electrochem)
        )
    if not checks:
        return {}
    return {
        "policy": "dissipation_reaction_identity_v55",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
    }
