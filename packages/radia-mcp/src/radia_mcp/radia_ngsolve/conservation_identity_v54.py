"""Piezoelectric-work and reacting-flow conservation identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping

from .dissipation_reaction_identity_v55 import validate_public_v55_identity


PIEZO = "piezoelectric_energy_reciprocity_voltage_charge_work_phase_owner_identity"
SPECIES = "reactingflow_species_massfraction_rate_flux_time_solution_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], *names: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _number(value: object, *, nonnegative: bool = False) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
        and (not nonnegative or float(value) >= 0.0)
    )


def _close(left: object, right: object) -> bool:
    return _number(left) and _number(right) and math.isclose(float(left), float(right), rel_tol=1.0e-10, abs_tol=1.0e-12)


def _result_identity(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _piezo_ok(row: Mapping[str, object]) -> bool:
    voltage = row.get("voltage_v")
    charge = row.get("charge_c")
    electric_work = row.get("electric_work_j")
    mechanical_work = row.get("mechanical_work_j")
    phase = row.get("harmonic_phase_rad")
    return (
        _generation(row, "electric_generation", "mechanical_generation", "phase_generation", "owner_generation", "result_generation")
        and all(_number(value) for value in (voltage, charge, electric_work, mechanical_work, phase))
        and -math.pi <= float(phase) <= math.pi
        and _close(electric_work, 0.5 * float(voltage) * float(charge))
        and _close(mechanical_work, electric_work)
        and all(row.get("result_" + name) == row.get(name) for name in ("voltage_v", "charge_c", "electric_work_j", "mechanical_work_j", "harmonic_phase_rad"))
        and str(row.get("solution_owner") or "").startswith("solution:")
        and row.get("result_solution_owner") == row.get("solution_owner")
        and _result_identity(row)
    )


def _species_ok(row: Mapping[str, object]) -> bool:
    fractions = row.get("species_mass_fraction")
    rates = row.get("stoichiometric_rate_mol_m3_s")
    mass_rate = row.get("total_species_mass_rate_kg_s")
    boundary_flux = row.get("boundary_mass_flux_kg_s")
    time_s = row.get("time_s")
    fractions_ok = (
        isinstance(fractions, Mapping)
        and len(fractions) >= 2
        and all(isinstance(name, str) and name and _number(value, nonnegative=True) for name, value in fractions.items())
        and _close(sum(float(value) for value in fractions.values()), 1.0)
    )
    rates_ok = (
        isinstance(rates, Mapping)
        and set(rates) == set(fractions or {})
        and all(_number(value) for value in rates.values())
        and any(float(value) < 0.0 for value in rates.values())
        and any(float(value) > 0.0 for value in rates.values())
    )
    return (
        _generation(row, "fraction_generation", "rate_generation", "flux_generation", "time_generation", "owner_generation", "result_generation")
        and fractions_ok and rates_ok
        and _number(mass_rate) and _number(boundary_flux) and _close(float(mass_rate) + float(boundary_flux), 0.0)
        and _number(time_s, nonnegative=True)
        and all(row.get("result_" + name) == row.get(name) for name in ("species_mass_fraction", "stoichiometric_rate_mol_m3_s", "total_species_mass_rate_kg_s", "boundary_mass_flux_kg_s", "time_s"))
        and str(row.get("solution_owner") or "").startswith("solution:")
        and row.get("result_solution_owner") == row.get("solution_owner")
        and _result_identity(row)
    )


def validate_public_v54_identity(payload: object) -> dict[str, object]:
    """Validate optional v54 conservation identities in a public result payload."""
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    v55 = validate_public_v55_identity(payload)
    if v55:
        checks.update(v55["checks"])
    piezo = payload.get(PIEZO)
    species = payload.get(SPECIES)
    if piezo is not None:
        checks["v54_piezoelectric_work_reciprocity_phase_owner"] = isinstance(piezo, Mapping) and _piezo_ok(piezo)
    if species is not None:
        checks["v54_reacting_species_fraction_rate_flux_time_owner"] = isinstance(species, Mapping) and _species_ok(species)
    if not checks:
        return {}
    return {
        "policy": "conservation_identity_v54",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
    }
