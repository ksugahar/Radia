"""Magnetostatic virtual-work and heat-balance artifact identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


FORCE = "magnetostatic_energy_coenergy_force_displacement_derivative_owner_identity"
HEAT = "heatflow_joulesource_temperature_flux_balance_region_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _number(value: object, *, positive: bool = False, nonnegative: bool = False) -> bool:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return False
    number = float(value)
    return math.isfinite(number) and (not positive or number > 0.0) and (not nonnegative or number >= 0.0)


def _close(left: object, right: object) -> bool:
    return _number(left) and _number(right) and math.isclose(
        float(left), float(right), rel_tol=1.0e-9, abs_tol=1.0e-12
    )


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _force_ok(row: Mapping[str, object]) -> bool:
    samples = row.get("coenergy_samples")
    if not isinstance(samples, Sequence) or isinstance(samples, (str, bytes)) or len(samples) != 2:
        return False
    parsed: list[tuple[float, float]] = []
    for sample in samples:
        if not isinstance(sample, Mapping) or set(sample) != {"displacement_m", "coenergy_j"} or not _number(sample["displacement_m"]) or not _number(sample["coenergy_j"], nonnegative=True):
            return False
        parsed.append((float(sample["displacement_m"]), float(sample["coenergy_j"])))
    parsed.sort()
    delta = parsed[1][0] - parsed[0][0]
    derivative = (parsed[1][1] - parsed[0][1]) / delta if delta > 0.0 else math.nan
    names = ("magnetic_energy_j", "coenergy_samples", "displacement_step_m", "coenergy_derivative_n", "force_n", "solution_owner")
    return (
        _generations(row, "energy_generation", "coenergy_generation", "displacement_generation", "derivative_generation", "force_generation", "owner_generation", "result_generation")
        and _number(row.get("magnetic_energy_j"), nonnegative=True)
        and _number(row.get("displacement_step_m"), positive=True)
        and _close(2.0 * float(row["displacement_step_m"]), delta)
        and _close(row.get("coenergy_derivative_n"), derivative)
        and _close(row.get("force_n"), derivative)
        and str(row.get("solution_owner") or "").startswith("solution:")
        and all(row.get("result_" + name) == row.get(name) for name in names)
        and _result(row)
    )


def _heat_ok(row: Mapping[str, object]) -> bool:
    temperatures = row.get("temperature_field_k")
    regions = row.get("region_joule_balance_w")
    temperatures_ok = (
        isinstance(temperatures, Mapping)
        and set(temperatures) == {"minimum_k", "maximum_k", "mean_k"}
        and all(_number(value, positive=True) for value in temperatures.values())
        and float(temperatures["minimum_k"]) <= float(temperatures["mean_k"]) <= float(temperatures["maximum_k"])
    )
    regions_ok = (
        isinstance(regions, Mapping)
        and bool(regions)
        and all(isinstance(name, str) and name.startswith("region:") and _number(value, nonnegative=True) for name, value in regions.items())
    )
    source = row.get("joule_source_w")
    flux = row.get("outward_boundary_flux_w")
    names = ("joule_source_w", "temperature_field_k", "outward_boundary_flux_w", "region_joule_balance_w", "balance_residual_w", "solution_owner")
    return (
        _generations(row, "source_generation", "temperature_generation", "flux_generation", "balance_generation", "region_generation", "owner_generation", "result_generation")
        and _number(source, nonnegative=True)
        and temperatures_ok
        and _number(flux, nonnegative=True)
        and regions_ok
        and _close(sum(float(value) for value in regions.values()), source)
        and _close(row.get("balance_residual_w"), float(source) - float(flux))
        and _close(source, flux)
        and str(row.get("solution_owner") or "").startswith("solution:")
        and all(row.get("result_" + name) == row.get(name) for name in names)
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    if identity.get(FORCE) is not None:
        checks["v56_magnetostatic_coenergy_derivative_force_owner"] = isinstance(identity[FORCE], Mapping) and _force_ok(identity[FORCE])
    if identity.get(HEAT) is not None:
        checks["v56_heat_joule_temperature_flux_region_owner"] = isinstance(identity[HEAT], Mapping) and _heat_ok(identity[HEAT])
    return checks
