"""FEMM-derived electromagnetic artifact identity checks for v49."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


NONLINEAR = "nonlinear_bh_branch_incremental_permeability_temperature_lamination_owner_identity"
HARMONIC = "harmonic_geometry_depth_frequency_phase_circuit_loss_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_vector(value: object, length: int | None = None) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and (length is None or len(value) == length)
        and bool(value)
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def _bh_rows(value: object) -> bool:
    if not isinstance(value, list) or len(value) < 3 or not all(_finite_vector(row, 2) for row in value):
        return False
    magnetic_flux_density = [float(row[0]) for row in value]
    magnetic_field = [float(row[1]) for row in value]
    return (
        magnetic_flux_density[0] == 0.0
        and magnetic_field[0] == 0.0
        and all(right > left for left, right in zip(magnetic_flux_density, magnetic_flux_density[1:]))
        and all(right > left for left, right in zip(magnetic_field, magnetic_field[1:]))
    )


def _positive_definite_2x2(value: object) -> bool:
    if not isinstance(value, list) or len(value) != 2 or not all(_finite_vector(row, 2) for row in value):
        return False
    a, b = (float(item) for item in value[0])
    c, d = (float(item) for item in value[1])
    return a > 0.0 and d > 0.0 and math.isclose(b, c, rel_tol=1.0e-12, abs_tol=1.0e-12) and a * d - b * c > 0.0


def _lamination(value: object) -> bool:
    if not isinstance(value, Mapping):
        return False
    fill = value.get("fill_factor")
    thickness = value.get("sheet_thickness_m")
    return (
        isinstance(fill, (int, float))
        and math.isfinite(float(fill))
        and 0.0 < float(fill) <= 1.0
        and value.get("direction") in {"in-plane", "normal"}
        and isinstance(thickness, (int, float))
        and math.isfinite(float(thickness))
        and float(thickness) > 0.0
    )


def _nonlinear_ok(row: Mapping[str, object]) -> bool:
    branch = row.get("bh_branch")
    bh_rows = row.get("bh_rows_t_a_per_m")
    incremental = row.get("incremental_permeability_h_per_m")
    temperature = row.get("temperature_c")
    lamination = row.get("lamination")
    owner = str(row.get("material_owner") or "")
    return (
        _generations(
            row,
            "branch_generation",
            "incremental_generation",
            "temperature_generation",
            "lamination_generation",
            "result_generation",
        )
        and branch == "ascending:first-quadrant"
        and row.get("result_bh_branch") == branch
        and _bh_rows(bh_rows)
        and row.get("result_bh_rows_t_a_per_m") == bh_rows
        and _positive_definite_2x2(incremental)
        and row.get("result_incremental_permeability_h_per_m") == incremental
        and isinstance(temperature, (int, float))
        and math.isfinite(float(temperature))
        and -273.15 < float(temperature) <= 1000.0
        and row.get("result_temperature_c") == temperature
        and _lamination(lamination)
        and row.get("result_lamination") == lamination
        and owner.startswith("material:")
        and row.get("result_material_owner") == owner
        and _result(row)
    )


def _circuit(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and bool(str(value.get("name") or ""))
        and _finite_vector(value.get("current_a"), 2)
        and isinstance(value.get("turns"), int)
        and int(value["turns"]) > 0
    )


def _losses(value: object) -> bool:
    return (
        isinstance(value, Mapping)
        and set(value) == {"joule_w", "core_w"}
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) and float(item) >= 0.0 for item in value.values())
    )


def _harmonic_ok(row: Mapping[str, object]) -> bool:
    geometry = row.get("geometry_type")
    depth = row.get("depth_m")
    frequency = row.get("frequency_hz")
    phase = row.get("phase_convention")
    circuit = row.get("circuit")
    losses = row.get("losses_w")
    owner = str(row.get("result_owner") or "")
    return (
        _generations(
            row,
            "geometry_generation",
            "frequency_generation",
            "circuit_generation",
            "loss_generation",
            "result_generation",
        )
        and geometry in {"planar", "axisymmetric"}
        and row.get("result_geometry_type") == geometry
        and isinstance(depth, (int, float))
        and math.isfinite(float(depth))
        and float(depth) > 0.0
        and row.get("result_depth_m") == depth
        and isinstance(frequency, (int, float))
        and math.isfinite(float(frequency))
        and float(frequency) > 0.0
        and row.get("result_frequency_hz") == frequency
        and phase in {"exp(+jwt)", "exp(-jwt)"}
        and row.get("result_phase_convention") == phase
        and _circuit(circuit)
        and row.get("result_circuit") == circuit
        and _losses(losses)
        and row.get("result_losses_w") == losses
        and owner.startswith("harmonic-result:")
        and row.get("accepted_result_owner") == owner
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    """Validate optional nonlinear-material and harmonic-result identity records."""
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    nonlinear = identity.get(NONLINEAR)
    harmonic = identity.get(HARMONIC)
    if nonlinear is not None:
        checks["v49_nonlinear_bh_incremental_temperature_lamination_owner"] = isinstance(nonlinear, Mapping) and _nonlinear_ok(nonlinear)
    if harmonic is not None:
        checks["v49_harmonic_geometry_depth_frequency_phase_circuit_loss_owner"] = isinstance(harmonic, Mapping) and _harmonic_ok(harmonic)
    return checks
