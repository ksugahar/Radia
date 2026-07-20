"""H-B work and cross-stiffness artifact checks for v55."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


DEMAG = "demag_energy_hb_curvebranch_volume_material_owner_identity"
BEARING = "magneticbearing_crossstiffness_matrix_symmetry_coordinate_owner_identity"


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _vector(value: object, length: int) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and len(value) == length and all(_finite(item) for item in value)


def _demag_ok(row: Mapping[str, object]) -> bool:
    samples = row.get("hb_samples")
    samples_ok = isinstance(samples, Sequence) and not isinstance(samples, (str, bytes)) and len(samples) >= 3
    if samples_ok:
        samples_ok = all(isinstance(point, Mapping) and set(point) == {"h_a_per_m", "b_t"} and _finite(point["h_a_per_m"]) and _finite(point["b_t"]) for point in samples)
    if samples_ok:
        samples_ok = all(float(left["h_a_per_m"]) > float(right["h_a_per_m"]) and float(left["b_t"]) > float(right["b_t"]) for left, right in zip(samples, samples[1:]))
    if not samples_ok:
        return False
    work_density = sum(
        0.5 * (float(left["h_a_per_m"]) + float(right["h_a_per_m"])) * (float(right["b_t"]) - float(left["b_t"]))
        for left, right in zip(samples, samples[1:])
    )
    volume = row.get("material_volume_m3")
    energy = row.get("demag_energy_j")
    return (
        _generations(row, "hb_generation", "branch_generation", "volume_generation", "energy_generation", "material_generation", "owner_generation", "result_generation")
        and row.get("curve_branch") == "descending_demag"
        and row.get("result_curve_branch") == row.get("curve_branch")
        and row.get("result_hb_samples") == samples
        and work_density > 0.0
        and _finite(volume) and float(volume) > 0.0
        and row.get("result_material_volume_m3") == volume
        and _finite(energy) and float(energy) > 0.0
        and math.isclose(float(energy), work_density * float(volume), rel_tol=1.0e-10, abs_tol=1.0e-12)
        and row.get("result_demag_energy_j") == energy
        and bool(str(row.get("material_revision") or ""))
        and row.get("result_material_revision") == row.get("material_revision")
        and str(row.get("material_owner") or "").startswith("material:")
        and row.get("result_material_owner") == row.get("material_owner")
        and str(row.get("solution_owner") or "").startswith("solution:")
        and row.get("result_solution_owner") == row.get("solution_owner")
        and _result(row)
    )


def _bearing_ok(row: Mapping[str, object]) -> bool:
    matrix = row.get("cross_stiffness_n_per_m")
    tolerance = row.get("reciprocity_tolerance")
    matrix_ok = isinstance(matrix, Sequence) and not isinstance(matrix, (str, bytes)) and len(matrix) == 2 and all(_vector(line, 2) for line in matrix)
    if matrix_ok:
        a, b = (float(value) for value in matrix[0]); c, d = (float(value) for value in matrix[1])
        matrix_ok = _finite(tolerance) and 0.0 <= float(tolerance) <= 1.0e-6 and math.isclose(b, c, rel_tol=0.0, abs_tol=float(tolerance)) and a > 0.0 and d > 0.0 and a * d - b * c > 0.0
    basis = row.get("coordinate_basis")
    basis_ok = isinstance(basis, Mapping) and set(basis) == {"x", "y"} and all(_vector(vector, 3) and math.isclose(sum(float(item) ** 2 for item in vector), 1.0, abs_tol=1.0e-12) for vector in basis.values())
    if basis_ok:
        basis_ok = math.isclose(sum(float(x) * float(y) for x, y in zip(basis["x"], basis["y"])), 0.0, abs_tol=1.0e-12)
    load = row.get("load_point_m")
    return (
        _generations(row, "matrix_generation", "basis_generation", "reciprocity_generation", "loadpoint_generation", "owner_generation", "result_generation")
        and matrix_ok and row.get("result_cross_stiffness_n_per_m") == matrix
        and basis_ok and row.get("result_coordinate_basis") == basis
        and row.get("result_reciprocity_tolerance") == tolerance
        and _vector(load, 3) and row.get("result_load_point_m") == load
        and str(row.get("body_owner") or "").startswith("body:") and row.get("result_body_owner") == row.get("body_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    if DEMAG in identity:
        checks["magnetic_force_v55_demag_hb_work_volume_material_owner"] = isinstance(identity[DEMAG], Mapping) and _demag_ok(identity[DEMAG])
    if BEARING in identity:
        checks["magnetic_force_v55_bearing_cross_stiffness_basis_reciprocity_owner"] = isinstance(identity[BEARING], Mapping) and _bearing_ok(identity[BEARING])
    return checks
