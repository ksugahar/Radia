"""Mass-property, sweep-history, STEP, and BREP identity checks for v51."""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence

from .selector_exchange_identity_v52 import (
    validate_public_identity as validate_public_v52_identity,
    validate_source_identity as validate_source_v52_identity,
)


MASS = "mass_properties_frame_inertia_parallel_axis_density_shape_owner_identity"
SWEEP = "sweep_profile_path_trihedron_transition_selfintersection_history_owner_identity"
STEP = "step_external_reference_occurrence_name_schema_unit_color_owner_identity"
BREP = "brep_occt_version_location_precision_triangulation_cache_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *names: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _positive(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) and float(value) > 0.0


def _vector(value: object, size: int) -> list[float] | None:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != size:
        return None
    try:
        values = [float(item) for item in value]
    except (TypeError, ValueError):
        return None
    return values if all(math.isfinite(item) for item in values) else None


def _matrix3(value: object) -> list[list[float]] | None:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or len(value) != 3:
        return None
    rows = [_vector(row, 3) for row in value]
    return rows if all(row is not None for row in rows) else None  # type: ignore[return-value]


def _parallel_axis_ok(mass: float, centroidal: list[list[float]], shift: list[float], shifted: list[list[float]]) -> bool:
    norm2 = sum(value * value for value in shift)
    expected = [[centroidal[i][j] + mass * ((norm2 if i == j else 0.0) - shift[i] * shift[j]) for j in range(3)] for i in range(3)]
    return all(math.isclose(shifted[i][j], expected[i][j], rel_tol=1e-10, abs_tol=1e-12) for i in range(3) for j in range(3))


def _mass_ok(row: Mapping[str, object]) -> bool:
    mass = row.get("mass_kg")
    density = row.get("density_kg_m3")
    centroidal = _matrix3(row.get("centroidal_inertia_kg_m2"))
    shift = _vector(row.get("parallel_axis_shift_m"), 3)
    shifted = _matrix3(row.get("shifted_inertia_kg_m2"))
    tensors_ok = centroidal is not None and shifted is not None and all(
        math.isclose(matrix[i][j], matrix[j][i], abs_tol=1e-12)
        for matrix in (centroidal, shifted)
        for i in range(3)
        for j in range(3)
    ) and all(centroidal[i][i] > 0.0 and shifted[i][i] > 0.0 for i in range(3))
    return (
        _generations(row, "frame_generation", "inertia_generation", "shift_generation", "density_generation", "owner_generation", "result_generation")
        and row.get("coordinate_frame") == "global_cartesian"
        and row.get("result_coordinate_frame") == row.get("coordinate_frame")
        and _positive(mass)
        and row.get("result_mass_kg") == mass
        and tensors_ok
        and row.get("result_centroidal_inertia_kg_m2") == row.get("centroidal_inertia_kg_m2")
        and shift is not None
        and row.get("result_parallel_axis_shift_m") == row.get("parallel_axis_shift_m")
        and row.get("result_shifted_inertia_kg_m2") == row.get("shifted_inertia_kg_m2")
        and _parallel_axis_ok(float(mass), centroidal, shift, shifted)  # type: ignore[arg-type]
        and _positive(density)
        and row.get("result_density_kg_m3") == density
        and str(row.get("shape_owner") or "").startswith("shape:")
        and row.get("result_shape_owner") == row.get("shape_owner")
        and _result(row)
    )


def _history(value: object) -> bool:
    return isinstance(value, Mapping) and bool(value) and all(
        isinstance(name, str)
        and name
        and isinstance(items, Sequence)
        and not isinstance(items, (str, bytes))
        and bool(items)
        and all(isinstance(item, str) and item for item in items)
        for name, items in value.items()
    )


def _sweep_ok(row: Mapping[str, object]) -> bool:
    intersections = row.get("self_intersections")
    history = row.get("topology_history")
    return (
        _generations(row, "profile_generation", "path_generation", "trihedron_generation", "transition_generation", "intersection_generation", "history_generation", "owner_generation", "result_generation")
        and str(row.get("profile_id") or "").startswith("profile:")
        and row.get("result_profile_id") == row.get("profile_id")
        and str(row.get("path_id") or "").startswith("path:")
        and row.get("result_path_id") == row.get("path_id")
        and row.get("trihedron") in {"corrected_frenet", "frenet", "fixed"}
        and row.get("result_trihedron") == row.get("trihedron")
        and row.get("transition") in {"transformed", "right_corner", "round_corner"}
        and row.get("result_transition") == row.get("transition")
        and intersections == []
        and row.get("result_self_intersections") == intersections
        and _history(history)
        and row.get("result_topology_history") == history
        and str(row.get("shape_owner") or "").startswith("shape:")
        and row.get("result_shape_owner") == row.get("shape_owner")
        and _result(row)
    )


def _public_rows(payload: Mapping[str, object]) -> list[Mapping[str, object]]:
    rows: list[Mapping[str, object]] = []
    if isinstance(payload.get("reference"), list):
        rows.extend(row for row in payload["reference"] if isinstance(row, Mapping))
    measured = payload.get("measured")
    if isinstance(measured, Mapping):
        for values in measured.values():
            if isinstance(values, list):
                rows.extend(row for row in values if isinstance(row, Mapping))
    return rows


def _string_map(value: object, prefix: str) -> bool:
    return isinstance(value, Mapping) and bool(value) and all(isinstance(name, str) and name.startswith(prefix) and isinstance(item, str) and item for name, item in value.items())


def _step_ok(row: Mapping[str, object]) -> bool:
    references = row.get("external_references")
    names = row.get("occurrence_names")
    colors = row.get("occurrence_colors")
    return (
        _generations(row, "reference_generation", "name_generation", "schema_generation", "unit_generation", "color_generation", "owner_generation", "result_generation")
        and isinstance(references, Mapping)
        and bool(references)
        and all(isinstance(name, str) and name.startswith("document:") and isinstance(items, Sequence) and not isinstance(items, (str, bytes)) and bool(items) and all(isinstance(item, str) and item.startswith("part:") for item in items) for name, items in references.items())
        and row.get("result_external_references") == references
        and _string_map(names, "occurrence:")
        and row.get("result_occurrence_names") == names
        and row.get("step_schema") in {"AP203", "AP214", "AP242"}
        and row.get("result_step_schema") == row.get("step_schema")
        and row.get("length_unit") == "m"
        and row.get("result_length_unit") == row.get("length_unit")
        and isinstance(colors, Mapping)
        and set(colors) == set(names)
        and all((vector := _vector(color, 3)) is not None and all(0.0 <= value <= 1.0 for value in vector) for color in colors.values())
        and row.get("result_occurrence_colors") == colors
        and str(row.get("document_owner") or "").startswith("document:")
        and row.get("result_document_owner") == row.get("document_owner")
        and _result(row)
    )


def _rigid_location(value: object) -> bool:
    values = _vector(value, 7)
    return values is not None and math.isclose(sum(item * item for item in values[3:]), 1.0, rel_tol=1e-10, abs_tol=1e-12)


def _brep_ok(row: Mapping[str, object]) -> bool:
    precision = row.get("model_precision_m")
    return (
        _generations(row, "version_generation", "location_generation", "precision_generation", "triangulation_generation", "owner_generation", "result_generation")
        and re.fullmatch(r"\d+\.\d+\.\d+", str(row.get("occt_version") or "")) is not None
        and row.get("result_occt_version") == row.get("occt_version")
        and _rigid_location(row.get("shape_location"))
        and row.get("result_shape_location") == row.get("shape_location")
        and _positive(precision)
        and float(precision) <= 1e-3
        and row.get("result_model_precision_m") == precision
        and _digest(row.get("triangulation_cache_sha256"))
        and row.get("result_triangulation_cache_sha256") == row.get("triangulation_cache_sha256")
        and str(row.get("shape_owner") or "").startswith("shape:")
        and row.get("result_shape_owner") == row.get("shape_owner")
        and _result(row)
    )


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {"policy": policy, "status": "ok" if all(checks.values()) else "needs_attention", "checks": checks, "issues": [name for name, accepted in checks.items() if not accepted]}


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    rows = _public_rows(payload)
    checks: dict[str, bool] = {}
    v52 = validate_public_v52_identity(payload)
    if v52:
        checks.update(v52["checks"])
    mass = [row.get(MASS) for row in rows if MASS in row]
    sweep = [row.get(SWEEP) for row in rows if SWEEP in row]
    if mass:
        checks["v51_mass_frame_inertia_parallel_axis_density_owner"] = len(mass) == len(rows) and all(isinstance(item, Mapping) and _mass_ok(item) for item in mass)
    if sweep:
        checks["v51_sweep_profile_path_trihedron_transition_history_owner"] = len(sweep) == len(rows) and all(isinstance(item, Mapping) and _sweep_ok(item) for item in sweep)
    return _report("build123d_v51_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping) or not isinstance(payload.get("replay_identity"), Mapping):
        return {}
    identity = payload["replay_identity"]
    checks: dict[str, bool] = {}
    v52 = validate_source_v52_identity(payload)
    if v52:
        checks.update(v52["checks"])
    step = identity.get(STEP)
    brep = identity.get(BREP)
    if step is not None:
        checks["v51_step_references_names_schema_units_colors_owner"] = isinstance(step, Mapping) and _step_ok(step)
    if brep is not None:
        checks["v51_brep_occt_location_precision_triangulation_owner"] = isinstance(brep, Mapping) and _brep_ok(brep)
    return _report("build123d_v51_source_identity_v1", checks) if checks else {}
