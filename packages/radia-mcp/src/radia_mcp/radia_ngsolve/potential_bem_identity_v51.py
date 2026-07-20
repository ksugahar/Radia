"""Coupled-potential and BEM-matrix artifact identity checks for v51."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


POTENTIAL = "scalar_vector_potential_gauge_domain_interface_trace_solution_owner_identity"
BEM = "bem_matrix_reciprocity_symmetry_panel_orientation_cache_revision_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _domains(value: object) -> bool:
    return isinstance(value, list) and bool(value) and all(isinstance(name, str) and bool(name) for name in value) and len(value) == len(set(value))


def _potential_ok(row: Mapping[str, object]) -> bool:
    scalar = row.get("scalar_potential_domains")
    vector = row.get("vector_potential_domains")
    traces = row.get("interface_traces")
    traces_ok = (
        isinstance(traces, list)
        and bool(traces)
        and all(
            isinstance(item, Mapping)
            and item.get("scalar_domain") in scalar
            and item.get("vector_domain") in vector
            and item.get("trace") == "tangential_continuity"
            for item in traces
        )
    ) if _domains(scalar) and _domains(vector) else False
    return (
        _generations(row, "gauge_generation", "domain_generation", "interface_generation", "trace_generation", "owner_generation", "result_generation")
        and row.get("gauge") == "coulomb"
        and row.get("result_gauge") == row.get("gauge")
        and _domains(scalar)
        and row.get("result_scalar_potential_domains") == scalar
        and _domains(vector)
        and set(scalar).isdisjoint(vector)
        and row.get("result_vector_potential_domains") == vector
        and traces_ok
        and row.get("result_interface_traces") == traces
        and _digest(row.get("interface_trace_sha256"))
        and row.get("result_interface_trace_sha256") == row.get("interface_trace_sha256")
        and str(row.get("solution_owner") or "").startswith("solution:")
        and row.get("result_solution_owner") == row.get("solution_owner")
        and _result(row)
    )


def _bem_ok(row: Mapping[str, object]) -> bool:
    shape = row.get("matrix_shape")
    reciprocity = row.get("reciprocity_relative_error")
    return (
        _generations(row, "reciprocity_generation", "symmetry_generation", "orientation_generation", "cache_generation", "owner_generation", "result_generation")
        and isinstance(shape, Sequence)
        and not isinstance(shape, (str, bytes))
        and len(shape) == 2
        and all(isinstance(value, int) and not isinstance(value, bool) and value > 0 for value in shape)
        and shape[0] == shape[1]
        and row.get("result_matrix_shape") == shape
        and isinstance(reciprocity, (int, float))
        and not isinstance(reciprocity, bool)
        and math.isfinite(float(reciprocity))
        and 0.0 <= float(reciprocity) <= 1.0e-10
        and row.get("result_reciprocity_relative_error") == reciprocity
        and row.get("symmetry_class") in {"symmetric", "hermitian"}
        and row.get("result_symmetry_class") == row.get("symmetry_class")
        and row.get("panel_orientation") == "outward"
        and row.get("result_panel_orientation") == row.get("panel_orientation")
        and _digest(row.get("panel_orientation_sha256"))
        and row.get("result_panel_orientation_sha256") == row.get("panel_orientation_sha256")
        and str(row.get("cache_revision") or "").startswith("cache:")
        and row.get("result_cache_revision") == row.get("cache_revision")
        and str(row.get("matrix_owner") or "").startswith("matrix:")
        and row.get("result_matrix_owner") == row.get("matrix_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    potential = identity.get(POTENTIAL)
    bem = identity.get(BEM)
    if potential is not None:
        checks["magnetic_force_v51_coupled_potential_gauge_domain_trace_owner"] = isinstance(potential, Mapping) and _potential_ok(potential)
    if bem is not None:
        checks["magnetic_force_v51_bem_reciprocity_symmetry_orientation_cache_owner"] = isinstance(bem, Mapping) and _bem_ok(bem)
    return checks
