"""BEM interaction quadrature and maglev equilibrium identity checks for v53."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .energy_derivative_identity_v54 import validate_public_identity as validate_public_v54_identity


QUADRATURE = "bem_singular_quadrature_self_near_far_panel_owner_identity"
MAGLEV = "maglev_equilibrium_force_stiffness_displacement_body_owner_identity"


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generations(row: Mapping[str, object], *fields: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(field) == generation for field in fields)


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _quadrature_ok(row: Mapping[str, object]) -> bool:
    interactions = row.get("panel_interactions")
    interactions_ok = isinstance(interactions, Sequence) and not isinstance(interactions, (str, bytes)) and bool(interactions)
    classifications: set[str] = set()
    seen: set[tuple[int, int]] = set()
    if interactions_ok:
        for item in interactions:
            if not isinstance(item, Mapping) or set(item) != {"source_panel", "target_panel", "separation_over_size", "classification", "rule"}:
                interactions_ok = False
                break
            source = item["source_panel"]; target = item["target_panel"]; ratio = item["separation_over_size"]; classification = item["classification"]
            if not all(isinstance(panel, int) and not isinstance(panel, bool) and panel > 0 for panel in (source, target)) or not _finite(ratio) or float(ratio) < 0.0 or (source, target) in seen:
                interactions_ok = False
                break
            expected = "self" if source == target and math.isclose(float(ratio), 0.0, abs_tol=1.0e-15) else "near" if source != target and 0.0 < float(ratio) <= 2.0 else "far" if source != target and float(ratio) > 2.0 else None
            expected_rule = {"self": "duffy_singular", "near": "adaptive_near", "far": "gauss_far"}.get(expected)
            if classification != expected or item["rule"] != expected_rule:
                interactions_ok = False
                break
            seen.add((source, target)); classifications.add(classification)
    return (
        _generations(row, "classification_generation", "quadrature_generation", "panel_generation", "owner_generation", "result_generation")
        and interactions_ok
        and classifications == {"self", "near", "far"}
        and row.get("result_panel_interactions") == interactions
        and str(row.get("panel_owner") or "").startswith("panel-set:")
        and row.get("result_panel_owner") == row.get("panel_owner")
        and _result(row)
    )


def _numeric_list(value: object) -> bool:
    return isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and bool(value) and all(_finite(item) for item in value)


def _maglev_ok(row: Mapping[str, object]) -> bool:
    displacement = row.get("displacement_path_m")
    force = row.get("force_path_n")
    index = row.get("equilibrium_index")
    paths_ok = _numeric_list(displacement) and _numeric_list(force) and len(displacement) == len(force) and len(displacement) >= 3 and len(displacement) % 2 == 1 and all(float(left) < float(right) for left, right in zip(displacement, displacement[1:]))
    index_ok = paths_ok and isinstance(index, int) and not isinstance(index, bool) and 0 < index < len(displacement) - 1
    derivative = (float(force[index + 1]) - float(force[index - 1])) / (float(displacement[index + 1]) - float(displacement[index - 1])) if index_ok else math.nan
    return (
        _generations(row, "equilibrium_generation", "force_generation", "stiffness_generation", "displacement_generation", "owner_generation", "result_generation")
        and index_ok
        and math.isclose(float(displacement[index]), 0.0, abs_tol=1.0e-15)
        and row.get("result_displacement_path_m") == displacement
        and row.get("result_force_path_n") == force
        and row.get("result_equilibrium_index") == index
        and _finite(row.get("equilibrium_force_n"))
        and math.isclose(float(row["equilibrium_force_n"]), float(force[index]), rel_tol=1.0e-12, abs_tol=1.0e-12)
        and math.isclose(float(row["equilibrium_force_n"]), 0.0, abs_tol=1.0e-12)
        and row.get("result_equilibrium_force_n") == row.get("equilibrium_force_n")
        and _finite(row.get("stiffness_n_per_m"))
        and float(row["stiffness_n_per_m"]) < 0.0
        and math.isclose(float(row["stiffness_n_per_m"]), derivative, rel_tol=1.0e-12, abs_tol=1.0e-12)
        and row.get("result_stiffness_n_per_m") == row.get("stiffness_n_per_m")
        and str(row.get("body_owner") or "").startswith("body:")
        and row.get("result_body_owner") == row.get("body_owner")
        and _result(row)
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks = validate_public_v54_identity(identity)
    quadrature = identity.get(QUADRATURE)
    maglev = identity.get(MAGLEV)
    if quadrature is not None:
        checks["magnetic_force_v53_bem_quadrature_classification_panel_owner"] = isinstance(quadrature, Mapping) and _quadrature_ok(quadrature)
    if maglev is not None:
        checks["magnetic_force_v53_maglev_equilibrium_stiffness_body_owner"] = isinstance(maglev, Mapping) and _maglev_ok(maglev)
    return checks
