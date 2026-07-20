"""Solver-neutral magnetic force and nonlinear material lineage checks."""

from __future__ import annotations

import math
from collections.abc import Mapping


FORCE = "force_method_body_owner_sign_displacement_pair_causal_identity"
BH = "nonlinear_bh_operating_point_row_hysteresis_branch_mapping_identity"


def _sha(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _generation(row: Mapping[str, object], names: tuple[str, ...]) -> bool:
    value = str(row.get("generation") or "")
    return bool(value) and all(row.get(name) == value for name in names)


def _force_ok(row: Mapping[str, object]) -> bool:
    displacement = row.get("displacement_pair_m")
    coenergy = row.get("coenergy_pair_j")
    return (
        _generation(row, ("force_method_generation", "body_owner_generation", "displacement_pair_generation", "result_generation"))
        and row.get("force_method") == row.get("result_force_method") == "weighted_stress_tensor"
        and str(row.get("body_owner") or "").startswith("group:")
        and row.get("result_body_owner") == row.get("body_owner")
        and row.get("force_sign_convention") == row.get("result_force_sign_convention") == "positive_displacement_direction"
        and isinstance(displacement, list) and len(displacement) == 2
        and all(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) for value in displacement)
        and float(displacement[0]) < float(displacement[1])
        and row.get("result_displacement_pair_m") == displacement
        and isinstance(coenergy, list) and len(coenergy) == 2
        and all(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value)) for value in coenergy)
        and row.get("result_coenergy_pair_j") == coenergy
        and _sha(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def _bh_ok(row: Mapping[str, object]) -> bool:
    keys = row.get("operating_point_row_keys")
    branches = row.get("hysteresis_branches")
    history = row.get("excitation_history_sha256")
    return (
        _generation(row, ("operating_point_generation", "branch_generation", "history_generation", "result_generation"))
        and isinstance(keys, list) and bool(keys) and len(set(keys)) == len(keys)
        and row.get("result_operating_point_row_keys") == keys
        and isinstance(branches, list) and len(branches) == len(keys)
        and all(branch in {"ascending", "descending", "initial"} for branch in branches)
        and row.get("result_hysteresis_branches") == branches
        and _sha(history) and row.get("result_excitation_history_sha256") == history
        and str(row.get("material_owner") or "").startswith("material:")
        and row.get("result_material_owner") == row.get("material_owner")
        and _sha(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def validate_public_identity(identity: object) -> dict[str, bool]:
    if not isinstance(identity, Mapping):
        return {}
    checks: dict[str, bool] = {}
    force = identity.get(FORCE)
    bh = identity.get(BH)
    if force is not None:
        checks["v47_force_body_sign_displacement_causality"] = isinstance(force, Mapping) and _force_ok(force)
    if bh is not None:
        checks["v47_bh_operating_row_branch_history"] = isinstance(bh, Mapping) and _bh_ok(bh)
    return checks
