"""Modal and nonlinear-continuation identity checks for v51 artifacts."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .adjoint_weakform_identity_v52 import validate_public_v52_identity


_EIGENMODE = "eigenmode_frequency_normalization_phase_subspace_mesh_owner_identity"
_CONTINUATION = "continuation_branch_predictor_corrector_loadpath_turningpoint_owner_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation_closed(row: Mapping[str, object], *names: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _result_identity_ok(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def _finite_sequence(value: object) -> list[float] | None:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)) or not value:
        return None
    try:
        values = [float(item) for item in value]
    except (TypeError, ValueError):
        return None
    return values if all(math.isfinite(item) for item in values) else None


def _eigenmode_ok(row: Mapping[str, object]) -> bool:
    frequencies = _finite_sequence(row.get("frequency_hz"))
    phase = row.get("phase_anchor")
    return (
        _generation_closed(
            row,
            "frequency_generation",
            "normalization_generation",
            "phase_generation",
            "subspace_generation",
            "mesh_generation",
            "owner_generation",
            "result_generation",
        )
        and frequencies is not None
        and len(frequencies) >= 2
        and all(value > 0.0 for value in frequencies)
        and max(frequencies) - min(frequencies) <= 1e-9 * max(frequencies)
        and row.get("result_frequency_hz") == row.get("frequency_hz")
        and row.get("normalization") == "unit_generalized_mass"
        and row.get("result_normalization") == row.get("normalization")
        and isinstance(phase, Mapping)
        and set(phase) == {"dof", "component", "sign"}
        and isinstance(phase.get("dof"), int)
        and not isinstance(phase.get("dof"), bool)
        and int(phase["dof"]) >= 0
        and phase.get("component") in {"real", "imag"}
        and phase.get("sign") in {"positive", "negative"}
        and row.get("result_phase_anchor") == phase
        and _digest(row.get("degenerate_subspace_basis_sha256"))
        and row.get("result_degenerate_subspace_basis_sha256") == row.get("degenerate_subspace_basis_sha256")
        and str(row.get("mesh_revision") or "").startswith("mesh:")
        and row.get("result_mesh_revision") == row.get("mesh_revision")
        and str(row.get("mode_owner") or "").startswith("mode-set:")
        and row.get("result_mode_owner") == row.get("mode_owner")
        and _result_identity_ok(row)
    )


def _continuation_ok(row: Mapping[str, object]) -> bool:
    states = row.get("predictor_corrector_states")
    load_path = _finite_sequence(row.get("load_path"))
    turning = row.get("turning_point_index")
    turning_valid = (
        load_path is not None
        and isinstance(turning, int)
        and not isinstance(turning, bool)
        and 0 < turning < len(load_path) - 1
        and load_path[turning] > load_path[turning - 1]
        and load_path[turning] > load_path[turning + 1]
    )
    return (
        _generation_closed(
            row,
            "branch_generation",
            "state_generation",
            "loadpath_generation",
            "turningpoint_generation",
            "owner_generation",
            "result_generation",
        )
        and str(row.get("branch_id") or "").startswith("branch:")
        and row.get("result_branch_id") == row.get("branch_id")
        and isinstance(states, Sequence)
        and not isinstance(states, (str, bytes))
        and len(states) == 2
        and str(states[0]).startswith("predictor:")
        and str(states[1]).startswith("corrector:")
        and row.get("result_predictor_corrector_states") == states
        and turning_valid
        and row.get("result_load_path") == row.get("load_path")
        and row.get("result_turning_point_index") == turning
        and str(row.get("solution_owner") or "").startswith("solution:")
        and row.get("result_solution_owner") == row.get("solution_owner")
        and _result_identity_ok(row)
    )


def validate_public_v51_identity(payload: object) -> dict[str, object]:
    """Validate optional degenerate-mode and continuation result records."""
    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    v52 = validate_public_v52_identity(payload)
    if v52:
        checks.update(v52["checks"])
    eigenmode = payload.get(_EIGENMODE)
    continuation = payload.get(_CONTINUATION)
    if eigenmode is not None:
        checks["v51_eigenmode_frequency_normalization_phase_subspace_mesh_owner"] = (
            isinstance(eigenmode, Mapping) and _eigenmode_ok(eigenmode)
        )
    if continuation is not None:
        checks["v51_continuation_branch_predictor_corrector_loadpath_turningpoint_owner"] = (
            isinstance(continuation, Mapping) and _continuation_ok(continuation)
        )
    if not checks:
        return {}
    return {
        "policy": "modal_continuation_identity_v51",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
    }
