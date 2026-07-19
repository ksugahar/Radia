"""Cross-artifact lineage checks for v47 public solver summaries."""

from __future__ import annotations

from collections.abc import Mapping


_ROWS = "force_torque_energy_parameter_row_key_identity"
_CACHE = "model_mesh_study_result_cache_owner_chain_identity"


def _digest(value: object) -> bool:
    text = str(value or "").lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _unique_strings(value: object) -> bool:
    return (
        isinstance(value, list)
        and bool(value)
        and all(isinstance(item, str) and bool(item.strip()) for item in value)
        and len(set(value)) == len(value)
    )


def _result_identity_ok(row: Mapping[str, object]) -> bool:
    return (
        bool(str(row.get("owner") or ""))
        and row.get("accepted_owner") == row.get("owner")
        and _digest(row.get("result_sha256"))
        and row.get("accepted_result_sha256") == row.get("result_sha256")
    )


def _rows_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation") or "")
    keys = row.get("parameter_row_keys")
    return (
        bool(generation)
        and all(row.get(name) == generation for name in (
            "force_generation", "torque_generation", "energy_generation", "result_generation"
        ))
        and _unique_strings(keys)
        and row.get("force_parameter_row_keys") == keys
        and row.get("torque_parameter_row_keys") == keys
        and row.get("energy_parameter_row_keys") == keys
        and _digest(row.get("parameter_row_order_sha256"))
        and row.get("result_parameter_row_order_sha256") == row.get("parameter_row_order_sha256")
        and _result_identity_ok(row)
    )


def _cache_ok(row: Mapping[str, object]) -> bool:
    generation = str(row.get("generation") or "")
    chain = row.get("owner_chain")
    return (
        bool(generation)
        and all(row.get(name) == generation for name in (
            "model_generation", "mesh_generation", "study_generation", "solution_generation",
            "result_generation", "cache_generation"
        ))
        and _unique_strings(chain)
        and len(chain) == 5
        and row.get("cached_result_owner_chain") == chain
        and _digest(row.get("model_mesh_study_result_sha256"))
        and row.get("cached_owner_chain_sha256") == row.get("model_mesh_study_result_sha256")
        and _result_identity_ok(row)
    )


def validate_public_identity(payload: object) -> dict[str, object]:
    """Validate optional v47 records without changing older artifact behavior."""

    if not isinstance(payload, Mapping):
        return {}
    checks: dict[str, bool] = {}
    rows = payload.get(_ROWS)
    cache = payload.get(_CACHE)
    if rows is not None:
        checks["v47_multi_output_parameter_rows"] = isinstance(rows, Mapping) and _rows_ok(rows)
    if cache is not None:
        checks["v47_model_mesh_study_result_cache_chain"] = isinstance(cache, Mapping) and _cache_ok(cache)
    if not checks:
        return {}
    return {
        "policy": "comsol_v47_public_identity_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
    }
