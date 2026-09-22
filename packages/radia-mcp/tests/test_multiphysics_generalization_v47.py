from __future__ import annotations

from radia_mcp.radia_ngsolve.cross_artifact_lineage_v47 import validate_public_identity


PROMOTED_CASE_IDS = {
    "v47_public_force_torque_energy_parameter_row_key_permutation_mismatch",
    "v47_public_model_mesh_study_result_cache_generation_owner_chain_mismatch",
}


def _records() -> dict[str, object]:
    rows_generation = "multi-output-rows-v47-901"
    cache_generation = "cache-chain-v47-901"
    keys = ["speed=1000|current=5", "speed=2000|current=5", "speed=3000|current=5"]
    chain = ["model:m1", "mesh:mesh1", "study:std1", "solution:sol1", "result:r1"]
    return {
        "force_torque_energy_parameter_row_key_identity": {
            "generation": rows_generation,
            "force_generation": rows_generation,
            "torque_generation": rows_generation,
            "energy_generation": rows_generation,
            "result_generation": rows_generation,
            "parameter_row_keys": keys,
            "force_parameter_row_keys": keys,
            "torque_parameter_row_keys": keys,
            "energy_parameter_row_keys": keys,
            "parameter_row_order_sha256": "1" * 64,
            "result_parameter_row_order_sha256": "1" * 64,
            "owner": "result/multi-output-v47-901",
            "accepted_owner": "result/multi-output-v47-901",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        },
        "model_mesh_study_result_cache_owner_chain_identity": {
            "generation": cache_generation,
            "model_generation": cache_generation,
            "mesh_generation": cache_generation,
            "study_generation": cache_generation,
            "solution_generation": cache_generation,
            "result_generation": cache_generation,
            "cache_generation": cache_generation,
            "owner_chain": chain,
            "cached_result_owner_chain": chain,
            "model_mesh_study_result_sha256": "3" * 64,
            "cached_owner_chain_sha256": "3" * 64,
            "owner": "cache/result-v47-901",
            "accepted_owner": "cache/result-v47-901",
            "result_sha256": "4" * 64,
            "accepted_result_sha256": "4" * 64,
        },
    }


def test_v47_positive_replays_are_accepted() -> None:
    result = validate_public_identity(_records())
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_v47_row_permutation_is_rejected() -> None:
    records = _records()
    row = records["force_torque_energy_parameter_row_key_identity"]
    row["torque_parameter_row_keys"] = list(reversed(row["parameter_row_keys"]))
    assert validate_public_identity(records)["status"] == "needs_attention"


def test_v47_stale_cache_chain_is_rejected() -> None:
    records = _records()
    row = records["model_mesh_study_result_cache_owner_chain_identity"]
    row["study_generation"] = "old"
    row["cached_result_owner_chain"] = ["model:m1", "mesh:mesh0", "study:std0", "solution:sol0", "result:r1"]
    assert validate_public_identity(records)["status"] == "needs_attention"
