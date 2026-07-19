from __future__ import annotations

from radia_mcp.radia_ngsolve.network_artifact_lineage_v47 import FIELD, SMATRIX, validate_public_v47_identity


PROMOTED_CASE_IDS = {
    "v47_public_smatrix_port_order_reference_plane_network_owner_mismatch",
    "v47_public_field_energy_loss_q_monitor_frequency_row_key_mismatch",
}


def _payload() -> dict[str, object]:
    smatrix_generation = "smatrix-v47"
    field_generation = "field-v47"
    planes = {"P1": 0.0, "P2": 0.0}
    keys = ["f=1.0GHz", "f=1.1GHz"]
    return {
        "runs": [
            {
                SMATRIX: {
                    "generation": smatrix_generation,
                    **{
                        key: smatrix_generation
                        for key in (
                            "port_generation",
                            "reference_plane_generation",
                            "normalization_generation",
                            "network_generation",
                            "result_generation",
                        )
                    },
                    "port_order": ["P1", "P2"],
                    "result_port_order": ["P1", "P2"],
                    "reference_plane_m": planes,
                    "result_reference_plane_m": planes,
                    "network_normalization_ohm": 50.0,
                    "result_network_normalization_ohm": 50.0,
                    "network_owner": "network:test",
                    "result_network_owner": "network:test",
                    "result_sha256": "1" * 64,
                    "accepted_result_sha256": "1" * 64,
                },
                FIELD: {
                    "generation": field_generation,
                    **{
                        key: field_generation
                        for key in (
                            "monitor_generation",
                            "frequency_generation",
                            "energy_generation",
                            "loss_generation",
                            "q_generation",
                            "result_generation",
                        )
                    },
                    "monitor_identity": "monitor:test",
                    "result_monitor_identity": "monitor:test",
                    "frequency_row_keys": keys,
                    "result_frequency_row_keys": keys,
                    "field_energy_j": [1.0, 2.0],
                    "result_field_energy_j": [1.0, 2.0],
                    "loss_w": [0.1, 0.2],
                    "result_loss_w": [0.1, 0.2],
                    "q_factor": [10.0, 10.0],
                    "result_q_factor": [10.0, 10.0],
                    "result_owner": "result:test",
                    "result_result_owner": "result:test",
                    "result_sha256": "2" * 64,
                    "accepted_result_sha256": "2" * 64,
                },
            }
        ]
    }


def test_v47_positive_network_artifacts_are_accepted() -> None:
    assert all(validate_public_v47_identity(_payload()).values())


def test_v47_smatrix_owner_mapping_mutation_is_rejected() -> None:
    payload = _payload()
    row = payload["runs"][0][SMATRIX]
    row["result_port_order"] = ["P2", "P1"]
    row["result_reference_plane_m"] = {"P1": 0.01, "P2": 0.0}
    row["result_network_normalization_ohm"] = 75.0
    row["result_network_owner"] = "network:other"
    assert not all(validate_public_v47_identity(payload).values())


def test_v47_field_monitor_frequency_rows_mutation_is_rejected() -> None:
    payload = _payload()
    row = payload["runs"][0][FIELD]
    row["result_monitor_identity"] = "monitor:other"
    row["result_frequency_row_keys"] = ["f=1.1GHz", "f=1.0GHz"]
    row["result_q_factor"] = [20.0, 5.0]
    assert not all(validate_public_v47_identity(payload).values())
