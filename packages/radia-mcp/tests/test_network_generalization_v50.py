from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.network_artifact_identity_v50 import FLOQUET, TD_PORT, validate_public_v50_identity


PROMOTED_CASE_IDS = {
    "v50_public_periodic_floquet_phase_lattice_vector_mode_normalization_owner_mismatch",
    "v50_public_time_domain_port_pulse_bandwidth_deembedding_reference_plane_owner_mismatch",
}


def _payload() -> dict[str, object]:
    floquet_generation = "floquet-v50"
    port_generation = "td-port-v50"
    lattice = [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0]]
    normalization = {"kind": "unit-power", "value_w": 1.0}
    bandwidth = [1.0e9, 10.0e9]
    return {"runs": [{
        FLOQUET: {
            "generation": floquet_generation,
            **{key: floquet_generation for key in ("phase_generation", "lattice_generation", "mode_generation", "boundary_generation", "result_generation")},
            "floquet_phase_deg": [15.0, -10.0], "result_floquet_phase_deg": [15.0, -10.0],
            "lattice_vectors_m": lattice, "result_lattice_vectors_m": lattice,
            "mode_normalization": normalization, "result_mode_normalization": normalization,
            "boundary_owner": "boundary:floquet-v50", "result_boundary_owner": "boundary:floquet-v50",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        TD_PORT: {
            "generation": port_generation,
            **{key: port_generation for key in ("pulse_generation", "deembedding_generation", "reference_generation", "port_generation", "result_generation")},
            "pulse_bandwidth_hz": bandwidth, "result_pulse_bandwidth_hz": bandwidth,
            "deembedding_distance_m": 0.002, "result_deembedding_distance_m": 0.002,
            "reference_plane": "port-plane:z0", "result_reference_plane": "port-plane:z0",
            "port_owner": "port:td-v50", "result_port_owner": "port:td-v50",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
    }]}


def test_v50_positive_periodic_and_time_port_artifacts_are_accepted() -> None:
    assert all(validate_public_v50_identity(_payload()).values())


def test_v50_floquet_phase_lattice_normalization_and_owner_drift_is_rejected() -> None:
    payload = deepcopy(_payload())
    payload["runs"][0][FLOQUET]["result_floquet_phase_deg"] = [-15.0, 10.0]
    payload["runs"][0][FLOQUET]["result_lattice_vectors_m"] = list(reversed(payload["runs"][0][FLOQUET]["lattice_vectors_m"]))
    payload["runs"][0][FLOQUET]["result_mode_normalization"] = {"kind": "unit-voltage", "value_v": 1.0}
    payload["runs"][0][FLOQUET]["result_boundary_owner"] = "boundary:foreign"
    assert not all(validate_public_v50_identity(payload).values())


def test_v50_time_port_bandwidth_deembedding_reference_and_owner_drift_is_rejected() -> None:
    payload = deepcopy(_payload())
    payload["runs"][0][TD_PORT]["result_pulse_bandwidth_hz"] = [2.0e9, 5.0e9]
    payload["runs"][0][TD_PORT]["result_deembedding_distance_m"] = -0.002
    payload["runs"][0][TD_PORT]["result_reference_plane"] = "port-plane:foreign"
    payload["runs"][0][TD_PORT]["result_port_owner"] = "port:foreign"
    assert not all(validate_public_v50_identity(payload).values())
