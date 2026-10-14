"""Network artifact-identity gates (v47..v50).

One positive closure per gate; every negative test keeps its own failure
signal.  Payload builders live in ``_network_generalization_payloads``.
"""

from __future__ import annotations

from copy import deepcopy

from _network_generalization_payloads import (
    FARFIELD,
    FIELD,
    FLOQUET,
    MODAL,
    SMATRIX,
    TD_PORT,
    THERMAL,
    WAKE,
    _payload_v47,
    _payload_v48,
    _payload_v49,
    _payload_v50,
)
from radia_mcp.radia_ngsolve.network_artifact_identity_v49 import validate_public_v49_identity
from radia_mcp.radia_ngsolve.network_artifact_identity_v50 import validate_public_v50_identity
from radia_mcp.radia_ngsolve.network_artifact_lineage_v47 import validate_public_v47_identity
from radia_mcp.radia_ngsolve.network_semantic_identity_v48 import validate_public_v48_identity


# --- network_artifact_lineage_v47 (v47) -------------------------------------


def test_v47_positive_network_artifacts_are_accepted() -> None:
    assert all(validate_public_v47_identity(_payload_v47()).values())


def test_v47_smatrix_owner_mapping_mutation_is_rejected() -> None:
    payload = _payload_v47()
    row = payload["runs"][0][SMATRIX]
    row["result_port_order"] = ["P2", "P1"]
    row["result_reference_plane_m"] = {"P1": 0.01, "P2": 0.0}
    row["result_network_normalization_ohm"] = 75.0
    row["result_network_owner"] = "network:other"
    assert not all(validate_public_v47_identity(payload).values())


def test_v47_field_monitor_frequency_rows_mutation_is_rejected() -> None:
    payload = _payload_v47()
    row = payload["runs"][0][FIELD]
    row["result_monitor_identity"] = "monitor:other"
    row["result_frequency_row_keys"] = ["f=1.1GHz", "f=1.0GHz"]
    row["result_q_factor"] = [20.0, 5.0]
    assert not all(validate_public_v47_identity(payload).values())


# --- network_semantic_identity_v48 (v48) ------------------------------------


def test_v48_positive_network_artifacts_are_accepted() -> None:
    assert all(validate_public_v48_identity(_payload_v48()).values())


def test_v48_farfield_identity_mutation_is_rejected() -> None:
    payload = _payload_v48()
    payload["runs"][0][FARFIELD].update({"result_polarization_basis": "cartesian_xy", "result_normalization": "peak_field", "result_theta_deg": [90.0, 45.0, 0.0], "result_monitor_owner": "monitor:old"})
    assert not all(validate_public_v48_identity(payload).values())


def test_v48_em_thermal_mapping_mutation_is_rejected() -> None:
    payload = _payload_v48()
    payload["runs"][0][THERMAL].update({"mapped_loss_component_ids": ["dielectric", "conductor"], "mapped_interpolation_method": "nearest", "mapped_time_average": "peak", "mapped_frequency_hz": 2.4e9})
    assert not all(validate_public_v48_identity(payload).values())


def test_v48_self_consistent_nonphysical_conventions_are_rejected() -> None:
    payload = deepcopy(_payload_v48())
    farfield = payload["runs"][0][FARFIELD]
    farfield["polarization_basis"] = farfield["result_polarization_basis"] = "cartesian_xy"
    thermal = payload["runs"][0][THERMAL]
    thermal["interpolation_method"] = thermal["mapped_interpolation_method"] = "nearest"
    assert not all(validate_public_v48_identity(payload).values())


# --- network_artifact_identity_v49 (v49) ------------------------------------


def test_v49_positive_network_artifacts_are_accepted() -> None:
    assert all(validate_public_v49_identity(_payload_v49()).values())


def test_v49_modal_identity_mutation_is_rejected() -> None:
    payload = _payload_v49()
    payload["runs"][0][MODAL].update({"result_mode_index": 2, "result_degeneracy_basis": [[0.0, 1.0], [1.0, 0.0]], "result_polarization_phase_deg": -90.0, "result_reference_impedance_ohm": 75.0, "result_mesh_sha256": "8" * 64, "result_port_owner": "port:old"})
    assert not all(validate_public_v49_identity(payload).values())


def test_v49_wake_identity_mutation_is_rejected() -> None:
    payload = _payload_v49()
    payload["runs"][0][WAKE].update({"result_bunch_charge_c": 2.0e-9, "result_time_reference": "simulation_start", "result_monitor_time_s": [0.0, 2.0e-12, 1.0e-12, 3.0e-12], "result_normalization": "absolute_voltage", "accepted_result_owner": "result:old"})
    assert not all(validate_public_v49_identity(payload).values())


def test_v49_self_consistent_nonphysical_network_artifacts_are_rejected() -> None:
    payload = deepcopy(_payload_v49())
    modal = payload["runs"][0][MODAL]
    modal["degeneracy_basis"] = modal["result_degeneracy_basis"] = [[1.0, 0.0], [1.0, 0.0]]
    wake = payload["runs"][0][WAKE]
    wake["normalization"] = wake["result_normalization"] = "absolute_voltage"
    assert not all(validate_public_v49_identity(payload).values())


# --- network_artifact_identity_v50 (v50) ------------------------------------


def test_v50_positive_periodic_and_time_port_artifacts_are_accepted() -> None:
    assert all(validate_public_v50_identity(_payload_v50()).values())


def test_v50_floquet_phase_lattice_normalization_and_owner_drift_is_rejected() -> None:
    payload = deepcopy(_payload_v50())
    payload["runs"][0][FLOQUET]["result_floquet_phase_deg"] = [-15.0, 10.0]
    payload["runs"][0][FLOQUET]["result_lattice_vectors_m"] = list(reversed(payload["runs"][0][FLOQUET]["lattice_vectors_m"]))
    payload["runs"][0][FLOQUET]["result_mode_normalization"] = {"kind": "unit-voltage", "value_v": 1.0}
    payload["runs"][0][FLOQUET]["result_boundary_owner"] = "boundary:foreign"
    assert not all(validate_public_v50_identity(payload).values())


def test_v50_time_port_bandwidth_deembedding_reference_and_owner_drift_is_rejected() -> None:
    payload = deepcopy(_payload_v50())
    payload["runs"][0][TD_PORT]["result_pulse_bandwidth_hz"] = [2.0e9, 5.0e9]
    payload["runs"][0][TD_PORT]["result_deembedding_distance_m"] = -0.002
    payload["runs"][0][TD_PORT]["result_reference_plane"] = "port-plane:foreign"
    payload["runs"][0][TD_PORT]["result_port_owner"] = "port:foreign"
    assert not all(validate_public_v50_identity(payload).values())
