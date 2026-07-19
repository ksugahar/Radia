from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.network_artifact_identity_v49 import MODAL, WAKE, validate_public_v49_identity


PROMOTED_CASE_IDS = {
    "v49_public_modal_port_mode_index_degeneracy_polarization_phase_impedance_mesh_owner_mismatch",
    "v49_public_particle_wakefield_bunch_charge_time_reference_monitor_normalization_owner_mismatch",
}


def _payload() -> dict[str, object]:
    modal_generation = "modal-port-v49"
    wake_generation = "wake-v49"
    basis = [[1.0, 0.0], [0.0, 1.0]]
    times = [0.0, 1.0e-12, 2.0e-12, 3.0e-12]
    wake = [0.0, 1.2e12, 0.8e12, 0.2e12]
    return {"runs": [{
        MODAL: {
            "generation": modal_generation,
            **{key: modal_generation for key in ("mode_generation", "degeneracy_generation", "polarization_generation", "impedance_generation", "mesh_generation", "port_generation", "result_generation")},
            "mode_index": 1, "result_mode_index": 1,
            "degeneracy_basis": basis, "result_degeneracy_basis": basis,
            "polarization_phase_deg": 90.0, "result_polarization_phase_deg": 90.0,
            "reference_impedance_ohm": 50.0, "result_reference_impedance_ohm": 50.0,
            "mesh_sha256": "1" * 64, "result_mesh_sha256": "1" * 64,
            "port_owner": "port:modal-v49", "result_port_owner": "port:modal-v49",
            "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
        },
        WAKE: {
            "generation": wake_generation,
            **{key: wake_generation for key in ("charge_generation", "time_generation", "monitor_generation", "normalization_generation", "result_generation")},
            "bunch_charge_c": 1.0e-9, "result_bunch_charge_c": 1.0e-9,
            "time_reference": "bunch_center", "result_time_reference": "bunch_center",
            "monitor_time_s": times, "result_monitor_time_s": times,
            "wake_v_per_c": wake, "result_wake_v_per_c": wake,
            "normalization": "per_coulomb", "result_normalization": "per_coulomb",
            "result_owner": "result:wake-v49", "accepted_result_owner": "result:wake-v49",
            "result_sha256": "3" * 64, "accepted_result_sha256": "3" * 64,
        },
    }]}


def test_v49_positive_network_artifacts_are_accepted() -> None:
    assert all(validate_public_v49_identity(_payload()).values())


def test_v49_modal_identity_mutation_is_rejected() -> None:
    payload = _payload()
    payload["runs"][0][MODAL].update({"result_mode_index": 2, "result_degeneracy_basis": [[0.0, 1.0], [1.0, 0.0]], "result_polarization_phase_deg": -90.0, "result_reference_impedance_ohm": 75.0, "result_mesh_sha256": "8" * 64, "result_port_owner": "port:old"})
    assert not all(validate_public_v49_identity(payload).values())


def test_v49_wake_identity_mutation_is_rejected() -> None:
    payload = _payload()
    payload["runs"][0][WAKE].update({"result_bunch_charge_c": 2.0e-9, "result_time_reference": "simulation_start", "result_monitor_time_s": [0.0, 2.0e-12, 1.0e-12, 3.0e-12], "result_normalization": "absolute_voltage", "accepted_result_owner": "result:old"})
    assert not all(validate_public_v49_identity(payload).values())


def test_v49_self_consistent_nonphysical_network_artifacts_are_rejected() -> None:
    payload = deepcopy(_payload())
    modal = payload["runs"][0][MODAL]
    modal["degeneracy_basis"] = modal["result_degeneracy_basis"] = [[1.0, 0.0], [1.0, 0.0]]
    wake = payload["runs"][0][WAKE]
    wake["normalization"] = wake["result_normalization"] = "absolute_voltage"
    assert not all(validate_public_v49_identity(payload).values())
