"""Payload builders for ``test_network_generalization.py``.

Not collected by pytest (leading underscore).  Each builder feeds its own
``validate_public_v<N>_identity`` gate.
"""

from __future__ import annotations

from radia_mcp.radia_ngsolve.network_artifact_identity_v49 import MODAL, WAKE
from radia_mcp.radia_ngsolve.network_artifact_identity_v50 import FLOQUET, TD_PORT
from radia_mcp.radia_ngsolve.network_artifact_lineage_v47 import FIELD, SMATRIX
from radia_mcp.radia_ngsolve.network_semantic_identity_v48 import FARFIELD, THERMAL


# --- network_artifact_lineage_v47 (v47) -------------------------------------


def _payload_v47() -> dict[str, object]:
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


# --- network_semantic_identity_v48 (v48) ------------------------------------


def _payload_v48() -> dict[str, object]:
    farfield_generation = "farfield-v48"
    thermal_generation = "em-thermal-v48"
    theta = [0.0, 45.0, 90.0]
    phi = [0.0, 90.0, 180.0]
    samples = [[1.0, 0.2], [0.8, 0.3], [0.4, 0.1]]
    components = ["conductor", "dielectric"]
    losses = [2.5, 0.5]
    return {"runs": [{
        FARFIELD: {
            "generation": farfield_generation,
            **{key: farfield_generation for key in ("polarization_generation", "normalization_generation", "power_generation", "angular_grid_generation", "monitor_generation", "result_generation")},
            "polarization_basis": "spherical_theta_phi", "result_polarization_basis": "spherical_theta_phi",
            "normalization": "accepted_radiated_power", "result_normalization": "accepted_radiated_power",
            "radiated_power_w": 4.0, "result_radiated_power_w": 4.0,
            "theta_deg": theta, "result_theta_deg": theta, "phi_deg": phi, "result_phi_deg": phi,
            "field_samples": samples, "result_field_samples": samples,
            "monitor_owner": "monitor:farfield-v48", "result_monitor_owner": "monitor:farfield-v48",
            "result_sha256": "1" * 64, "accepted_result_sha256": "1" * 64,
        },
        THERMAL: {
            "generation": thermal_generation,
            **{key: thermal_generation for key in ("loss_mapping_generation", "source_mesh_generation", "target_mesh_generation", "interpolation_generation", "time_average_generation", "frequency_generation", "task_generation", "result_generation")},
            "loss_component_ids": components, "mapped_loss_component_ids": components,
            "loss_w": losses, "mapped_loss_w": losses,
            "source_mesh_sha256": "2" * 64, "mapped_source_mesh_sha256": "2" * 64,
            "target_mesh_sha256": "3" * 64, "mapped_target_mesh_sha256": "3" * 64,
            "interpolation_method": "conservative_nodal", "mapped_interpolation_method": "conservative_nodal",
            "time_average": "cycle_average", "mapped_time_average": "cycle_average",
            "frequency_hz": 2.45e9, "mapped_frequency_hz": 2.45e9,
            "task_owner": "task:em-thermal-v48", "mapped_task_owner": "task:em-thermal-v48",
            "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
        },
    }]}


# --- network_artifact_identity_v49 (v49) ------------------------------------


def _payload_v49() -> dict[str, object]:
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


# --- network_artifact_identity_v50 (v50) ------------------------------------


def _payload_v50() -> dict[str, object]:
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
