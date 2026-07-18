from __future__ import annotations

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v30 import _summary_v30

_PROMOTED_CASE_IDS = (
    "v31_public_wave_port_modal_power_impedance_deembed_phase_reference_balance_mismatch",
    "v31_public_farfield_spherical_basis_handedness_polarization_phase_radiated_power_mismatch",
)


def _summary_v31():
    summary = _summary_v30()
    for index, row in enumerate(summary["runs"]):
        generation = f"wave-port-reference-{361 + index}"
        row["wave_port_modal_power_impedance_deembed_phase_balance_result_identity"] = {
            "port_generation": generation,
            "mode_port_generation": generation,
            "power_port_generation": generation,
            "impedance_port_generation": generation,
            "deembed_port_generation": generation,
            "phase_port_generation": generation,
            "owner_port_generation": generation,
            "balance_port_generation": generation,
            "result_port_generation": generation,
            "mode_ids": ["port1:TE10", "port2:TE10"],
            "result_mode_ids": ["port1:TE10", "port2:TE10"],
            "modal_power_normalization_w": [1.0, 1.0],
            "result_modal_power_normalization_w": [1.0, 1.0],
            "reference_impedance_ohm": [50.0, 50.0],
            "result_reference_impedance_ohm": [50.0, 50.0],
            "deembed_plane_m": [0.0, 0.1],
            "result_deembed_plane_m": [0.0, 0.1],
            "phase_reference_rad": [0.0, 0.0],
            "result_phase_reference_rad": [0.0, 0.0],
            "port_mode_owner_ids": ["project-361:port-1:TE10", "project-361:port-2:TE10"],
            "result_port_mode_owner_ids": ["project-361:port-1:TE10", "project-361:port-2:TE10"],
            "incident_power_w": 1.0,
            "reflected_power_w": 0.04,
            "transmitted_power_w": 0.94,
            "dissipated_power_w": 0.02,
            "result_power_balance_w": 1.0,
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        }
        generation = f"farfield-basis-{361 + index}"
        row["farfield_spherical_basis_handedness_polarization_phase_power_result_identity"] = {
            "farfield_generation": generation,
            "basis_farfield_generation": generation,
            "handedness_farfield_generation": generation,
            "order_farfield_generation": generation,
            "polarization_farfield_generation": generation,
            "phase_farfield_generation": generation,
            "weights_farfield_generation": generation,
            "power_farfield_generation": generation,
            "owner_farfield_generation": generation,
            "result_farfield_generation": generation,
            "spherical_basis": "e_theta_e_phi",
            "result_spherical_basis": "e_theta_e_phi",
            "coordinate_handedness": "right_handed",
            "result_coordinate_handedness": "right_handed",
            "angular_order": "theta_major_phi_minor",
            "result_angular_order": "theta_major_phi_minor",
            "polarization_phase_convention": "exp_plus_j_phase",
            "result_polarization_phase_convention": "exp_plus_j_phase",
            "theta_weights": [0.25, 0.5, 0.25],
            "result_theta_weights": [0.25, 0.5, 0.25],
            "phi_weights": [0.5, 0.5],
            "result_phi_weights": [0.5, 0.5],
            "radiated_power_w": 0.94,
            "integrated_radiated_power_w": 0.94,
            "farfield_owner_id": "project-361:monitor-ff-1",
            "accepted_farfield_owner_id": "project-361:monitor-ff-1",
            "result_sha256": "2" * 64,
            "accepted_result_sha256": "2" * 64,
        }
    return summary


def test_v31_public_positive_wave_port_and_farfield_identities():
    assert nonlinear_inductance_sweep_gate(_summary_v31())["status"] == "ok"


def test_v31_public_wave_port_modal_power_impedance_deembed_phase_reference_balance_mismatch():
    summary = _summary_v31()
    identity = summary["runs"][0][
        "wave_port_modal_power_impedance_deembed_phase_balance_result_identity"
    ]
    identity.update(
        {
            "phase_port_generation": "wave-port-reference-360",
            "balance_port_generation": "wave-port-reference-359",
            "result_mode_ids": ["port1:TM01", "port2:TE10"],
            "result_modal_power_normalization_w": [0.5, 2.0],
            "result_reference_impedance_ohm": [75.0, 50.0],
            "result_deembed_plane_m": [0.02, 0.08],
            "result_phase_reference_rad": [3.141592653589793, 0.0],
            "result_port_mode_owner_ids": ["old:port-1", "old:port-2"],
            "result_power_balance_w": 1.2,
            "accepted_result_sha256": "9" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "wave_ports_use_current_modal_power_impedance_deembed_phase_owner_balance_and_result"
    ]


def test_v31_public_farfield_spherical_basis_handedness_polarization_phase_radiated_power_mismatch():
    summary = _summary_v31()
    identity = summary["runs"][0][
        "farfield_spherical_basis_handedness_polarization_phase_power_result_identity"
    ]
    identity.update(
        {
            "basis_farfield_generation": "farfield-basis-360",
            "power_farfield_generation": "farfield-basis-359",
            "result_spherical_basis": "e_phi_e_theta",
            "result_coordinate_handedness": "left_handed",
            "result_angular_order": "phi_major_theta_minor",
            "result_polarization_phase_convention": "exp_minus_j_phase",
            "result_theta_weights": [1.0, -1.0, 1.0],
            "result_phi_weights": [1.0],
            "integrated_radiated_power_w": 1.25,
            "accepted_farfield_owner_id": "old-project:old-monitor",
            "accepted_result_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "farfields_use_current_spherical_basis_handedness_order_polarization_weights_power_owner_and_result"
    ]
