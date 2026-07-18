from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v34 import _summary_v34


_PROMOTED_CASE_IDS = (
    "v35_public_waveguide_port_mode_power_orthogonality_impedance_deembed_cutoff_mismatch",
    "v35_public_nearfar_directivity_gain_efficiency_polarization_sphere_power_mismatch",
)


def _summary_v35():
    summary = _summary_v34()
    c0 = 299_792_458.0
    cutoff = 6.5e9
    frequencies = [8.0e9, 10.0e9, 12.0e9]
    impedance = [377.0 / math.sqrt(1.0 - (cutoff / frequency) ** 2) for frequency in frequencies]
    beta = [2.0 * math.pi / c0 * math.sqrt(frequency**2 - cutoff**2) for frequency in frequencies]
    for index, row in enumerate(summary["runs"]):
        generation = f"waveguide-port-{411 + index}"
        row[
            "waveguide_port_mode_power_orthogonality_impedance_deembed_cutoff_frequency_owner_result_identity"
        ] = {
            "waveguide_port_generation": generation,
            **{
                key: generation
                for key in (
                    "mode_generation", "power_generation", "orthogonality_generation",
                    "impedance_generation", "deembed_generation", "cutoff_generation",
                    "frequency_generation", "mesh_generation", "owner_generation",
                    "result_generation",
                )
            },
            "mode_name": "TE10", "result_mode_name": "TE10",
            "normalization": "accepted_power_1w", "result_normalization": "accepted_power_1w",
            "modal_power_w": [1.0, 1.0, 1.0], "result_modal_power_w": [1.0, 1.0, 1.0],
            "mode_gram_real": [[1.0, 0.0], [0.0, 1.0]],
            "result_mode_gram_real": [[1.0, 0.0], [0.0, 1.0]],
            "mode_gram_imag": [[0.0, 0.0], [0.0, 0.0]],
            "result_mode_gram_imag": [[0.0, 0.0], [0.0, 0.0]],
            "impedance_definition": "te_wave_impedance",
            "result_impedance_definition": "te_wave_impedance",
            "modal_impedance_ohm": impedance, "result_modal_impedance_ohm": list(impedance),
            "frequency_grid_hz": frequencies, "result_frequency_grid_hz": list(frequencies),
            "cutoff_frequency_hz": cutoff, "result_cutoff_frequency_hz": cutoff,
            "propagation_constant_rad_m": beta, "result_propagation_constant_rad_m": list(beta),
            "reference_plane_m": 0.0, "result_reference_plane_m": 0.0,
            "deembedded_reference_plane_m": 0.01, "result_deembedded_reference_plane_m": 0.01,
            "deembed_phase_rad": [-item * 0.01 for item in beta],
            "result_deembed_phase_rad": [-item * 0.01 for item in beta],
            "port_mesh_sha256": "1" * 64, "result_port_mesh_sha256": "1" * 64,
            "port_owner": "waveguide-port/case-411",
            "accepted_port_owner": "waveguide-port/case-411",
            "port_result_sha256": "2" * 64, "accepted_port_result_sha256": "2" * 64,
        }
        generation = f"nearfar-{411 + index}"
        weights = [math.pi / 2.0] * 8
        intensity = [8.0 / math.pi, 8.0 / math.pi] + [0.0] * 6
        row[
            "nearfar_sphere_power_directivity_gain_efficiency_polarization_quadrature_mesh_owner_result_identity"
        ] = {
            "nearfar_generation": generation,
            **{
                key: generation
                for key in (
                    "sphere_generation", "power_generation", "directivity_generation",
                    "gain_generation", "efficiency_generation", "polarization_generation",
                    "quadrature_generation", "mesh_generation", "owner_generation",
                    "result_generation",
                )
            },
            "frequency_hz": 10.0e9, "result_frequency_hz": 10.0e9,
            "accepted_power_w": 10.0, "result_accepted_power_w": 10.0,
            "enclosing_sphere_power_w": 8.0, "result_enclosing_sphere_power_w": 8.0,
            "radiated_power_w": 8.0, "result_radiated_power_w": 8.0,
            "radiation_efficiency": 0.8, "result_radiation_efficiency": 0.8,
            "maximum_directivity_linear": 4.0, "result_maximum_directivity_linear": 4.0,
            "realized_gain_linear": 3.2, "result_realized_gain_linear": 3.2,
            "polarization_basis": "theta_phi_right_handed",
            "result_polarization_basis": "theta_phi_right_handed",
            "copolar_definition": "ludwig3", "result_copolar_definition": "ludwig3",
            "angular_quadrature_weights_sr": weights,
            "result_angular_quadrature_weights_sr": list(weights),
            "radiation_intensity_w_sr": intensity,
            "result_radiation_intensity_w_sr": list(intensity),
            "farfield_mesh_sha256": "3" * 64, "result_farfield_mesh_sha256": "3" * 64,
            "nearfar_owner": "nearfar/case-411", "accepted_nearfar_owner": "nearfar/case-411",
            "nearfar_result_sha256": "4" * 64, "accepted_nearfar_result_sha256": "4" * 64,
        }
    return summary


def test_v35_public_positive_waveguide_port_and_nearfar_closure():
    assert nonlinear_inductance_sweep_gate(_summary_v35())["status"] == "ok"


def test_v35_public_waveguide_port_mode_power_orthogonality_impedance_deembed_cutoff_mismatch():
    summary = _summary_v35()
    record = summary["runs"][0][
        "waveguide_port_mode_power_orthogonality_impedance_deembed_cutoff_frequency_owner_result_identity"
    ]
    record.update(
        {
            "mode_generation": "waveguide-port-410", "deembed_generation": "waveguide-port-409",
            "result_generation": "waveguide-port-408", "result_mode_name": "TM01",
            "result_normalization": "peak_field", "result_modal_power_w": [-1.0, 2.0],
            "result_mode_gram_real": [[1.0, 0.5], [0.5, 0.1]],
            "result_mode_gram_imag": [[0.0, 1.0], [-1.0, 0.0]],
            "result_impedance_definition": "lumped_50_ohm",
            "result_modal_impedance_ohm": [50.0, -50.0],
            "result_frequency_grid_hz": [12.0e9, 6.0e9],
            "result_cutoff_frequency_hz": 15.0e9,
            "result_propagation_constant_rad_m": [-1.0, -2.0],
            "result_deembedded_reference_plane_m": -0.01,
            "result_deembed_phase_rad": [1.0, 2.0],
            "result_port_mesh_sha256": "a" * 64,
            "accepted_port_owner": "waveguide-port/old",
            "accepted_port_result_sha256": "b" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "waveguide_port_modes_use_current_power_orthogonality_impedance_deembed_cutoff_frequency_mesh_owner_and_result"
    ]


def test_v35_public_nearfar_directivity_gain_efficiency_polarization_sphere_power_mismatch():
    summary = _summary_v35()
    record = summary["runs"][0][
        "nearfar_sphere_power_directivity_gain_efficiency_polarization_quadrature_mesh_owner_result_identity"
    ]
    record.update(
        {
            "power_generation": "nearfar-410", "polarization_generation": "nearfar-409",
            "result_generation": "nearfar-408", "result_frequency_hz": 9.0e9,
            "result_accepted_power_w": -10.0, "result_enclosing_sphere_power_w": 20.0,
            "result_radiated_power_w": 4.0, "result_radiation_efficiency": 1.5,
            "result_maximum_directivity_linear": -4.0, "result_realized_gain_linear": 9.0,
            "result_polarization_basis": "left_handed_xy", "result_copolar_definition": "unknown",
            "result_angular_quadrature_weights_sr": [1.0, -1.0],
            "result_radiation_intensity_w_sr": [100.0],
            "result_farfield_mesh_sha256": "c" * 64,
            "accepted_nearfar_owner": "nearfar/old",
            "accepted_nearfar_result_sha256": "d" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "nearfar_results_use_current_sphere_power_directivity_gain_efficiency_polarization_quadrature_mesh_owner_and_result"
    ]


def test_v35_public_rejects_self_consistent_nearfar_gain_over_efficiency():
    summary = _summary_v35()
    for row in summary["runs"]:
        record = row[
            "nearfar_sphere_power_directivity_gain_efficiency_polarization_quadrature_mesh_owner_result_identity"
        ]
        record["realized_gain_linear"] = 5.0
        record["result_realized_gain_linear"] = 5.0
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"
