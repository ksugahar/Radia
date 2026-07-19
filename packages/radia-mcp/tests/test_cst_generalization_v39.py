from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v38 import _summary_v38


_PROMOTED_CASE_IDS = (
    "v39_public_waveguide_mode_cutoff_impedance_power_orthogonality_deembed_mismatch",
    "v39_public_antenna_nearfar_directivity_gain_efficiency_polarization_power_mismatch",
)

_WAVEGUIDE_KEY = (
    "waveguide_mode_cutoff_impedance_power_orthogonality_propagation_deembed_"
    "mesh_owner_result_identity"
)
_ANTENNA_KEY = (
    "antenna_nearfar_directivity_gain_efficiency_polarization_power_mesh_owner_"
    "result_identity"
)
C0 = 299_792_458.0
ETA0 = 376.730313668


def _waveguide(index: int) -> dict:
    generation = f"waveguide-mode-{715 + index}"
    width, height, frequency = 0.02, 0.01, 12.0e9
    cutoff = C0 / (2.0 * width)
    factor = math.sqrt(1.0 - (cutoff / frequency) ** 2)
    impedance = ETA0 / factor
    beta = 2.0 * math.pi * frequency * factor / C0
    distance, raw_phase = 0.012, -1.2
    mirrored = {
        "waveguide_width_m": width,
        "waveguide_height_m": height,
        "mode_name": "TE10",
        "mode_index": [1, 0],
        "cutoff_frequency_hz": cutoff,
        "frequency_hz": frequency,
        "modal_impedance_ohm": impedance,
        "propagation_constant_rad_m": beta,
        "normalized_forward_power_w": 1.0,
        "mode_overlap_real": [[1.0, 0.0], [0.0, 1.0]],
        "mode_overlap_imag": [[0.0, 0.0], [0.0, 0.0]],
        "orthogonality_tolerance": 1.0e-9,
        "port_plane_m": 0.0,
        "reference_plane_m": distance,
        "deembed_distance_m": distance,
        "deembed_convention": "port_to_reference_add_beta_l",
        "raw_s21_phase_rad": raw_phase,
        "deembedded_s21_phase_rad": raw_phase + beta * distance,
        "waveguide_mesh_sha256": "1" * 64,
    }
    return {
        "waveguide_mode_generation": generation,
        **{
            key: generation
            for key in (
                "cutoff_generation",
                "impedance_generation",
                "power_generation",
                "orthogonality_generation",
                "propagation_generation",
                "deembed_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        },
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "waveguide_mode_owner": f"waveguide/mode-{715 + index}",
        "accepted_waveguide_mode_owner": f"waveguide/mode-{715 + index}",
        "waveguide_mode_result_sha256": "2" * 64,
        "accepted_waveguide_mode_result_sha256": "2" * 64,
    }


def _antenna(index: int) -> dict:
    generation = f"antenna-nearfar-{715 + index}"
    directivity, efficiency, mismatch = 4.0, 0.8, 0.9
    gain = directivity * efficiency
    realized = gain * mismatch
    mirrored = {
        "frequency_hz": 2.4e9,
        "nearfield_surface_closed": True,
        "near_to_far_transform": "equivalence_surface_stratton_chu",
        "polarization_basis": "ieee_theta_phi",
        "co_polar_component": "theta",
        "cross_polar_component": "phi",
        "directivity_linear": directivity,
        "directivity_dbi": 10.0 * math.log10(directivity),
        "radiation_efficiency": efficiency,
        "gain_linear": gain,
        "mismatch_efficiency": mismatch,
        "realized_gain_linear": realized,
        "realized_gain_dbi": 10.0 * math.log10(realized),
        "accepted_power_w": 1.0,
        "radiated_power_w": efficiency,
        "loss_power_w": 1.0 - efficiency,
        "power_balance_residual_w": 0.0,
        "farfield_sphere_samples": 2592,
        "antenna_mesh_sha256": "3" * 64,
    }
    return {
        "antenna_generation": generation,
        **{
            key: generation
            for key in (
                "nearfield_generation",
                "farfield_generation",
                "directivity_generation",
                "gain_generation",
                "efficiency_generation",
                "polarization_generation",
                "power_generation",
                "mesh_generation",
                "owner_generation",
                "result_generation",
            )
        },
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "antenna_owner": f"antenna/nearfar-{715 + index}",
        "accepted_antenna_owner": f"antenna/nearfar-{715 + index}",
        "antenna_result_sha256": "4" * 64,
        "accepted_antenna_result_sha256": "4" * 64,
    }


def _summary_v39() -> dict:
    summary = _summary_v38()
    for index, row in enumerate(summary["runs"]):
        row[_WAVEGUIDE_KEY] = _waveguide(index)
        row[_ANTENNA_KEY] = _antenna(index)
    return summary


def test_v39_public_positive_waveguide_and_antenna_closure() -> None:
    assert nonlinear_inductance_sweep_gate(_summary_v39())["status"] == "ok"


def test_v39_public_waveguide_mode_cutoff_impedance_power_orthogonality_deembed_mismatch() -> None:
    summary = _summary_v39()
    row = summary["runs"][0][_WAVEGUIDE_KEY]
    row.update(
        {
            "cutoff_generation": "waveguide-mode-714",
            "deembed_generation": "waveguide-mode-713",
            "result_generation": "waveguide-mode-712",
            "result_cutoff_frequency_hz": -1.0,
            "result_modal_impedance_ohm": -50.0,
            "result_normalized_forward_power_w": -1.0,
            "result_mode_overlap_real": [[1.0, 1.0], [1.0, 1.0]],
            "result_propagation_constant_rad_m": -1.0,
            "result_deembed_distance_m": -0.01,
            "result_deembedded_s21_phase_rad": 9.0,
            "accepted_waveguide_mode_owner": "waveguide/old",
            "accepted_waveguide_mode_result_sha256": "a" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "waveguide_modes_use_current_cutoff_impedance_power_orthogonality_propagation_deembed_mesh_owner_and_result"
    ]


def test_v39_public_antenna_nearfar_directivity_gain_efficiency_polarization_power_mismatch() -> None:
    summary = _summary_v39()
    row = summary["runs"][0][_ANTENNA_KEY]
    row.update(
        {
            "farfield_generation": "antenna-nearfar-714",
            "power_generation": "antenna-nearfar-713",
            "result_generation": "antenna-nearfar-712",
            "result_nearfield_surface_closed": False,
            "result_near_to_far_transform": "unknown",
            "result_directivity_linear": -1.0,
            "result_realized_gain_linear": 9.0,
            "result_radiation_efficiency": 2.0,
            "result_polarization_basis": "left_handed_local",
            "result_accepted_power_w": -1.0,
            "result_radiated_power_w": 2.0,
            "result_power_balance_residual_w": 3.0,
            "accepted_antenna_owner": "antenna/old",
            "accepted_antenna_result_sha256": "b" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "antenna_nearfar_uses_current_transform_directivity_gain_efficiency_polarization_power_mesh_owner_and_result"
    ]


def test_v39_public_rejects_self_consistent_wrong_waveguide_impedance() -> None:
    summary = _summary_v39()
    for run in summary["runs"]:
        row = run[_WAVEGUIDE_KEY]
        row["modal_impedance_ohm"] *= 2.0
        row["result_modal_impedance_ohm"] = row["modal_impedance_ohm"]
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v39_public_rejects_self_consistent_unknown_nearfar_transform() -> None:
    summary = _summary_v39()
    for run in summary["runs"]:
        row = run[_ANTENNA_KEY]
        row["near_to_far_transform"] = "unknown"
        row["result_near_to_far_transform"] = "unknown"
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"
