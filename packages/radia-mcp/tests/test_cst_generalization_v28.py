from __future__ import annotations

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v27 import _summary_v27


_PROMOTED_CASE_IDS = (
    "v28_public_waveguide_port_mode_cutoff_normalization_reference_plane_mesh_field_result_mismatch",
    "v28_public_wake_impedance_bunch_profile_time_grid_frequency_transform_normalization_result_mismatch",
)


def _summary_v28():
    summary = _summary_v27()
    for index, row in enumerate(summary["runs"]):
        generation = f"waveguide-port-{321 + index}"
        row[
            "waveguide_port_mode_cutoff_normalization_reference_plane_mesh_field_result_generation_identity"
        ] = {
            "port_generation": generation,
            "mode_port_generation": generation,
            "cutoff_port_generation": generation,
            "normalization_port_generation": generation,
            "reference_plane_port_generation": generation,
            "mesh_port_generation": generation,
            "field_port_generation": generation,
            "result_port_generation": generation,
            "port_id": "P1",
            "result_port_id": "P1",
            "mode_id": "TE10",
            "result_mode_id": "TE10",
            "cutoff_frequency_hz": 6.557e9,
            "result_cutoff_frequency_hz": 6.557e9,
            "evaluation_frequency_hz": 10.0e9,
            "result_evaluation_frequency_hz": 10.0e9,
            "normalization": "unit-power-wave",
            "result_normalization": "unit-power-wave",
            "reference_plane_m": 0.005,
            "result_reference_plane_m": 0.005,
            "port_mesh_sha256": "1" * 64,
            "result_port_mesh_sha256": "1" * 64,
            "field_eigenvector_sha256": "2" * 64,
            "result_field_eigenvector_sha256": "2" * 64,
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        }
        generation = f"wake-impedance-{321 + index}"
        row[
            "wake_impedance_bunch_profile_time_grid_frequency_transform_normalization_mesh_result_generation_identity"
        ] = {
            "wake_generation": generation,
            "bunch_wake_generation": generation,
            "time_wake_generation": generation,
            "transform_wake_generation": generation,
            "frequency_wake_generation": generation,
            "normalization_wake_generation": generation,
            "mesh_wake_generation": generation,
            "result_wake_generation": generation,
            "bunch_profile": "gaussian",
            "result_bunch_profile": "gaussian",
            "bunch_sigma_s": 1.0e-12,
            "result_bunch_sigma_s": 1.0e-12,
            "bunch_charge_c": 1.0e-9,
            "result_bunch_charge_c": 1.0e-9,
            "time_grid_s": [0.0, 1.0e-12, 2.0e-12, 3.0e-12],
            "result_time_grid_s": [0.0, 1.0e-12, 2.0e-12, 3.0e-12],
            "wake_potential_v_c": [0.0, 2.0e12, 1.0e12, 0.0],
            "result_wake_potential_v_c": [0.0, 2.0e12, 1.0e12, 0.0],
            "fft_convention": "exp-minus-i-omega-t",
            "result_fft_convention": "exp-minus-i-omega-t",
            "frequency_grid_hz": [0.0, 1.0e9, 2.0e9],
            "result_frequency_grid_hz": [0.0, 1.0e9, 2.0e9],
            "impedance_normalization": "longitudinal-v-per-coulomb",
            "result_impedance_normalization": "longitudinal-v-per-coulomb",
            "mesh_sha256": "4" * 64,
            "result_mesh_sha256": "4" * 64,
            "result_sha256": "5" * 64,
            "accepted_result_sha256": "5" * 64,
        }
    return summary


def test_v28_public_positive_waveguide_port_and_wake_identities() -> None:
    assert nonlinear_inductance_sweep_gate(_summary_v28())["status"] == "ok"


def test_v28_public_rejects_waveguide_port_identity_mismatch() -> None:
    summary = _summary_v28()
    identity = summary["runs"][0][
        "waveguide_port_mode_cutoff_normalization_reference_plane_mesh_field_result_generation_identity"
    ]
    identity.update(
        {
            "mode_port_generation": "waveguide-port-320",
            "mesh_port_generation": "waveguide-port-319",
            "result_port_id": "P2",
            "result_mode_id": "TM11",
            "result_cutoff_frequency_hz": 12.0e9,
            "result_evaluation_frequency_hz": 9.0e9,
            "result_normalization": "unit-voltage-wave",
            "result_reference_plane_m": 0.0,
            "result_port_mesh_sha256": "a" * 64,
            "result_field_eigenvector_sha256": "b" * 64,
            "accepted_result_sha256": "c" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "waveguide_ports_use_current_mode_cutoff_normalization_plane_mesh_field_and_result"
    ]


def test_v28_public_rejects_wake_impedance_identity_mismatch() -> None:
    summary = _summary_v28()
    identity = summary["runs"][0][
        "wake_impedance_bunch_profile_time_grid_frequency_transform_normalization_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "bunch_wake_generation": "wake-impedance-320",
            "frequency_wake_generation": "wake-impedance-319",
            "result_bunch_profile": "rectangular",
            "result_bunch_sigma_s": 2.0e-12,
            "result_bunch_charge_c": 2.0e-9,
            "result_time_grid_s": [0.0, 2.0e-12, 5.0e-12],
            "result_wake_potential_v_c": [0.0, -1.0e12, 0.0],
            "result_fft_convention": "exp-plus-i-omega-t",
            "result_frequency_grid_hz": [0.0, 1.4e9, 3.0e9],
            "result_impedance_normalization": "transverse-ohm-per-metre",
            "result_mesh_sha256": "d" * 64,
            "accepted_result_sha256": "e" * 64,
        }
    )
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "wake_impedance_uses_current_bunch_time_grid_transform_frequency_normalization_mesh_and_result"
    ]
