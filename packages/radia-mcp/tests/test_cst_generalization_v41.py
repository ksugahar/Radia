from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import nonlinear_inductance_sweep_gate
from test_cst_generalization_v40 import _summary_v40


_WAVEGUIDE = "waveguide_cutoff_propagation_impedance_mode_orthogonality_sparameter_power_owner_result_identity"
_ANTENNA = "antenna_farfield_directivity_gain_efficiency_power_polarization_mesh_owner_result_identity"
_PROMOTED_CASE_IDS = (
    "v41_public_waveguide_cutoff_propagation_impedance_mode_orthogonality_sparameter_power_mismatch",
    "v41_public_antenna_farfield_directivity_gain_efficiency_power_polarization_mismatch",
)
C0 = 299_792_458.0
MU0 = 4.0e-7 * math.pi


def _summary_v41() -> dict:
    summary = _summary_v40()
    for index, run in enumerate(summary["runs"]):
        generation = f"waveguide-te10-{724 + index}"
        width, height, frequency = 22.86e-3, 10.16e-3, 10.0e9
        omega = 2.0 * math.pi * frequency
        cutoff = C0 / (2.0 * width)
        propagation = math.sqrt((omega / C0) ** 2 - (math.pi / width) ** 2)
        s11, s21 = [0.1, 0.0], [math.sqrt(0.97), 0.0]
        incident = 1.0
        reflected = incident * sum(item * item for item in s11)
        transmitted = incident * sum(item * item for item in s21)
        values = {
            "waveguide_width_m": width, "waveguide_height_m": height,
            "mode_name": "TE10", "frequency_hz": frequency,
            "cutoff_frequency_hz": cutoff,
            "propagation_constant_rad_per_m": propagation,
            "guide_wavelength_m": 2.0 * math.pi / propagation,
            "guide_impedance_ohm": omega * MU0 / propagation,
            "modal_overlap_matrix": [[1.0, 0.0], [0.0, 1.0]],
            "s11_complex": s11, "s21_complex": s21,
            "incident_power_w": incident, "reflected_power_w": reflected,
            "transmitted_power_w": transmitted,
            "wall_loss_w": incident - reflected - transmitted,
            "power_balance_residual_w": 0.0,
        }
        run[_WAVEGUIDE] = {
            "waveguide_generation": generation,
            **{key: generation for key in ("geometry_generation", "mode_generation", "cutoff_generation", "propagation_generation", "impedance_generation", "orthogonality_generation", "sparameter_generation", "power_generation", "owner_generation", "result_generation")},
            **values, **{f"result_{key}": value for key, value in values.items()},
            "waveguide_owner": f"waveguide/te10-{724 + index}",
            "accepted_waveguide_owner": f"waveguide/te10-{724 + index}",
            "waveguide_result_sha256": "5" * 64,
            "accepted_waveguide_result_sha256": "5" * 64,
        }

        generation = f"antenna-farfield-{724 + index}"
        incident, reflected = 10.0, 1.0
        accepted, conductor_loss, dielectric_loss = 9.0, 0.5, 0.5
        radiated = accepted - conductor_loss - dielectric_loss
        efficiency, directivity = radiated / accepted, 6.0
        values = {
            "incident_power_w": incident, "reflected_power_w": reflected,
            "accepted_power_w": accepted, "radiated_power_w": radiated,
            "conductor_loss_w": conductor_loss, "dielectric_loss_w": dielectric_loss,
            "radiation_efficiency": efficiency, "directivity_linear": directivity,
            "gain_linear": directivity * efficiency,
            "maximum_radiation_intensity_w_per_sr": directivity * radiated / (4.0 * math.pi),
            "polarization_basis": "linear_xy", "co_polar_fraction": 1.0,
            "cross_polar_fraction": 0.0, "power_balance_residual_w": 0.0,
        }
        run[_ANTENNA] = {
            "antenna_generation": generation,
            **{key: generation for key in ("excitation_generation", "farfield_generation", "directivity_generation", "gain_generation", "efficiency_generation", "power_generation", "polarization_generation", "mesh_generation", "owner_generation", "result_generation")},
            **values, **{f"result_{key}": value for key, value in values.items()},
            "antenna_owner": f"antenna/farfield-{724 + index}",
            "accepted_antenna_owner": f"antenna/farfield-{724 + index}",
            "antenna_result_sha256": "6" * 64,
            "accepted_antenna_result_sha256": "6" * 64,
        }
    return summary


def test_v41_public_positive_waveguide_and_antenna_closure() -> None:
    assert nonlinear_inductance_sweep_gate(_summary_v41())["status"] == "ok"


def test_v41_public_waveguide_cutoff_propagation_impedance_mode_orthogonality_sparameter_power_mismatch() -> None:
    summary = _summary_v41()
    summary["runs"][0][_WAVEGUIDE].update({"cutoff_generation": "waveguide-te10-723", "result_cutoff_frequency_hz": -1.0, "result_modal_overlap_matrix": [[1.0, 1.0], [1.0, 1.0]], "accepted_waveguide_owner": "waveguide/old"})
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v41_public_antenna_farfield_directivity_gain_efficiency_power_polarization_mismatch() -> None:
    summary = _summary_v41()
    summary["runs"][0][_ANTENNA].update({"farfield_generation": "antenna-farfield-723", "result_radiation_efficiency": 2.0, "result_gain_linear": 20.0, "result_polarization_basis": "unknown", "accepted_antenna_owner": "antenna/old"})
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_wrong_waveguide_impedance() -> None:
    summary = _summary_v41()
    for run in summary["runs"]:
        row = run[_WAVEGUIDE]
        row["guide_impedance_ohm"] = row["result_guide_impedance_ohm"] = 50.0
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_wrong_antenna_gain() -> None:
    summary = _summary_v41()
    for run in summary["runs"]:
        row = run[_ANTENNA]
        row["gain_linear"] = row["result_gain_linear"] = row["directivity_linear"]
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"
