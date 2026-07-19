from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import nonlinear_inductance_sweep_gate
from test_cst_generalization_v39 import _summary_v39


_VIA = "pcb_via_stub_geometry_impedance_resonance_sparameter_current_loss_power_mesh_owner_result_identity"
_CAVITY = "cavity_degenerate_mode_frequency_orthogonality_energy_symmetry_quality_mesh_owner_result_identity"
_PROMOTED_CASE_IDS = (
    "v40_public_pcb_via_stub_impedance_resonance_sparameter_current_loss_power_mismatch",
    "v40_public_cavity_degenerate_modes_frequency_orthogonality_energy_symmetry_q_mismatch",
)
C0 = 299_792_458.0


def _via(index: int) -> dict:
    generation = f"pcb-via-{724 + index}"
    board, diameter, stub, epsilon_r = 1.6e-3, 0.30e-3, 0.80e-3, 4.0
    impedance = 60.0 / math.sqrt(epsilon_r) * math.log(4.0 * board / diameter)
    resonance = C0 / (4.0 * stub * math.sqrt(epsilon_r))
    s11, s21 = [0.2, -0.1], [0.85, -0.05]
    incident = 1.0
    reflected = incident * sum(item * item for item in s11)
    accepted = incident - reflected
    transmitted = incident * sum(item * item for item in s21)
    conductor_loss = 0.12
    values = {
        "board_thickness_m": board, "via_diameter_m": diameter,
        "stub_length_m": stub, "relative_permittivity": epsilon_r,
        "characteristic_impedance_ohm": impedance,
        "quarterwave_resonance_hz": resonance, "s11_complex": s11,
        "s21_complex": s21, "incident_power_w": incident,
        "reflected_power_w": reflected, "accepted_power_w": accepted,
        "transmitted_power_w": transmitted,
        "barrel_current_rms_a": math.sqrt(accepted / impedance),
        "conductor_loss_w": conductor_loss,
        "dielectric_loss_w": accepted - transmitted - conductor_loss,
        "power_balance_residual_w": 0.0, "pcb_mesh_sha256": "1" * 64,
    }
    return {
        "via_generation": generation,
        **{key: generation for key in ("geometry_generation", "impedance_generation", "resonance_generation", "sparameter_generation", "current_generation", "loss_generation", "power_generation", "mesh_generation", "owner_generation", "result_generation")},
        **values, **{f"result_{key}": value for key, value in values.items()},
        "via_owner": f"pcb/via-{724 + index}",
        "accepted_via_owner": f"pcb/via-{724 + index}",
        "via_result_sha256": "2" * 64, "accepted_via_result_sha256": "2" * 64,
    }


def _cavity(index: int) -> dict:
    generation = f"cavity-degenerate-{724 + index}"
    side = 0.03
    frequency = C0 / 2.0 * math.sqrt(2.0 / side**2)
    values = {
        "cavity_dimensions_m": [side, side, side],
        "mode_names": ["TE101", "TE011"],
        "mode_frequencies_hz": [frequency, frequency],
        "degeneracy_tolerance_relative": 1.0e-9,
        "modal_overlap_matrix": [[1.0, 0.0], [0.0, 1.0]],
        "electric_energy_j": [0.5, 0.5], "magnetic_energy_j": [0.5, 0.5],
        "symmetry_classes": ["even_x_odd_y", "odd_x_even_y"],
        "quality_factors": [5000.0, 5000.0], "cavity_mesh_sha256": "3" * 64,
    }
    return {
        "cavity_generation": generation,
        **{key: generation for key in ("frequency_generation", "orthogonality_generation", "energy_generation", "symmetry_generation", "quality_generation", "mesh_generation", "owner_generation", "result_generation")},
        **values, **{f"result_{key}": value for key, value in values.items()},
        "cavity_owner": f"cavity/degenerate-{724 + index}",
        "accepted_cavity_owner": f"cavity/degenerate-{724 + index}",
        "cavity_result_sha256": "4" * 64, "accepted_cavity_result_sha256": "4" * 64,
    }


def _summary_v40() -> dict:
    summary = _summary_v39()
    for index, run in enumerate(summary["runs"]):
        run[_VIA] = _via(index)
        run[_CAVITY] = _cavity(index)
    return summary


def test_v40_public_positive_via_and_cavity_closure() -> None:
    assert nonlinear_inductance_sweep_gate(_summary_v40())["status"] == "ok"


def test_v40_public_pcb_via_stub_impedance_resonance_sparameter_current_loss_power_mismatch() -> None:
    summary = _summary_v40()
    summary["runs"][0][_VIA].update({"geometry_generation": "pcb-via-723", "result_characteristic_impedance_ohm": -1.0, "result_s11_complex": [2.0, 0.0], "result_accepted_power_w": -1.0, "accepted_via_owner": "pcb/old"})
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v40_public_cavity_degenerate_modes_frequency_orthogonality_energy_symmetry_q_mismatch() -> None:
    summary = _summary_v40()
    summary["runs"][0][_CAVITY].update({"frequency_generation": "cavity-degenerate-723", "result_mode_frequencies_hz": [1.0, 2.0], "result_modal_overlap_matrix": [[1.0, 1.0], [1.0, 1.0]], "result_quality_factors": [-1.0, -1.0], "accepted_cavity_owner": "cavity/old"})
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_via_power_gap() -> None:
    summary = _summary_v40()
    for run in summary["runs"]:
        row = run[_VIA]
        row["dielectric_loss_w"] = row["result_dielectric_loss_w"] = 0.0
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_cavity_energy_imbalance() -> None:
    summary = _summary_v40()
    for run in summary["runs"]:
        row = run[_CAVITY]
        row["magnetic_energy_j"] = row["result_magnetic_energy_j"] = [0.4, 0.4]
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"
