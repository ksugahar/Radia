from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v41 import _summary_v41


_CAVITY = "cavity_eigenmode_frequency_qfactor_fieldenergy_orthogonality_modevolume_mesh_result_identity"
_TDR = "tdr_time_distance_impedance_reflection_deembed_risetime_loss_waveform_energy_result_identity"
_PROMOTED_CASE_IDS = (
    "v42_public_cavity_eigenmode_frequency_qfactor_fieldenergy_orthogonality_modevolume_mismatch",
    "v42_public_tdr_impedance_reflection_distance_deembed_risetime_loss_energy_mismatch",
)
C0 = 299_792_458.0


def _summary_v42() -> dict:
    summary = _summary_v41()
    for index, run in enumerate(summary["runs"]):
        generation = f"cavity-mode-842-{index}"
        dimensions = [0.1, 0.08, 0.05]
        modes = [[1, 0, 1], [0, 1, 1]]
        frequencies = [0.5 * C0 * math.sqrt(sum((mode[i] / dimensions[i]) ** 2 for i in range(3))) for mode in modes]
        electric, magnetic = [0.5, 0.4], [0.5, 0.4]
        values = {
            "cavity_dimensions_m": dimensions, "mode_indices": modes,
            "eigenfrequency_hz": frequencies, "unloaded_q": [10000.0, 8000.0],
            "electric_energy_j": electric, "magnetic_energy_j": magnetic,
            "total_field_energy_j": [a + b for a, b in zip(electric, magnetic)],
            "mode_overlap_matrix": [[1.0, 0.0], [0.0, 1.0]],
            "mode_volume_m3": [2.0e-4, 1.5e-4],
        }
        run[_CAVITY] = {
            "cavity_mode_generation": generation,
            **{key: generation for key in ("geometry_generation", "frequency_generation", "qfactor_generation", "energy_generation", "orthogonality_generation", "modevolume_generation", "mesh_generation", "result_generation")},
            **values, **{f"result_{key}": value for key, value in values.items()},
            "mesh_owner": "mesh:cavity-mode-842", "accepted_mesh_owner": "mesh:cavity-mode-842",
            "cavity_result_sha256": "1" * 64, "accepted_cavity_result_sha256": "1" * 64,
        }

        generation = f"tdr-line-842-{index}"
        permittivity, velocity = 4.0, C0 / 2.0
        times, deembed = [0.0, 1.0e-9, 2.0e-9], 0.2e-9
        impedance, reference = [50.0, 75.0, 50.0], 50.0
        reflection = [(value - reference) / (value + reference) for value in impedance]
        losses = [0.0, 0.1, 0.2]
        waveform = [(1.0 + gamma) * 10.0 ** (-loss / 20.0) for gamma, loss in zip(reflection, losses)]
        energy = sum(0.5 * (waveform[i] ** 2 + waveform[i + 1] ** 2) / reference * (times[i + 1] - times[i]) for i in range(2))
        values = {
            "time_s": times, "relative_permittivity": permittivity,
            "propagation_velocity_m_per_s": velocity, "deembed_time_s": deembed,
            "distance_m": [0.5 * velocity * max(0.0, time - deembed) for time in times],
            "reference_impedance_ohm": reference, "impedance_ohm": impedance,
            "reflection_coefficient": reflection, "source_rise_time_s": 0.1e-9,
            "line_loss_db": losses, "tdr_waveform_v": waveform, "waveform_energy_j": energy,
        }
        run[_TDR] = {
            "tdr_generation": generation,
            **{key: generation for key in ("time_generation", "distance_generation", "impedance_generation", "reflection_generation", "deembed_generation", "risetime_generation", "loss_generation", "energy_generation", "result_generation")},
            **values, **{f"result_{key}": value for key, value in values.items()},
            "tdr_result_sha256": "2" * 64, "accepted_tdr_result_sha256": "2" * 64,
        }
    return summary


def test_v42_public_positive_cavity_and_tdr_closure() -> None:
    assert nonlinear_inductance_sweep_gate(_summary_v42())["status"] == "ok"


def test_v42_public_cavity_mismatch() -> None:
    summary = _summary_v42()
    summary["runs"][0][_CAVITY].update({"frequency_generation": "cavity-mode-841", "result_eigenfrequency_hz": [-1.0], "result_mode_overlap_matrix": [[1.0, 1.0]], "accepted_mesh_owner": "stale:mesh"})
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v42_public_tdr_mismatch() -> None:
    summary = _summary_v42()
    summary["runs"][0][_TDR].update({"distance_generation": "tdr-line-841", "result_distance_m": [-1.0], "result_reflection_coefficient": [2.0], "result_waveform_energy_j": -1.0})
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v42_public_rejects_self_consistent_wrong_cavity_frequency() -> None:
    summary = _summary_v42()
    for run in summary["runs"]:
        row = run[_CAVITY]
        row["eigenfrequency_hz"] = row["result_eigenfrequency_hz"] = [1.0, 2.0]
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v42_public_rejects_self_consistent_wrong_tdr_reflection() -> None:
    summary = _summary_v42()
    for run in summary["runs"]:
        row = run[_TDR]
        row["reflection_coefficient"] = row["result_reflection_coefficient"] = [0.0, 0.0, 0.0]
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"
