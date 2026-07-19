from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from test_magnetic_force_generalization_v41_elf import _summary_v41


_DEMAG = (
    "demagnetizing_tensor_symmetry_trace_eigenvalue_reciprocity_energy_mesh_"
    "magnetization_result_identity"
)
_LINEAR = (
    "linear_motor_cogging_force_position_periodicity_work_coenergy_phase_"
    "thrust_mesh_result_identity"
)
_PROMOTED_CASE_IDS = (
    "v42_public_demagtensor_symmetry_trace_eigenvalue_reciprocity_energy_mesh_mismatch",
    "v42_public_linearmotor_cogging_force_position_periodicity_workenergy_phase_mismatch",
)
_BEARING = "magneticbearing_stiffnessmatrix_symmetry_crosscoupling_force_energy_stability_mesh_result_identity"
_HYSTERESIS = "hysteresis_minorloop_remanence_coercivity_loss_path_energy_material_mesh_result_identity"
_PROMOTED_V43_CASE_IDS = (
    "v43_public_magneticbearing_stiffnessmatrix_symmetry_crosscoupling_force_energy_stability_mismatch",
    "v43_public_hysteresis_minorloop_remanence_coercivity_loss_path_energy_mesh_mismatch",
)


def _matvec(matrix: list[list[float]], vector: list[float]) -> list[float]:
    return [sum(left * right for left, right in zip(row, vector)) for row in matrix]


def _dot(left: list[float], right: list[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def _summary_v42() -> dict:
    summary = _summary_v41()
    identity = summary["artifact_identity"]

    generation = "demag-tensor-842"
    tensor = [[0.2, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.5]]
    magnetization = [1.0e5, 2.0e5, 3.0e5]
    probe = [-2.0e5, 1.0e5, 0.5e5]
    tensor_m = _matvec(tensor, magnetization)
    tensor_probe = _matvec(tensor, probe)
    volume = 1.0e-6
    values = {
        "demag_tensor": tensor,
        "tensor_trace": 1.0,
        "tensor_eigenvalues": [0.2, 0.3, 0.5],
        "magnetization_a_per_m": magnetization,
        "probe_magnetization_a_per_m": probe,
        "demag_field_a_per_m": [-value for value in tensor_m],
        "reciprocity_left": _dot(magnetization, tensor_probe),
        "reciprocity_right": _dot(probe, tensor_m),
        "magnet_volume_m3": volume,
        "demag_energy_j": (
            0.5
            * 4.0e-7
            * math.pi
            * volume
            * _dot(magnetization, tensor_m)
        ),
    }
    identity[_DEMAG] = {
        "demag_tensor_generation": generation,
        **{
            key: generation
            for key in (
                "symmetry_generation", "trace_generation", "eigenvalue_generation",
                "reciprocity_generation", "energy_generation", "mesh_generation",
                "magnetization_generation", "result_generation",
            )
        },
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "mesh_owner": "mesh:demag-tensor-842",
        "accepted_mesh_owner": "mesh:demag-tensor-842",
        "magnetization_owner": "magnetization:demag-tensor-842",
        "accepted_magnetization_owner": "magnetization:demag-tensor-842",
        "demag_result_sha256": "1" * 64,
        "accepted_demag_result_sha256": "1" * 64,
    }

    generation = "linear-motor-period-842"
    period = 0.04
    positions = [0.0, 0.01, 0.02, 0.03, 0.04]
    wave_number = 2.0 * math.pi / period
    amplitude = 10.0 / wave_number
    current_peak = 5.0
    force = [
        amplitude * wave_number * math.sin(wave_number * position)
        for position in positions
    ]
    phase_currents = [
        [
            current_peak * math.sin(wave_number * position),
            current_peak * math.sin(wave_number * position - 2.0 * math.pi / 3.0),
            current_peak * math.sin(wave_number * position + 2.0 * math.pi / 3.0),
        ]
        for position in positions
    ]
    values = {
        "position_m": positions,
        "period_m": period,
        "phase_order": ["U", "V", "W"],
        "phase_currents_a": phase_currents,
        "current_peak_a": current_peak,
        "coenergy_amplitude_j": amplitude,
        "coenergy_j": [-amplitude * math.cos(wave_number * position) for position in positions],
        "cogging_force_n": force,
        "base_thrust_n": 100.0,
        "thrust_n": [100.0 + value for value in force],
        "periodic_work_j": 0.0,
    }
    identity[_LINEAR] = {
        "linear_motor_generation": generation,
        **{
            key: generation
            for key in (
                "position_generation", "periodicity_generation", "phase_generation",
                "force_generation", "thrust_generation", "work_generation",
                "coenergy_generation", "mesh_generation", "result_generation",
            )
        },
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "mesh_owner": "mesh:linear-motor-period-842",
        "accepted_mesh_owner": "mesh:linear-motor-period-842",
        "linear_motor_result_sha256": "2" * 64,
        "accepted_linear_motor_result_sha256": "2" * 64,
    }
    return summary


def test_v42_public_positive_demag_tensor_and_linear_motor_closure() -> None:
    assert magnetic_force_method_profile_gate(_summary_v42())["status"] == "ok"


def test_v42_public_demag_tensor_mismatch() -> None:
    summary = _summary_v42()
    summary["artifact_identity"][_DEMAG].update(
        {
            "symmetry_generation": "demag-tensor-841",
            "result_tensor_trace": 1.2,
            "result_tensor_eigenvalues": [-0.1, 0.3, 1.0],
            "result_reciprocity_right": -1.0,
            "result_demag_energy_j": -1.0,
            "accepted_mesh_owner": "stale:mesh",
        }
    )
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v42_public_linear_motor_mismatch() -> None:
    summary = _summary_v42()
    summary["artifact_identity"][_LINEAR].update(
        {
            "position_generation": "linear-motor-period-841",
            "result_position_m": [0.04, 0.03, 0.02, 0.01, 0.0],
            "result_phase_currents_a": [[5.0, 5.0, 5.0]],
            "result_periodic_work_j": 1.0,
            "accepted_mesh_owner": "stale:mesh",
        }
    )
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v42_public_rejects_self_consistent_wrong_demag_trace() -> None:
    summary = _summary_v42()
    row = summary["artifact_identity"][_DEMAG]
    tensor = [[0.2, 0.0, 0.0], [0.0, 0.3, 0.0], [0.0, 0.0, 0.6]]
    row["demag_tensor"] = row["result_demag_tensor"] = tensor
    row["tensor_trace"] = row["result_tensor_trace"] = 1.1
    row["tensor_eigenvalues"] = row["result_tensor_eigenvalues"] = [0.2, 0.3, 0.6]
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v42_public_rejects_self_consistent_wrong_cogging_force() -> None:
    summary = _summary_v42()
    row = summary["artifact_identity"][_LINEAR]
    row["cogging_force_n"] = row["result_cogging_force_n"] = [0.0] * 5
    row["thrust_n"] = row["result_thrust_n"] = [100.0] * 5
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def _summary_v43() -> dict:
    summary = _summary_v42()
    generation = "bearing-stiffness-843"
    summary["artifact_identity"][_BEARING] = {
        "bearing_generation": generation,
        **{key: generation for key in ("equilibrium_generation", "stiffness_generation", "symmetry_generation", "energy_generation", "stability_generation", "gap_generation", "mesh_generation", "result_generation")},
        "force_equilibrium_n": [0.0, 0.0], "result_force_equilibrium_n": [0.0, 0.0],
        "stiffness_matrix_n_per_m": [[100.0, -5.0], [-5.0, 80.0]],
        "result_stiffness_matrix_n_per_m": [[100.0, -5.0], [-5.0, 80.0]],
        "energy_curvature_n_per_m": 79.0, "result_energy_curvature_n_per_m": 79.0,
        "stability_sign": "stable", "result_stability_sign": "stable",
        "gap_m": 1.0e-3, "result_gap_m": 1.0e-3,
        "mesh_owner": "mesh:bearing-stiffness-843", "result_mesh_owner": "mesh:bearing-stiffness-843",
        "bearing_result_sha256": "5" * 64, "accepted_bearing_result_sha256": "5" * 64,
    }
    generation = "hyst-minor-843"
    summary["artifact_identity"][_HYSTERESIS] = {
        "hysteresis_generation": generation,
        **{key: generation for key in ("path_generation", "branch_generation", "remanence_generation", "coercivity_generation", "loss_generation", "energy_generation", "material_generation", "mesh_generation", "result_generation")},
        "field_path_a_per_m": [-1.0, 1.0, -1.0], "result_field_path_a_per_m": [-1.0, 1.0, -1.0],
        "magnetization_a_per_m": [-0.8, 1.0, -0.8], "result_magnetization_a_per_m": [-0.8, 1.0, -0.8],
        "remanence_a_per_m": 0.8, "result_remanence_a_per_m": 0.8,
        "coercivity_a_per_m": 0.4, "result_coercivity_a_per_m": 0.4,
        "loop_area_loss_j_per_m3": 0.16, "result_loop_area_loss_j_per_m3": 0.16,
        "cycle_energy_j": 0.16, "result_cycle_energy_j": 0.16,
        "material_owner": "material:hyst-minor-843", "result_material_owner": "material:hyst-minor-843",
        "mesh_owner": "mesh:hyst-minor-843", "result_mesh_owner": "mesh:hyst-minor-843",
        "hysteresis_result_sha256": "6" * 64, "accepted_hysteresis_result_sha256": "6" * 64,
    }
    return summary


def test_v43_public_positive_bearing_and_hysteresis_closure() -> None:
    assert magnetic_force_method_profile_gate(_summary_v43())["status"] == "ok"


def test_v43_public_bearing_mismatch() -> None:
    summary = _summary_v43()
    summary["artifact_identity"][_BEARING]["result_stability_sign"] = "unstable"
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["magnetic_bearings_close_stiffness_symmetry_crosscoupling_force_energy_stability_gap_mesh_and_result"]


def test_v43_public_hysteresis_mismatch() -> None:
    summary = _summary_v43()
    summary["artifact_identity"][_HYSTERESIS]["result_remanence_a_per_m"] = -0.8
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["hysteresis_minorloops_close_field_path_branch_remanence_coercivity_loss_energy_material_mesh_and_result"]
