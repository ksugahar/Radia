from __future__ import annotations

import math

from test_magnetic_force_generalization_v31 import _gate
from test_magnetic_force_generalization_v40 import _identity_v40


_MAGNET = (
    "permanent_magnet_recoil_loadline_operating_demag_energy_virtualwork_"
    "force_mesh_result_generation_identity"
)
_CAPACITANCE = (
    "electrostatic_capacitance_matrix_symmetry_psd_charge_energy_reciprocity_"
    "mesh_result_generation_identity"
)
_PROMOTED_CASE_IDS = (
    "v41_public_permanentmagnet_recoilline_loadline_operatingpoint_demag_energy_force_mismatch",
    "v41_public_electrostatic_capacitancematrix_symmetry_psd_charge_energy_reciprocity_mismatch",
)


def _identity_v41():
    identity = _identity_v40()
    generation = "pm-loadline-724"
    mu0 = 4.0e-7 * math.pi
    recoil_mu_r = 1.05
    remanence = 1.2
    coercive_field = remanence / (mu0 * recoil_mu_r)
    permeance = 2.0
    operating_h = -remanence / (mu0 * (recoil_mu_r + permeance))
    operating_b = remanence + mu0 * recoil_mu_r * operating_h
    demag_knee_h = -8.0e5
    demag_margin = abs(demag_knee_h) - abs(operating_h)
    volume = 1.0e-6
    energy = 0.5 * operating_b * operating_b * volume / mu0
    displacement = 1.0e-4
    energy_minus = energy - 1.0e-3
    energy_plus = energy + 1.0e-3
    force = -(energy_plus - energy_minus) / (2.0 * displacement)
    identity[_MAGNET] = {
        "magnet_generation": generation,
        **{
            key: generation
            for key in (
                "recoil_generation",
                "loadline_generation",
                "operating_generation",
                "demag_generation",
                "energy_generation",
                "force_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "recoil_relative_permeability": recoil_mu_r,
        "result_recoil_relative_permeability": recoil_mu_r,
        "remanence_t": remanence,
        "result_remanence_t": remanence,
        "coercive_field_a_per_m": coercive_field,
        "result_coercive_field_a_per_m": coercive_field,
        "loadline_permeance_coefficient": permeance,
        "result_loadline_permeance_coefficient": permeance,
        "operating_h_a_per_m": operating_h,
        "result_operating_h_a_per_m": operating_h,
        "operating_b_t": operating_b,
        "result_operating_b_t": operating_b,
        "demag_knee_h_a_per_m": demag_knee_h,
        "result_demag_knee_h_a_per_m": demag_knee_h,
        "demag_margin_a_per_m": demag_margin,
        "result_demag_margin_a_per_m": demag_margin,
        "magnet_volume_m3": volume,
        "result_magnet_volume_m3": volume,
        "field_energy_j": energy,
        "result_field_energy_j": energy,
        "virtual_work_displacement_m": displacement,
        "result_virtual_work_displacement_m": displacement,
        "energy_minus_j": energy_minus,
        "result_energy_minus_j": energy_minus,
        "energy_plus_j": energy_plus,
        "result_energy_plus_j": energy_plus,
        "virtual_work_force_n": force,
        "result_virtual_work_force_n": force,
        "mesh_owner": "magnetics:pm-mesh-724",
        "accepted_mesh_owner": "magnetics:pm-mesh-724",
        "magnet_result_sha256": "5" * 64,
        "accepted_magnet_result_sha256": "5" * 64,
    }

    generation = "capacitance-matrix-724"
    matrix = [[2.0e-9, -0.5e-9], [-0.5e-9, 1.5e-9]]
    voltages = [10.0, 0.0]
    charges = [
        sum(matrix[row][column] * voltages[column] for column in range(2))
        for row in range(2)
    ]
    energy = 0.5 * sum(voltages[index] * charges[index] for index in range(2))
    identity[_CAPACITANCE] = {
        "capacitance_generation": generation,
        **{
            key: generation
            for key in (
                "matrix_generation",
                "symmetry_generation",
                "psd_generation",
                "charge_generation",
                "energy_generation",
                "reciprocity_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "conductor_names": ["electrode_1", "electrode_2"],
        "result_conductor_names": ["electrode_1", "electrode_2"],
        "capacitance_matrix_f": matrix,
        "result_capacitance_matrix_f": matrix,
        "drive_voltage_v": voltages,
        "result_drive_voltage_v": voltages,
        "conductor_charge_c": charges,
        "result_conductor_charge_c": charges,
        "stored_energy_j": energy,
        "result_stored_energy_j": energy,
        "symmetry_residual_f": 0.0,
        "result_symmetry_residual_f": 0.0,
        "reciprocity_residual_f": 0.0,
        "result_reciprocity_residual_f": 0.0,
        "mesh_owner": "electrostatics:cap-mesh-724",
        "accepted_mesh_owner": "electrostatics:cap-mesh-724",
        "capacitance_result_sha256": "6" * 64,
        "accepted_capacitance_result_sha256": "6" * 64,
    }
    return identity


def test_v41_public_positive_magnet_and_capacitance_closure():
    assert _gate(_identity_v41())["status"] == "ok"


def test_v41_public_permanent_magnet_mismatch():
    identity = _identity_v41()
    identity[_MAGNET].update(
        {
            "recoil_generation": "pm-loadline-723",
            "force_generation": "pm-loadline-722",
            "result_generation": "pm-loadline-721",
            "result_recoil_relative_permeability": -1.05,
            "result_loadline_permeance_coefficient": -2.0,
            "result_operating_h_a_per_m": 1.0e6,
            "result_operating_b_t": -1.0,
            "result_demag_margin_a_per_m": -1.0,
            "result_field_energy_j": -1.0,
            "result_virtual_work_force_n": 10.0,
            "accepted_mesh_owner": "stale:mesh",
            "accepted_magnet_result_sha256": "9" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "permanent_magnets_close_recoil_loadline_operating_demag_energy_virtualwork_mesh_and_result"
    ]


def test_v41_public_capacitance_matrix_mismatch():
    identity = _identity_v41()
    identity[_CAPACITANCE].update(
        {
            "matrix_generation": "capacitance-matrix-723",
            "energy_generation": "capacitance-matrix-722",
            "result_generation": "capacitance-matrix-721",
            "result_conductor_names": ["electrode_2", "electrode_1"],
            "result_capacitance_matrix_f": [[-2.0e-9, 1.0e-9], [0.0, 1.5e-9]],
            "result_conductor_charge_c": [0.0, 0.0],
            "result_stored_energy_j": -1.0,
            "result_symmetry_residual_f": 1.0,
            "result_reciprocity_residual_f": 1.0,
            "accepted_mesh_owner": "stale:mesh",
            "accepted_capacitance_result_sha256": "a" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "electrostatic_capacitance_matrices_close_symmetry_psd_charge_energy_reciprocity_mesh_and_result"
    ]


def test_v41_public_rejects_self_consistent_wrong_loadline_point():
    identity = _identity_v41()
    identity[_MAGNET]["operating_b_t"] = 0.6
    identity[_MAGNET]["result_operating_b_t"] = 0.6
    assert _gate(identity)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_indefinite_capacitance_matrix():
    identity = _identity_v41()
    matrix = [[1.0e-9, 2.0e-9], [2.0e-9, 1.0e-9]]
    charges = [10.0e-9, 20.0e-9]
    identity[_CAPACITANCE]["capacitance_matrix_f"] = matrix
    identity[_CAPACITANCE]["result_capacitance_matrix_f"] = matrix
    identity[_CAPACITANCE]["conductor_charge_c"] = charges
    identity[_CAPACITANCE]["result_conductor_charge_c"] = charges
    identity[_CAPACITANCE]["stored_energy_j"] = 0.5e-7
    identity[_CAPACITANCE]["result_stored_energy_j"] = 0.5e-7
    assert _gate(identity)["status"] == "needs_attention"
