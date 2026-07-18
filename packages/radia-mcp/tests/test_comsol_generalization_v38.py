from __future__ import annotations

from test_comsol_generalization_v36 import _summary, gate
from test_comsol_generalization_v37 import (
    _with_v37_capacitance_and_thermoelastic_identity,
)

_PROMOTED_CASE_IDS = (
    "v38_public_thermoviscous_pressure_acoustics_interface_velocity_traction_dissipation_power_mismatch",
    "v38_public_piezoelectric_charge_strain_reciprocity_electromechanical_energy_terminal_mismatch",
)


def _with_v38_thermoviscous_and_piezoelectric_identity(summary: dict) -> dict:
    summary = _with_v37_capacitance_and_thermoelastic_identity(summary)
    generation = "thermoviscous-interface-258"
    summary[
        "thermoviscous_pressure_interface_velocity_traction_dissipation_power_normal_mesh_result_generation_identity"
    ] = {
        "interface_generation": generation,
        **{
            key: generation
            for key in (
                "velocity_generation",
                "traction_generation",
                "viscous_generation",
                "thermal_generation",
                "power_generation",
                "normal_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "frequency_hz": 1.0e4,
        "result_frequency_hz": 1.0e4,
        "interface_area_m2": 1.0e-2,
        "result_interface_area_m2": 1.0e-2,
        "normal_velocity_complex_m_per_s": [1.0e-2, 2.0e-3],
        "result_normal_velocity_complex_m_per_s": [1.0e-2, 2.0e-3],
        "pressure_complex_pa": [200.0, 40.0],
        "result_pressure_complex_pa": [200.0, 40.0],
        "traction_sign": "minus_pressure_times_outward_normal",
        "result_traction_sign": "minus_pressure_times_outward_normal",
        "normal_orientation": "thermoviscous_to_pressure_acoustics",
        "result_normal_orientation": "thermoviscous_to_pressure_acoustics",
        "interface_power_w": 1.04e-2,
        "result_interface_power_w": 1.04e-2,
        "viscous_loss_w": 3.0e-3,
        "result_viscous_loss_w": 3.0e-3,
        "thermal_loss_w": 2.0e-3,
        "result_thermal_loss_w": 2.0e-3,
        "outgoing_acoustic_power_w": 5.4e-3,
        "result_outgoing_acoustic_power_w": 5.4e-3,
        "mesh_owner": "comp1/mesh1:thermoviscous-258",
        "accepted_mesh_owner": "comp1/mesh1:thermoviscous-258",
        "mesh_sha256": "1" * 64,
        "accepted_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }

    generation = "piezoelectric-reciprocity-258"
    summary[
        "piezoelectric_charge_strain_reciprocity_electromechanical_energy_polarization_mesh_result_generation_identity"
    ] = {
        "piezo_generation": generation,
        **{
            key: generation
            for key in (
                "charge_generation",
                "strain_generation",
                "reciprocity_generation",
                "electrical_energy_generation",
                "elastic_energy_generation",
                "coupling_energy_generation",
                "polarization_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "direct_coefficient_c_per_n": 2.0e-10,
        "result_direct_coefficient_c_per_n": 2.0e-10,
        "converse_coefficient_m_per_v": 2.0e-10,
        "result_converse_coefficient_m_per_v": 2.0e-10,
        "electric_field_v_per_m": 1.0e5,
        "result_electric_field_v_per_m": 1.0e5,
        "mechanical_stress_pa": 1.0e6,
        "result_mechanical_stress_pa": 1.0e6,
        "induced_strain": 2.0e-5,
        "result_induced_strain": 2.0e-5,
        "induced_charge_density_c_per_m2": 2.0e-4,
        "result_induced_charge_density_c_per_m2": 2.0e-4,
        "terminal_charge_c": 4.0e-6,
        "result_terminal_charge_c": 4.0e-6,
        "electrical_work_j": 3.0e-2,
        "result_electrical_work_j": 3.0e-2,
        "elastic_energy_j": 2.0e-2,
        "result_elastic_energy_j": 2.0e-2,
        "coupling_energy_j": 1.0e-2,
        "result_coupling_energy_j": 1.0e-2,
        "total_stored_energy_j": 4.0e-2,
        "result_total_stored_energy_j": 4.0e-2,
        "polarization_frame": "material_axis_3",
        "result_polarization_frame": "material_axis_3",
        "mesh_owner": "comp1/mesh1:piezo-258",
        "accepted_mesh_owner": "comp1/mesh1:piezo-258",
        "mesh_sha256": "3" * 64,
        "accepted_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def test_v38_public_positive_thermoviscous_and_piezoelectric_contracts() -> None:
    result = gate(_with_v38_thermoviscous_and_piezoelectric_identity(_summary()))
    assert result["status"] == "ok"


def test_v38_public_thermoviscous_pressure_acoustics_interface_velocity_traction_dissipation_power_mismatch() -> None:
    summary = _with_v38_thermoviscous_and_piezoelectric_identity(_summary())
    identity = summary[
        "thermoviscous_pressure_interface_velocity_traction_dissipation_power_normal_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "velocity_generation": "thermoviscous-interface-257",
            "power_generation": "thermoviscous-interface-256",
            "result_generation": "thermoviscous-interface-255",
            "result_frequency_hz": -1.0e4,
            "result_normal_velocity_complex_m_per_s": [-1.0e-2, 2.0e-3],
            "result_pressure_complex_pa": [200.0, -40.0],
            "result_traction_sign": "plus_pressure_times_normal",
            "result_normal_orientation": "pressure_to_thermoviscous",
            "result_interface_power_w": -1.04e-2,
            "result_viscous_loss_w": -3.0e-3,
            "result_outgoing_acoustic_power_w": 1.0e-1,
            "accepted_mesh_owner": "comp1/mesh0:old",
            "accepted_result_sha256": "9" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "thermoviscous_interfaces_use_current_velocity_traction_dissipation_power_normal_mesh_and_result"
    ]


def test_v38_public_piezoelectric_charge_strain_reciprocity_electromechanical_energy_terminal_mismatch() -> None:
    summary = _with_v38_thermoviscous_and_piezoelectric_identity(_summary())
    identity = summary[
        "piezoelectric_charge_strain_reciprocity_electromechanical_energy_polarization_mesh_result_generation_identity"
    ]
    identity.update(
        {
            "reciprocity_generation": "piezoelectric-reciprocity-257",
            "coupling_energy_generation": "piezoelectric-reciprocity-256",
            "result_generation": "piezoelectric-reciprocity-255",
            "result_direct_coefficient_c_per_n": -2.0e-10,
            "result_converse_coefficient_m_per_v": 4.0e-10,
            "result_induced_strain": -2.0e-5,
            "result_induced_charge_density_c_per_m2": 1.0e-3,
            "result_terminal_charge_c": -4.0e-6,
            "result_coupling_energy_j": -1.0e-2,
            "result_total_stored_energy_j": 6.0e-2,
            "result_polarization_frame": "global_z",
            "accepted_mesh_sha256": "a" * 64,
            "accepted_result_sha256": "b" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "piezoelectric_results_use_current_charge_strain_reciprocity_energy_polarization_mesh_and_result"
    ]


def test_v38_public_rejects_self_consistent_negative_thermal_loss() -> None:
    summary = _with_v38_thermoviscous_and_piezoelectric_identity(_summary())
    identity = summary[
        "thermoviscous_pressure_interface_velocity_traction_dissipation_power_normal_mesh_result_generation_identity"
    ]
    identity["thermal_loss_w"] = identity["result_thermal_loss_w"] = -2.0e-3
    identity["outgoing_acoustic_power_w"] = identity[
        "result_outgoing_acoustic_power_w"
    ] = 9.4e-3
    assert gate(summary)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_nonreciprocal_piezo_coefficients() -> None:
    summary = _with_v38_thermoviscous_and_piezoelectric_identity(_summary())
    identity = summary[
        "piezoelectric_charge_strain_reciprocity_electromechanical_energy_polarization_mesh_result_generation_identity"
    ]
    identity["converse_coefficient_m_per_v"] = identity[
        "result_converse_coefficient_m_per_v"
    ] = 4.0e-10
    identity["induced_strain"] = identity["result_induced_strain"] = 4.0e-5
    assert gate(summary)["status"] == "needs_attention"
