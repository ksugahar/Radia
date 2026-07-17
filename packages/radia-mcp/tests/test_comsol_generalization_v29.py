from __future__ import annotations

from test_rotational_eddy_brake_energy_gate import (
    _summary,
    _with_v28_thermoelastic_and_field_circuit_identity,
    gate,
)


_PROMOTED_CASE_IDS = (
    "v29_public_rotating_sliding_interface_sector_pitch_azimuth_interpolation_torque_periodicity_mismatch",
    "v29_public_acoustic_radiation_impedance_modal_projection_reference_area_power_flux_mismatch",
)


def _with_v29_sliding_and_radiation_identity(summary: dict) -> dict:
    summary = _with_v28_thermoelastic_and_field_circuit_identity(summary)
    generation = "sliding-interface-161"
    summary[
        "rotating_sliding_interface_sector_pitch_azimuth_interpolation_frame_periodicity_mesh_torque_generation_identity"
    ] = {
        "sliding_generation": generation,
        "sector_sliding_generation": generation,
        "azimuth_sliding_generation": generation,
        "interpolation_sliding_generation": generation,
        "frame_sliding_generation": generation,
        "periodicity_sliding_generation": generation,
        "mesh_sliding_generation": generation,
        "result_sliding_generation": generation,
        "sector_pitch_deg": 30.0,
        "result_sector_pitch_deg": 30.0,
        "sector_count": 12,
        "result_sector_count": 12,
        "azimuth_origin_deg": 0.0,
        "result_azimuth_origin_deg": 0.0,
        "source_interface_tag": "rotor_if",
        "result_source_interface_tag": "rotor_if",
        "target_interface_tag": "stator_if",
        "result_target_interface_tag": "stator_if",
        "interpolation": "conservative_mortar_azimuth",
        "result_interpolation": "conservative_mortar_azimuth",
        "rotor_frame": "rotor_cylindrical",
        "result_rotor_frame": "rotor_cylindrical",
        "periodic_phase_deg": 0.0,
        "result_periodic_phase_deg": 0.0,
        "azimuth_samples_deg": [0.0, 7.5, 15.0, 22.5, 30.0],
        "result_azimuth_samples_deg": [0.0, 7.5, 15.0, 22.5, 30.0],
        "torque_nm": [10.0, 10.5, 10.0, 9.5, 10.0],
        "result_torque_nm": [10.0, 10.5, 10.0, 9.5, 10.0],
        "sliding_mesh_sha256": "1" * 64,
        "result_sliding_mesh_sha256": "1" * 64,
        "torque_result_sha256": "2" * 64,
        "accepted_torque_result_sha256": "2" * 64,
    }
    generation = "acoustic-radiation-161"
    summary[
        "acoustic_radiation_impedance_modal_trace_reference_area_pressure_velocity_power_frequency_result_generation_identity"
    ] = {
        "radiation_generation": generation,
        "modal_radiation_generation": generation,
        "trace_radiation_generation": generation,
        "area_radiation_generation": generation,
        "convention_radiation_generation": generation,
        "power_radiation_generation": generation,
        "frequency_radiation_generation": generation,
        "result_radiation_generation": generation,
        "modal_basis_id": "p1_interface_modes_mass_normalized",
        "result_modal_basis_id": "p1_interface_modes_mass_normalized",
        "mode_indices": [1, 2, 3],
        "result_mode_indices": [1, 2, 3],
        "trace_projection": "l2_p1_boundary",
        "result_trace_projection": "l2_p1_boundary",
        "reference_area_m2": 0.1,
        "result_reference_area_m2": 0.1,
        "pressure_velocity_convention": "outward_positive_velocity",
        "result_pressure_velocity_convention": "outward_positive_velocity",
        "frequency_grid_hz": [100.0, 200.0, 500.0],
        "result_frequency_grid_hz": [100.0, 200.0, 500.0],
        "radiation_impedance_ri": [[20.0, 5.0], [30.0, 8.0], [50.0, 12.0]],
        "result_radiation_impedance_ri": [[20.0, 5.0], [30.0, 8.0], [50.0, 12.0]],
        "outward_power_flux_w": [1.0, 1.5, 2.0],
        "result_outward_power_flux_w": [1.0, 1.5, 2.0],
        "radiation_mesh_sha256": "3" * 64,
        "result_radiation_mesh_sha256": "3" * 64,
        "radiation_result_sha256": "4" * 64,
        "accepted_radiation_result_sha256": "4" * 64,
    }
    return summary


def test_v29_public_positive_sliding_interface_and_acoustic_radiation() -> None:
    assert gate(_with_v29_sliding_and_radiation_identity(_summary()))["status"] == "ok"


def test_v29_public_rotating_sliding_interface_sector_pitch_azimuth_interpolation_torque_periodicity_mismatch() -> None:
    summary = _with_v29_sliding_and_radiation_identity(_summary())
    summary[
        "rotating_sliding_interface_sector_pitch_azimuth_interpolation_frame_periodicity_mesh_torque_generation_identity"
    ].update(
        {
            "sector_sliding_generation": "sliding-interface-160",
            "interpolation_sliding_generation": "sliding-interface-159",
            "result_sliding_generation": "sliding-interface-158",
            "result_sector_pitch_deg": 60.0,
            "result_sector_count": 10,
            "result_azimuth_origin_deg": 15.0,
            "result_source_interface_tag": "stator_if",
            "result_target_interface_tag": "rotor_if",
            "result_interpolation": "nearest_neighbor",
            "result_rotor_frame": "global_cartesian",
            "result_periodic_phase_deg": 180.0,
            "result_azimuth_samples_deg": [0.0, 10.0, 20.0],
            "result_torque_nm": [10.0, 8.0, 12.0],
            "result_sliding_mesh_sha256": "8" * 64,
            "accepted_torque_result_sha256": "9" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "rotating_sliding_interface_uses_current_sector_azimuth_interpolation_frame_periodicity_mesh_and_torque"
    ]


def test_v29_public_acoustic_radiation_impedance_modal_projection_reference_area_power_flux_mismatch() -> None:
    summary = _with_v29_sliding_and_radiation_identity(_summary())
    summary[
        "acoustic_radiation_impedance_modal_trace_reference_area_pressure_velocity_power_frequency_result_generation_identity"
    ].update(
        {
            "modal_radiation_generation": "acoustic-radiation-160",
            "convention_radiation_generation": "acoustic-radiation-159",
            "result_radiation_generation": "acoustic-radiation-158",
            "result_modal_basis_id": "p0_unnormalized",
            "result_mode_indices": [3, 2, 1],
            "result_trace_projection": "point_sample",
            "result_reference_area_m2": 1.0,
            "result_pressure_velocity_convention": "inward_positive_velocity",
            "result_frequency_grid_hz": [100.0, 300.0],
            "result_radiation_impedance_ri": [[-20.0, 5.0]],
            "result_outward_power_flux_w": [-1.0],
            "result_radiation_mesh_sha256": "a" * 64,
            "accepted_radiation_result_sha256": "b" * 64,
        }
    )
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "acoustic_radiation_uses_current_modes_trace_area_convention_frequency_power_and_result"
    ]
