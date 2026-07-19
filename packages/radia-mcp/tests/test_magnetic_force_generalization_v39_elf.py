from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import (
    magnetic_force_method_profile_gate,
)
from test_magnetic_force_generalization_v38_elf import _summary_v38


_PROMOTED_CASE_IDS = (
    "v39_public_thin_conductor_eddy_surface_impedance_skin_current_complex_power_mismatch",
    "v39_public_magnetic_gear_harmonic_polepair_phase_torque_actionreaction_mismatch",
)

_THIN_KEY = (
    "thin_conductor_surface_impedance_skin_sheetcurrent_fieldjump_complexpower_"
    "surface_owner_result_identity"
)
_GEAR_KEY = (
    "magnetic_gear_harmonic_polepair_modulation_phase_ratio_torque_"
    "actionreaction_power_owner_result_identity"
)


def _summary_v39() -> dict:
    summary = _summary_v38()
    identity = summary["artifact_identity"]

    generation = "thin-conductor-271"
    frequency = 100_000.0
    conductivity = 5.8e7
    permeability = 4.0e-7 * math.pi
    skin_depth = math.sqrt(
        2.0 / (2.0 * math.pi * frequency * permeability * conductivity)
    )
    surface_resistance = math.sqrt(
        math.pi * frequency * permeability / conductivity
    )
    sheet_current = [100.0, 0.0]
    area = 2.0e-2
    loss = 0.5 * surface_resistance * sheet_current[0] ** 2 * area
    identity[_THIN_KEY] = {
        "thin_generation": generation,
        **{
            key: generation
            for key in (
                "impedance_generation",
                "skin_generation",
                "current_generation",
                "field_generation",
                "power_generation",
                "surface_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "frequency_hz": frequency,
        "result_frequency_hz": frequency,
        "conductivity_s_m": conductivity,
        "result_conductivity_s_m": conductivity,
        "relative_permeability": 1.0,
        "result_relative_permeability": 1.0,
        "skin_depth_m": skin_depth,
        "result_skin_depth_m": skin_depth,
        "surface_impedance_ohm": [surface_resistance, surface_resistance],
        "result_surface_impedance_ohm": [surface_resistance, surface_resistance],
        "sheet_current_peak_a_m": sheet_current,
        "result_sheet_current_peak_a_m": sheet_current,
        "tangential_field_jump_peak_a_m": sheet_current,
        "result_tangential_field_jump_peak_a_m": sheet_current,
        "surface_area_m2": area,
        "result_surface_area_m2": area,
        "joule_loss_w": loss,
        "result_joule_loss_w": loss,
        "reactive_power_var": loss,
        "result_reactive_power_var": loss,
        "surface_owner": "surface:thin-conductor-271",
        "accepted_surface_owner": "surface:thin-conductor-271",
        "thin_result_sha256": "1" * 64,
        "accepted_thin_result_sha256": "1" * 64,
    }

    generation = "magnetic-gear-271"
    high_speed = 10.0
    low_speed = -20.0
    high_torque = 20.0
    low_torque = 10.0
    identity[_GEAR_KEY] = {
        "gear_generation": generation,
        **{
            key: generation
            for key in (
                "harmonic_generation",
                "pole_generation",
                "phase_generation",
                "ratio_generation",
                "torque_generation",
                "reaction_generation",
                "power_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "high_speed_pole_pairs": 4,
        "result_high_speed_pole_pairs": 4,
        "low_speed_pole_pairs": 2,
        "result_low_speed_pole_pairs": 2,
        "modulator_segment_count": 6,
        "result_modulator_segment_count": 6,
        "working_harmonic_order": 6,
        "result_working_harmonic_order": 6,
        "mechanical_phase_rad": math.pi / 12.0,
        "result_mechanical_phase_rad": math.pi / 12.0,
        "high_speed_rad_s": high_speed,
        "result_high_speed_rad_s": high_speed,
        "low_speed_rad_s": low_speed,
        "result_low_speed_rad_s": low_speed,
        "gear_ratio": -2.0,
        "result_gear_ratio": -2.0,
        "high_speed_torque_nm": high_torque,
        "result_high_speed_torque_nm": high_torque,
        "low_speed_torque_nm": low_torque,
        "result_low_speed_torque_nm": low_torque,
        "modulator_reaction_torque_nm": -(high_torque + low_torque),
        "result_modulator_reaction_torque_nm": -(high_torque + low_torque),
        "power_balance_residual_w": high_torque * high_speed + low_torque * low_speed,
        "result_power_balance_residual_w": high_torque * high_speed + low_torque * low_speed,
        "model_owner": "gear:magnetic-271",
        "accepted_model_owner": "gear:magnetic-271",
        "gear_result_sha256": "2" * 64,
        "accepted_gear_result_sha256": "2" * 64,
    }
    return summary


def test_v39_public_positive_thin_conductor_and_magnetic_gear_closure() -> None:
    assert magnetic_force_method_profile_gate(_summary_v39())["status"] == "ok"


def test_v39_public_thin_conductor_eddy_surface_impedance_skin_current_complex_power_mismatch() -> None:
    summary = _summary_v39()
    row = summary["artifact_identity"][_THIN_KEY]
    row.update(
        {
            "impedance_generation": "thin-conductor-270",
            "power_generation": "thin-conductor-269",
            "result_generation": "thin-conductor-268",
            "result_skin_depth_m": -1.0,
            "result_surface_impedance_ohm": [-1.0, 1.0],
            "result_sheet_current_peak_a_m": [-100.0, 0.0],
            "result_tangential_field_jump_peak_a_m": [0.0, 0.0],
            "result_joule_loss_w": -1.0,
            "result_reactive_power_var": -1.0,
            "accepted_surface_owner": "stale:surface",
            "accepted_thin_result_sha256": "a" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "thin_conductors_close_surface_impedance_skin_current_field_jump_complex_power_owner_and_result"
    ]


def test_v39_public_magnetic_gear_harmonic_polepair_phase_torque_actionreaction_mismatch() -> None:
    summary = _summary_v39()
    row = summary["artifact_identity"][_GEAR_KEY]
    row.update(
        {
            "pole_generation": "magnetic-gear-270",
            "reaction_generation": "magnetic-gear-269",
            "result_generation": "magnetic-gear-268",
            "result_modulator_segment_count": 5,
            "result_working_harmonic_order": 3,
            "result_mechanical_phase_rad": -1.0,
            "result_gear_ratio": 2.0,
            "result_low_speed_rad_s": 20.0,
            "result_low_speed_torque_nm": -10.0,
            "result_modulator_reaction_torque_nm": 0.0,
            "result_power_balance_residual_w": 400.0,
            "accepted_model_owner": "stale:gear",
            "accepted_gear_result_sha256": "b" * 64,
        }
    )
    result = magnetic_force_method_profile_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "magnetic_gears_close_harmonics_poles_ratio_torque_reaction_power_owner_and_result"
    ]


def test_v39_public_rejects_self_consistent_wrong_surface_impedance() -> None:
    summary = _summary_v39()
    row = summary["artifact_identity"][_THIN_KEY]
    row["surface_impedance_ohm"] = [2.0 * value for value in row["surface_impedance_ohm"]]
    row["result_surface_impedance_ohm"] = list(row["surface_impedance_ohm"])
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v39_public_rejects_self_consistent_wrong_gear_direction() -> None:
    summary = _summary_v39()
    row = summary["artifact_identity"][_GEAR_KEY]
    row["low_speed_rad_s"] = 20.0
    row["result_low_speed_rad_s"] = 20.0
    row["gear_ratio"] = 2.0
    row["result_gear_ratio"] = 2.0
    row["power_balance_residual_w"] = 400.0
    row["result_power_balance_residual_w"] = 400.0
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"
