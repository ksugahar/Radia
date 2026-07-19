from __future__ import annotations

import math

from test_comsol_generalization_v36 import _summary, gate
from test_comsol_generalization_v40 import _with_v40_thermoacoustic_and_battery_identity


_PROMOTED_CASE_IDS = (
    "v41_public_piezoelectric_admittance_resonance_antiresonance_coupling_energy_phase_mismatch",
    "v41_public_fluidfilm_bearing_reynolds_pressure_load_friction_temperature_power_mismatch",
)
_PIEZO_KEY = (
    "piezoelectric_admittance_resonance_antiresonance_coupling_energy_phase_"
    "mesh_result_generation_identity"
)
_BEARING_KEY = (
    "fluidfilm_bearing_reynolds_pressure_load_friction_temperature_power_mesh_"
    "result_generation_identity"
)


def _with_v41_piezoelectric_and_bearing_identity(summary: dict) -> dict:
    summary = _with_v40_thermoacoustic_and_battery_identity(summary)
    generation = "piezoelectric-724"
    resonance, antiresonance = 100_000.0, 105_000.0
    voltage, current, phase_deg = 10.0, 0.02, -30.0
    apparent_power = voltage * current
    real_power = apparent_power * math.cos(math.radians(phase_deg))
    mechanical_power = 0.15
    mirrored = {
        "resonance_frequency_hz": resonance,
        "antiresonance_frequency_hz": antiresonance,
        "electromechanical_coupling_squared": 1.0 - (resonance / antiresonance) ** 2,
        "voltage_rms_v": voltage,
        "current_rms_a": current,
        "admittance_magnitude_s": current / voltage,
        "admittance_phase_deg": phase_deg,
        "real_electrical_power_w": real_power,
        "reactive_electrical_power_var": apparent_power * math.sin(math.radians(phase_deg)),
        "mechanical_output_power_w": mechanical_power,
        "dielectric_loss_w": real_power - mechanical_power,
        "mechanical_stored_energy_j": 2.0e-6,
        "electric_stored_energy_j": 3.0e-6,
        "power_balance_residual_w": 0.0,
    }
    summary[_PIEZO_KEY] = {
        "piezoelectric_generation": generation,
        **{key: generation for key in (
            "admittance_generation", "resonance_generation", "coupling_generation",
            "phase_generation", "energy_generation", "power_generation",
            "mesh_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "comp1/mesh1:piezoelectric-724",
        "accepted_mesh_owner": "comp1/mesh1:piezoelectric-724",
        "mesh_sha256": "1" * 64,
        "accepted_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }

    generation = "fluidfilm-bearing-724"
    clearance, eccentricity = 50.0e-6, 0.6
    angular_speed, friction_torque = 100.0, 2.0
    shaft_power, viscous_power = angular_speed * friction_torque, 150.0
    mirrored = {
        "journal_radius_m": 0.025,
        "bearing_length_m": 0.05,
        "radial_clearance_m": clearance,
        "eccentricity_ratio": eccentricity,
        "minimum_film_thickness_m": clearance * (1.0 - eccentricity),
        "maximum_pressure_pa": 4.0e6,
        "integrated_load_n": 1000.0,
        "attitude_angle_deg": 55.0,
        "angular_speed_rad_per_s": angular_speed,
        "friction_torque_nm": friction_torque,
        "shaft_power_w": shaft_power,
        "viscous_dissipation_w": viscous_power,
        "removed_heat_w": shaft_power - viscous_power,
        "inlet_temperature_k": 313.15,
        "maximum_temperature_k": 333.15,
        "power_balance_residual_w": 0.0,
    }
    summary[_BEARING_KEY] = {
        "fluidfilm_generation": generation,
        **{key: generation for key in (
            "film_generation", "pressure_generation", "load_generation",
            "friction_generation", "temperature_generation", "power_generation",
            "mesh_generation", "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "comp1/mesh1:fluidfilm-724",
        "accepted_mesh_owner": "comp1/mesh1:fluidfilm-724",
        "mesh_sha256": "3" * 64,
        "accepted_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def test_v41_public_positive_piezoelectric_and_bearing_contracts() -> None:
    assert gate(_with_v41_piezoelectric_and_bearing_identity(_summary()))["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 2


def test_v41_public_piezoelectric_admittance_resonance_antiresonance_coupling_energy_phase_mismatch() -> None:
    summary = _with_v41_piezoelectric_and_bearing_identity(_summary())
    summary[_PIEZO_KEY].update({
        "coupling_generation": "piezoelectric-723",
        "result_generation": "piezoelectric-722",
        "result_antiresonance_frequency_hz": 90_000.0,
        "result_electromechanical_coupling_squared": -1.0,
        "result_admittance_phase_deg": 150.0,
        "result_real_electrical_power_w": -2.0,
        "result_power_balance_residual_w": 9.0,
        "accepted_mesh_owner": "comp1/mesh0:old",
        "accepted_result_sha256": "a" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["piezoelectric_admittance_results_use_current_resonance_coupling_phase_energy_power_mesh_and_result"]


def test_v41_public_fluidfilm_bearing_reynolds_pressure_load_friction_temperature_power_mismatch() -> None:
    summary = _with_v41_piezoelectric_and_bearing_identity(_summary())
    summary[_BEARING_KEY].update({
        "pressure_generation": "fluidfilm-bearing-723",
        "result_generation": "fluidfilm-bearing-722",
        "result_minimum_film_thickness_m": -1.0,
        "result_maximum_pressure_pa": -1.0,
        "result_integrated_load_n": -1.0,
        "result_friction_torque_nm": -2.0,
        "result_maximum_temperature_k": 250.0,
        "result_power_balance_residual_w": 9.0,
        "accepted_mesh_owner": "comp1/mesh0:old",
        "accepted_result_sha256": "b" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["fluidfilm_bearing_results_use_current_film_pressure_load_friction_temperature_power_mesh_and_result"]


def test_v41_public_rejects_self_consistent_antiresonance_below_resonance() -> None:
    summary = _with_v41_piezoelectric_and_bearing_identity(_summary())
    identity = summary[_PIEZO_KEY]
    antiresonance = 90_000.0
    coupling = 1.0 - (identity["resonance_frequency_hz"] / antiresonance) ** 2
    identity["antiresonance_frequency_hz"] = identity["result_antiresonance_frequency_hz"] = antiresonance
    identity["electromechanical_coupling_squared"] = identity["result_electromechanical_coupling_squared"] = coupling
    assert gate(summary)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_overclosed_bearing_film() -> None:
    summary = _with_v41_piezoelectric_and_bearing_identity(_summary())
    identity = summary[_BEARING_KEY]
    identity["eccentricity_ratio"] = identity["result_eccentricity_ratio"] = 1.2
    thickness = identity["radial_clearance_m"] * (1.0 - 1.2)
    identity["minimum_film_thickness_m"] = identity["result_minimum_film_thickness_m"] = thickness
    assert gate(summary)["status"] == "needs_attention"
