from __future__ import annotations

import math

from test_comsol_generalization_v36 import _summary, gate
from test_comsol_generalization_v42 import _with_v42_induction_and_species_identity


_PROMOTED_CASE_IDS = (
    "v43_public_microwaveheating_sparameter_absorbedpower_jouleheat_temperature_energy_mismatch",
    "v43_public_poroelastic_wave_pressure_displacement_flux_dissipation_mass_energy_mismatch",
)
_MICROWAVE_KEY = (
    "microwaveheating_sparameter_absorbedpower_jouleheat_temperature_energy_"
    "mesh_result_generation_identity"
)
_POROELASTIC_KEY = (
    "poroelastic_wave_pressure_displacement_flux_dissipation_mass_energy_"
    "mesh_result_generation_identity"
)


def _with_v43_microwave_and_poroelastic_identity(summary: dict) -> dict:
    summary = _with_v42_induction_and_species_identity(summary)
    generation = "microwave-heating-726"
    incident, reflected, transmitted = 100.0, 10.0, 5.0
    absorbed = incident - reflected - transmitted
    mirrored = {
        "frequency_hz": 2.45e9,
        "reference_impedance_ohm": 50.0,
        "s11_magnitude": math.sqrt(reflected / incident),
        "s21_magnitude": math.sqrt(transmitted / incident),
        "incident_power_w": incident,
        "reflected_power_w": reflected,
        "transmitted_power_w": transmitted,
        "absorbed_power_w": absorbed,
        "joule_heat_w": 80.0,
        "dielectric_heat_w": 5.0,
        "electromagnetic_power_residual_w": 0.0,
        "outward_thermal_flux_w": absorbed,
        "ambient_temperature_k": 293.15,
        "maximum_temperature_k": 335.65,
        "temperature_rise_k": 42.5,
        "thermal_power_residual_w": 0.0,
    }
    summary[_MICROWAVE_KEY] = {
        "microwave_generation": generation,
        **{key: generation for key in (
            "sparameter_generation", "power_generation", "heat_generation",
            "temperature_generation", "energy_generation", "mesh_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "component/mesh:microwave-726",
        "accepted_mesh_owner": "component/mesh:microwave-726",
        "mesh_sha256": "1" * 64,
        "accepted_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }

    generation = "poroelastic-wave-726"
    mirrored = {
        "frequency_hz": 125.0,
        "porosity": 0.32,
        "solid_displacement_amplitude_m": 2.5e-6,
        "pore_pressure_amplitude_pa": 1250.0,
        "darcy_flux_amplitude_m_per_s": 3.0e-5,
        "pressure_displacement_phase_deg": -35.0,
        "fluid_mass_kg": 0.032,
        "fluid_mass_rate_kg_per_s": 0.004,
        "net_inward_mass_flux_kg_per_s": 0.004,
        "solid_energy_j": 1.2,
        "fluid_energy_j": 0.8,
        "dissipated_power_w": 0.4,
        "input_power_w": 0.4,
        "mass_balance_residual_kg_per_s": 0.0,
        "energy_balance_residual_w": 0.0,
    }
    summary[_POROELASTIC_KEY] = {
        "poroelastic_generation": generation,
        **{key: generation for key in (
            "pressure_generation", "displacement_generation", "flux_generation",
            "mass_generation", "energy_generation", "mesh_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "mesh_owner": "component/mesh:poroelastic-726",
        "accepted_mesh_owner": "component/mesh:poroelastic-726",
        "mesh_sha256": "3" * 64,
        "accepted_mesh_sha256": "3" * 64,
        "result_sha256": "4" * 64,
        "accepted_result_sha256": "4" * 64,
    }
    return summary


def test_v43_public_positive_microwave_and_poroelastic_contracts() -> None:
    assert gate(_with_v43_microwave_and_poroelastic_identity(_summary()))["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 2


def test_v43_public_microwaveheating_sparameter_absorbedpower_jouleheat_temperature_energy_mismatch() -> None:
    summary = _with_v43_microwave_and_poroelastic_identity(_summary())
    summary[_MICROWAVE_KEY].update({
        "power_generation": "microwave-heating-725",
        "result_generation": "microwave-heating-724",
        "result_s11_magnitude": 1.5,
        "result_absorbed_power_w": -10.0,
        "result_electromagnetic_power_residual_w": 95.0,
        "result_thermal_power_residual_w": 85.0,
        "accepted_mesh_owner": "component/mesh:old",
        "accepted_result_sha256": "a" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["microwave_heating_results_use_current_sparameters_power_heat_temperature_energy_mesh_and_result"]


def test_v43_public_poroelastic_wave_pressure_displacement_flux_dissipation_mass_energy_mismatch() -> None:
    summary = _with_v43_microwave_and_poroelastic_identity(_summary())
    summary[_POROELASTIC_KEY].update({
        "flux_generation": "poroelastic-wave-725",
        "result_generation": "poroelastic-wave-724",
        "result_porosity": 1.4,
        "result_mass_balance_residual_kg_per_s": 1.0,
        "result_energy_balance_residual_w": 2.0,
        "accepted_mesh_owner": "component/mesh:old",
        "accepted_result_sha256": "b" * 64,
    })
    result = gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"]["poroelastic_wave_results_use_current_pressure_displacement_flux_mass_energy_mesh_and_result"]


def test_v43_public_rejects_self_consistent_nonpassive_sparameters() -> None:
    summary = _with_v43_microwave_and_poroelastic_identity(_summary())
    identity = summary[_MICROWAVE_KEY]
    identity["s11_magnitude"] = identity["result_s11_magnitude"] = 1.1
    assert gate(summary)["status"] == "needs_attention"


def test_v43_public_rejects_self_consistent_invalid_poroelastic_phase() -> None:
    summary = _with_v43_microwave_and_poroelastic_identity(_summary())
    identity = summary[_POROELASTIC_KEY]
    identity["pressure_displacement_phase_deg"] = 220.0
    identity["result_pressure_displacement_phase_deg"] = 220.0
    assert gate(summary)["status"] == "needs_attention"
