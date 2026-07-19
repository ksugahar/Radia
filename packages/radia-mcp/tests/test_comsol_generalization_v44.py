from __future__ import annotations

from test_comsol_generalization_v36 import _summary, gate
from test_comsol_generalization_v43 import _with_v43_microwave_and_poroelastic_identity


_PROMOTED_CASE_IDS = (
    "v44_public_microwave_boundaryport_sparameter_power_normalization_temperature_coupling_restart_owner_mismatch",
    "v44_public_acoustics_poroelastic_impedance_phase_flux_energy_timewindow_dataset_owner_mismatch",
)


def _with_v44(summary: dict) -> dict:
    summary = _with_v43_microwave_and_poroelastic_identity(summary)
    generation = "microwave-port-744"
    values = {
        "frequency_hz": 2.45e9, "port_normalization_ohm": 50.0,
        "incident_power_w": 100.0, "reflected_power_w": 10.0,
        "transmitted_power_w": 5.0, "absorbed_power_w": 85.0,
        "s11_power_fraction": 0.10, "s21_power_fraction": 0.05,
        "thermal_coupling_power_w": 85.0, "temperature_rise_k": 42.5,
    }
    key = "microwave_boundaryport_sparameter_power_normalization_temperature_coupling_restart_owner_result_identity"
    summary[key] = {
        "microwave_port_generation": generation, "port_generation": generation,
        "power_generation": generation, "thermal_generation": generation,
        "restart_generation": generation, "owner_generation": generation,
        "result_generation": generation, **values,
        **{f"result_{name}": value for name, value in values.items()},
        "restart_checkpoint_id": "chk-microwave-744", "result_restart_checkpoint_id": "chk-microwave-744",
        "mesh_owner": "component/mesh:microwave-744", "accepted_mesh_owner": "component/mesh:microwave-744",
        "mesh_sha256": "7" * 64, "accepted_mesh_sha256": "7" * 64,
        "result_sha256": "8" * 64, "accepted_result_sha256": "8" * 64,
    }
    generation = "acoustic-poroelastic-744"
    values = {
        "frequency_hz": 125.0, "impedance_magnitude_pa_s_per_m": 2.5e5,
        "impedance_phase_deg": -35.0, "pressure_flux_w": 0.4,
        "mechanical_energy_j": 1.2, "fluid_energy_j": 0.8,
        "dissipated_power_w": 0.4, "time_window_s": 0.008,
        "energy_balance_residual_w": 0.0,
    }
    key = "acoustics_poroelastic_impedance_phase_flux_energy_timewindow_dataset_owner_result_identity"
    summary[key] = {
        "acoustic_generation": generation, "impedance_generation": generation,
        "phase_generation": generation, "flux_generation": generation,
        "energy_generation": generation, "dataset_generation": generation,
        "owner_generation": generation, "result_generation": generation,
        **values, **{f"result_{name}": value for name, value in values.items()},
        "dataset_tag": "dset-acoustic-744", "result_dataset_tag": "dset-acoustic-744",
        "boundary_owner": "component/boundary:acoustic-744", "accepted_boundary_owner": "component/boundary:acoustic-744",
        "boundary_sha256": "9" * 64, "accepted_boundary_sha256": "9" * 64,
        "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
    }
    return summary


def test_v44_public_positive_port_and_acoustic_closure() -> None:
    assert gate(_with_v44(_summary()))["status"] == "ok"
    assert len(_PROMOTED_CASE_IDS) == 2


def test_v44_public_rejects_port_power_mismatch() -> None:
    summary = _with_v44(_summary())
    summary["microwave_boundaryport_sparameter_power_normalization_temperature_coupling_restart_owner_result_identity"]["result_absorbed_power_w"] = -1.0
    assert gate(summary)["status"] == "needs_attention"


def test_v44_public_rejects_acoustic_phase_mismatch() -> None:
    summary = _with_v44(_summary())
    summary["acoustics_poroelastic_impedance_phase_flux_energy_timewindow_dataset_owner_result_identity"]["result_impedance_phase_deg"] = 45.0
    assert gate(summary)["status"] == "needs_attention"
