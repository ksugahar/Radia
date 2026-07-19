from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import magnetic_force_method_profile_gate
from test_magnetic_force_generalization_v40_elf import _summary_v40


_MAGLEV = "maglev_equilibrium_airgap_position_force_gradient_stiffness_potential_energy_stability_geometry_result_identity"
_EDDY = "eddy_shield_frequency_conductivity_permeability_skin_depth_thickness_phase_attenuation_loss_energy_geometry_result_identity"
_PROMOTED_CASE_IDS = (
    "v41_public_maglev_equilibrium_forcegradient_stiffness_energy_stability_gap_mismatch",
    "v41_public_eddyshield_skin_depth_phase_lag_loss_frequency_energy_mismatch",
)


def _summary_v41() -> dict:
    summary = _summary_v40()
    identity = summary["artifact_identity"]
    generation = "maglev-equilibrium-724"
    positions = [-0.01, 0.0, 0.01]
    stiffness = 200.0
    values = {
        "air_gap_m": 0.02,
        "sample_position_m": positions,
        "force_n": [-stiffness * position for position in positions],
        "equilibrium_position_m": 0.0,
        "equilibrium_force_n": 0.0,
        "force_gradient_n_per_m": -stiffness,
        "stiffness_n_per_m": stiffness,
        "potential_energy_j": [0.5 * stiffness * position**2 for position in positions],
        "energy_curvature_j_per_m2": stiffness,
        "stability": "stable",
    }
    identity[_MAGLEV] = {
        "maglev_generation": generation,
        **{key: generation for key in ("gap_generation", "position_generation", "force_generation", "gradient_generation", "stiffness_generation", "energy_generation", "stability_generation", "geometry_generation", "result_generation")},
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "geometry_owner": "geometry:maglev-equilibrium-724",
        "accepted_geometry_owner": "geometry:maglev-equilibrium-724",
        "maglev_result_sha256": "5" * 64,
        "accepted_maglev_result_sha256": "5" * 64,
    }

    generation = "eddy-shield-724"
    mu0 = 4.0e-7 * math.pi
    frequency = 1000.0
    omega = 2.0 * math.pi * frequency
    mu_r = 1.0
    conductivity = 5.8e7
    permeability = mu0 * mu_r
    skin_depth = math.sqrt(2.0 / (omega * permeability * conductivity))
    thickness = 2.0e-3
    attenuation = math.exp(-thickness / skin_depth)
    incident = 1.0e-3
    area = 0.05
    surface_resistance = 1.0 / (conductivity * skin_depth)
    transmitted = incident * attenuation
    volume = area * thickness
    values = {
        "frequency_hz": frequency,
        "angular_frequency_rad_s": omega,
        "relative_permeability": mu_r,
        "conductivity_s_per_m": conductivity,
        "skin_depth_m": skin_depth,
        "shield_thickness_m": thickness,
        "attenuation_factor": attenuation,
        "phase_lag_rad": -thickness / skin_depth,
        "incident_field_t": incident,
        "transmitted_field_t": transmitted,
        "surface_resistance_ohm": surface_resistance,
        "shield_area_m2": area,
        "eddy_loss_w": 0.5 * surface_resistance * (incident / permeability) ** 2 * area,
        "shield_volume_m3": volume,
        "stored_energy_j": transmitted**2 * volume / (2.0 * permeability),
    }
    identity[_EDDY] = {
        "eddy_shield_generation": generation,
        **{key: generation for key in ("frequency_generation", "material_generation", "skin_depth_generation", "thickness_generation", "phase_generation", "field_generation", "loss_generation", "energy_generation", "geometry_generation", "result_generation")},
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "geometry_owner": "geometry:eddy-shield-724",
        "accepted_geometry_owner": "geometry:eddy-shield-724",
        "eddy_shield_result_sha256": "6" * 64,
        "accepted_eddy_shield_result_sha256": "6" * 64,
    }
    return summary


def test_v41_public_positive_maglev_and_eddy_shield_closure() -> None:
    assert magnetic_force_method_profile_gate(_summary_v41())["status"] == "ok"


def test_v41_public_maglev_equilibrium_forcegradient_stiffness_energy_stability_gap_mismatch() -> None:
    summary = _summary_v41()
    summary["artifact_identity"][_MAGLEV].update({"force_generation": "maglev-equilibrium-723", "result_air_gap_m": -0.01, "result_stiffness_n_per_m": -200.0, "result_stability": "unstable", "accepted_geometry_owner": "stale:geometry"})
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v41_public_eddyshield_skin_depth_phase_lag_loss_frequency_energy_mismatch() -> None:
    summary = _summary_v41()
    summary["artifact_identity"][_EDDY].update({"frequency_generation": "eddy-shield-723", "result_skin_depth_m": -1.0, "result_attenuation_factor": 2.0, "result_eddy_loss_w": -1.0, "accepted_geometry_owner": "stale:geometry"})
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_wrong_maglev_energy() -> None:
    summary = _summary_v41()
    row = summary["artifact_identity"][_MAGLEV]
    row["potential_energy_j"] = row["result_potential_energy_j"] = [0.0, 1.0, 0.0]
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v41_public_rejects_self_consistent_wrong_eddy_skin_depth() -> None:
    summary = _summary_v41()
    row = summary["artifact_identity"][_EDDY]
    row["skin_depth_m"] = row["result_skin_depth_m"] = 1.0
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"
