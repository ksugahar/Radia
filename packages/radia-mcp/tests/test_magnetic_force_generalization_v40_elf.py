from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import magnetic_force_method_profile_gate
from test_magnetic_force_generalization_v39_elf import _summary_v39


_SHIELD = "multilayer_magnetic_shield_permeability_thickness_radius_interface_flux_attenuation_leakage_energy_geometry_result_identity"
_TRANSFORMER = "transformer_leakage_mutual_inductance_fluxlinkage_reciprocity_psd_coenergy_force_winding_result_identity"
_PROMOTED_CASE_IDS = (
    "v40_public_multilayer_magnetic_shield_permeability_thickness_attenuation_flux_energy_mismatch",
    "v40_public_transformer_leakage_mutual_inductance_fluxlinkage_reciprocity_coenergy_force_mismatch",
)


def _summary_v40() -> dict:
    summary = _summary_v39()
    identity = summary["artifact_identity"]
    generation = "multilayer-shield-280"
    permeability = [2000.0, 5000.0]
    thickness = [1.0e-3, 0.8e-3]
    radii = [0.12, 0.10]
    factors = [1.0 + (mu_r - 1.0) * layer / (2.0 * radius) for mu_r, layer, radius in zip(permeability, thickness, radii)]
    attenuation = math.prod(factors)
    external = 1.0e-3
    leakage = external / attenuation
    volume = 4.0e-3
    energy = leakage**2 * volume / (2.0 * 4.0e-7 * math.pi)
    values = {
        "relative_permeability": permeability, "layer_thickness_m": thickness,
        "layer_mean_radius_m": radii, "layer_shielding_factor": factors,
        "interface_normal_flux_t": [leakage, leakage, leakage],
        "external_field_t": external, "attenuation_factor": attenuation,
        "leakage_field_t": leakage, "cavity_volume_m3": volume,
        "stored_energy_j": energy,
    }
    identity[_SHIELD] = {
        "shield_generation": generation,
        **{key: generation for key in ("material_generation", "thickness_generation", "geometry_generation", "flux_generation", "attenuation_generation", "field_generation", "energy_generation", "owner_generation", "result_generation")},
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "geometry_owner": "geometry:multilayer-shield-280",
        "accepted_geometry_owner": "geometry:multilayer-shield-280",
        "shield_result_sha256": "1" * 64,
        "accepted_shield_result_sha256": "1" * 64,
    }

    generation = "transformer-coupling-280"
    primary = 12.0e-3
    secondary = 3.0e-3
    mutual = 5.4e-3
    currents = [4.0, -6.0]
    linkages = [primary * currents[0] + mutual * currents[1], mutual * currents[0] + secondary * currents[1]]
    coenergy = 0.5 * sum(current * linkage for current, linkage in zip(currents, linkages))
    gradient = -2.0e-2
    values = {
        "inductance_matrix_h": [[primary, mutual], [mutual, secondary]],
        "primary_leakage_inductance_h": primary - mutual**2 / secondary,
        "winding_currents_a": currents, "flux_linkages_wb_turn": linkages,
        "reciprocity_residual_h": 0.0, "coenergy_j": coenergy,
        "mutual_inductance_gradient_h_per_m": gradient,
        "force_n": gradient * currents[0] * currents[1],
    }
    identity[_TRANSFORMER] = {
        "transformer_generation": generation,
        **{key: generation for key in ("inductance_generation", "leakage_generation", "flux_generation", "reciprocity_generation", "energy_generation", "force_generation", "winding_generation", "result_generation")},
        **values,
        **{f"result_{key}": value for key, value in values.items()},
        "winding_owner": "winding:transformer-280",
        "accepted_winding_owner": "winding:transformer-280",
        "transformer_result_sha256": "2" * 64,
        "accepted_transformer_result_sha256": "2" * 64,
    }
    return summary


def test_v40_public_positive_shield_and_transformer_closure() -> None:
    assert magnetic_force_method_profile_gate(_summary_v40())["status"] == "ok"


def test_v40_public_multilayer_magnetic_shield_permeability_thickness_attenuation_flux_energy_mismatch() -> None:
    summary = _summary_v40()
    summary["artifact_identity"][_SHIELD].update({"material_generation": "multilayer-shield-279", "result_attenuation_factor": 0.5, "result_interface_normal_flux_t": [1.0, 2.0], "result_stored_energy_j": -1.0, "accepted_geometry_owner": "stale:geometry"})
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v40_public_transformer_leakage_mutual_inductance_fluxlinkage_reciprocity_coenergy_force_mismatch() -> None:
    summary = _summary_v40()
    summary["artifact_identity"][_TRANSFORMER].update({"inductance_generation": "transformer-coupling-279", "result_inductance_matrix_h": [[-1.0, 2.0], [3.0, -1.0]], "result_coenergy_j": -1.0, "result_force_n": -1.0, "accepted_winding_owner": "stale:winding"})
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_wrong_shield_attenuation() -> None:
    summary = _summary_v40()
    row = summary["artifact_identity"][_SHIELD]
    row["attenuation_factor"] = row["result_attenuation_factor"] = 2.0
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"


def test_v40_public_rejects_self_consistent_nonreciprocal_transformer() -> None:
    summary = _summary_v40()
    row = summary["artifact_identity"][_TRANSFORMER]
    row["inductance_matrix_h"] = row["result_inductance_matrix_h"] = [[12.0e-3, 5.4e-3], [4.0e-3, 3.0e-3]]
    assert magnetic_force_method_profile_gate(summary)["status"] == "needs_attention"
