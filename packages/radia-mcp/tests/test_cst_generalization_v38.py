from __future__ import annotations

import math

from radia_mcp.radia_ngsolve.nonlinear_inductance_sweep_gate import (
    nonlinear_inductance_sweep_gate,
)
from test_cst_generalization_v37 import _summary_v37


_PROMOTED_CASE_IDS = (
    "v38_public_microstrip_quasitem_impedance_effective_permittivity_delay_loss_sparameter_mismatch",
    "v38_public_shielding_aperture_incident_transmitted_field_power_se_frequency_mesh_mismatch",
)
C0 = 299_792_458.0


def _microstrip(index: int):
    generation = f"microstrip-quasitem-{614 + index}"
    width, height, copper, relative_permittivity, length = (
        2.0e-3, 1.0e-3, 35.0e-6, 4.0, 0.1
    )
    ratio = width / height
    effective_permittivity = (
        (relative_permittivity + 1.0) / 2.0
        + (relative_permittivity - 1.0)
        / (2.0 * math.sqrt(1.0 + 12.0 / ratio))
    )
    impedance = 120.0 * math.pi / (
        math.sqrt(effective_permittivity)
        * (ratio + 1.393 + 0.667 * math.log(ratio + 1.444))
    )
    conductor_loss, dielectric_loss, s11 = 0.02, 0.01, 0.1
    s21 = math.sqrt(1.0 - s11**2 - conductor_loss - dielectric_loss)
    mirrored = {
        "trace_width_m": width,
        "substrate_height_m": height,
        "copper_thickness_m": copper,
        "relative_permittivity": relative_permittivity,
        "line_length_m": length,
        "quasitem_model": "hammerstad_zero_thickness_core",
        "effective_permittivity": effective_permittivity,
        "characteristic_impedance_ohm": impedance,
        "propagation_delay_s": length * math.sqrt(effective_permittivity) / C0,
        "frequency_hz": [1.0e9, 2.0e9, 3.0e9],
        "conductor_loss_fraction": conductor_loss,
        "dielectric_loss_fraction": dielectric_loss,
        "s11_magnitude": s11,
        "s21_magnitude": s21,
        "incident_power_w": 1.0,
        "mesh_dof": [20_000, 40_000, 80_000],
        "mesh_impedance_ohm": [impedance * 1.02, impedance * 1.004, impedance],
        "maximum_final_relative_change": 0.01,
        "microstrip_mesh_sha256": "1" * 64,
    }
    return {
        "microstrip_generation": generation,
        **{key: generation for key in (
            "geometry_generation", "material_generation", "impedance_generation",
            "delay_generation", "loss_generation", "sparameter_generation",
            "frequency_generation", "mesh_generation", "owner_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "microstrip_owner": f"pcb/microstrip-{614 + index}",
        "accepted_microstrip_owner": f"pcb/microstrip-{614 + index}",
        "microstrip_result_sha256": "2" * 64,
        "accepted_microstrip_result_sha256": "2" * 64,
    }


def _shielding(index: int):
    generation = f"shield-aperture-{614 + index}"
    incident = [1.0, 1.0, 1.0]
    transmitted = [0.1, 0.05, 0.02]
    incident_power = [item**2 for item in incident]
    transmitted_power = [item**2 for item in transmitted]
    field_se = [20.0 * math.log10(a / b) for a, b in zip(incident, transmitted)]
    power_se = [
        10.0 * math.log10(a / b)
        for a, b in zip(incident_power, transmitted_power)
    ]
    mirrored = {
        "aperture_plane": "xy_normal_positive_z",
        "incident_polarization": "x_linear",
        "frequency_hz": [1.0e9, 2.0e9, 3.0e9],
        "incident_field_v_m": incident,
        "transmitted_field_v_m": transmitted,
        "incident_power_density_normalized": incident_power,
        "transmitted_power_density_normalized": transmitted_power,
        "shielding_effectiveness_field_db": field_se,
        "shielding_effectiveness_power_db": power_se,
        "probe_frame": "global_cartesian_xyz_m",
        "mesh_dof": [30_000, 60_000, 120_000],
        "mesh_selected_se_db": [power_se[1] - 1.0, power_se[1] - 0.2, power_se[1]],
        "maximum_final_se_change_db": 0.5,
        "shield_mesh_sha256": "3" * 64,
    }
    return {
        "shield_generation": generation,
        **{key: generation for key in (
            "aperture_generation", "polarization_generation", "field_generation",
            "power_generation", "se_generation", "frequency_generation",
            "probe_generation", "mesh_generation", "owner_generation",
            "result_generation",
        )},
        **mirrored,
        **{f"result_{key}": value for key, value in mirrored.items()},
        "shield_owner": f"emc/shield-{614 + index}",
        "accepted_shield_owner": f"emc/shield-{614 + index}",
        "shield_result_sha256": "4" * 64,
        "accepted_shield_result_sha256": "4" * 64,
    }


def _summary_v38():
    summary = _summary_v37()
    for index, row in enumerate(summary["runs"]):
        row[
            "microstrip_quasitem_geometry_permittivity_impedance_delay_loss_sparameter_frequency_mesh_owner_result_identity"
        ] = _microstrip(index)
        row[
            "shielding_aperture_orientation_polarization_field_power_se_frequency_probe_mesh_owner_result_identity"
        ] = _shielding(index)
    return summary


def test_v38_public_positive_microstrip_and_shielding_closure():
    assert nonlinear_inductance_sweep_gate(_summary_v38())["status"] == "ok"


def test_v38_public_microstrip_quasitem_impedance_effective_permittivity_delay_loss_sparameter_mismatch():
    summary = _summary_v38()
    row = summary["runs"][0][
        "microstrip_quasitem_geometry_permittivity_impedance_delay_loss_sparameter_frequency_mesh_owner_result_identity"
    ]
    row.update({
        "impedance_generation": "microstrip-quasitem-613",
        "loss_generation": "microstrip-quasitem-612",
        "result_generation": "microstrip-quasitem-611",
        "result_effective_permittivity": -1.0,
        "result_characteristic_impedance_ohm": -50.0,
        "result_propagation_delay_s": -1.0,
        "result_conductor_loss_fraction": -0.2,
        "result_dielectric_loss_fraction": 2.0,
        "result_s11_magnitude": 2.0,
        "result_s21_magnitude": 2.0,
        "result_mesh_impedance_ohm": [10.0, 100.0],
        "accepted_microstrip_owner": "pcb/old",
        "accepted_microstrip_result_sha256": "a" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "microstrip_results_use_current_quasitem_impedance_permittivity_delay_loss_passivity_mesh_owner_and_result"
    ]


def test_v38_public_shielding_aperture_incident_transmitted_field_power_se_frequency_mesh_mismatch():
    summary = _summary_v38()
    row = summary["runs"][0][
        "shielding_aperture_orientation_polarization_field_power_se_frequency_probe_mesh_owner_result_identity"
    ]
    row.update({
        "aperture_generation": "shield-aperture-613",
        "se_generation": "shield-aperture-612",
        "result_generation": "shield-aperture-611",
        "result_aperture_plane": "unknown",
        "result_incident_polarization": "left_circular",
        "result_frequency_hz": [3.0e9, 1.0e9],
        "result_incident_field_v_m": [-1.0],
        "result_transmitted_field_v_m": [2.0],
        "result_transmitted_power_density_normalized": [4.0],
        "result_shielding_effectiveness_field_db": [-6.0],
        "result_shielding_effectiveness_power_db": [6.0],
        "result_probe_frame": "local_spherical",
        "result_mesh_selected_se_db": [40.0, 20.0],
        "accepted_shield_owner": "emc/old",
        "accepted_shield_result_sha256": "b" * 64,
    })
    result = nonlinear_inductance_sweep_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["runs"][0]["checks"][
        "shielding_results_use_current_aperture_polarization_field_power_se_frequency_probe_mesh_owner_and_result"
    ]


def test_v38_public_rejects_self_consistent_wrong_microstrip_impedance():
    summary = _summary_v38()
    for row in summary["runs"]:
        record = row[
            "microstrip_quasitem_geometry_permittivity_impedance_delay_loss_sparameter_frequency_mesh_owner_result_identity"
        ]
        record["characteristic_impedance_ohm"] *= 2.0
        record["result_characteristic_impedance_ohm"] = record["characteristic_impedance_ohm"]
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"


def test_v38_public_rejects_self_consistent_field_power_se_sign_error():
    summary = _summary_v38()
    for row in summary["runs"]:
        record = row[
            "shielding_aperture_orientation_polarization_field_power_se_frequency_probe_mesh_owner_result_identity"
        ]
        wrong = [-item for item in record["shielding_effectiveness_power_db"]]
        record["shielding_effectiveness_power_db"] = wrong
        record["result_shielding_effectiveness_power_db"] = wrong
    assert nonlinear_inductance_sweep_gate(summary)["status"] == "needs_attention"
