from __future__ import annotations

import math

from test_magnetic_force_generalization_v31 import _gate
from test_magnetic_force_generalization_v36 import _identity_v36

_PROMOTED_CASE_IDS = (
    "v37_public_axisymmetric_to_3d_force_revolution_volume_energy_direction_owner_mismatch",
    "v37_public_harmonic_circuit_impedance_voltage_current_complex_power_loss_owner_mismatch",
)


def _identity_v37():
    identity = _identity_v36()
    generation = "axisym-revolution-246"
    factor = 2.0 * math.pi
    area = 0.01
    radius = 0.05
    identity[
        "axisymmetric_3d_force_revolution_volume_energy_coenergy_direction_displacement_field_mesh_result_identity"
    ] = {
        "revolution_generation": generation,
        **{
            key: generation
            for key in (
                "factor_generation",
                "volume_generation",
                "energy_generation",
                "coenergy_generation",
                "direction_generation",
                "displacement_generation",
                "field_generation",
                "mesh_generation",
                "result_generation",
            )
        },
        "revolution_factor": factor,
        "result_revolution_factor": factor,
        "meridional_area_m2": area,
        "result_meridional_area_m2": area,
        "centroid_radius_m": radius,
        "result_centroid_radius_m": radius,
        "swept_volume_m3": factor * radius * area,
        "result_swept_volume_m3": factor * radius * area,
        "axisymmetric_energy_j_per_rad": 2.0,
        "result_axisymmetric_energy_j_per_rad": 2.0,
        "revolved_energy_j": factor * 2.0,
        "result_revolved_energy_j": factor * 2.0,
        "axisymmetric_coenergy_j_per_rad": 2.1,
        "result_axisymmetric_coenergy_j_per_rad": 2.1,
        "revolved_coenergy_j": factor * 2.1,
        "result_revolved_coenergy_j": factor * 2.1,
        "virtual_displacement_m": [0.0, 0.0, 1.0e-4],
        "result_virtual_displacement_m": [0.0, 0.0, 1.0e-4],
        "force_direction_unit": [0.0, 0.0, 1.0],
        "result_force_direction_unit": [0.0, 0.0, 1.0],
        "force_n": [0.0, 0.0, 12.0],
        "result_force_n": [0.0, 0.0, 12.0],
        "field_owner": "axisym/field-246",
        "accepted_field_owner": "axisym/field-246",
        "mesh_sha256": "1" * 64,
        "accepted_mesh_sha256": "1" * 64,
        "force_result_sha256": "2" * 64,
        "accepted_force_result_sha256": "2" * 64,
    }
    generation = "harmonic-circuit-246"
    voltage = [10.0, 2.0]
    current = [2.0, -1.0]
    denominator = current[0] ** 2 + current[1] ** 2
    impedance = [
        (voltage[0] * current[0] + voltage[1] * current[1]) / denominator,
        (voltage[1] * current[0] - voltage[0] * current[1]) / denominator,
    ]
    power = [
        0.5 * (voltage[0] * current[0] + voltage[1] * current[1]),
        0.5 * (voltage[1] * current[0] - voltage[0] * current[1]),
    ]
    identity[
        "harmonic_circuit_voltage_current_impedance_complex_power_copper_field_loss_rms_owner_result_identity"
    ] = {
        "circuit_generation": generation,
        **{
            key: generation
            for key in (
                "voltage_generation",
                "current_generation",
                "impedance_generation",
                "power_generation",
                "loss_generation",
                "rms_generation",
                "owner_generation",
                "result_generation",
            )
        },
        "voltage_peak_phasor_v": voltage,
        "result_voltage_peak_phasor_v": voltage,
        "current_peak_phasor_a": current,
        "result_current_peak_phasor_a": current,
        "impedance_ohm": impedance,
        "result_impedance_ohm": impedance,
        "complex_power_va": power,
        "result_complex_power_va": power,
        "copper_loss_w": 7.0,
        "result_copper_loss_w": 7.0,
        "field_loss_w": 2.0,
        "result_field_loss_w": 2.0,
        "phasor_convention": "peak_cosine",
        "result_phasor_convention": "peak_cosine",
        "circuit_owner": "circuit:winding-246",
        "accepted_circuit_owner": "circuit:winding-246",
        "circuit_result_sha256": "3" * 64,
        "accepted_circuit_result_sha256": "3" * 64,
    }
    return identity


def test_v37_public_positive_axisymmetric_revolution_and_harmonic_circuit_closure():
    assert _gate(_identity_v37())["status"] == "ok"


def test_v37_public_axisymmetric_to_3d_force_revolution_volume_energy_direction_owner_mismatch():
    identity = _identity_v37()
    row = identity[
        "axisymmetric_3d_force_revolution_volume_energy_coenergy_direction_displacement_field_mesh_result_identity"
    ]
    row.update(
        {
            "factor_generation": "axisym-revolution-245",
            "direction_generation": "axisym-revolution-244",
            "result_generation": "axisym-revolution-243",
            "result_revolution_factor": math.pi,
            "result_swept_volume_m3": -1.0,
            "result_revolved_energy_j": 2.0,
            "result_revolved_coenergy_j": 2.1,
            "result_virtual_displacement_m": [1.0e-4, 0.0, 0.0],
            "result_force_direction_unit": [-1.0, 0.0, 0.0],
            "result_force_n": [-12.0, 0.0, 0.0],
            "accepted_field_owner": "stale:field",
            "accepted_mesh_sha256": "a" * 64,
            "accepted_force_result_sha256": "b" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "axisymmetric_3d_forces_use_current_revolution_volume_energy_coenergy_direction_displacement_field_mesh_and_result"
    ]


def test_v37_public_harmonic_circuit_impedance_voltage_current_complex_power_loss_owner_mismatch():
    identity = _identity_v37()
    row = identity[
        "harmonic_circuit_voltage_current_impedance_complex_power_copper_field_loss_rms_owner_result_identity"
    ]
    row.update(
        {
            "voltage_generation": "harmonic-circuit-245",
            "power_generation": "harmonic-circuit-244",
            "result_generation": "harmonic-circuit-243",
            "result_voltage_peak_phasor_v": [10.0, -2.0],
            "result_current_peak_phasor_a": [-2.0, 1.0],
            "result_impedance_ohm": [-4.0, 1.0],
            "result_complex_power_va": [-9.0, -7.0],
            "result_copper_loss_w": -7.0,
            "result_field_loss_w": 20.0,
            "result_phasor_convention": "rms_sine",
            "accepted_circuit_owner": "stale:circuit",
            "accepted_circuit_result_sha256": "c" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "harmonic_circuits_use_current_voltage_current_impedance_power_losses_rms_owner_and_result"
    ]


def test_v37_public_rejects_self_consistent_wrong_revolution_factor():
    identity = _identity_v37()
    row = identity[
        "axisymmetric_3d_force_revolution_volume_energy_coenergy_direction_displacement_field_mesh_result_identity"
    ]
    row["revolution_factor"] = row["result_revolution_factor"] = math.pi
    row["swept_volume_m3"] = row["result_swept_volume_m3"] = (
        math.pi * row["centroid_radius_m"] * row["meridional_area_m2"]
    )
    row["revolved_energy_j"] = row["result_revolved_energy_j"] = math.pi * 2.0
    row["revolved_coenergy_j"] = row["result_revolved_coenergy_j"] = math.pi * 2.1
    assert _gate(identity)["status"] == "needs_attention"


def test_v37_public_rejects_self_consistent_harmonic_power_loss_imbalance():
    identity = _identity_v37()
    row = identity[
        "harmonic_circuit_voltage_current_impedance_complex_power_copper_field_loss_rms_owner_result_identity"
    ]
    row["copper_loss_w"] = row["result_copper_loss_w"] = 1.0
    row["field_loss_w"] = row["result_field_loss_w"] = 1.0
    assert _gate(identity)["status"] == "needs_attention"
