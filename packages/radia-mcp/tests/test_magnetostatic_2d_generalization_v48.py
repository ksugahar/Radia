from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.electromagnetic_semantic_identity_v48 import (
    ELECTROSTATIC,
    INCREMENTAL,
    validate_public_identity,
)


PROMOTED_CASE_IDS = {
    "v48_public_incremental_permeability_frozen_bias_harmonic_phasor_operating_point_owner_mismatch",
    "v48_public_electrostatic_capacitance_charge_energy_voltage_sweep_identity_mismatch",
}


def _identity() -> dict[str, object]:
    incremental_generation = "incremental-bias-v48-901"
    electrostatic_generation = "capacitance-sweep-v48-901"
    tangent = {"material:core": [[0.0020, 0.0001], [0.0001, 0.0018]]}
    voltage = [0.0, 1.0, 2.0, 3.0]
    capacitance = 2.0e-9
    charge = [capacitance * value for value in voltage]
    energy = [0.5 * capacitance * value * value for value in voltage]
    owners = ["conductor:electrode-a" for _ in voltage]
    return {
        INCREMENTAL: {
            "generation": incremental_generation,
            "bias_generation": incremental_generation,
            "phasor_generation": incremental_generation,
            "tangent_generation": incremental_generation,
            "operating_point_generation": incremental_generation,
            "result_generation": incremental_generation,
            "frozen_bias_sha256": "6" * 64,
            "result_frozen_bias_sha256": "6" * 64,
            "harmonic_phasor_a": [1.5, -0.25],
            "result_harmonic_phasor_a": [1.5, -0.25],
            "material_tangent_h_per_m": tangent,
            "result_material_tangent_h_per_m": tangent,
            "operating_point_owner": "operating-point:incremental-v48-901",
            "result_operating_point_owner": "operating-point:incremental-v48-901",
            "result_sha256": "7" * 64,
            "accepted_result_sha256": "7" * 64,
        },
        ELECTROSTATIC: {
            "generation": electrostatic_generation,
            "voltage_generation": electrostatic_generation,
            "charge_generation": electrostatic_generation,
            "energy_generation": electrostatic_generation,
            "conductor_generation": electrostatic_generation,
            "result_generation": electrostatic_generation,
            "voltage_v": voltage,
            "result_voltage_v": voltage,
            "charge_c": charge,
            "result_charge_c": charge,
            "field_energy_j": energy,
            "result_field_energy_j": energy,
            "capacitance_f": capacitance,
            "result_capacitance_f": capacitance,
            "conductor_owner_rows": owners,
            "result_conductor_owner_rows": owners,
            "result_sha256": "8" * 64,
            "accepted_result_sha256": "8" * 64,
        },
    }


def test_v48_positive_incremental_and_electrostatic_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v48_incremental_owner_and_phasor_mutations_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[INCREMENTAL]["result_harmonic_phasor_a"] = [1.5, 0.25]
    identity[INCREMENTAL]["result_operating_point_owner"] = "operating-point:old"
    checks = validate_public_identity(identity)
    assert checks["v48_incremental_bias_phasor_tangent_owner"] is False


def test_v48_electrostatic_row_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[ELECTROSTATIC]["result_charge_c"] = [0.0, 2.0e-9, 6.0e-9, 4.0e-9]
    identity[ELECTROSTATIC]["result_conductor_owner_rows"][2] = "conductor:electrode-b"
    checks = validate_public_identity(identity)
    assert checks["v48_electrostatic_charge_energy_voltage_owner_closure"] is False
