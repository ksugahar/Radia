from copy import deepcopy
import math

from radia_mcp.radia_ngsolve.energy_derivative_identity_v52 import MAGNET_TORQUE, VIRTUAL_FORCE, validate_public_identity


PROMOTED_CASE_IDS = {
    "v52_public_virtualdisplacement_force_meshperturbation_energy_owner_mismatch",
    "v52_public_magnet_torque_angularenergy_periodicunwrap_owner_mismatch",
}


def _generations(generation: str, fields: tuple[str, ...]) -> dict[str, str]:
    return {"generation": generation, **{field: generation for field in fields}}


def _identity() -> dict[str, object]:
    result = {"result_sha256": "f" * 64, "accepted_result_sha256": "f" * 64}
    wrapped = [350.0, 355.0, 0.0, 5.0, 10.0]
    unwrapped = [350.0, 355.0, 360.0, 365.0, 370.0]
    energy = [1.00, 0.95, 0.90, 0.85, 0.80]
    torque = -((energy[3] - energy[1]) / math.radians(unwrapped[3] - unwrapped[1]))
    return {
        VIRTUAL_FORCE: {
            **_generations("virtual-force-v52", ("mesh_generation", "energy_generation", "displacement_generation", "force_generation", "owner_generation", "result_generation")),
            "minus_mesh_sha256": "a" * 64, "result_minus_mesh_sha256": "a" * 64,
            "plus_mesh_sha256": "b" * 64, "result_plus_mesh_sha256": "b" * 64,
            "displacement_axis": [1.0, 0.0, 0.0], "result_displacement_axis": [1.0, 0.0, 0.0],
            "displacement_step_m": 1.0e-4, "result_displacement_step_m": 1.0e-4,
            "energy_minus_j": 1.002, "result_energy_minus_j": 1.002,
            "energy_plus_j": 0.998, "result_energy_plus_j": 0.998,
            "force_n": [20.0, 0.0, 0.0], "result_force_n": [20.0, 0.0, 0.0],
            "force_sign_convention": "negative_energy_gradient", "result_force_sign_convention": "negative_energy_gradient",
            "solution_owner": "solution:virtual-force-v52", "result_solution_owner": "solution:virtual-force-v52", **result,
        },
        MAGNET_TORQUE: {
            **_generations("magnet-torque-v52", ("angle_generation", "unwrap_generation", "energy_generation", "derivative_generation", "owner_generation", "result_generation")),
            "angles_wrapped_deg": wrapped, "result_angles_wrapped_deg": wrapped,
            "angles_unwrapped_deg": unwrapped, "result_angles_unwrapped_deg": unwrapped,
            "angular_energy_j": energy, "result_angular_energy_j": energy,
            "torque_at_center_nm": torque, "result_torque_at_center_nm": torque,
            "derivative_sign_convention": "negative_energy_gradient", "result_derivative_sign_convention": "negative_energy_gradient",
            "magnet_owner": "magnet:torque-v52", "result_magnet_owner": "magnet:torque-v52", **result,
        },
    }


def test_v52_positive_public_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v52_frozen_counterfactuals_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[VIRTUAL_FORCE]["result_force_n"] = [-20.0, 0.0, 0.0]
    identity[MAGNET_TORQUE]["result_angles_unwrapped_deg"] = [350.0, 355.0, 0.0, 5.0, 10.0]
    assert not all(validate_public_identity(identity).values())


def test_v52_self_consistent_wrong_derivative_sign_is_rejected() -> None:
    identity = deepcopy(_identity())
    identity[VIRTUAL_FORCE]["force_sign_convention"] = identity[VIRTUAL_FORCE]["result_force_sign_convention"] = "positive_energy_gradient"
    identity[MAGNET_TORQUE]["torque_at_center_nm"] = identity[MAGNET_TORQUE]["result_torque_at_center_nm"] = -identity[MAGNET_TORQUE]["torque_at_center_nm"]
    assert not all(validate_public_identity(identity).values())
