from __future__ import annotations

from copy import deepcopy

from radia_mcp.radia_ngsolve.demag_virtual_work_identity_v49 import DEMAG, VIRTUAL_WORK, validate_public_identity


PROMOTED_CASE_IDS = {
    "v49_public_nonlinear_demag_recoil_branch_temperature_loadstep_result_owner_mismatch",
    "v49_public_virtual_work_displacement_frame_mesh_state_energy_owner_mismatch",
}


def _identity() -> dict[str, object]:
    demag_generation = "demag-recoil-v49-901"
    virtual_work_generation = "virtual-work-v49-901"
    recoil = [[0.9, -120000.0], [1.0, -80000.0], [1.1, -40000.0]]
    displacement = [-0.0001, 0.0, 0.0001]
    energy = [1.00012, 1.00000, 0.99988]
    return {
        DEMAG: {
            "generation": demag_generation,
            "branch_generation": demag_generation,
            "temperature_generation": demag_generation,
            "loadstep_generation": demag_generation,
            "result_generation": demag_generation,
            "recoil_branch_t_a_per_m": recoil,
            "result_recoil_branch_t_a_per_m": recoil,
            "temperature_c": 140.0,
            "result_temperature_c": 140.0,
            "load_step": 12,
            "result_load_step": 12,
            "magnet_owner": "magnet:rotor-v49-901",
            "result_magnet_owner": "magnet:rotor-v49-901",
            "result_sha256": "1" * 64,
            "accepted_result_sha256": "1" * 64,
        },
        VIRTUAL_WORK: {
            "generation": virtual_work_generation,
            "displacement_generation": virtual_work_generation,
            "frame_generation": virtual_work_generation,
            "mesh_generation": virtual_work_generation,
            "energy_generation": virtual_work_generation,
            "result_generation": virtual_work_generation,
            "displacement_m": displacement,
            "result_displacement_m": displacement,
            "displacement_frame": "frame:global-x",
            "result_displacement_frame": "frame:global-x",
            "mesh_state_sha256": "2" * 64,
            "result_mesh_state_sha256": "2" * 64,
            "energy_j": energy,
            "result_energy_j": energy,
            "force_owner": "force:body-v49-901",
            "result_force_owner": "force:body-v49-901",
            "result_sha256": "3" * 64,
            "accepted_result_sha256": "3" * 64,
        },
    }


def test_v49_positive_demag_and_virtual_work_artifacts_are_accepted() -> None:
    assert all(validate_public_identity(_identity()).values())


def test_v49_demag_branch_temperature_step_digest_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[DEMAG]["result_recoil_branch_t_a_per_m"] = list(reversed(identity[DEMAG]["recoil_branch_t_a_per_m"]))
    identity[DEMAG]["result_temperature_c"] = 20.0
    identity[DEMAG]["accepted_result_sha256"] = "8" * 64
    identity[DEMAG]["result_magnet_owner"] = "magnet:old"
    assert validate_public_identity(identity)["demag_v49_recoil_temperature_loadstep_digest_owner"] is False


def test_v49_virtual_work_frame_mesh_energy_and_owner_mutations_are_rejected() -> None:
    identity = deepcopy(_identity())
    identity[VIRTUAL_WORK]["result_displacement_frame"] = "frame:local-y"
    identity[VIRTUAL_WORK]["result_mesh_state_sha256"] = "9" * 64
    identity[VIRTUAL_WORK]["result_energy_j"] = list(reversed(identity[VIRTUAL_WORK]["energy_j"]))
    identity[VIRTUAL_WORK]["result_force_owner"] = "force:old"
    assert validate_public_identity(identity)["force_v49_virtual_work_displacement_frame_mesh_energy_owner"] is False
