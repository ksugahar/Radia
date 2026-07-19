from __future__ import annotations

from radia_mcp.radia_ngsolve.magnetic_force_v44_identity import validate_public_identity


_PROMOTED_CASE_IDS = (
    "v45_public_magnetic_bearing_dynamic_stiffness_phase_damping_force_power_stability_owner_mismatch",
    "v45_public_demagnetization_minor_loop_field_path_remanence_coercivity_loss_temperature_mismatch",
)


def _identity():
    generation = "test-845"
    return {
        "v45_public_magnetic_bearing_dynamic_stiffness_phase_damping_force_power_stability_owner_mismatch": {
            "generation": generation, **{key: generation for key in ("stiffness_generation", "phase_generation", "damping_generation", "force_generation", "power_generation", "stability_generation", "mesh_generation", "result_generation")},
            "frequency_hz": [100.0, 200.0], "result_frequency_hz": [100.0, 200.0], "dynamic_stiffness_n_per_m": [100.0, 110.0], "result_dynamic_stiffness_n_per_m": [100.0, 110.0], "phase_deg": [0.0, 10.0], "result_phase_deg": [0.0, 10.0], "damping_n_s_per_m": [2.0, 2.5], "result_damping_n_s_per_m": [2.0, 2.5], "force_n": [10.0, 11.0], "result_force_n": [10.0, 11.0], "power_w": [1.0, 1.2], "result_power_w": [1.0, 1.2], "stability_sign": "stable", "result_stability_sign": "stable", "mesh_owner": "mesh:test", "result_mesh_owner": "mesh:test", "result_sha256": "a" * 64, "accepted_result_sha256": "a" * 64,
        },
        "v45_public_demagnetization_minor_loop_field_path_remanence_coercivity_loss_temperature_mismatch": {
            "generation": generation, **{key: generation for key in ("field_path_generation", "remanence_generation", "coercivity_generation", "loss_generation", "temperature_generation", "material_generation", "result_generation")}, "field_path_a_per_m": [-1.0, 0.5, -0.5, 1.0], "result_field_path_a_per_m": [-1.0, 0.5, -0.5, 1.0], "remanence_a_per_m": 0.8, "result_remanence_a_per_m": 0.8, "coercivity_a_per_m": 0.4, "result_coercivity_a_per_m": 0.4, "loss_energy_j": 0.16, "result_loss_energy_j": 0.16, "temperature_k": 293.15, "result_temperature_k": 293.15, "material_owner": "material:test", "result_material_owner": "material:test", "result_sha256": "b" * 64, "accepted_result_sha256": "b" * 64,
        },
    }


def test_v45_elf_public_identity_accepts_closed_artifacts():
    checks = validate_public_identity(_identity())
    assert checks and all(checks.values())


def test_v45_elf_public_identity_rejects_dynamic_phase_mutation():
    identity = _identity()
    identity["v45_public_magnetic_bearing_dynamic_stiffness_phase_damping_force_power_stability_owner_mismatch"]["result_phase_deg"] = [0.0, -10.0]
    checks = validate_public_identity(identity)
    assert checks and not all(checks.values())
