from __future__ import annotations

from radia_mcp.radia_ngsolve.magnetic_force_v44_identity import validate_public_identity


_DYNAMIC = "magneticbearing_dynamicstiffness_phase_damping_force_power_stability_mesh_result_identity"
_DEMAG = "demag_minorloop_fieldpath_remanence_loss_energy_temperature_material_mesh_result_identity"
_PROMOTED_CASE_IDS = (
    "v44_public_magneticbearing_dynamicstiffness_phase_damping_force_power_stability_result_mismatch",
    "v44_public_demag_minorloop_fieldpath_remanence_loss_energy_temperature_mesh_mismatch",
)


def _identity() -> dict:
    return {
        _DYNAMIC: {
            "bearing_dynamic_generation": "bearing-dynamic-844",
            **{key: "bearing-dynamic-844" for key in ("dynamic_stiffness_generation", "phase_generation", "damping_generation", "force_generation", "power_generation", "stability_generation", "mesh_generation", "result_generation")},
            "frequency_hz": [100.0, 200.0, 300.0], "result_frequency_hz": [100.0, 200.0, 300.0],
            "dynamic_stiffness_n_per_m": [100.0, 110.0, 120.0], "result_dynamic_stiffness_n_per_m": [100.0, 110.0, 120.0],
            "phase_deg": [0.0, 10.0, 20.0], "result_phase_deg": [0.0, 10.0, 20.0],
            "damping_n_s_per_m": [2.0, 2.5, 3.0], "result_damping_n_s_per_m": [2.0, 2.5, 3.0],
            "force_n": [10.0, 11.0, 12.0], "result_force_n": [10.0, 11.0, 12.0],
            "power_w": [1.0, 1.2, 1.4], "result_power_w": [1.0, 1.2, 1.4],
            "stability_sign": "stable", "result_stability_sign": "stable",
            "mesh_owner": "mesh:bearing-dynamic-844", "result_mesh_owner": "mesh:bearing-dynamic-844",
            "bearing_dynamic_result_sha256": "9" * 64, "accepted_bearing_dynamic_result_sha256": "9" * 64,
        },
        _DEMAG: {
            "demag_minor_generation": "demag-minor-844",
            **{key: "demag-minor-844" for key in ("fieldpath_generation", "remanence_generation", "branch_generation", "loss_generation", "energy_generation", "temperature_generation", "material_generation", "mesh_generation", "result_generation")},
            "field_path_a_per_m": [-1.0, 0.5, -0.5, 1.0], "result_field_path_a_per_m": [-1.0, 0.5, -0.5, 1.0],
            "magnetization_a_per_m": [-0.8, 0.4, -0.3, 0.9], "result_magnetization_a_per_m": [-0.8, 0.4, -0.3, 0.9],
            "remanence_a_per_m": 0.8, "result_remanence_a_per_m": 0.8,
            "coercivity_a_per_m": 0.4, "result_coercivity_a_per_m": 0.4,
            "loss_energy_j": 0.16, "result_loss_energy_j": 0.16,
            "temperature_k": 293.15, "result_temperature_k": 293.15,
            "material_owner": "material:demag-minor-844", "result_material_owner": "material:demag-minor-844",
            "mesh_owner": "mesh:demag-minor-844", "result_mesh_owner": "mesh:demag-minor-844",
            "demag_minor_result_sha256": "a" * 64, "accepted_demag_minor_result_sha256": "a" * 64,
        },
    }


def test_v44_public_dynamic_and_demag_identity_positive() -> None:
    result = validate_public_identity(_identity())
    assert result == {
        "magnetic_force_v44_dynamic_bearing_identity": True,
        "magnetic_force_v44_demag_minorloop_identity": True,
    }


def test_v44_public_identity_rejects_phase_and_temperature_mutations() -> None:
    identity = _identity()
    identity[_DYNAMIC]["result_phase_deg"] = [0.0, -10.0, -20.0]
    identity[_DEMAG]["result_temperature_k"] = 350.0
    result = validate_public_identity(identity)
    assert result["magnetic_force_v44_dynamic_bearing_identity"] is False
    assert result["magnetic_force_v44_demag_minorloop_identity"] is False
