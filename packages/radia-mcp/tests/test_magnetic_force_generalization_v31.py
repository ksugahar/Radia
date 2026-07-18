from __future__ import annotations

from test_magnetic_force_generalization_v30 import _gate, _identity_v30


_PROMOTED_CASE_IDS = (
    "v31_public_nonlinear_coenergy_incremental_inductance_path_derivative_current_state_mismatch",
    "v31_public_lamination_anisotropy_fill_factor_orientation_frequency_loss_balance_mismatch",
)


def _identity_v31():
    identity = _identity_v30()
    generation = "incremental-inductance-181"
    identity[
        "nonlinear_coenergy_incremental_inductance_path_derivative_current_state_mesh_result_identity"
    ] = {
        "analysis_generation": generation,
        **{key: generation for key in (
            "current_path_generation", "magnetic_state_generation",
            "perturbation_generation", "derivative_generation",
            "circuit_generation", "mesh_generation", "result_generation",
        )},
        "nonlinear_material": True,
        "result_nonlinear_material": True,
        "current_path_a": [8.0, 9.0, 10.0, 11.0, 12.0],
        "result_current_path_a": [8.0, 9.0, 10.0, 11.0, 12.0],
        "current_state_index": 2,
        "result_current_state_index": 2,
        "nominal_current_a": 10.0,
        "result_nominal_current_a": 10.0,
        "perturbation_a": 0.1,
        "result_perturbation_a": 0.1,
        "current_samples_a": [9.9, 10.0, 10.1],
        "result_current_samples_a": [9.9, 10.0, 10.1],
        "flux_linkage_wb_turn": [0.792, 0.8, 0.808],
        "result_flux_linkage_wb_turn": [0.792, 0.8, 0.808],
        "derivative_rule": "symmetric_central_difference",
        "result_derivative_rule": "symmetric_central_difference",
        "incremental_inductance_h": 0.08,
        "result_incremental_inductance_h": 0.08,
        "circuit_name": "coil_a",
        "result_circuit_name": "coil_a",
        "magnetic_state_sha256": "1" * 64,
        "result_magnetic_state_sha256": "1" * 64,
        "mesh_sha256": "2" * 64,
        "result_mesh_sha256": "2" * 64,
        "result_sha256": "3" * 64,
        "accepted_result_sha256": "3" * 64,
    }
    generation = "lamination-loss-181"
    components = {"hysteresis_w": 3.2, "classical_eddy_w": 1.1, "excess_w": 0.4}
    coefficients = {"kh": 0.021, "kc": 0.00031, "ke": 0.0012}
    identity[
        "lamination_anisotropy_fill_orientation_frequency_loss_volume_balance_result_identity"
    ] = {
        "analysis_generation": generation,
        **{key: generation for key in (
            "anisotropy_generation", "fill_generation", "orientation_generation",
            "frequency_generation", "loss_generation", "volume_generation",
            "result_generation",
        )},
        "mu_axis_order": ["rolling", "transverse"],
        "result_mu_axis_order": ["rolling", "transverse"],
        "relative_permeability_axes": [1500.0, 40.0],
        "result_relative_permeability_axes": [1500.0, 40.0],
        "lamination_fill_factor": 0.95,
        "result_lamination_fill_factor": 0.95,
        "stacking_direction": "global_z",
        "result_stacking_direction": "global_z",
        "material_frame": "local_xy",
        "result_material_frame": "local_xy",
        "frequency_hz": 400.0,
        "result_frequency_hz": 400.0,
        "loss_coefficient_basis": "bertotti_three_term",
        "result_loss_coefficient_basis": "bertotti_three_term",
        "loss_coefficients": coefficients,
        "result_loss_coefficients": dict(coefficients),
        "gross_volume_m3": 0.001,
        "result_gross_volume_m3": 0.001,
        "active_iron_volume_m3": 0.00095,
        "result_active_iron_volume_m3": 0.00095,
        "loss_components_w": components,
        "result_loss_components_w": dict(components),
        "total_core_loss_w": 4.7,
        "result_total_core_loss_w": 4.7,
        "mesh_sha256": "4" * 64,
        "result_mesh_sha256": "4" * 64,
        "result_sha256": "5" * 64,
        "accepted_result_sha256": "5" * 64,
    }
    return identity


def test_v31_public_positive_incremental_inductance_and_lamination_loss():
    assert _gate(_identity_v31())["status"] == "ok"


def test_v31_public_nonlinear_coenergy_incremental_inductance_path_derivative_current_state_mismatch():
    identity = _identity_v31()
    record = identity[
        "nonlinear_coenergy_incremental_inductance_path_derivative_current_state_mesh_result_identity"
    ]
    record.update({
        "current_path_generation": "incremental-inductance-180",
        "result_current_state_index": 3,
        "result_current_samples_a": [9.9, 10.0, 10.2],
        "result_derivative_rule": "forward_difference",
        "result_incremental_inductance_h": 0.19,
        "result_magnetic_state_sha256": "a" * 64,
        "accepted_result_sha256": "b" * 64,
    })
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "nonlinear_incremental_inductance_uses_current_path_state_symmetric_derivative_circuit_mesh_and_result"
    ]


def test_v31_public_lamination_anisotropy_fill_factor_orientation_frequency_loss_balance_mismatch():
    identity = _identity_v31()
    record = identity[
        "lamination_anisotropy_fill_orientation_frequency_loss_volume_balance_result_identity"
    ]
    record.update({
        "anisotropy_generation": "lamination-loss-180",
        "result_mu_axis_order": ["transverse", "rolling"],
        "result_lamination_fill_factor": 1.0,
        "result_frequency_hz": 50.0,
        "result_active_iron_volume_m3": 0.001,
        "result_total_core_loss_w": 9.0,
        "accepted_result_sha256": "c" * 64,
    })
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "laminated_loss_uses_current_anisotropy_fill_orientation_frequency_volume_balance_and_result"
    ]
