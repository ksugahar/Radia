from __future__ import annotations

from test_femm_generalization_v19 import _gate
from test_femm_generalization_v20 import _identity_v20
from test_force_coenergy_gate import _quadratic_case


def _identity_v21(sample_count):
    identity = _identity_v20(sample_count)
    identity["weighted_stress_force_mask_material_mesh_generation_identity"] = {
        "solve_generation": "magnetostatic-solve-31",
        "weighted_stress_mask_solve_generation": "magnetostatic-solve-31",
        "material_label_solve_generation": "magnetostatic-solve-31",
        "mesh_field_solve_generation": "magnetostatic-solve-31",
        "force_integral_solve_generation": "magnetostatic-solve-31",
        "body_group_ids": [4, 5],
        "mask_body_group_ids": [4, 5],
        "material_labels": ["steel", "magnet"],
        "resolved_material_labels": ["steel", "magnet"],
        "weighted_force_n": [12.5, -3.2],
        "reported_force_n": [12.5, -3.2],
        "mask_table_sha256": "1" * 64,
        "force_mask_table_sha256": "1" * 64,
        "mesh_field_sha256": "2" * 64,
        "force_mesh_field_sha256": "2" * 64,
    }
    identity["harmonic_loss_phase_frequency_lamination_generation_identity"] = {
        "analysis_generation": "harmonic-loss-31",
        "phase_convention_analysis_generation": "harmonic-loss-31",
        "frequency_analysis_generation": "harmonic-loss-31",
        "lamination_analysis_generation": "harmonic-loss-31",
        "material_loss_analysis_generation": "harmonic-loss-31",
        "loss_result_analysis_generation": "harmonic-loss-31",
        "phase_convention": "exp(+jwt)",
        "loss_phase_convention": "exp(+jwt)",
        "frequency_hz": 400.0,
        "loss_frequency_hz": 400.0,
        "lamination_orientations": ["in_plane", "stacking_z"],
        "loss_lamination_orientations": ["in_plane", "stacking_z"],
        "material_loss_coefficients": [[0.02, 1.6], [0.01, 2.0]],
        "loss_material_coefficients": [[0.02, 1.6], [0.01, 2.0]],
        "material_loss_table_sha256": "3" * 64,
        "loss_material_table_sha256": "3" * 64,
    }
    return identity


def test_v21_public_positive_force_and_harmonic_loss_identity():
    positions, _, _ = _quadratic_case()
    assert _gate(_identity_v21(len(positions)))["status"] == "ok"


def test_v21_public_weighted_stress_force_mask_material_mesh_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v21(len(positions))
    identity["weighted_stress_force_mask_material_mesh_generation_identity"].update(
        {
            "weighted_stress_mask_solve_generation": "magnetostatic-solve-30",
            "material_label_solve_generation": "magnetostatic-solve-29",
            "mesh_field_solve_generation": "magnetostatic-solve-28",
            "mask_body_group_ids": [4, 6],
            "resolved_material_labels": ["air", "magnet"],
            "reported_force_n": [10.1, -1.4],
            "force_mask_table_sha256": "a" * 64,
            "force_mesh_field_sha256": "b" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "weighted_stress_force_uses_current_mask_materials_and_mesh_field"
    ]


def test_v21_public_harmonic_loss_phase_frequency_lamination_generation_mismatch():
    positions, _, _ = _quadratic_case()
    identity = _identity_v21(len(positions))
    identity["harmonic_loss_phase_frequency_lamination_generation_identity"].update(
        {
            "phase_convention_analysis_generation": "harmonic-loss-30",
            "frequency_analysis_generation": "harmonic-loss-29",
            "lamination_analysis_generation": "harmonic-loss-28",
            "material_loss_analysis_generation": "harmonic-loss-27",
            "loss_phase_convention": "exp(-jwt)",
            "loss_frequency_hz": 60.0,
            "loss_lamination_orientations": ["stacking_z", "in_plane"],
            "loss_material_coefficients": [[0.01, 2.0], [0.02, 1.6]],
            "loss_material_table_sha256": "c" * 64,
        }
    )
    result = _gate(identity)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "harmonic_loss_uses_current_phase_frequency_lamination_and_material_data"
    ]
