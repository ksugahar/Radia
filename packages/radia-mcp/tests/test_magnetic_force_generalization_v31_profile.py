from __future__ import annotations

from radia_mcp.radia_ngsolve.magnetic_force_method_profile_gate import magnetic_force_method_profile_gate
from test_magnetic_force_generalization_v30_profile import _summary_v30

_PROMOTED_CASE_IDS = (
    "v31_public_bem_near_singular_gap_quadrature_normal_force_reciprocity_mismatch",
    "v31_public_hysteresis_minor_loop_state_remanence_return_point_energy_dissipation_mismatch",
)

def _summary_v31():
    summary = _summary_v30(); identity = summary["artifact_identity"]
    g = "near-gap-bem-351"
    identity["bem_near_singular_gap_quadrature_normal_order_force_reciprocity_geometry_result_identity"] = {
        "bem_generation": g, **{key: g for key in ("gap_bem_generation", "quadrature_bem_generation", "normal_bem_generation", "order_bem_generation", "force_bem_generation", "reciprocity_bem_generation", "geometry_bem_generation", "result_bem_generation")},
        "gap_m": 2e-5, "result_gap_m": 2e-5, "panel_size_m": 1e-3, "result_panel_size_m": 1e-3,
        "quadrature_policy": "gap_adaptive_duffy", "result_quadrature_policy": "gap_adaptive_duffy", "quadrature_order": 12, "result_quadrature_order": 12,
        "source_normal": [0., 0., 1.], "result_source_normal": [0., 0., 1.], "target_normal": [0., 0., -1.], "result_target_normal": [0., 0., -1.],
        "source_target_order": ["body_a", "body_b"], "result_source_target_order": ["body_a", "body_b"],
        "force_on_source_n": [0., 0., 5.], "result_force_on_source_n": [0., 0., 5.], "force_on_target_n": [0., 0., -5.], "result_force_on_target_n": [0., 0., -5.],
        "action_reaction_residual_n": 0., "result_action_reaction_residual_n": 0., "geometry_sha256": "1" * 64, "result_geometry_sha256": "1" * 64, "result_sha256": "2" * 64, "accepted_result_sha256": "2" * 64,
    }
    g = "minor-loop-351"
    identity["hysteresis_minor_loop_state_reversal_return_memory_remanence_energy_time_material_identity"] = {
        "loop_generation": g, **{key: g for key in ("state_loop_generation", "reversal_loop_generation", "memory_loop_generation", "remanence_loop_generation", "energy_loop_generation", "time_loop_generation", "material_loop_generation", "result_loop_generation")},
        "initial_state_sha256": "3" * 64, "result_initial_state_sha256": "3" * 64, "time_s": [0., .1, .2, .3, .4], "result_time_s": [0., .1, .2, .3, .4],
        "drive_h_a_per_m": [0., 100., 20., 100., 0.], "result_drive_h_a_per_m": [0., 100., 20., 100., 0.], "reversal_indices": [1, 2, 3], "result_reversal_indices": [1, 2, 3],
        "return_point_memory_closed": True, "result_return_point_memory_closed": True, "remanence_t": .35, "result_remanence_t": .35, "loop_energy_j_per_m3": 42., "result_loop_energy_j_per_m3": 42.,
        "material_owner": "steel-a", "result_material_owner": "steel-a", "result_sha256": "4" * 64, "accepted_result_sha256": "4" * 64,
    }
    return summary

def test_v31_public_positive_bem_and_minor_loop(): assert magnetic_force_method_profile_gate(_summary_v31())["status"] == "ok"

def test_v31_public_bem_near_singular_gap_quadrature_normal_force_reciprocity_mismatch():
    summary = _summary_v31(); record = summary["artifact_identity"]["bem_near_singular_gap_quadrature_normal_order_force_reciprocity_geometry_result_identity"]
    record.update({"gap_bem_generation": "old", "result_quadrature_policy": "fixed_gauss", "result_source_normal": [0., 0., -1.], "result_force_on_target_n": [0., 0., 4.], "result_action_reaction_residual_n": 9., "accepted_result_sha256": "8" * 64})
    result = magnetic_force_method_profile_gate(summary); assert result["status"] == "needs_attention"; assert not result["checks"]["bem_near_contact_force_uses_gap_adaptive_quadrature_normals_order_reciprocity_geometry_and_result"]

def test_v31_public_hysteresis_minor_loop_state_remanence_return_point_energy_dissipation_mismatch():
    summary = _summary_v31(); record = summary["artifact_identity"]["hysteresis_minor_loop_state_reversal_return_memory_remanence_energy_time_material_identity"]
    record.update({"state_loop_generation": "old", "result_time_s": [0., .2, .1, .3, .4], "result_return_point_memory_closed": False, "result_remanence_t": -.2, "result_loop_energy_j_per_m3": -42., "result_material_owner": "old", "accepted_result_sha256": "a" * 64})
    result = magnetic_force_method_profile_gate(summary); assert result["status"] == "needs_attention"; assert not result["checks"]["hysteresis_minor_loop_uses_initial_state_reversals_return_memory_remanence_energy_time_and_material"]
