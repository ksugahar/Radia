from __future__ import annotations

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v25 import _summary_v25


_PROMOTED_CASE_IDS = (
    "v26_public_cq_contour_radius_timestep_laplace_branch_transfer_operator_inverse_transform_mismatch",
    "v26_public_fembem_trace_map_normal_material_wavenumber_coupling_matrix_mesh_generation_mismatch",
)


def _summary_v26():
    summary = _summary_v25()
    generation = "cq-transfer-301"
    summary[
        "cq_contour_radius_timestep_laplace_branch_transfer_operator_inverse_transform_generation_identity"
    ] = {
        "cq_generation": generation,
        "contour_cq_generation": generation,
        "timestep_cq_generation": generation,
        "laplace_branch_cq_generation": generation,
        "transfer_operator_cq_generation": generation,
        "inverse_transform_cq_generation": generation,
        "result_cq_generation": generation,
        "contour_radius": 0.95,
        "result_contour_radius": 0.95,
        "time_step_s": 1.0e-4,
        "result_time_step_s": 1.0e-4,
        "laplace_branch": "principal_sqrt_outgoing",
        "result_laplace_branch": "principal_sqrt_outgoing",
        "laplace_points_ri": [[100.0, 0.0], [80.0, 20.0], [60.0, 35.0], [80.0, -20.0]],
        "result_laplace_points_ri": [[100.0, 0.0], [80.0, 20.0], [60.0, 35.0], [80.0, -20.0]],
        "transfer_operator_id": "helmholtz_calderon_p1",
        "result_transfer_operator_id": "helmholtz_calderon_p1",
        "transfer_operator_sha256": "1" * 64,
        "result_transfer_operator_sha256": "1" * 64,
        "inverse_transform": "fft_conjugate_symmetric",
        "result_inverse_transform": "fft_conjugate_symmetric",
        "time_history_sha256": "2" * 64,
        "reported_time_history_sha256": "2" * 64,
    }
    generation = "fembem-coupling-301"
    summary[
        "fembem_trace_map_normal_material_wavenumber_coupling_matrix_mesh_generation_identity"
    ] = {
        "coupling_generation": generation,
        "trace_map_coupling_generation": generation,
        "normal_coupling_generation": generation,
        "material_coupling_generation": generation,
        "wavenumber_coupling_generation": generation,
        "matrix_coupling_generation": generation,
        "mesh_coupling_generation": generation,
        "result_coupling_generation": generation,
        "trace_map_sha256": "3" * 64,
        "result_trace_map_sha256": "3" * 64,
        "normal_orientation": "volume_outward",
        "result_normal_orientation": "volume_outward",
        "normal_field_sha256": "4" * 64,
        "result_normal_field_sha256": "4" * 64,
        "fluid_density_kg_m3": 1.2,
        "result_fluid_density_kg_m3": 1.2,
        "sound_speed_m_s": 343.0,
        "result_sound_speed_m_s": 343.0,
        "wavenumber_ri_m_inv": [18.318324511, 0.02],
        "result_wavenumber_ri_m_inv": [18.318324511, 0.02],
        "coupling_matrix_sha256": "5" * 64,
        "result_coupling_matrix_sha256": "5" * 64,
        "volume_mesh_sha256": "6" * 64,
        "result_volume_mesh_sha256": "6" * 64,
        "boundary_mesh_sha256": "7" * 64,
        "result_boundary_mesh_sha256": "7" * 64,
        "coupled_result_sha256": "8" * 64,
        "reported_coupled_result_sha256": "8" * 64,
    }
    return summary


def test_v26_public_positive_cq_transfer_and_fembem_coupling_identity() -> None:
    assert regularized_trace_inverse_path_gate(_summary_v26())["status"] == "ok"


def test_v26_public_rejects_cq_transfer_identity_mismatch() -> None:
    summary = _summary_v26()
    identity = summary[
        "cq_contour_radius_timestep_laplace_branch_transfer_operator_inverse_transform_generation_identity"
    ]
    identity.update(
        {
            "contour_cq_generation": "cq-transfer-300",
            "result_contour_radius": 0.8,
            "result_time_step_s": 2.0e-4,
            "result_laplace_branch": "negative_sqrt_incoming",
            "result_laplace_points_ri": [[-100.0, 0.0]],
            "result_transfer_operator_id": "laplace_single_layer_p0",
            "result_transfer_operator_sha256": "1" * 63 + "2",
            "result_inverse_transform": "direct_dft_unsigned",
            "reported_time_history_sha256": "2" * 63 + "3",
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "cq_history_uses_current_contour_timestep_branch_transfer_and_inverse_transform"
    ]


def test_v26_public_rejects_fembem_coupling_identity_mismatch() -> None:
    summary = _summary_v26()
    identity = summary[
        "fembem_trace_map_normal_material_wavenumber_coupling_matrix_mesh_generation_identity"
    ]
    identity.update(
        {
            "trace_map_coupling_generation": "fembem-coupling-300",
            "result_trace_map_sha256": "3" * 63 + "4",
            "result_normal_orientation": "volume_inward",
            "result_fluid_density_kg_m3": 1000.0,
            "result_sound_speed_m_s": 1480.0,
            "result_wavenumber_ri_m_inv": [4.2, -0.02],
            "result_coupling_matrix_sha256": "5" * 63 + "6",
            "result_volume_mesh_sha256": "6" * 63 + "7",
            "result_boundary_mesh_sha256": "7" * 63 + "8",
            "reported_coupled_result_sha256": "8" * 63 + "9",
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "fembem_coupling_uses_current_trace_normals_material_wavenumber_matrices_and_mesh"
    ]
