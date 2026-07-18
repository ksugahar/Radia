from __future__ import annotations

from radia_mcp.radia_ngsolve.regularized_trace_inverse_gate import (
    regularized_trace_inverse_path_gate,
)
from test_matlab_generalization_v28 import _summary_v28


_PROMOTED_CASE_IDS = (
    "v29_public_bem_near_singular_quadrature_distance_element_size_adaptive_order_reference_mismatch",
    "v29_public_fembem_energy_flux_reciprocity_interface_trace_orientation_frequency_mismatch",
)


def _summary_v29():
    summary = _summary_v28()
    generation = "near-singular-331"
    summary[
        "bem_near_singular_quadrature_distance_element_size_adaptive_order_reference_result_generation_identity"
    ] = {
        "quadrature_generation": generation,
        "target_quadrature_generation": generation,
        "geometry_quadrature_generation": generation,
        "order_quadrature_generation": generation,
        "map_quadrature_generation": generation,
        "kernel_quadrature_generation": generation,
        "reference_quadrature_generation": generation,
        "mesh_quadrature_generation": generation,
        "result_quadrature_generation": generation,
        "target_distance_m": 2.0e-4,
        "result_target_distance_m": 2.0e-4,
        "element_size_m": 1.0e-3,
        "result_element_size_m": 1.0e-3,
        "distance_size_ratio": 0.2,
        "result_distance_size_ratio": 0.2,
        "adaptive_order": 16,
        "result_adaptive_order": 16,
        "quadrature_rule": "adaptive-duffy-p1",
        "result_quadrature_rule": "adaptive-duffy-p1",
        "coordinate_map": "target-aligned-barycentric",
        "result_coordinate_map": "target-aligned-barycentric",
        "kernel": "helmholtz-single-layer-p1",
        "result_kernel": "helmholtz-single-layer-p1",
        "reference_integral_ri": [0.125, -0.03125],
        "computed_integral_ri": [0.1250000001, -0.0312499999],
        "relative_error": 1.1e-9,
        "relative_tolerance": 1.0e-7,
        "element_mesh_sha256": "1" * 64,
        "result_element_mesh_sha256": "1" * 64,
        "result_sha256": "2" * 64,
        "accepted_result_sha256": "2" * 64,
    }
    generation = "fembem-energy-331"
    summary[
        "fembem_energy_flux_reciprocity_interface_trace_orientation_frequency_incident_result_generation_identity"
    ] = {
        "coupling_generation": generation,
        "trace_coupling_generation": generation,
        "normal_coupling_generation": generation,
        "frequency_coupling_generation": generation,
        "incident_coupling_generation": generation,
        "reciprocity_coupling_generation": generation,
        "energy_coupling_generation": generation,
        "result_coupling_generation": generation,
        "interface_trace_basis": "p1-nodal-boundary-trace",
        "result_interface_trace_basis": "p1-nodal-boundary-trace",
        "interface_trace_shape": [48, 120],
        "result_interface_trace_shape": [48, 120],
        "normal_orientation": "volume-outward",
        "result_normal_orientation": "volume-outward",
        "frequency_hz": 800.0,
        "result_frequency_hz": 800.0,
        "incident_field_sha256": "3" * 64,
        "result_incident_field_sha256": "3" * 64,
        "reciprocity_pair_ids": ["source-a/receiver-b", "source-b/receiver-a"],
        "result_reciprocity_pair_ids": ["source-a/receiver-b", "source-b/receiver-a"],
        "reciprocity_values_ri": [[0.5, -0.1], [0.5000000002, -0.1000000001]],
        "result_reciprocity_values_ri": [[0.5, -0.1], [0.5000000002, -0.1000000001]],
        "reciprocity_relative_error": 4.4e-10,
        "reciprocity_relative_tolerance": 1.0e-7,
        "fem_outward_power_w": 1.25,
        "bem_radiated_power_w": 1.249999999,
        "energy_flux_relative_error": 8.0e-10,
        "energy_flux_relative_tolerance": 1.0e-7,
        "coupled_result_sha256": "4" * 64,
        "accepted_coupled_result_sha256": "4" * 64,
    }
    return summary


def test_v29_public_positive_near_singular_and_energy_reciprocity_identities() -> None:
    assert regularized_trace_inverse_path_gate(_summary_v29())["status"] == "ok"


def test_v29_public_bem_near_singular_quadrature_distance_element_size_adaptive_order_reference_mismatch() -> None:
    summary = _summary_v29()
    identity = summary[
        "bem_near_singular_quadrature_distance_element_size_adaptive_order_reference_result_generation_identity"
    ]
    identity.update(
        {
            "target_quadrature_generation": "near-singular-330",
            "mesh_quadrature_generation": "near-singular-329",
            "result_target_distance_m": 1.0e-3,
            "result_element_size_m": 2.0e-3,
            "result_distance_size_ratio": 0.5,
            "result_adaptive_order": 4,
            "result_quadrature_rule": "triangle-centroid",
            "result_coordinate_map": "global-cartesian",
            "result_kernel": "laplace-p0",
            "computed_integral_ri": [0.2, 0.0],
            "relative_error": 0.5,
            "result_element_mesh_sha256": "a" * 64,
            "accepted_result_sha256": "b" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "bem_near_singular_quadrature_uses_current_distance_size_order_map_kernel_reference_mesh_and_result"
    ]


def test_v29_public_fembem_energy_flux_reciprocity_interface_trace_orientation_frequency_mismatch() -> None:
    summary = _summary_v29()
    identity = summary[
        "fembem_energy_flux_reciprocity_interface_trace_orientation_frequency_incident_result_generation_identity"
    ]
    identity.update(
        {
            "trace_coupling_generation": "fembem-energy-330",
            "frequency_coupling_generation": "fembem-energy-329",
            "result_interface_trace_basis": "p0-cell",
            "result_interface_trace_shape": [47, 120],
            "result_normal_orientation": "volume-inward",
            "result_frequency_hz": 1000.0,
            "result_incident_field_sha256": "c" * 64,
            "result_reciprocity_pair_ids": ["source-b/receiver-a", "source-a/receiver-b"],
            "result_reciprocity_values_ri": [[0.4, 0.0], [0.7, 0.0]],
            "reciprocity_relative_error": 0.3,
            "bem_radiated_power_w": 0.5,
            "energy_flux_relative_error": 0.6,
            "accepted_coupled_result_sha256": "d" * 64,
        }
    )
    result = regularized_trace_inverse_path_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "fembem_energy_flux_and_reciprocity_use_current_trace_normal_frequency_incident_field_and_result"
    ]
