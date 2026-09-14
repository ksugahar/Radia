from __future__ import annotations

import asyncio
import json
import os
import sys
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from radia_mcp.radia_ngsolve.field_profile_gate import (
    dual_formulation_symmetric_field_profile_gate,
    nonlinear_constitutive_response_parity_gate,
    nonlinear_magnetic_field_energy_parity_gate,
    nonlinear_magnetic_refinement_energy_gate,
    nonlinear_magnetic_spatial_evidence_gate,
)
from radia_mcp.radia_ngsolve.server import (
    nonlinear_magnetic_field_energy_parity_gate as mcp_nonlinear_parity_gate,
    nonlinear_constitutive_response_parity_gate as mcp_constitutive_parity_gate,
    nonlinear_magnetic_refinement_energy_gate as mcp_nonlinear_refinement_gate,
    nonlinear_magnetic_spatial_evidence_gate as mcp_nonlinear_gate,
)


def _summary():
    return {
        "component_id": "B_axis",
        "field_unit": "T",
        "axis_unit": "m",
        "sample_count": 201,
        "axis_min": -1.0,
        "axis_max": 1.0,
        "center_axis": 0.0,
        "center_value_a": -5.61756e-9,
        "center_value_b": -5.61627e-9,
        "profile_relative_l2_difference": 9.42e-4,
        "center_relative_difference": 2.31e-4,
        "symmetry_relative_a": 2.51e-4,
        "symmetry_relative_b": 2.63e-4,
    }


def test_dual_formulation_profile_gate_accepts_agreement_and_symmetry():
    result = dual_formulation_symmetric_field_profile_gate(_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())


def test_dual_formulation_profile_gate_rejects_single_point_only_or_asymmetry():
    bad = _summary()
    bad["sample_count"] = 1
    bad["profile_relative_l2_difference"] = 0.2
    bad["symmetry_relative_b"] = 0.3
    result = dual_formulation_symmetric_field_profile_gate(bad)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "sample_count_sufficient",
        "profile_formulations_agree",
        "formulation_b_is_symmetric",
    }


def test_dual_formulation_profile_gate_rejects_trivial_zero_field():
    bad = _summary()
    bad["center_value_a"] = 0.0
    bad["center_value_b"] = 0.0
    assert dual_formulation_symmetric_field_profile_gate(bad)["checks"]["center_field_nonzero"] is False


def _nonlinear_volume_summary():
    identity = {
        "observable_id": "iron_volume_B",
        "field_unit": "T",
        "coordinate_system": "right-handed Cartesian",
    }
    return {
        **identity,
        "nonlinear": True,
        "solver_converged": True,
        "linear_reference_only": False,
        "response_order": 1,
        "material_update_order": 0,
        "spatial_observable": "volume_integral",
        "integration_order": 8,
        "sample_count": 0,
        "volume_m3": 1.0e-6,
        "average_field_T": [0.0, 0.0, -1.50],
        "rms_magnitude_T": 1.62,
        "reference": {
            **identity,
            "average_field_T": [0.0, 0.0, -1.48],
            "rms_magnitude_T": 1.60,
        },
    }


def test_nonlinear_magnetic_spatial_gate_accepts_matched_volume_evidence():
    result = nonlinear_magnetic_spatial_evidence_gate(_nonlinear_volume_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())
    assert result["metrics"]["average_vector_relative_difference"] < 0.02


def test_nonlinear_magnetic_spatial_gate_rejects_point_only_and_order_mismatch():
    bad = _nonlinear_volume_summary()
    bad.update(
        response_order=2,
        material_update_order=0,
        spatial_observable="point",
        integration_order=0,
        sample_count=1,
    )
    result = nonlinear_magnetic_spatial_evidence_gate(bad)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "response_material_orders_compatible",
        "spatial_observable_supported",
        "spatial_sampling_sufficient",
    }


def test_nonlinear_magnetic_spatial_gate_rejects_unconverged_or_linear_only():
    bad = _nonlinear_volume_summary()
    bad["solver_converged"] = False
    bad["linear_reference_only"] = True
    result = nonlinear_magnetic_spatial_evidence_gate(bad)
    assert set(result["issues"]) >= {
        "solver_converged",
        "not_linear_only_evidence",
    }


def test_nonlinear_magnetic_spatial_gate_rejects_rms_disagreement_and_wraps_mcp():
    bad = _nonlinear_volume_summary()
    bad["rms_magnitude_T"] = 2.1
    result = nonlinear_magnetic_spatial_evidence_gate(bad)
    assert result["checks"]["rms_magnitude_matches_reference"] is False
    wrapped = json.loads(mcp_nonlinear_gate(json.dumps(bad)))
    assert wrapped["status"] == "needs_attention"


def _nonlinear_refinement_summary():
    def level(mesh_size, average_z, rms, energy):
        identity = {
            "material_domain": "air|iron",
            "coordinate_system": "right-handed Cartesian",
            "unit_system": "SI",
            "nonlinear_state_id": "bh-state-17",
            "mesh_topology_geometry_sha256": "1" * 64,
            "bh_table_sha256": "2" * 64,
            "material_state_identity_sha256": "3" * 64,
            "solution_sha256": "4" * 64,
            "refinement_parent_identity_sha256": "5" * 64,
        }
        return {
            "mesh_size_m": mesh_size,
            "solver_converged": True,
            "response_order": 2,
            "material_update_order": 1,
            "physical_relative_permeability_bounds": [1.0, 2100.0],
            "volume_m3": 1.0e-6,
            "average_field_T": [0.0, 0.0, average_z],
            "rms_magnitude_T": rms,
            "magnetic_energy_J": energy,
            "magnetic_coenergy_J": 1.1 * energy,
            "h_dot_b_integral_J": 2.1 * energy,
            "legendre_residual_relative": 0.0,
            "field_identity": {**identity, "unit": "T"},
            "energy_identity": {**identity, "unit": "J"},
        }

    return {
        "levels": [
            level(8.0e-3, -1.42, 1.66, 1.20e-3),
            level(4.0e-3, -1.48, 1.62, 1.24e-3),
            level(2.0e-3, -1.50, 1.61, 1.25e-3),
        ]
    }


def test_nonlinear_refinement_energy_gate_accepts_contracting_matched_ladder():
    result = nonlinear_magnetic_refinement_energy_gate(_nonlinear_refinement_summary())
    assert result["status"] == "ok"
    assert all(result["checks"].values())
    assert result["pair_relative_changes"][1]["combined"] < result["pair_relative_changes"][0]["combined"]


def test_nonlinear_refinement_energy_gate_rejects_nonpositive_material_state():
    bad = _nonlinear_refinement_summary()
    bad["levels"][1]["physical_relative_permeability_bounds"] = [0.0, 2100.0]
    result = nonlinear_magnetic_refinement_energy_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_levels_valid"] is False
    assert result["level_checks"][1]["physical_permeability_positive"] is False


def test_nonlinear_refinement_energy_gate_rejects_noncontracting_ladder():
    bad = _nonlinear_refinement_summary()
    bad["levels"][2]["average_field_T"] = [0.0, 0.0, -1.25]
    result = nonlinear_magnetic_refinement_energy_gate(bad)
    assert result["checks"]["field_rms_energy_changes_contract"] is False


def test_nonlinear_refinement_energy_gate_rejects_field_energy_identity_mismatch_and_wraps_mcp():
    bad = _nonlinear_refinement_summary()
    bad["levels"][2]["energy_identity"]["nonlinear_state_id"] = "stale-state"
    result = nonlinear_magnetic_refinement_energy_gate(bad)
    assert result["checks"]["all_levels_valid"] is False
    assert result["level_checks"][2]["field_energy_identity_matches"] is False
    wrapped = json.loads(mcp_nonlinear_refinement_gate(json.dumps(bad)))
    assert wrapped["status"] == "needs_attention"


def test_nonlinear_refinement_energy_gate_rejects_stale_solution_digest():
    bad = _nonlinear_refinement_summary()
    bad["levels"][1]["energy_identity"]["solution_sha256"] = "5" * 64
    result = nonlinear_magnetic_refinement_energy_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["level_checks"][1]["field_energy_identity_matches"] is False


def test_nonlinear_refinement_energy_gate_rejects_material_domain_mismatch():
    bad = _nonlinear_refinement_summary()
    bad["levels"][1]["energy_identity"]["material_domain"] = "iron"
    result = nonlinear_magnetic_refinement_energy_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_levels_valid"] is False
    assert result["level_checks"][1]["field_energy_identity_matches"] is False


def test_nonlinear_refinement_energy_gate_rejects_stale_parent_lineage():
    bad = _nonlinear_refinement_summary()
    for identity_name in ("field_identity", "energy_identity"):
        bad["levels"][1][identity_name][
            "refinement_parent_identity_sha256"
        ] = "6" * 64
    result = nonlinear_magnetic_refinement_energy_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["refinement_parent_identity_consistent"] is False
    assert result["level_checks"][1]["field_energy_identity_matches"] is True


def test_nonlinear_refinement_energy_gate_rejects_legendre_residual():
    bad = _nonlinear_refinement_summary()
    bad["levels"][2]["legendre_residual_relative"] = 2.0e-4
    result = nonlinear_magnetic_refinement_energy_gate(bad)
    assert result["checks"]["all_levels_valid"] is False
    assert result["level_checks"][2]["legendre_energy_identity_satisfied"] is False


def _nonlinear_parity_summary():
    identity = {
        "geometry_identity_sha256": "1" * 64,
        "material_identity_sha256": "2" * 64,
        "excitation_identity_sha256": "3" * 64,
        "bh_table_sha256": "4" * 64,
        "constitutive_interpolation": "monotone_pchip",
        "constitutive_extrapolation": "vacuum_slope",
        "magnetic_anisotropy": "isotropic",
        "region_labels": ["iron"],
        "coordinate_system": "right-handed Cartesian",
        "unit_system": "SI",
        "analysis_kind": "magnetostatic",
        "case_index": 1,
        "time_semantics": "static",
    }
    return {
        "candidate": {
            "identity": dict(identity),
            "average_field_T": [0.0, 0.0, 1.02],
            "rms_magnitude_T": 1.12,
            "magnetic_energy_J": 2.04,
            "magnetic_coenergy_J": 2.97,
        },
        "reference": {
            "identity": dict(identity),
            "average_field_T": [0.0, 0.0, 1.0],
            "rms_magnitude_T": 1.1,
            "magnetic_energy_J": 2.0,
            "magnetic_coenergy_J": 3.0,
        },
    }


def test_nonlinear_field_energy_parity_gate_accepts_identity_matched_observables():
    summary = _nonlinear_parity_summary()
    result = nonlinear_magnetic_field_energy_parity_gate(summary)
    assert result["status"] == "ok"
    assert result["comparison_performed"] is True
    wrapped = json.loads(mcp_nonlinear_parity_gate(json.dumps(summary)))
    assert wrapped["status"] == "ok"


def test_nonlinear_field_energy_parity_gate_rejects_material_law_before_numbers():
    bad = _nonlinear_parity_summary()
    bad["reference"]["identity"]["constitutive_interpolation"] = "piecewise_linear"
    result = nonlinear_magnetic_field_energy_parity_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["physical_identity_matches"] is False
    assert result["identity_checks"]["constitutive_interpolation"] is False
    assert result["comparison_performed"] is False
    assert all(value is None for value in result["relative_differences"].values())


def _constitutive_response_summary():
    h_values = [0.0, 100.0, 500.0, 1000.0, 5000.0, 10000.0]
    b_values = [0.0, 0.4, 1.0, 1.3, 1.6, 1.7]
    coenergy = [0.0, 20.0, 300.0, 875.0, 6675.0, 14925.0]
    energy = [h * b - ws for h, b, ws in zip(h_values, b_values, coenergy)]
    grid_sha = "ab5a4869477578a2f3d7244d52af91efb58ec84f39d5ab1454f307452bed9700"
    identity = {
        "bh_table_sha256": "4" * 64,
        "response_grid_sha256": grid_sha,
        "material_model": "single_valued_isotropic_soft_magnetic",
        "magnetic_anisotropy": "isotropic",
        "H_unit": "A/m",
        "B_unit": "T",
        "differential_permeability_unit": "H/m",
        "energy_density_unit": "J/m^3",
    }
    lane = {
        "identity": identity,
        "H_A_per_m": h_values,
        "B_T": b_values,
        "differential_permeability_H_per_m": [0.004, 0.003, 0.001, 1.0e-4, 2.0e-5, 1.256637e-6],
        "energy_density_J_per_m3": energy,
        "coenergy_density_J_per_m3": coenergy,
        "d_energy_d_B_A_per_m": h_values,
        "d_coenergy_d_H_T": b_values,
    }
    return {"candidate": json.loads(json.dumps(lane)), "reference": lane}


def test_constitutive_response_gate_accepts_realized_common_grid_and_wraps_mcp():
    result = nonlinear_constitutive_response_parity_gate(_constitutive_response_summary())
    assert result["status"] == "ok"
    assert result["comparison_performed"] is True
    wrapped = json.loads(mcp_constitutive_parity_gate(json.dumps(_constitutive_response_summary())))
    assert wrapped["status"] == "ok"


def test_constitutive_response_gate_rejects_same_table_digest_but_different_response():
    bad = _constitutive_response_summary()
    bad["candidate"]["B_T"][3] *= 1.02
    result = nonlinear_constitutive_response_parity_gate(bad)
    assert result["status"] == "needs_attention"
    assert result["checks"]["B_T_matches"] is False


def test_constitutive_response_gate_rejects_saturation_transition_mismatch():
    bad = _constitutive_response_summary()
    bad["candidate"]["differential_permeability_H_per_m"][-1] = 1.0e-4
    result = nonlinear_constitutive_response_parity_gate(bad)
    assert result["checks"]["differential_permeability_H_per_m_matches"] is False


def test_constitutive_response_gate_rejects_non_si_grid_units_before_comparison():
    bad = _constitutive_response_summary()
    bad["candidate"]["identity"]["H_unit"] = "Oe"
    bad["reference"]["identity"]["H_unit"] = "Oe"
    result = nonlinear_constitutive_response_parity_gate(bad)
    assert result["comparison_performed"] is False
    assert result["lane_checks"]["candidate"]["SI_units_explicit"] is False


def test_constitutive_response_gate_rejects_broken_energy_derivative_identity():
    bad = _constitutive_response_summary()
    bad["candidate"]["d_energy_d_B_A_per_m"][2] += 20.0
    result = nonlinear_constitutive_response_parity_gate(bad)
    assert result["checks"]["candidate_response_valid"] is False
    assert result["lane_checks"]["candidate"]["energy_derivative_identity_satisfied"] is False


async def _probe_constitutive_gate_stdio():
    repo = Path(__file__).resolve().parents[3]
    environment = os.environ.copy()
    environment["PYTHONPATH"] = os.pathsep.join(
        [
            str(repo / "packages" / "radia-mcp" / "src"),
            str(repo / "src"),
            environment.get("PYTHONPATH", ""),
        ]
    ).rstrip(os.pathsep)
    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "radia_mcp.radia_ngsolve.server"],
        cwd=str(repo),
        env=environment,
    )
    async with stdio_client(params) as (read_stream, write_stream):
        async with ClientSession(read_stream, write_stream) as session:
            initialized = await session.initialize()
            tools = {tool.name for tool in (await session.list_tools()).tools}
            called = await session.call_tool(
                "radia_ngsolve_validation_run",
                {
                    "name": "nonlinear_constitutive_response_parity_gate",
                    "arguments": {
                        "summary_json": json.dumps(_constitutive_response_summary())
                    },
                },
            )
            return initialized.serverInfo.name, tools, json.loads(called.content[0].text)


def test_constitutive_response_gate_passes_real_stdio_protocol():
    server_name, tools, result = asyncio.run(
        asyncio.wait_for(_probe_constitutive_gate_stdio(), timeout=45)
    )
    assert server_name == "mcp-server-radia-ngsolve"
    assert "radia_ngsolve_validation_run" in tools
    assert result["status"] == "ok"
