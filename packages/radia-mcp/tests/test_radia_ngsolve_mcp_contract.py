import asyncio
import json
import os
import sys
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from cae_mcp_core.common.mcp_contract import SCHEMA
from radia_mcp.radia_ngsolve.server import mcp


def test_mixed_omega_knowledge_preserves_candidate_and_evidence_boundaries():
    from radia_mcp.radia_ngsolve.knowledge.kelvin import KELVIN_SOURCE_IN_OMEGA_FORM

    text = KELVIN_SOURCE_IN_OMEGA_FORM
    assert "development solver" in text
    assert "not a capability enabled by upgrading radia-mcp alone" in text
    assert "material_update_order" in text
    assert "v5-v11" in text
    assert "unsmoothed element-local" in text
    assert "Historical C-type" in text
    assert "FFAG" in text and "remains pending" in text


def test_radia_ngsolve_server_has_mathworks_style_runtime_contract():
    low_level = mcp._mcp_server
    assert low_level.version
    assert low_level.instructions

    tools = mcp._tool_manager._tools
    assert tools
    for name, tool in tools.items():
        assert tool.title, name
        assert tool.annotations is not None, name
        assert tool.annotations.readOnlyHint is not None, name
        assert tool.annotations.destructiveHint is not None, name
        assert tool.annotations.idempotentHint is not None, name
        assert tool.annotations.openWorldHint is not None, name
        assert tool.meta["caeai.contract"] == SCHEMA


def test_hdiv_vim_knowledge_tool_is_read_only_and_closed_world():
    tool = mcp._tool_manager._tools["hdiv_vim"]

    assert tool.annotations.readOnlyHint is True
    assert tool.annotations.destructiveHint is False
    assert tool.annotations.idempotentHint is True
    assert tool.annotations.openWorldHint is False


async def _probe_hdiv_vim_stdio() -> dict[str, object]:
    package_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(package_root / "src"), env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "radia_mcp.radia_ngsolve.server"],
        cwd=str(package_root),
        env=env,
    )
    async with stdio_client(params) as (read_stream, write_stream):
        async with ClientSession(read_stream, write_stream) as session:
            initialized = await session.initialize()
            listed = await session.list_tools()
            called = await session.call_tool("hdiv_vim", {"topic": "eddy_bubble"})
            identity = {
                "observable_id": "iron_volume_B",
                "field_unit": "T",
                "coordinate_system": "right-handed Cartesian",
            }
            summary = {
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
            catalog = await session.call_tool(
                "radia_ngsolve_validation_catalog",
                {"query": "nonlinear_magnetic_spatial_evidence_gate"},
            )
            catalog_payload = catalog.structuredContent or json.loads(catalog.content[0].text)
            nonlinear_gate = await session.call_tool(
                "radia_ngsolve_validation_run",
                {
                    "name": "nonlinear_magnetic_spatial_evidence_gate",
                    "arguments": {"summary_json": json.dumps(summary)},
                },
            )
            nonlinear_payload = json.loads(nonlinear_gate.content[0].text)
            refinement_catalog = await session.call_tool(
                "radia_ngsolve_validation_catalog",
                {"query": "nonlinear_magnetic_refinement_energy_gate"},
            )
            refinement_catalog_payload = (
                refinement_catalog.structuredContent
                or json.loads(refinement_catalog.content[0].text)
            )
            refinement_identity = {
                "material_domain": "iron",
                "coordinate_system": "right-handed Cartesian",
                "unit_system": "SI",
                "nonlinear_state_id": "state-1",
                "mesh_topology_geometry_sha256": "1" * 64,
                "bh_table_sha256": "2" * 64,
                "material_state_identity_sha256": "3" * 64,
                "solution_sha256": "4" * 64,
                "refinement_parent_identity_sha256": "5" * 64,
            }
            refinement_levels = []
            for mesh_size, average_z, rms in (
                (8.0e-3, -1.42, 1.66),
                (4.0e-3, -1.48, 1.62),
                (2.0e-3, -1.20, 1.61),
            ):
                refinement_levels.append(
                    {
                        "mesh_size_m": mesh_size,
                        "solver_converged": True,
                        "response_order": 2,
                        "material_update_order": 1,
                        "physical_relative_permeability_bounds": [1.0, 2100.0],
                        "volume_m3": 1.0e-6,
                        "average_field_T": [0.0, 0.0, average_z],
                        "rms_magnitude_T": rms,
                        "magnetic_energy_J": 1.2e-3,
                        "magnetic_coenergy_J": 1.3e-3,
                        "h_dot_b_integral_J": 2.5e-3,
                        "legendre_residual_relative": 0.0,
                        "field_identity": {**refinement_identity, "unit": "T"},
                        "energy_identity": {**refinement_identity, "unit": "J"},
                    }
                )
            refinement_gate = await session.call_tool(
                "radia_ngsolve_validation_run",
                {
                    "name": "nonlinear_magnetic_refinement_energy_gate",
                    "arguments": {
                        "summary_json": json.dumps({"levels": refinement_levels})
                    },
                },
            )
            refinement_payload = json.loads(refinement_gate.content[0].text)
            parity_identity = {
                "geometry_identity_sha256": "6" * 64,
                "material_identity_sha256": "7" * 64,
                "excitation_identity_sha256": "8" * 64,
                "bh_table_sha256": "9" * 64,
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
            parity_summary = {
                "candidate": {
                    "identity": parity_identity,
                    "average_field_T": [0.0, 0.0, 1.01],
                    "rms_magnitude_T": 1.11,
                    "magnetic_energy_J": 2.01,
                    "magnetic_coenergy_J": 2.99,
                },
                "reference": {
                    "identity": dict(parity_identity),
                    "average_field_T": [0.0, 0.0, 1.0],
                    "rms_magnitude_T": 1.1,
                    "magnetic_energy_J": 2.0,
                    "magnetic_coenergy_J": 3.0,
                },
            }
            parity_catalog = await session.call_tool(
                "radia_ngsolve_validation_catalog",
                {"query": "nonlinear_magnetic_field_energy_parity_gate"},
            )
            parity_catalog_payload = (
                parity_catalog.structuredContent
                or json.loads(parity_catalog.content[0].text)
            )
            parity_gate = await session.call_tool(
                "radia_ngsolve_validation_run",
                {
                    "name": "nonlinear_magnetic_field_energy_parity_gate",
                    "arguments": {"summary_json": json.dumps(parity_summary)},
                },
            )
            parity_payload = json.loads(parity_gate.content[0].text)
            text = called.content[0].text
            normalized = " ".join(text.split())
            return {
                "server_name": initialized.serverInfo.name,
                "protocol_version": initialized.protocolVersion,
                "tool_count": len(listed.tools),
                "listed": any(tool.name == "hdiv_vim" for tool in listed.tools),
                "is_error": bool(called.isError),
                "nonlinear_gate_discovered": any(
                    operation["name"] == "nonlinear_magnetic_spatial_evidence_gate"
                    for operation in catalog_payload["operations"]
                ),
                "nonlinear_gate_is_error": bool(nonlinear_gate.isError),
                "nonlinear_gate_status": nonlinear_payload["status"],
                "refinement_gate_discovered": any(
                    operation["name"] == "nonlinear_magnetic_refinement_energy_gate"
                    for operation in refinement_catalog_payload["operations"]
                ),
                "refinement_gate_is_error": bool(refinement_gate.isError),
                "refinement_gate_status": refinement_payload["status"],
                "refinement_gate_issues": refinement_payload["issues"],
                "parity_gate_discovered": any(
                    operation["name"]
                    == "nonlinear_magnetic_field_energy_parity_gate"
                    for operation in parity_catalog_payload["operations"]
                ),
                "parity_gate_is_error": bool(parity_gate.isError),
                "parity_gate_status": parity_payload["status"],
                "teaches_direct_q2": "direct-Q2" in text,
                "bounds_h_convergence": (
                    "for this thin magnetic-conductor disk lane" in normalized
                    and "not for every geometry" in normalized
                ),
                "teaches_mapped_bdm2_gate": (
                    "mapped/non-affine pure-HEX BDM2 primal material and field lane"
                    in normalized
                    and "production C++ composite operator" in normalized
                    and "complete-host tensor rules" in normalized
                    and "same 207 active mapped-body BDM2 DoFs" in normalized
                    and "reflection-invariant whole-host Duffy rules" in normalized
                    and "not an open-boundary accuracy oracle" in normalized
                ),
                "teaches_h1_hodge_mixed_gate": (
                    "H1HodgeDemagOperator" in normalized
                    and "hdiv.FreeDofs()" in normalized
                    and "snapshot residual below `3e-12`" in normalized
                    and "independently of the repaired open-boundary operator"
                    in normalized
                    and "not an accuracy oracle" in normalized
                ),
                "teaches_h1_hodge_accuracy_ladder": (
                    "1.0916977441" in normalized
                    and "4.65%" in normalized
                    and "0.98%" in normalized
                    and "strict h/p ladder" in normalized
                    and "not a universal open-boundary" in normalized
                ),
                "teaches_coupled_local_esim_boundary": (
                    "solve_frequency_local_esim" in normalized
                    and "50/1000/5000 A/m" in normalized
                    and "fixed-Gram replay to machine precision" in normalized
                    and "bulk nonlinear B-H operator" in normalized
                    and "is not yet implemented" in normalized
                ),
                "teaches_local_hdiv_response_gate": (
                    "HDiv local-response completeness gate" in normalized
                    and "NgsolveHDivLocalPolynomialTrainingPorts" in normalized
                    and "30 normalized vector-polynomial probes" in normalized
                    and "rejects them as physical `external_fields`" in normalized
                    and "not a universal accuracy claim" in normalized
                ),
            }


def test_hdiv_vim_passes_real_stdio_initialize_list_call():
    result = asyncio.run(asyncio.wait_for(_probe_hdiv_vim_stdio(), timeout=45))

    assert result["server_name"] == "mcp-server-radia-ngsolve"
    assert result["protocol_version"]
    assert 40 <= result["tool_count"] < 100
    assert result["listed"] is True
    assert result["is_error"] is False
    assert result["nonlinear_gate_discovered"] is True
    assert result["nonlinear_gate_is_error"] is False
    assert result["nonlinear_gate_status"] == "ok"
    assert result["refinement_gate_discovered"] is True
    assert result["refinement_gate_is_error"] is False
    assert result["refinement_gate_status"] == "needs_attention"
    assert "field_rms_energy_changes_contract" in result["refinement_gate_issues"]
    assert result["parity_gate_discovered"] is True
    assert result["parity_gate_is_error"] is False
    assert result["parity_gate_status"] == "ok"
    assert result["teaches_direct_q2"] is True
    assert result["bounds_h_convergence"] is True
    assert result["teaches_mapped_bdm2_gate"] is True
    assert result["teaches_h1_hodge_mixed_gate"] is True
    assert result["teaches_h1_hodge_accuracy_ladder"] is True
    assert result["teaches_coupled_local_esim_boundary"] is True
    assert result["teaches_local_hdiv_response_gate"] is True
