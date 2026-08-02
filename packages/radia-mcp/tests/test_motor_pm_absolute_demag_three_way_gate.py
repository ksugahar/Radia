from __future__ import annotations

import asyncio
import copy
import json
import os
import sys
from pathlib import Path

import pytest

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from radia_mcp.motor.pm_armature_reaction_gate import (
    pm_absolute_demag_three_way_gate,
)
from radia_mcp.motor.server import motor_pm_absolute_demag_three_way_gate


def _fixture() -> dict[str, object]:
    identity = "b" * 64
    return {
        "schema": "radia.pm-ring-absolute-demag-attribution-evidence.v1",
        "physics_identity_sha256": identity,
        "reference": {"executed": True, "physics_identity_sha256": identity},
        "hdiv": {
            "cell_family": "HEX",
            "project_lanes": ["BDM1", "BDM2"],
            "bdm1_source_nrms": [0.1802, 0.1104, 0.07244],
            "bdm2_source_nrms": [0.1681, 0.0949, 0.05926],
            "bdm2_final_vs_fem_nrms": 0.04004,
            "bdm2_final_relative_step": 4.1e-7,
        },
        "independent_fem": {
            "formulation": "continuous_H1_scalar_potential",
            "cell_family": "TET",
            "equal_sector_volumes": True,
            "source_nrms_outer_p2": [0.02232, 0.02042, 0.02026],
            "source_nrms_p3": 0.01957,
            "p2_p3_same_mesh_nrms": 0.00270,
        },
        "research_lab_retirement_ready": False,
        "product_or_market_retirement_ready": False,
    }


def test_attributes_remaining_gap_without_validating_absolute_demag() -> None:
    result = pm_absolute_demag_three_way_gate(_fixture())
    assert result["status"] == "attributed_pending"
    assert result["independent_fem_corroborates_reference"] is True
    assert result["absolute_self_demagnetizing_field_validated"] is False
    assert result["attribution"] == "hdiv_resolution_or_cross_body_operator_gap"
    assert result["metrics"]["final_bdm1_source_nrms"] == 0.07244


def test_validates_only_when_bdm2_and_independent_fem_close() -> None:
    artifact = _fixture()
    artifact["hdiv"]["bdm2_source_nrms"] = [0.08, 0.04, 0.02]
    artifact["hdiv"]["bdm2_final_vs_fem_nrms"] = 0.01
    result = pm_absolute_demag_three_way_gate(artifact)
    assert result["status"] == "validated"
    assert result["absolute_self_demagnetizing_field_validated"] is True


def test_rejects_unstable_fem_or_retirement_overclaim() -> None:
    artifact = _fixture()
    artifact["independent_fem"]["p2_p3_same_mesh_nrms"] = 0.02
    artifact["research_lab_retirement_ready"] = True
    result = pm_absolute_demag_three_way_gate(artifact)
    assert result["status"] == "needs_attention"
    assert result["checks"]["fem_order_is_stable"] is False
    assert result["checks"]["scope_does_not_overclaim_retirement"] is False


def test_mcp_tool_returns_attribution_and_handles_invalid_input() -> None:
    response = json.loads(motor_pm_absolute_demag_three_way_gate(json.dumps(_fixture())))
    assert response["status"] == "attributed_pending"

    broken = copy.deepcopy(_fixture())
    broken["hdiv"]["bdm2_source_nrms"] = [0.1]
    response = json.loads(motor_pm_absolute_demag_three_way_gate(json.dumps(broken)))
    assert response["status"] == "invalid_input"


async def _probe_stdio() -> dict[str, object]:
    package_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([
        str(package_root / "src"), env.get("PYTHONPATH", "")
    ]).rstrip(os.pathsep)
    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "radia_mcp.motor.server"],
        cwd=str(package_root),
        env=env,
    )
    async with stdio_client(params) as (read_stream, write_stream):
        async with ClientSession(read_stream, write_stream) as session:
            initialized = await session.initialize()
            listed = await session.list_tools()
            called = await session.call_tool(
                "motor_pm_absolute_demag_three_way_gate",
                {"summary_json": json.dumps(_fixture())},
            )
            return {
                "server_name": initialized.serverInfo.name,
                "listed": any(
                    tool.name == "motor_pm_absolute_demag_three_way_gate"
                    for tool in listed.tools
                ),
                "is_error": bool(called.isError),
                "status": json.loads(called.content[0].text)["status"],
            }


def test_changed_tool_passes_real_stdio_protocol() -> None:
    result = asyncio.run(asyncio.wait_for(_probe_stdio(), timeout=45))
    assert result == {
        "server_name": "mcp-server-motor",
        "listed": True,
        "is_error": False,
        "status": "attributed_pending",
    }


@pytest.mark.parametrize(
    ("section", "field"),
    [
        ("hdiv", "bdm2_final_vs_fem_nrms"),
        ("hdiv", "bdm2_final_relative_step"),
        ("independent_fem", "source_nrms_p3"),
    ],
)
def test_rejects_negative_error_or_iteration_metrics(section, field) -> None:
    artifact = _fixture()
    artifact[section][field] = -1.0e-6

    with pytest.raises(ValueError, match="nonnegative"):
        pm_absolute_demag_three_way_gate(artifact)
