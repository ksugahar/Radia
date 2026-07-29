from __future__ import annotations

import asyncio
import json
import os
import sys
from pathlib import Path

import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

pytest.importorskip("ngsolve")

from radia_mcp.radia_ngsolve.vol2d_circuit import write_structured_rect_vol


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
PACKAGE_ROOT = REPOSITORY_ROOT / "packages" / "radia-mcp"


def _transient_request(tmp_path: Path) -> dict:
    path = tmp_path / "mcp_transient_heat.vol"
    write_structured_rect_vol(
        path,
        x0=0.0,
        x1=1.0,
        y0=0.0,
        y1=1.0,
        nx=2,
        ny=2,
        material="domain",
    )
    return {
        "operation": "transient_heat",
        "physics": "transient_heat",
        "vol_text": path.read_text(encoding="utf-8"),
        "source_name": "generated_mcp_transient_heat.vol",
        "element_family": "P1",
        "formulation": "planar",
        "model_depth_m": 0.25,
        "dirichlet_values": {},
        "robin_boundaries": {},
        "materials": {
            "domain": {
                "coefficient_si": 2.0,
                "volumetric_source_si": 8.0,
                "volumetric_heat_capacity_j_per_m3_k": 4.0,
            }
        },
        "initial_temperature_k": 300.0,
        "time_s": [0.0, 0.1],
        "theta": 1.0,
        "export_basename": "mcp_transient_heat",
    }


async def _probe_stdio(request: dict) -> dict:
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [
            str(PACKAGE_ROOT / "src"),
            str(REPOSITORY_ROOT / "src"),
            env.get("PYTHONPATH", ""),
        ]
    ).rstrip(os.pathsep)
    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "radia_mcp.fem.server"],
        cwd=str(PACKAGE_ROOT),
        env=env,
    )
    async with stdio_client(params) as (read_stream, write_stream):
        async with ClientSession(read_stream, write_stream) as session:
            initialized = await session.initialize()
            listed = await session.list_tools()
            called = await session.call_tool(
                "fem_vol2d_scalar_analysis",
                {"analysis_json": json.dumps(request)},
            )
            payload = json.loads(called.content[0].text)
            return {
                "server_name": initialized.serverInfo.name,
                "listed": any(
                    tool.name == "fem_vol2d_scalar_analysis" for tool in listed.tools
                ),
                "is_error": bool(called.isError),
                "payload": payload,
            }


def test_transient_heat_passes_real_mcp_stdio_and_owned_worker(tmp_path: Path) -> None:
    result = asyncio.run(
        asyncio.wait_for(_probe_stdio(_transient_request(tmp_path)), timeout=45)
    )

    assert result["server_name"] == "mcp-server-fem"
    assert result["listed"] is True
    assert result["is_error"] is False
    payload = result["payload"]
    assert payload["status"] == "solved"
    contract = payload["result_contract"]
    assert contract["maximum_temperature_history_k"][-1] == pytest.approx(300.2)
    assert contract["maximum_step_residual_inf"] < 1.0e-10
