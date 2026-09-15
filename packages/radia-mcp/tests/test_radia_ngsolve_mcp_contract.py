import asyncio
import os
import sys
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from radia_mcp.common.mcp_contract import SCHEMA
from radia_mcp.radia_ngsolve.server import mcp


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
            text = called.content[0].text
            normalized = " ".join(text.split())
            return {
                "server_name": initialized.serverInfo.name,
                "protocol_version": initialized.protocolVersion,
                "tool_count": len(listed.tools),
                "listed": any(tool.name == "hdiv_vim" for tool in listed.tools),
                "is_error": bool(called.isError),
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
    assert result["teaches_direct_q2"] is True
    assert result["bounds_h_convergence"] is True
    assert result["teaches_mapped_bdm2_gate"] is True
    assert result["teaches_h1_hodge_mixed_gate"] is True
    assert result["teaches_h1_hodge_accuracy_ladder"] is True
    assert result["teaches_coupled_local_esim_boundary"] is True
    assert result["teaches_local_hdiv_response_gate"] is True
