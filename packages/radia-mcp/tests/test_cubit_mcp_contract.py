import asyncio
import json
import os
import sys
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client


TWO_DOMAIN_VOL = """mesh3d
dimension
3

surfaceelements
2
1 1 1 0 3 1 2 3
2 2 2 0 3 5 6 7

volumeelements
2
1 4 1 2 3 4
2 4 5 6 7 8

points
8
0 0 0
1 0 0
0 1 0
0 0 1
2 0 0
3 0 0
2 1 0
2 0 1

materials
2
1 left
2 right

endmesh
"""

TWO_DOMAIN_INTERFACE_VOL = TWO_DOMAIN_VOL.replace(
    "surfaceelements\n2",
    "surfaceelements\n3",
).replace(
    "\nvolumeelements",
    "\n3 3 1 2 3 2 3 4\n\nvolumeelements",
)


async def _probe_cubit_inventory_stdio() -> dict[str, object]:
    package_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(package_root / "src"), env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "radia_mcp.cubit.server"],
        cwd=str(package_root),
        env=env,
    )
    stale = TWO_DOMAIN_VOL.replace("2 2 2 0 3 5 6 7", "2 2 1 0 3 5 6 7")
    duplicate_interface = TWO_DOMAIN_INTERFACE_VOL.replace(
        "surfaceelements\n3",
        "surfaceelements\n4",
    ).replace(
        "\nvolumeelements",
        "\n4 4 2 1 3 2 3 4\n\nvolumeelements",
    )
    async with stdio_client(params) as (read_stream, write_stream):
        async with ClientSession(read_stream, write_stream) as session:
            initialized = await session.initialize()
            listed = await session.list_tools()
            positive = await session.call_tool(
                "cubit_vol_inventory", {"text": TWO_DOMAIN_VOL}
            )
            negative = await session.call_tool(
                "cubit_vol_inventory", {"text": stale}
            )
            interface = await session.call_tool(
                "cubit_vol_inventory", {"text": TWO_DOMAIN_INTERFACE_VOL}
            )
            duplicate = await session.call_tool(
                "cubit_vol_inventory", {"text": duplicate_interface}
            )
            return {
                "server_name": initialized.serverInfo.name,
                "protocol_version": initialized.protocolVersion,
                "listed": any(
                    tool.name == "cubit_vol_inventory" for tool in listed.tools
                ),
                "positive": json.loads(positive.content[0].text),
                "positive_error": bool(positive.isError),
                "negative": json.loads(negative.content[0].text),
                "negative_error": bool(negative.isError),
                "interface": json.loads(interface.content[0].text),
                "interface_error": bool(interface.isError),
                "duplicate": json.loads(duplicate.content[0].text),
                "duplicate_error": bool(duplicate.isError),
            }


def test_cubit_inventory_passes_real_stdio_initialize_list_and_positive_negative_calls():
    result = asyncio.run(
        asyncio.wait_for(_probe_cubit_inventory_stdio(), timeout=45)
    )

    assert "cubit" in str(result["server_name"]).lower()
    assert result["protocol_version"]
    assert result["listed"] is True
    assert result["positive_error"] is False
    assert result["positive"]["status"] == "ok"
    assert result["positive"]["boundary_domain_ownership"]["passed"] is True
    assert result["negative_error"] is False
    assert result["negative"]["status"] == "needs_attention"
    assert result["negative"]["boundary_domain_ownership"]["passed"] is False
    assert result["negative"]["boundary_domain_ownership"][
        "unreferenced_volume_domains"
    ] == [2]
    assert result["interface_error"] is False
    assert result["interface"]["status"] == "ok"
    assert result["interface"]["boundary_domain_ownership"][
        "internal_interface_domain_pair_counts"
    ] == {"1<->2": 1}
    assert result["duplicate_error"] is False
    assert result["duplicate"]["status"] == "needs_attention"
    assert result["duplicate"]["boundary_domain_ownership"][
        "duplicate_surface_connectivity_rows"
    ] == [4]
