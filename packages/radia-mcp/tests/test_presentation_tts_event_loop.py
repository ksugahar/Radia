"""The TTS embedding tool synthesises audio with asyncio; called from the MCP
server it runs inside an event loop, where asyncio.run() raises.  The helper
must work both with and without a running loop."""
import asyncio

from radia_mcp.presentation.plans.T31 import _run_blocking


async def _double(x):
    await asyncio.sleep(0)
    return 2 * x


def test_run_blocking_without_loop():
    assert _run_blocking(_double(21)) == 42


def test_run_blocking_inside_running_loop():
    async def handler():
        # this is what the MCP tool handler looks like: sync code awaited by
        # the server loop, which then needs the coroutine run to completion
        return _run_blocking(_double(5))

    assert asyncio.run(handler()) == 10
