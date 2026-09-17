"""Grouped validation calls must not block the MCP event loop."""

import asyncio
import threading

from mcp.server.fastmcp import FastMCP

from cubit_mesh_export.mcp._support.tool_group import CoarseToolRegistry


def test_sync_grouped_operation_runs_off_event_loop():
    mcp = FastMCP("grouped-offload-test")
    registry = CoarseToolRegistry(mcp, namespace="cubit", profile="core",
                                  min_group_size=1)

    @registry.tool()
    def identify_thread() -> int:
        """Return the executing thread identity."""
        return threading.get_ident()

    registry.install()
    runner = mcp._tool_manager._tools["cubit_validation_run"].fn

    async def check():
        event_loop_thread = threading.get_ident()
        worker_thread = await runner("identify_thread")
        return event_loop_thread, worker_thread

    event_loop_thread, worker_thread = asyncio.run(check())
    assert worker_thread != event_loop_thread
