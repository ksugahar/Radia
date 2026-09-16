"""Security and server-responsiveness regression contracts."""
import asyncio
import json
import time

from mcp.server.fastmcp import FastMCP

from cubit_mesh_export.mcp import server


def test_race_history_rejects_path_like_identifiers(monkeypatch, tmp_path):
    monkeypatch.setattr(server, "_race_history_dir", lambda: tmp_path)
    outside = tmp_path.parent / "outside.json"
    outside.write_text(
        json.dumps({"all_results": [{"name": "bad", "is_ai": True,
                                     "recipe": ["exit"]}]}),
        encoding="utf-8",
    )

    status = json.loads(server.cubit_mesh_race_status("../outside"))
    applied = json.loads(server.cubit_mesh_apply_choice("../outside", "bad"))

    assert status["status"] == "error"
    assert applied["status"] == "error"
    assert "opaque identifier" in status["error"]
    assert outside.exists()


def test_registered_tools_are_async_and_do_not_block_event_loop():
    probe = FastMCP("thread-offload-probe")

    @probe.tool()
    def slow_probe() -> str:
        time.sleep(0.08)
        return "done"

    assert server._offload_sync_tool_functions(probe) == 1
    tool = probe._tool_manager._tools["slow_probe"]
    assert tool.is_async is True

    async def exercise():
        started = time.perf_counter()
        slow = asyncio.create_task(tool.run({}))
        await asyncio.sleep(0.01)
        heartbeat_at = time.perf_counter() - started
        result = await slow
        return heartbeat_at, result

    heartbeat_at, result = asyncio.run(exercise())
    assert heartbeat_at < 0.05
    assert result == "done"
    assert server._OFFLOADED_TOOL_COUNT > 0
    assert all(tool.is_async
               for tool in server.mcp._tool_manager._tools.values())
