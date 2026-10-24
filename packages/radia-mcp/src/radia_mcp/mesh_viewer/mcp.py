"""Thin asynchronous adapter: the server never owns a GUI event loop."""
from __future__ import annotations

import asyncio
import json
import os
from pathlib import Path
import subprocess
import sys

from .launcher import require_interactive_desktop, validate_mesh_path

MODULE = "radia_mcp.mesh_viewer"
SOURCE_ROOT = Path(__file__).resolve().parents[2]


async def open_mesh_viewer(filepath: str) -> dict:
    """Open an existing .vol for the current interactive Windows user.

    Returns launch acceptance, not proof that rendering succeeded. No file
    writes, solvers, associations, shell, DDE, COM or runtime installation.
    RADIA_MESH_VIEWER_PYTHON is an operator-selected interpreter, not a tool arg.
    """
    path = await asyncio.to_thread(validate_mesh_path, filepath)
    desktop = require_interactive_desktop()
    python = Path(os.environ.get("RADIA_MESH_VIEWER_PYTHON", sys.executable)).resolve(strict=True)
    if not python.is_file():
        raise ValueError("Viewer interpreter must be an existing executable")
    # Use the SAME reviewed package as the server with the selected dependency
    # interpreter. This changes only the child's environment, never installs.
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "PYTHONNOUSERSITE": "1",
           "PYTHONPATH": str(SOURCE_ROOT)}
    argv = [str(python), "-B", "-m", MODULE, str(path)]
    probe = await asyncio.create_subprocess_exec(
        *argv, "--check", env=env, stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE, creationflags=subprocess.CREATE_NO_WINDOW)
    try:
        stdout, stderr = await asyncio.wait_for(probe.communicate(), timeout=30)
    except (TimeoutError, asyncio.CancelledError):
        probe.kill()
        await probe.wait()
        raise
    if probe.returncode:
        raise RuntimeError("Viewer preflight failed: " + stderr.decode("utf-8", errors="replace").strip())
    identity = json.loads(stdout.decode("utf-8"))
    child = await asyncio.to_thread(subprocess.Popen,
        argv, env=env, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL, close_fds=True,
        creationflags=subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP)
    return {"status": "launched", "pid": child.pid, "argv": argv,
            "runtime_identity": identity, "desktop": desktop}


def register_mesh_viewer_tool(mcp) -> None:
    """Register only this adapter on an existing MCP server."""
    from mcp.types import ToolAnnotations
    mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False,
                                       idempotentHint=False, openWorldHint=False))(open_mesh_viewer)
