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


async def open_mesh_viewer(filepath: str, post_output: str | None = None,
                           checked_report: str | None = None,
                           source_sha256: str | None = None,
                           group_scope: str | None = None) -> dict:
    """Open VOL/strict MSH4.2 for the current interactive Windows user.

    Returns launch acceptance, not proof that rendering succeeded. No file
    writes unless an explicit VOL post_output is supplied. No solvers,
    associations, shell, DDE, COM or runtime installation.
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
    options = {'post-output': post_output, 'checked-report': checked_report,
               'source-sha256': source_sha256, 'group-scope': group_scope}
    if post_output is None and any(value is not None for key, value in options.items() if key != 'post-output'):
        raise ValueError('Post-view options require post_output')
    if post_output is not None:
        if any(value is None for value in options.values()):
            raise ValueError('Post-view requires checked_report, source_sha256 and group_scope')
        for key, value in options.items():
            argv.extend(['--' + key, value])
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
