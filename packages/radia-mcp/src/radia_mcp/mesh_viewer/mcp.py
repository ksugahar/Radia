"""Asynchronous conversion transaction, followed by detached human viewing."""
from __future__ import annotations
import asyncio
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from .launcher import require_interactive_desktop, resolve_mesh_path
from ._worker_process import run_worker, transaction_deadline
from .post_view import prepare_post_view_async

MODULE = 'radia_mcp.mesh_viewer'
SOURCE_ROOT = Path(__file__).resolve().parents[2]

async def open_mesh_viewer(filepath: str, post_output: str | None = None,
                           checked_report: str | None = None,
                           source_sha256: str | None = None,
                           group_scope: str | None = None) -> dict:
    """Await verified conversion before launch; launch is not rendering proof."""
    deadline = transaction_deadline()
    path = await asyncio.to_thread(resolve_mesh_path, filepath)
    desktop = require_interactive_desktop()
    python = Path(os.environ.get('RADIA_MESH_VIEWER_PYTHON', sys.executable)).resolve(strict=True)
    if not python.is_file(): raise ValueError('Viewer interpreter must be an existing executable')
    env = {**os.environ, 'PYTHONDONTWRITEBYTECODE': '1', 'PYTHONNOUSERSITE': '1',
           'PYTHONIOENCODING': 'utf-8', 'PYTHONUTF8': '1', 'PYTHONPATH': str(SOURCE_ROOT)}
    options = (checked_report, source_sha256, group_scope)
    if post_output is None and any(value is not None for value in options):
        raise ValueError('Post-view options require post_output')
    conversion = None
    if post_output is not None:
        if any(value is None for value in options):
            raise ValueError('Post-view requires checked_report, source_sha256 and group_scope')
        conversion = await prepare_post_view_async(path, post_output, checked_report,
            source_sha256, group_scope, python=python, env=env, deadline=deadline)
        path = await asyncio.to_thread(resolve_mesh_path, conversion['output'])
    argv = [str(python), '-B', '-m', MODULE, str(path)]
    if conversion: argv.extend(['--expected-sha256', conversion['sha256']])
    probe = await run_worker([*argv, '--check'], env=env, deadline=deadline,
                             phase=f'Viewer preflight {path}')
    identity = json.loads(probe.stdout.decode('utf-8'))
    if time.monotonic() >= deadline:
        raise RuntimeError(f'GUI launch {path}: transaction deadline exhausted; verified artifact retained')
    try:
        child = await asyncio.to_thread(subprocess.Popen, argv, env=env,
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            close_fds=True, creationflags=subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP)
    except Exception as error:
        retained = f'; verified artifact retained: {conversion}' if conversion else ''
        raise RuntimeError(f'GUI launch failed: {error}{retained}') from error
    result = {'status': 'launched', 'pid': child.pid, 'argv': argv,
              'runtime_identity': identity, 'desktop': desktop}
    if conversion: result['post_view'] = conversion
    return result

def register_mesh_viewer_tool(mcp) -> None:
    from mcp.types import ToolAnnotations
    mcp.tool(annotations=ToolAnnotations(readOnlyHint=False, destructiveHint=False,
                                       idempotentHint=False, openWorldHint=False))(open_mesh_viewer)
