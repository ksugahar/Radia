# Cubit-owned implementation; intentionally maintained independently of Radia MCP.
# SPDX-License-Identifier: Apache-2.0
"""Shared subprocess runner for gmsh-API-backed checks and rendering.

gmsh keeps process-global state and can hard-crash on corrupt input, so
every gmsh API call made on behalf of an MCP tool runs in a child
process.  Results come back through a temp JSON file with stdout/stderr
redirected to a log file and ``stdin=DEVNULL`` (the stdio-MCP pipe-hang
pattern: a lingering grandchild must never hold the server's pipes).
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from importlib import util as _importlib_util
from pathlib import Path
from typing import Any


def gmsh_available() -> bool:
    """True when the gmsh Python package is importable in this env."""
    return _importlib_util.find_spec("gmsh") is not None


def run_gmsh_json_subprocess(script: str, args: list[str], *,
                             timeout_s: float,
                             prefix: str = "cubit_mcp_gmsh_") -> dict[str, Any]:
    """Run ``python <script.py> <args...> <out_json>`` and read the JSON back.

    The script MUST take the output JSON path as its LAST argv entry and
    write a JSON object there even on failure (wrap its body in
    try/except and record ``error``).  The script is written to a temp
    file rather than passed via ``-c`` -- long scripts overflow the
    Windows 32k command-line limit (WinError 206).  Returns
    ``{"ok": False, "ran": False, "error": ...}`` when gmsh is missing,
    the child dies without writing a result, or the timeout expires.
    """
    if not gmsh_available():
        return {"ok": False, "ran": False,
                "error": "Gmsh 5 missing; install Gmsh alone with --pre -i https://gmsh.info/python-packages-dev gmsh"}

    with tempfile.TemporaryDirectory(prefix=prefix) as work:
        out_json = Path(work) / "result.json"
        log_path = Path(work) / "run.log"
        script_path = Path(work) / "script.py"
        # Every reference operation must pass the same native/format boundary.
        # Keep native state in this child, including the readiness roundtrip.
        wrapped = f'''import json, sys
from pathlib import Path
sys.path.insert(0, {str(Path(__file__).resolve().parents[3])!r})
from cubit_mesh_export._gmsh_runtime import probe_gmsh5, require_msh42
try:
    runtime = probe_gmsh5()
    for argument in sys.argv[1:-1]:
        candidate = Path(argument)
        if candidate.is_file():
            if candidate.name.lower().endswith('.msh.gz'):
                raise ValueError('Compressed MSH input is unsupported; use MSH 4.2 directly')
            with candidate.open('rb') as source:
                marker = source.read(11).strip()
            if candidate.suffix.lower() == '.msh' or marker == b'$MeshFormat':
                require_msh42(candidate)
    exec(compile({script!r}, '<gmsh-operation>', 'exec'))
    output = Path(sys.argv[-1])
    result = json.loads(output.read_text(encoding='utf-8'))
    result['gmsh_runtime'] = runtime
except Exception as exc:
    result = {{'ok': False, 'ran': False,
              'error': f'{{type(exc).__name__}}: {{exc}}'}}
Path(sys.argv[-1]).write_text(json.dumps(result), encoding='utf-8')
'''
        script_path.write_text(wrapped, encoding="utf-8")
        cmd = [sys.executable, str(script_path), *args, str(out_json)]
        try:
            with open(log_path, "w", encoding="utf-8") as log:
                proc = subprocess.run(
                    cmd, stdin=subprocess.DEVNULL, stdout=log,
                    stderr=subprocess.STDOUT, timeout=timeout_s, check=False)
        except subprocess.TimeoutExpired:
            return {"ok": False, "ran": False,
                    "error": f"gmsh subprocess timed out after {timeout_s}s"}
        if not out_json.is_file():
            tail = log_path.read_text(encoding="utf-8", errors="replace")[-1500:]
            return {"ok": False, "ran": False,
                    "error": f"gmsh subprocess exited {proc.returncode} "
                             f"without result: {tail}"}
        result = json.loads(out_json.read_text(encoding="utf-8"))
        if proc.returncode != 0:
            return {"ok": False, "ran": False,
                    "error": f"Gmsh operation exited {proc.returncode}"}
        return result
