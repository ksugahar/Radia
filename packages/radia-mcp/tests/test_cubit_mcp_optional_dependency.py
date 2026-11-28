"""cubit-mcp works on the base radia-mcp install.

Guards: the Cubit usage tools (knowledge, documentation, lint) load and run
without the ``[cubit]`` extra, so another MCP server can use cubit-mcp
without pulling in radia, the NGSolve fork or numpy.
Assumes: the ``.vol`` export and check-vol keep their heavy imports inside
the functions that need them.
"""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

# Modules that belong to the [cubit] extra or to Cubit itself.
HEAVY = ("radia", "ngsolve", "netgen", "numpy", "scipy", "cubit")


def test_cubit_mcp_runs_without_the_cubit_extra(tmp_path):
    package_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(package_root / "src"), env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    journal = tmp_path / "probe.py"
    journal.write_text('cubit.cmd("brick x 1")\n', encoding="utf-8")
    script = rf"""
import asyncio
import sys

HEAVY = {HEAVY!r}
for name in HEAVY:
    sys.modules[name] = None  # any import of these now raises ImportError

from radia_mcp.cubit import server

tools = {{tool.name for tool in asyncio.run(server.mcp.list_tools())}}
for name in ("cubit_docs", "lint_cubit_script", "cubit_exec", "cubit_export_vol"):
    assert name in tools, (name, sorted(tools))

assert server.cubit_docs("all").strip(), "empty documentation overview"
assert server.lint_cubit_script({str(journal)!r}).strip(), "empty lint report"

loaded = sorted(n for n in HEAVY if sys.modules.get(n) is not None)
assert not loaded, loaded
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=package_root,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
