"""LTspice discovery stays usable without deploying radia-gpu to the client."""
import os
from pathlib import Path
import subprocess
import sys


def test_ltspice_cold_discovery_and_missing_backend_error():
    source = Path(__file__).resolve().parents[1] / "src"
    code = """
import sys
from radia_mcp.ltspice import server
from radia_mcp.meta.catalog import get
assert len(server.mcp._tool_manager._tools) == 32
assert get('mcp-server-radia-ltspice')['subpackage'] == 'radia_mcp.ltspice'
assert not any(n == 'radia_gpu' or n.startswith('radia_gpu.') for n in sys.modules)
try:
    server.netlist_to_asc('R1 a 0 1k\\n.end\\n')
except RuntimeError as exc:
    assert 'independent radia-gpu' in str(exc)
else:
    raise AssertionError('Absent backend must fail explicitly')
"""
    env = dict(os.environ, PYTHONPATH=str(source))
    subprocess.run([sys.executable, "-c", code], env=env, check=True,
                   capture_output=True, text=True)
