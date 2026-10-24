"""Guards: catalogued LTspice tools remain discoverable without their backend.

Assumes: a fresh child using this package source and an explicitly unavailable
radia-simulink backend; discovery and fail-loud behavior, not circuit accuracy.
"""
import os
from pathlib import Path
import subprocess
import sys


def test_ltspice_cold_discovery_and_missing_backend_error():
    source = Path(__file__).resolve().parents[1] / "src"
    code = """
import asyncio
import importlib.abc
import sys
from pathlib import Path

class AbsentBackend(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'radia_simulink' or fullname.startswith('radia_simulink.'):
            raise ModuleNotFoundError('Backend unavailable in this contract', name=fullname)

sys.meta_path.insert(0, AbsentBackend())
from radia_mcp.ltspice import server
from radia_mcp.meta.catalog import get
assert Path(server.__file__).resolve() == Path(sys.argv[1]) / 'radia_mcp/ltspice/server.py'
catalog = get('mcp-server-radia-ltspice')
tools = asyncio.run(server.mcp.list_tools())
names = {tool.name for tool in tools}
required = set(catalog['primary_tools'])
assert required and required <= names, {'missing': sorted(required - names), 'live': sorted(names)}
assert len(names) == len(tools), 'Discovery must not duplicate tool names'
assert catalog['subpackage'] == 'radia_mcp.ltspice'
sys.argv.append('--selftest')
assert server.main() == 0
assert not any(n == 'radia_simulink' or n.startswith('radia_simulink.') for n in sys.modules)
try:
    server.netlist_to_asc('R1 a 0 1k\\n.end\\n')
except RuntimeError as exc:
    assert 'independent radia-simulink' in str(exc)
else:
    raise AssertionError('Absent backend must fail explicitly')
"""
    env = dict(os.environ, PYTHONPATH=str(source), PYTHONUTF8="1", PYTHONIOENCODING="utf-8")
    result = subprocess.run([sys.executable, "-c", code, str(source.resolve())], env=env,
                            capture_output=True, encoding="utf-8", errors="replace", timeout=60)
    assert result.returncode == 0, (result.stdout, result.stderr)
