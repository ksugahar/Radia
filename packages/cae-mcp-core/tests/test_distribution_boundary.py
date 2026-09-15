"""Independent foundation acceptance, without either consuming product."""
import ast
import json
from pathlib import Path
import subprocess
import sys

import cae_mcp_core


def test_core_has_no_product_imports():
    root = Path(cae_mcp_core.__file__).parent
    for path in root.rglob('*.py'):
        for node in ast.walk(ast.parse(path.read_text(encoding='utf-8'))):
            imports = ([a.name for a in node.names] if isinstance(node, ast.Import)
                       else [node.module or ''] if isinstance(node, ast.ImportFrom) else [])
            assert not any(n.split('.')[0] in {'radia', 'radia_mcp', 'cubit_mesh_export'}
                           for n in imports), path


def test_core_registers_without_optional_native_or_product_imports():
    code = '''
import sys, json, importlib.abc
class ProductBlocker(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'radia','radia_mcp','cubit_mesh_export','netgen','ngsolve','gmsh','PySide6'}:
            raise AssertionError('Forbidden import: ' + fullname)
sys.meta_path.insert(0, ProductBlocker())
from mcp.server.fastmcp import FastMCP
from cae_mcp_core.common.status import register_status_tool
server = FastMCP('foundation-test')
register_status_tool(server, 'mcp-server-foundation', 'test', 'cae_mcp_core')
payload = server._tool_manager._tools['foundation_status'].fn()
assert payload['runtime_provenance']['distribution']['name'] == 'cae-mcp-core'
assert payload['runtime_contract']['complete']
print(json.dumps({'passed': True}))
'''
    run = subprocess.run([sys.executable, '-I', '-c', code], capture_output=True,
                         text=True, timeout=30)
    assert run.returncode == 0, run.stderr
    assert json.loads(run.stdout)['passed']
