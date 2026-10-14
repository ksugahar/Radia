"""Cubit-only installation owns its API, MCP entry point and runtime dependency."""
from importlib import metadata
import ast
from pathlib import Path

from packaging.requirements import Requirement


def test_exporter_owns_mcp_without_radia_dependency():
    distribution = metadata.distribution("cubit-mesh-export")
    requires = {Requirement(item).name for item in distribution.requires or []}
    assert "mcp" in requires
    assert not {"radia", "radia-mcp", "cae-mcp-core"}.intersection(requires)
    entry = next(e for e in distribution.entry_points if e.name == "mcp-server-cubit")
    assert entry.value == "cubit_mesh_export.mcp.server:main"
    assert {'mesh-quality', 'sculpt', 'youtube'} <= set(distribution.metadata.get_all('Provides-Extra', []))
    from cubit_mesh_export.mcp.api_reference import get_api_reference
    assert get_api_reference("all")


def test_cubit_runtime_has_no_other_product_imports():
    import cubit_mesh_export
    root = Path(cubit_mesh_export.__file__).parent
    violations = []
    for path in root.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8-sig"))
        for node in ast.walk(tree):
            modules = ([node.module or ""] if isinstance(node, ast.ImportFrom)
                       else [alias.name for alias in node.names] if isinstance(node, ast.Import)
                       else [])
            if any(name.split('.')[0] in {"radia", "radia_mcp", "cae_mcp_core"} for name in modules):
                violations.append(f"{path.relative_to(root)}:{node.lineno}")
    assert not violations, violations


def test_cubit_support_contains_only_cubit_example_providers():
    from cubit_mesh_export.mcp._support import examples
    assert set(examples.FAMILIES) == {'cubit'}
    assert set(examples.REFRESH_FUNCS) == set(examples.FAMILIES['cubit'])
    assert not hasattr(examples, 'refresh_build123d_examples')


def test_cubit_state_does_not_use_radia_settings(monkeypatch, tmp_path):
    from cubit_mesh_export.mcp._support.failure_log import state_dir
    monkeypatch.setenv('RADIA_MCP_STATE_DIR', str(tmp_path / 'radia'))
    monkeypatch.delenv('CUBIT_MCP_STATE_DIR', raising=False)
    assert state_dir().name in {'cubit-mesh-export', '.cubit-mesh-export'}
    monkeypatch.setenv('CUBIT_MCP_STATE_DIR', str(tmp_path / 'cubit'))
    assert state_dir() == tmp_path / 'cubit'


def test_optional_example_roots_survive_runtime_relocation(monkeypatch, tmp_path):
    from cubit_mesh_export.mcp._support import examples
    (tmp_path / '.git').mkdir()
    docs = tmp_path / 'docs'
    docs.mkdir()
    monkeypatch.setattr(examples, '__file__', str(tmp_path / 'arbitrary/deep/runtime/examples.py'))
    monkeypatch.chdir(tmp_path.parent)
    assert examples._resolve_local_root('repo:/docs') == docs


def test_support_ownership_notices_and_product_identity():
    from cubit_mesh_export.mcp import _support
    from cubit_mesh_export.mcp._support import examples, web_docs
    root = Path(_support.__file__).parent
    assert 'BSD 3-Clause License' in (root / 'LICENSE-BSD-3-Clause.txt').read_text()
    assert 'intentionally independent implementations' in (root / 'OWNERSHIP.md').read_text()
    for path in root.glob('*.py'):
        text = path.read_text(encoding='utf-8')
        assert 'intentionally maintained independently' in text, path
        assert 'RADIA_MCP_' not in text, path
        assert '"schema": "radia-mcp.' not in text, path
    assert 'cubit-mesh-export' in web_docs._USER_AGENT
    assert 'radia-mcp' not in web_docs._USER_AGENT
    assert 'bd_warehouse' not in examples.search_examples.__doc__
    assert 'build123d' not in examples.search_examples.__doc__


def test_status_helper_import_does_not_load_numerical_or_gui_runtimes():
    import subprocess
    import sys
    from cubit_mesh_export.mcp._support import status
    code = '''
import sys, importlib.util
class RejectHeavyImports:
    def find_spec(self, fullname, *args):
        if fullname.split('.')[0] in {'netgen', 'ngsolve', 'gmsh', 'PySide6'}:
            raise AssertionError('Unexpected runtime import: ' + fullname)
sys.meta_path.insert(0, RejectHeavyImports())
spec = importlib.util.spec_from_file_location('isolated_status_helper', sys.argv[1])
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
'''
    # Isolate the helper: the parent package intentionally initializes Netgen's
    # DLL search path for its bundled native extension, a separate contract.
    subprocess.run([sys.executable, '-c', code, status.__file__], check=True)
