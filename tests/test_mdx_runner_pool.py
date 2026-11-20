"""Guards: retired routes have no actions and probes stay solver-only.

Assumes: synthetic imports and no remote commands.
"""
import importlib.util
import importlib.metadata
import base64
import builtins
import sys
from types import SimpleNamespace
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_tool(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'tools' / f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_retired_compute_routes_are_deleted():
    """Guards: retired hosts have no command, target or shim left to call."""
    tool = load_tool("release_quad")
    for name in ("_deploy_mdx", "cmd_phase8e", "cmd_temp_shadows", "SIMULINK_TARGETS", "SSH_100"):
        assert not hasattr(tool, name), name
    assert tuple(tool.SOLVER_TARGETS) == ("validation",), tuple(tool.SOLVER_TARGETS)
    assert not (ROOT / "tools" / "release_temp_shadows.py").exists()


def test_compute_probe_is_solver_only_and_keeps_row_contract(monkeypatch, capsys, tmp_path):
    tool = load_tool('release_quad')
    radia = SimpleNamespace(__file__=str(tmp_path / '__init__.py'), __version__='1.0')
    monkeypatch.setitem(sys.modules, 'radia', radia)
    solver_versions = {'ngsolve': '6.2.2608', 'netgen-mesher': '6.2.2608'}
    monkeypatch.setattr(
        importlib.metadata, 'version', lambda name: solver_versions[name])
    original_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        assert name != 'cubit_mesh_export', 'compute probe must not load Cubit package'
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', guarded_import)
    exec(tool.CROSS_MACHINE_PROBE, {})
    compute_rows = capsys.readouterr().out.splitlines()
    rows = {" ".join(row.partition('=')[0].split()): row.partition('=')[2].strip()
            for row in compute_rows}
    assert len(compute_rows) == len(rows)  # duplicate rows must not disappear
    assert set(rows) == set(tool._PHASE9_FIELDS)
    assert all(value != 'N/A' for value in rows.values())
    assert compute_rows[0].endswith('= 1.0')
    assert 'cubit_mesh_export' not in tool.CROSS_MACHINE_PROBE
    assert 'radia-mcp' not in tool.CROSS_MACHINE_PROBE
