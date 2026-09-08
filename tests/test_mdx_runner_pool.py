"""The release requires both compute hosts; preflight selects idle capacity."""
import importlib.util
import base64
import builtins
import sys
from types import SimpleNamespace
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load_tool(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'tools' / f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def runner(host, *, busy=False, status='online'):
    return dict(labels=[dict(name='mdx'), dict(name=host)], busy=busy, status=status)


def test_preflight_selects_free_runner_and_rejects_unavailable_pool():
    tool = load_tool('ci_preflight_mdx')
    assert tool.select_idle_host([runner('mdx1', busy=True), runner('mdx2')]) == 'mdx2'
    assert tool.select_idle_host([runner('mdx1'), runner('mdx2', status='offline')]) == 'mdx1'
    with pytest.raises(RuntimeError, match='No idle'):
        tool.select_idle_host([runner('mdx1', busy=True), runner('mdx2', status='offline')])


def test_quad_deploys_both_compute_hosts_without_mcp(monkeypatch):
    tool = load_tool('release_quad')
    assert set(tool.SIMULINK_TARGETS) == {'lab', '100', 'mdx1', 'mdx2'}
    calls = []
    monkeypatch.setattr(tool, '_deploy_pypi', lambda host, label, **kwargs: calls.append((host, kwargs)) or 0)
    assert tool.cmd_phase8e(None) == 0
    assert [host for host, _ in calls] == ['mdx1', 'mdx2']
    assert all(options['include_mcp'] is False for _, options in calls)


def test_quad_stops_when_first_compute_deployment_fails(monkeypatch):
    tool = load_tool('release_quad')
    calls = []
    monkeypatch.setattr(tool, '_deploy_mdx', lambda host: calls.append(host) or 3)
    assert tool.cmd_phase8e(None) == 3
    assert calls == ['mdx1']


@pytest.mark.parametrize('host', ['mdx1', 'mdx2'])
def test_compute_deployment_never_installs_or_runs_cubit(monkeypatch, host):
    tool = load_tool('release_quad')
    calls = []
    monkeypatch.setattr(tool, '_read_repo_versions', lambda: {
        'radia': '1.0', 'radia-mcp': '2.0', 'cubit-mesh-export': '3.0'})

    def run(command, **kwargs):
        calls.append(command)
        return SimpleNamespace(stdout='radia (1.0)', returncode=0)

    monkeypatch.setattr(tool, 'run', run)
    assert tool._deploy_mdx(host) == 0
    assert [cmd[-1] for cmd in calls if cmd[:4] == ['python', '-m', 'pip', 'index']] == ['radia']
    command = next(cmd for cmd in calls if cmd[0] == 'ssh')
    assert command[1] == host
    script = base64.b64decode(command[-1]).decode('utf-16le')
    executable = '\n'.join(line for line in script.splitlines() if not line.lstrip().startswith('#'))
    assert '"radia==1.0"' in executable
    for forbidden in ('radia[cubit]', 'cubit-mesh-export', 'cubit-plugin-install',
                      'cubit-smoke-test', 'coreform_cubit.exe', 'cubit.exe'):
        assert forbidden not in executable


def test_compute_probe_does_not_import_cubit_and_keeps_row_contract(monkeypatch, capsys, tmp_path):
    tool = load_tool('release_quad')
    radia = SimpleNamespace(__file__=str(tmp_path / '__init__.py'), __version__='1.0')
    monkeypatch.setitem(sys.modules, 'radia', radia)
    original_import = builtins.__import__

    def guarded_import(name, *args, **kwargs):
        assert name != 'cubit_mesh_export', 'compute probe must not load Cubit package'
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', guarded_import)
    exec(tool.CROSS_MACHINE_PROBE_NO_MCP, {})
    compute_rows = capsys.readouterr().out.splitlines()
    assert len(compute_rows) == 11
    assert sum(row.endswith('= N/A') for row in compute_rows) == 4
    assert compute_rows[0].endswith('= 1.0')
    assert 'import radia, cubit_mesh_export' in tool.CROSS_MACHINE_PROBE_LAB
