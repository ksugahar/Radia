"""The release requires both compute hosts; preflight selects idle capacity."""
import importlib.util
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
