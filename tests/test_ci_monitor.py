"""CI monitoring must work without trusting the caller's Git worktree."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest


@pytest.mark.parametrize('repository', [None, 'owner/another-repository'])
def test_monitor_passes_explicit_repository_to_gh(monkeypatch, repository):
    path = Path(__file__).resolve().parents[1] / '.agents/skills/ci-monitor/monitor.py'
    spec = importlib.util.spec_from_file_location('ci_monitor_under_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.delenv('GH_REPO', raising=False)
    if repository:
        monkeypatch.setenv('GH_REPO', repository)
    calls = []

    def capture(command, **kwargs):
        calls.append(kwargs['env']['GH_REPO'])
        return SimpleNamespace(returncode=0, stdout='{}', stderr='')

    monkeypatch.setattr(module.subprocess, 'run', capture)
    assert module._gh(['run', 'view', '123'])[0] == 0
    assert calls == [repository or 'ksugahar/Radia']
