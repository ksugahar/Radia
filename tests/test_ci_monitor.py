"""CI monitoring must work without trusting the caller's Git worktree."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest


@pytest.fixture
def monitor():
    path = Path(__file__).resolve().parents[1] / '.agents/skills/ci-monitor/monitor.py'
    spec = importlib.util.spec_from_file_location('ci_monitor_under_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize('repository', [None, 'owner/another-repository'])
def test_monitor_passes_explicit_repository_to_gh(monkeypatch, repository, monitor):
    module = monitor
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


@pytest.mark.parametrize('conclusion', ['cancelled', 'failure', 'timed_out', 'action_required'])
def test_unsuccessful_runs_cannot_report_all_green(monkeypatch, capsys, monitor, conclusion):
    monkeypatch.setattr(monitor, 'get_run_state', lambda _: dict(
        status='completed', conclusion=conclusion, name='CI', headBranch='main'))
    monkeypatch.setattr(monitor, 'fetch_failing_log_tail', lambda *_: 'diagnostic')
    assert monitor.watch_runs(['123'], 0, 10) == 1
    assert 'ALL GREEN' not in capsys.readouterr().out


def test_empty_monitor_is_not_success(capsys, monitor):
    assert monitor.watch_runs([], 0, 10) == 2
    captured = capsys.readouterr()
    assert 'ALL GREEN' not in captured.out
    assert 'no CI runs' in captured.err


@pytest.mark.parametrize('response', [
    (124, '', 'timeout'),
    (1, '', 'permission denied'),
    (0, 'invalid JSON', ''),
    (0, '[]', ''),
    (0, '{}', ''),
    (0, '{"status": "completed", "conclusion": [], "name": "CI", "headBranch": "main"}', ''),
])
def test_monitor_errors_are_not_ci_failures(monkeypatch, capsys, monitor, response):
    monkeypatch.setattr(monitor, '_gh', lambda *_args, **_kwargs: response)
    def unexpected_log(*_):
        pytest.fail('No failed-job log is available for a monitor error')
    monkeypatch.setattr(monitor, 'fetch_failing_log_tail', unexpected_log)
    assert monitor.watch_runs(['123'], 0, 10) == 2
    output = capsys.readouterr().out
    assert 'ALL GREEN' not in output
    assert 'monitor error' in output


def test_completed_success_after_pending_state(monkeypatch, capsys, monitor):
    states = iter([
        (0, '{"status":"in_progress","conclusion":null,"name":"CI","headBranch":"main"}', ''),
        (0, '{"status":"completed","conclusion":"success","name":"CI","headBranch":"main"}', ''),
    ])
    monkeypatch.setattr(monitor, '_gh', lambda *_args, **_kwargs: next(states))
    assert monitor.watch_runs(['123'], 0, 10) == 0
    assert 'ALL GREEN (1 success' in capsys.readouterr().out


@pytest.mark.parametrize('args,expected', [
    ([], (None, 3)),
    (['--branch', 'release/v5'], ('release/v5', 3)),
    (['--branch', 'main', '--auto', '2'], ('main', 2)),
    (['--auto', '1'], (None, 1)),
])
def test_cli_preserves_discovery_scope(monkeypatch, monitor, args, expected):
    discovered = []
    watched = []
    def discover(branch, count):
        discovered.append((branch, count))
        return ['123']
    monkeypatch.setattr(monitor.sys, 'argv', ['monitor.py', *args])
    monkeypatch.setattr(monitor, 'discover_runs', discover)
    monkeypatch.setattr(monitor, 'watch_runs',
                        lambda ids, *_: watched.append(ids) or 0)
    assert monitor.main() == 0
    assert discovered == [expected]
    assert watched == [['123']]


def test_explicit_ids_never_trigger_discovery(monkeypatch, monitor):
    monkeypatch.setattr(monitor.sys, 'argv',
                        ['monitor.py', '123', '456', '--branch', 'main'])
    def unexpected(*_):
        pytest.fail('Explicit run IDs must not be replaced by discovery')
    monkeypatch.setattr(monitor, 'discover_runs', unexpected)
    watched = []
    monkeypatch.setattr(monitor, 'watch_runs',
                        lambda ids, *_: watched.append(ids) or 0)
    assert monitor.main() == 0
    assert watched == [['123', '456']]


@pytest.mark.parametrize('args', [
    ['--auto', '-1'], ['--poll', '0'], ['--poll', '-1'],
    ['--tail', '0'], ['--tail', '-1'],
])
def test_invalid_cli_limits_fail_before_discovery(monkeypatch, monitor, args):
    monkeypatch.setattr(monitor.sys, 'argv', ['monitor.py', *args])
    def unexpected(*_):
        pytest.fail('Invalid limits must fail before any network calls')
    monkeypatch.setattr(monitor, 'discover_runs', unexpected)
    with pytest.raises(SystemExit) as exc:
        monitor.main()
    assert exc.value.code == 2
