"""Owned headless transport, fail-closed startup and no ambiguous replay."""
import concurrent.futures
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import pytest

from cubit_mesh_export.mcp import session


def test_constructor_honors_explicit_install_override(tmp_path, monkeypatch):
    monkeypatch.setattr(session, '_OVERRIDE_BIN_DIR', tmp_path)
    assert session.CubitSession()._bin_dir == tmp_path


def test_install_directory_environment_is_honored(tmp_path, monkeypatch):
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    (bin_dir / 'cubit.py').touch()
    monkeypatch.delenv('CUBIT_BIN_DIR', raising=False)
    monkeypatch.setenv('CUBIT_INSTALL_DIR', str(tmp_path))
    assert session.find_cubit_install() == bin_dir.resolve()


def test_real_constructor_initializes_native_journal_identity(tmp_path, monkeypatch):
    monkeypatch.setenv('CUBIT_MCP_TEMP', str(tmp_path))
    sess = session.CubitSession(tmp_path)
    def accept(req, **kw):
        return {'id': req['id'], 'ok': True, 'result': [{'ok': True}]}
    monkeypatch.setattr(sess, '_call_via_stdio', accept)
    sess._start_native_journal_locked(1)
    assert sess._client_id in sess._native_journal_path.name
    assert 'cubit-mesh-export' in sess._native_journal_path.parts


@pytest.mark.parametrize('mode', ['gui', 'auto', 'existing', 'new'])
def test_constructor_rejects_retired_modes_before_discovery(mode, monkeypatch):
    monkeypatch.setattr(session, 'find_cubit_install', lambda: pytest.fail('must not discover'))
    with pytest.raises(ValueError, match='Only headless'):
        session.CubitSession(mode=mode)


def test_failed_command_is_not_replayed(tmp_path, monkeypatch):
    sess = session.CubitSession(tmp_path)
    calls = []
    monkeypatch.setattr(sess, 'ensure_started', lambda: {})
    monkeypatch.setattr(sess, '_start_native_journal_locked', lambda **kw: None)
    def fail(req, **kw):
        calls.append(req)
        raise session.CubitSessionError('response lost after execution')
    monkeypatch.setattr(sess, '_call_via_stdio', fail)
    with pytest.raises(session.CubitSessionError, match='No automatic replay'):
        sess.call('cmd', ['create brick x 1'])
    assert len(calls) == 1
    assert sess._proc is None
    assert not sess._owned


def test_dead_child_reports_state_loss_before_new_commands(tmp_path, monkeypatch):
    from types import SimpleNamespace
    sess = session.CubitSession(tmp_path)
    sess._proc = SimpleNamespace(poll=lambda: 1, wait=lambda **kw: 1,
                                 stdin=None, stdout=None, stderr=None)
    monkeypatch.setattr(sess, '_start_stdio_daemon', lambda: pytest.fail('silent restart'))
    with pytest.raises(session.CubitSessionError, match='geometry state was lost'):
        sess.ensure_started()
    assert sess._proc is None


def test_concurrent_startup_creates_one_child(tmp_path, monkeypatch):
    from types import SimpleNamespace
    sess = session.CubitSession(tmp_path)
    starts = []
    def start():
        starts.append(1)
        time.sleep(0.01)
        sess._proc = SimpleNamespace(poll=lambda: None)
        sess._ready_info = {'ready': True}
        return sess._ready_info
    monkeypatch.setattr(sess, '_start_stdio_daemon', start)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        assert all(pool.map(lambda _: sess.ensure_started()['ready'], range(16)))
    assert len(starts) == 1


def test_concurrent_calls_are_serialized(tmp_path, monkeypatch):
    sess = session.CubitSession(tmp_path)
    monkeypatch.setattr(sess, 'ensure_started', lambda: {})
    active = []
    def reply(req, **kw):
        assert not active
        active.append(req['id'])
        time.sleep(0.005)
        active.pop()
        return {'id': req['id'], 'ok': True, 'result': req['args']}
    monkeypatch.setattr(sess, '_call_via_stdio', reply)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(lambda i: sess.call('probe', [i]), range(24)))
    assert len({r['id'] for r in results}) == 24
    assert [r['result'] for r in results] == [[i] for i in range(24)]


@pytest.mark.parametrize('response', [b'not json', b'[]', b'{"ok":true}', b'{"id":99}'])
def test_invalid_response_fails_closed(tmp_path, monkeypatch, response):
    sess = session.CubitSession(tmp_path)
    monkeypatch.setattr(sess, '_send_stdio', lambda req: None)
    monkeypatch.setattr(sess, '_read_stdio_line', lambda **kw: response)
    with pytest.raises(session.CubitSessionError):
        sess._call_via_stdio({'id': 1}, 1)


@pytest.mark.parametrize('startup', ['silent', 'malformed', 'wrong_protocol'])
def test_failed_startup_reaps_child(tmp_path, monkeypatch, startup):
    sess = session.CubitSession(tmp_path)
    real_popen = subprocess.Popen
    children = []
    code = 'import time; time.sleep(60)'
    if startup != 'silent':
        line = 'not json' if startup == 'malformed' else json.dumps({'ready': True, 'protocol_version': 2})
        code = f'print({line!r}, flush=True); ' + code
    def fake_daemon(argv, **kwargs):
        proc = real_popen([sys.executable, '-c', code], **kwargs)
        children.append(proc)
        return proc
    monkeypatch.setattr(session, '_cubit_python_exe', lambda root: Path(sys.executable))
    monkeypatch.setattr(session.subprocess, 'Popen', fake_daemon)
    monkeypatch.setattr(session, 'CUBIT_READY_TIMEOUT_S', 0.2)
    try:
        with pytest.raises(session.CubitSessionError):
            sess.ensure_started()
        assert sess._proc is None
        assert children[0].poll() is not None
        assert all(p.closed for p in (children[0].stdin, children[0].stdout, children[0].stderr))
    finally:
        for child in children:
            if child.poll() is None:
                child.kill()
                child.wait(timeout=5)


def test_daemon_rejects_gui_before_importing_cubit():
    daemon = Path(session.__file__).with_name('daemon.py')
    result = subprocess.run([sys.executable, str(daemon)], capture_output=True,
                            env={**os.environ, 'CUBIT_DAEMON_MODE': 'gui'}, timeout=10)
    assert result.returncode == 2
    reply = json.loads(result.stdout)
    assert reply['ready'] is False
    assert 'Only headless' in reply['error']
