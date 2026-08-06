"""Tests for Wave-3 cubit MCP features: session-mode triad, --setup CLI
mode, session-status enrichment, and the version-robust lint golden
((rule, line) pairs, not message text)."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from radia_mcp.cubit import server as cubit_server
from radia_mcp.cubit import session as cubit_session
from radia_mcp.cubit.server import _lint_file, cubit_session_status


@pytest.mark.skipif(sys.platform != "win32", reason="Windows Job Object")
def test_kill_on_close_job_terminates_private_process():
    proc = subprocess.Popen(
        [sys.executable, "-c", "import time; time.sleep(60)"],
        creationflags=getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0),
    )
    handle = None
    try:
        handle = cubit_session._assign_kill_on_close_job(proc.pid)
        assert handle
        assert proc.poll() is None
        cubit_session._close_windows_handle(handle)
        handle = None
        assert proc.wait(timeout=5) is not None
    finally:
        if handle is not None:
            cubit_session._close_windows_handle(handle)
        if proc.poll() is None:
            proc.kill()
            proc.wait(timeout=5)


def test_private_drop_cleanup_retries_transient_windows_locks(
        monkeypatch, tmp_path):
    drop = tmp_path / "cubit-session-private"
    drop.mkdir()
    (drop / "cubit_stdout.log").write_text("log", encoding="utf-8")
    real_rmtree = cubit_session.shutil.rmtree
    calls = 0

    def temporarily_locked(path, ignore_errors=False):
        nonlocal calls
        calls += 1
        if calls >= 3:
            real_rmtree(path, ignore_errors=ignore_errors)

    monkeypatch.setattr(cubit_session.shutil, "rmtree", temporarily_locked)

    assert cubit_session._remove_tree_with_retry(drop, timeout_s=1.0)
    assert calls == 3
    assert not drop.exists()


def test_close_process_streams_closes_batch_pipes():
    proc = subprocess.Popen(
        [sys.executable, "-c", "pass"], stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait(timeout=10.0)

    cubit_session._close_process_streams(proc)

    assert proc.stdin.closed
    assert proc.stdout.closed
    assert proc.stderr.closed


# ---------------------------------------------------------------------------
# Session-mode triad
# ---------------------------------------------------------------------------

def _bare_session(tmp_path):
    sess = cubit_session.CubitSession.__new__(cubit_session.CubitSession)
    sess._bin_dir = tmp_path
    sess._mode = "gui"
    sess._proc = None
    sess._next_id = 1
    sess._ready_info = None
    sess._drop_dir = None
    sess._outbox = None
    sess._owned = False
    sess._last_license_warmup = {}
    sess._command_history = []
    sess._command_history_max = 10
    return sess


def test_invalid_session_mode_fails_loud(monkeypatch, tmp_path):
    monkeypatch.setenv("RADIA_CUBIT_SESSION_MODE", "sometimes")
    sess = _bare_session(tmp_path)
    with pytest.raises(cubit_session.CubitSessionError) as ei:
        sess._start_gui_bootstrap()
    assert "auto | new | existing" in str(ei.value)


def test_existing_mode_without_daemon_fails_loud(monkeypatch, tmp_path):
    monkeypatch.setenv("RADIA_CUBIT_SESSION_MODE", "existing")
    monkeypatch.setattr(cubit_session, "_user_daemon_dir",
                        lambda: tmp_path / "cubit-session")
    sess = _bare_session(tmp_path)
    with pytest.raises(cubit_session.CubitSessionError) as ei:
        sess._start_gui_bootstrap()
    assert "no live shared" in str(ei.value)


def test_new_mode_uses_private_drop_dir(monkeypatch, tmp_path):
    """Mode 'new' must ignore a live shared daemon and pick a private
    per-process drop dir (never clobber the shared one)."""
    import os

    shared = tmp_path / "cubit-session"
    shared.mkdir()
    (shared / "out").mkdir()
    (shared / "pid.lock").write_text(str(os.getpid()), encoding="utf-8")
    (shared / "ready").write_text(json.dumps(
        {"ready": True, "protocol_version": 2, "pid": os.getpid()}),
        encoding="utf-8")

    monkeypatch.setenv("RADIA_CUBIT_SESSION_MODE", "new")
    monkeypatch.setattr(cubit_session, "_user_daemon_dir", lambda: shared)
    sess = _bare_session(tmp_path)
    # Spawn path fails at launcher discovery (no Cubit under tmp_path) --
    # but by then the private drop dir must already exist and the shared
    # daemon's markers must be untouched.
    with pytest.raises((cubit_session.CubitSessionError, FileNotFoundError)):
        sess._start_gui_bootstrap()
    private = shared.parent / f"cubit-session-{os.getpid()}"
    assert private.is_dir()
    assert (shared / "pid.lock").is_file()          # shared daemon untouched
    assert (shared / "ready").is_file()


def test_attached_clients_use_distinct_request_filenames(tmp_path):
    first = _bare_session(tmp_path)
    second = _bare_session(tmp_path)

    first_stem = first._request_stem(1)
    second_stem = second._request_stem(1)

    assert first_stem != second_stem
    assert first_stem.endswith("-00000001")
    assert second_stem.endswith("-00000001")


# ---------------------------------------------------------------------------
# --setup CLI mode + status enrichment
# ---------------------------------------------------------------------------

def test_setup_mode_runs_and_reports(monkeypatch, capsys, tmp_path):
    monkeypatch.setattr(cubit_server._cs, "find_cubit_install",
                        lambda: tmp_path)
    monkeypatch.setattr(
        cubit_server,
        "cubit_doctor",
        lambda: json.dumps({"status": "ok", "problems": []}),
    )
    # Avoid a real license warmup: stub it.
    import radia_mcp.cubit.license_warmup as lw
    monkeypatch.setattr(lw, "warmup_license",
                        lambda *a, **k: {"status": "skipped",
                                         "reason": "test stub"})
    rc = cubit_server._setup_mode()
    out = capsys.readouterr().out
    assert "Doctor report" in out
    assert "SETUP OK" in out
    assert rc == 0


def test_setup_mode_fails_loud_when_cubit_is_missing(monkeypatch, capsys):
    monkeypatch.setattr(cubit_server._cs, "find_cubit_install", lambda: None)

    rc = cubit_server._setup_mode()
    out = capsys.readouterr().out

    assert rc == 1
    assert "Coreform Cubit not found" in out
    assert "Doctor report" not in out


def test_session_status_reports_mode_and_journal(monkeypatch):
    monkeypatch.setenv("RADIA_CUBIT_SESSION_MODE", "auto")
    out = json.loads(cubit_session_status())
    assert out["session_mode"] == "auto"
    # keys present whenever a singleton exists are optional here; the
    # baseline contract is bin_dir/alive/session_mode
    assert "alive" in out and "bin_dir" in out


# ---------------------------------------------------------------------------
# Version-robust lint golden: assert (rule, line) pairs, not messages
# ---------------------------------------------------------------------------

_FIXTURE = Path(__file__).resolve().parents[2].parent \
    / "tests" / "mcp_server" / "fixtures" / "bad_cubit_script.py"

# Locked 2026-08-05 against the committed fixture.  Assert only
# (rule, line) membership so rule WORDING can change freely; a rule
# deletion or a fixture edit must consciously update this set.
_EXPECTED_RULE_LINES = {
    ("deleted-api-usage", 11),
    ("deleted-api-usage", 46),
    ("missing-mesh-command", 50),
    ("geometry-block-2nd-order", 32),
    ("element-type-before-add", 35),
    ("wrong-connectivity-2nd-order", 39),
    ("nodeset-sideset-usage", 42),
    ("missing-boundary-block", 1),
    ("hardcoded-absolute-path", 23),
    ("wrong-file-extension", 50),
    ("missing-block-names", 1),
}


@pytest.mark.skipif(not _FIXTURE.is_file(),
                    reason="repo-level fixture not present")
def test_lint_golden_rule_line_pairs():
    findings = _lint_file(str(_FIXTURE))
    got = {(f["rule"], f["line"]) for f in findings}
    missing = _EXPECTED_RULE_LINES - got
    unexpected = got - _EXPECTED_RULE_LINES
    assert not missing, f"lint findings disappeared: {sorted(missing)}"
    assert not unexpected, (
        f"new lint findings on the golden fixture: {sorted(unexpected)} "
        "-- if intentional, update _EXPECTED_RULE_LINES")


@pytest.mark.skipif(
    not (_FIXTURE.parent / "clean_cubit_script.py").is_file(),
    reason="repo-level fixture not present")
def test_lint_clean_fixture_stays_clean():
    findings = _lint_file(str(_FIXTURE.parent / "clean_cubit_script.py"))
    assert findings == []
