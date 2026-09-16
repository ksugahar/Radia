"""Owned process cleanup, setup CLI and current headless status."""

import json
import subprocess
import sys

import pytest

from cubit_mesh_export.mcp import server as cubit_server
from cubit_mesh_export.mcp import session as cubit_session
from cubit_mesh_export.mcp.server import cubit_session_status


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


def test_close_process_streams_closes_batch_pipes():
    proc = subprocess.Popen(
        [sys.executable, "-c", "pass"], stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc.wait(timeout=10.0)

    cubit_session._close_proc_streams(proc)

    assert proc.stdin.closed
    assert proc.stdout.closed
    assert proc.stderr.closed


# ---------------------------------------------------------------------------
# Setup CLI and status
# ---------------------------------------------------------------------------


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


def test_session_status_reports_headless_mode_and_journal():
    out = json.loads(cubit_session_status())
    assert out["execution_mode"] == "batch"
    assert "session_mode" not in out
    assert out["gui_started"] is False
    # keys present whenever a singleton exists are optional here; the
    # baseline contract is bin_dir/alive/execution_mode
    assert "alive" in out and "bin_dir" in out
