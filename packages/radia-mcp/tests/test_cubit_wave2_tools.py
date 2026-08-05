"""Tests for the Wave-2 cubit MCP tools: cubit_doctor,
cubit_session_journal, the all-calls JSONL log, and the session
command-history recorder.  No Cubit license needed."""

import asyncio
import json

from radia_mcp.cubit import server as cubit_server
from radia_mcp.cubit import session as cubit_session
from radia_mcp.cubit.server import cubit_doctor, cubit_session_journal


def test_doctor_returns_per_check_statuses():
    out = json.loads(cubit_doctor())
    assert out["status"] in ("ok", "problems_found")
    checks = out["checks"]
    for key in ("install", "license", "plugin", "daemon", "drop_dir",
                "check_vol_deps"):
        assert key in checks, key
        assert checks[key].get("status") in ("ok", "warn", "error",
                                             "skipped"), checks[key]
    # problems entries must correspond to non-ok checks
    if out["status"] == "ok":
        assert out["problems"] == []


def test_session_journal_empty_history():
    out = json.loads(cubit_session_journal())
    assert out["status"] == "ok"
    assert out["n_commands"] == 0


def test_record_cmd_history_and_journal_format(monkeypatch, tmp_path):
    sess = cubit_session.CubitSession.__new__(cubit_session.CubitSession)
    sess._command_history = []
    sess._command_history_max = 10
    resp = {"ok": True, "result": [
        {"line": "brick x 1", "ok": True, "rc": 1},
        {"line": "mesh volume 99", "ok": False, "rc": 1},
    ]}
    sess._record_cmd_history(resp)
    assert [e["ok"] for e in sess._command_history] == [True, False]

    # bounded history
    for i in range(20):
        sess._record_cmd_history(
            {"result": [{"line": f"cmd {i}", "ok": True}]})
    assert len(sess._command_history) == 10

    # journal formatting through the tool (patch the singleton)
    sess._command_history = [
        {"ts": 0.0, "line": "brick x 1", "ok": True},
        {"ts": 0.0, "line": "mesh volume 99", "ok": False},
    ]
    monkeypatch.setattr(cubit_session, "_SINGLETON", sess)
    out_file = tmp_path / "session.jou"
    out = json.loads(cubit_session_journal(out_path=str(out_file)))
    assert out["n_commands"] == 2 and out["n_failed"] == 1
    text = out_file.read_text(encoding="utf-8")
    assert "brick x 1" in text
    assert "# FAILED: mesh volume 99" in text


def test_attach_ping_detects_unresponsive_daemon(monkeypatch, tmp_path):
    """A pid.lock pointing at a LIVE pid (this test process) with a ready
    marker but nobody answering the file-drop ping must fail fast with a
    clear hung-daemon message -- and must NOT kill the pid."""
    import os
    import pytest

    drop = tmp_path / "cubit-session"
    drop.mkdir()
    (drop / "out").mkdir()
    (drop / "pid.lock").write_text(str(os.getpid()), encoding="utf-8")
    (drop / "ready").write_text(json.dumps(
        {"ready": True, "protocol_version": 2, "pid": os.getpid(),
         "drop": str(drop)}), encoding="utf-8")

    monkeypatch.setattr(cubit_session, "_user_daemon_dir", lambda: drop)
    monkeypatch.setattr(cubit_session, "ATTACH_PING_TIMEOUT_S", 0.4)
    sess = cubit_session.CubitSession.__new__(cubit_session.CubitSession)
    sess._bin_dir = tmp_path          # never used before the attach branch
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

    with pytest.raises(cubit_session.CubitSessionError) as ei:
        sess._start_gui_bootstrap()
    msg = str(ei.value)
    assert "did not answer a ping" in msg
    assert "cubit_session_shutdown" in msg
    # non-destructive: markers left for the daemon's other clients
    assert (drop / "pid.lock").is_file()
    assert (drop / "ready").is_file()


def test_call_log_records_tool_calls(monkeypatch, tmp_path):
    # Redirect the state dir, then re-install the choke point so the
    # log lands in tmp_path.
    monkeypatch.setenv("RADIA_MCP_STATE_DIR", str(tmp_path))
    monkeypatch.setattr(cubit_server.mcp._tool_manager, "call_tool",
                        cubit_server.mcp._tool_manager.call_tool)
    cubit_server._install_call_log()
    asyncio.run(cubit_server.mcp.call_tool(
        "cubit_session_journal", {}))
    log = tmp_path / "logs" / "cubit_tool_calls.jsonl"
    assert log.is_file()
    lines = [json.loads(l) for l in
             log.read_text(encoding="utf-8").splitlines()]
    assert any(r["tool"] == "cubit_session_journal" and r["ok"]
               and "ms" in r for r in lines)
