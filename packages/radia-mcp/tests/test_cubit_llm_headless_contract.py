import inspect
import json
import threading
from pathlib import Path

import pytest

from radia_mcp.build123d import server as build123d_server
from radia_mcp.cubit import server, session


def test_cubit_session_defaults_to_batch(tmp_path, monkeypatch):
    monkeypatch.setattr(session, "find_cubit_install", lambda: tmp_path)

    sess = session.CubitSession()

    assert sess._mode == "batch"


def test_singleton_refuses_execution_mode_drift(tmp_path, monkeypatch):
    monkeypatch.setattr(session, "find_cubit_install", lambda: tmp_path)
    monkeypatch.setattr(session, "_SINGLETON", None)

    headless = session.CubitSession.get(mode="batch")
    assert headless._mode == "batch"
    with pytest.raises(session.CubitSessionError, match="refusing requested mode='gui'"):
        session.CubitSession.get(mode="gui")


def test_mcp_session_helper_requests_batch_explicitly(monkeypatch):
    requested = []

    class FakeSession:
        @classmethod
        def get(cls, mode):
            requested.append(mode)
            return cls()

        def ensure_started(self):
            return {"ready": True, "mode": "batch"}

    monkeypatch.setattr(server._cs, "CubitSession", FakeSession)

    sess, error = server._cubit_session_or_error()

    assert error is None
    assert isinstance(sess, FakeSession)
    assert requested == ["batch"]


def test_batch_rpc_response_reports_that_no_gui_started(monkeypatch):
    sess = session.CubitSession.__new__(session.CubitSession)
    sess._mode = "batch"
    sess._lock = threading.Lock()
    sess._next_id = 1
    sess._command_history = []
    sess._command_history_max = 10
    monkeypatch.setattr(sess, "ensure_started", lambda: {"ready": True})
    monkeypatch.setattr(
        sess, "_call_via_stdio",
        lambda request, timeout_s: {
            "id": request["id"], "ok": True, "result": "pong"
        },
    )

    response = sess.call("ping")

    assert response["execution_mode"] == "batch"
    assert response["gui_started"] is False


def test_cubit_load_uses_only_the_headless_session(tmp_path, monkeypatch):
    step = tmp_path / "part.step"
    step.write_text("ISO-10303-21;", encoding="ascii")

    class FakeSession:
        def call(self, op, args, timeout_s=None):
            if op == "cmd":
                return {"ok": True, "result": [
                    {"line": args[0], "ok": True, "rc": 1}
                ]}
            if op == "probe":
                return {"ok": True, "result": {"volumes": 1}}
            raise AssertionError(op)

    monkeypatch.setattr(
        server, "_cubit_session_or_error", lambda: (FakeSession(), None)
    )

    result = json.loads(server.cubit_load(path=str(step)))

    assert result["status"] == "ok"
    assert result["mode"] == "headless_persistent"
    assert result["gui_started"] is False
    assert result["summary"] == {"volumes": 1}
    source = inspect.getsource(server.cubit_load)
    assert "_cubit_gui_exe" not in source
    assert ".Popen(" not in source


def test_public_mesh_tools_expose_only_headless_session_language():
    public_tools = (
        server.cubit_load,
        server.cubit_mesh_auto,
        server.cubit_mesh_race,
        server.cubit_mesh_race_review,
        server.cubit_mesh_race_review_async,
        server.cubit_mesh_race_smart,
        server.cubit_mesh_race_smart_async,
        build123d_server.cadquery_to_cubit_hex,
        build123d_server.build123d_to_cubit_hex,
    )

    for tool in public_tools:
        parameters = inspect.signature(tool).parameters
        assert "commit_to_gui" not in parameters
        assert "cubit_exe" not in parameters
        assert "detach" not in parameters
        assert "include_human" not in parameters

    assert "apply_to_session" in inspect.signature(
        server.cubit_mesh_auto
    ).parameters
    assert "apply_to_session" in inspect.signature(
        server.cubit_mesh_race
    ).parameters


def test_mcp_surface_does_not_publish_gui_execution_tools():
    tool_names = set(server.mcp._tool_manager._tools)

    assert "cubit_load" in tool_names
    assert "open_in_cubit" not in tool_names
    assert "cubit_mesh_race_with_human" not in tool_names


def test_snapshot_fails_without_starting_a_gui():
    snapshot = json.loads(server.cubit_snapshot("unused.png"))

    assert snapshot["kind"] == "policy"
    assert snapshot["gui_started"] is False


def test_headless_journal_never_falls_back_to_gui_exe(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "coreform_cubit.exe").write_bytes(b"")
    monkeypatch.setattr(session, "get_cubit_bin_dir", lambda: bin_dir)

    def unexpected_run(*args, **kwargs):
        raise AssertionError("GUI fallback must not spawn a process")

    monkeypatch.setattr(session.subprocess, "run", unexpected_run)

    result = session.run_headless_journal(["reset"])

    assert result["status"] == "error"
    assert result["kind"] == "environment"
    assert "Refusing to fall back to the GUI launcher" in result["error"]


def test_all_llm_bridges_name_the_batch_session():
    cubit_source = inspect.getsource(server._cubit_session_or_error)
    build123d_source = inspect.getsource(build123d_server)
    daemon_source = Path(session.__file__).with_name("daemon.py").read_text(
        encoding="utf-8"
    )

    assert 'CubitSession.get(mode="batch")' in cubit_source
    assert 'CubitSession.get(mode="batch")' in build123d_source
    assert 'CubitSession.get(mode="gui")' not in cubit_source
    assert 'CubitSession.get(mode="gui")' not in build123d_source
    assert 'os.environ.get("CUBIT_DAEMON_MODE", "batch")' in daemon_source
    assert "Never launch or attach to the Cubit GUI" in server._SERVER_INSTRUCTIONS
