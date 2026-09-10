import inspect
import json
import threading
from pathlib import Path

import pytest

from radia_mcp.cubit import server, session


@pytest.fixture
def build123d_server():
    return pytest.importorskip("radia_mcp.build123d.server")


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
    for tool in server.mcp._tool_manager._tools.values():
        parameters = tool.parameters.get("properties", {})
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
    assert "cubit_stage" in tool_names
    assert "cubit_import_journal" in tool_names
    assert "cubit_show" not in tool_names
    assert "open_in_cubit" not in tool_names
    assert "cubit_mesh_race_with_human" not in tool_names


def test_snapshot_fails_without_starting_a_gui():
    snapshot = json.loads(server.cubit_snapshot("unused.png"))

    assert snapshot["kind"] == "policy"
    assert snapshot["gui_started"] is False


@pytest.mark.parametrize("platform", ["win32", "linux", "darwin"])
def test_headless_journal_never_falls_back_to_gui_exe(tmp_path, monkeypatch, platform):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "coreform_cubit.exe").write_bytes(b"")
    monkeypatch.setattr(session, "get_cubit_bin_dir", lambda: bin_dir)
    monkeypatch.setattr(session.sys, "platform", platform)

    def unexpected_run(*args, **kwargs):
        raise AssertionError("GUI fallback must not spawn a process")

    monkeypatch.setattr(session.subprocess, "run", unexpected_run)

    result = session.run_headless_journal(["reset"])

    assert result["status"] == "error"
    assert result["kind"] == "environment"
    assert "Refusing to fall back to the GUI launcher" in result["error"]


def test_all_llm_bridges_name_the_batch_session(build123d_server):
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


def test_import_journal_reads_without_starting_or_executing(tmp_path, monkeypatch):
    from types import SimpleNamespace

    journal = tmp_path / "human.jou"
    content = '# saved\ncreate brick x 2\ncreate brick x 2\nlist volume all\n{size=3}\n'
    journal.write_bytes(content.encode("utf-8-sig"))
    existing = SimpleNamespace(
        _lock=threading.Lock(),
        _command_history=[{"line": "create brick x 2", "ok": True}],
    )
    monkeypatch.setattr(session, "_SINGLETON", existing)
    monkeypatch.setattr(server, "_cubit_session_or_error",
                        lambda: pytest.fail("Import must not start or attach Cubit"))

    result = json.loads(server.cubit_import_journal(str(journal)))

    assert result["executed"] is False
    assert result["gui_started"] is False
    assert result["journal"] == content
    assert [row["line_number"] for row in result["candidate_commands"]] == [3, 5]
    assert result["excluded"][1]["reason"] == "matches_ai_history"
    assert len(result["sha256"]) == 64
    assert len(existing._command_history) == 1


def test_import_journal_without_session_and_bad_input(tmp_path, monkeypatch):
    monkeypatch.setattr(session, "_SINGLETON", None)
    journal = tmp_path / "human.jou"
    journal.write_text('reset\n', encoding="utf-8")
    result = json.loads(server.cubit_import_journal(str(journal)))
    assert result["history_entries"] == 0
    assert result["candidate_commands"] == [{"line_number": 1, "line": "reset"}]
    assert session._SINGLETON is None
    assert json.loads(server.cubit_import_journal(str(tmp_path / 'missing.jou')))["status"] == "error"
    journal.write_bytes(b'\xff\xff')
    assert json.loads(server.cubit_import_journal(str(journal)))["status"] == "error"
