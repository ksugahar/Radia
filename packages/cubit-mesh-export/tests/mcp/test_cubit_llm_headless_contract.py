import inspect
import json
import threading
from pathlib import Path

import pytest
from cubit_mesh_export.mcp import server, session


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


def test_cubit_llm_bridge_names_the_batch_session():
    cubit_source = inspect.getsource(server._cubit_session_or_error)
    daemon_source = Path(session.__file__).with_name("daemon.py").read_text(
        encoding="utf-8"
    )

    assert 'CubitSession.get(mode="batch")' in cubit_source
    assert 'CubitSession.get(mode="gui")' not in cubit_source
    assert 'os.environ.get("CUBIT_DAEMON_MODE", "batch")' in daemon_source
    assert "Never launch or attach to the Cubit GUI" in server._SERVER_INSTRUCTIONS


def test_import_journal_reads_without_starting_or_executing(tmp_path, monkeypatch):
    from types import SimpleNamespace

    journal = tmp_path / "human.jou"
    content = '# saved\n#{size=3}\ncreate brick x {size}\ncreate brick x {size}\n'
    journal.write_bytes(content.encode("utf-8-sig"))
    ai_journal = '#{size=3}\ncreate brick x {size}\n'
    existing = SimpleNamespace(
        _lock=threading.Lock(),
        native_journal_snapshot=lambda: {
            "journal": ai_journal,
            "paths": [str(tmp_path / "ai.jou")],
            "generation_count": 1,
            "errors": [],
            "recording_error": None,
        },
    )
    monkeypatch.setattr(session, "_SINGLETON", existing)
    monkeypatch.setattr(server, "_cubit_session_or_error",
                        lambda: pytest.fail("Import must not start or attach Cubit"))

    result = json.loads(server.cubit_import_journal(str(journal)))

    assert result["executed"] is False
    assert result["gui_started"] is False
    assert result["journal"] == content
    assert [row["line_number"] for row in result["candidate_commands"]] == [4]
    assert result["excluded"][1]["reason"] == "matches_cubit_recorded_ai_journal"
    assert result["ai_journal"] == ai_journal
    assert result["attribution"] == "cubit_native_record_exact_match"
    assert len(result["sha256"]) == 64


@pytest.mark.parametrize("platform", ["win32", "linux", "darwin"])
def test_native_record_uses_supported_file_syntax_and_preserves_aprepro(
        tmp_path, monkeypatch, platform):
    monkeypatch.setenv("RADIA_MCP_TEMP", str(tmp_path))
    monkeypatch.setattr(session.sys, "platform", platform)
    sess = session.CubitSession.__new__(session.CubitSession)
    sess._mode = "batch"
    sess._client_id = "test-client"
    sess._next_id = 1
    sess._native_journal_path = None
    sess._native_journal_paths = []
    sess._native_journal_error = None
    requests = []

    def call(request, timeout_s):
        requests.append(request)
        return {"ok": True, "result": [
            {"line": request["args"][0], "ok": True, "rc": 1},
        ]}

    monkeypatch.setattr(sess, "_call_via_stdio", call)
    sess._start_native_journal_locked(timeout_s=5)
    assert sess._native_journal_path.is_relative_to(tmp_path)

    command = requests[0]["args"][0]
    assert command.startswith('record "')
    assert "record journal" not in command
    assert "overwrite" not in command
    sess._native_journal_path.write_text(
        "#{mesh_size = 0.25}\ncreate brick x {mesh_size}\n",
        encoding="utf-8",
    )
    snapshot = sess.native_journal_snapshot()
    assert "#{mesh_size = 0.25}" in snapshot["journal"]
    assert server._normalized_cubit_journal_commands(snapshot["journal"]) == [
        "#{mesh_size = 0.25}", "create brick x {mesh_size}",
    ]


@pytest.mark.parametrize("platform", ["win32", "linux", "darwin"])
def test_cubit_temp_root_platform_default(tmp_path, monkeypatch, platform):
    monkeypatch.delenv("RADIA_MCP_TEMP", raising=False)
    monkeypatch.setattr(session.sys, "platform", platform)
    monkeypatch.setattr(session.tempfile, "gettempdir", lambda: str(tmp_path))
    expected = Path("C:/temp") if platform == "win32" else tmp_path
    assert session._cubit_temp_root() == expected


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
