"""Tests for the persistent gmsh session (matlab-mcp-core-server style).

Needs the gmsh Python package; each test leaves no session behind.
"""

import importlib.util

import pytest

from radia_mcp.gmsh.session import (
    GmshSession,
    GmshSessionError,
    session_exec,
    session_shutdown,
    session_status,
)

_GMSH_AVAILABLE = importlib.util.find_spec("gmsh") is not None

pytestmark = pytest.mark.skipif(not _GMSH_AVAILABLE,
                                reason="gmsh package not installed")

_TINY_MSH = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
3 1 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
1 1 1 1
3 1 4 1
1 1 2 3 4
$EndElements
"""


@pytest.fixture(autouse=True)
def _clean_session():
    yield
    session_shutdown()


def test_exec_state_persists_across_calls():
    out = session_exec("a = 21\nprint('hello from gmsh session')")
    assert out["ok"] is True
    assert "hello from gmsh session" in out["stdout"]
    assert "result" not in out

    out2 = session_exec("result = a * 2")
    assert out2["ok"] is True
    assert out2["result"] == 42
    assert out2["session_pid"] == out["session_pid"]  # same worker


def test_exec_opens_model_and_reports_status(tmp_path):
    msh = tmp_path / "tiny.msh"
    msh.write_text(_TINY_MSH, encoding="utf-8")
    path_literal = repr(str(msh))

    out = session_exec(f"gmsh.open({path_literal})\n"
                       "result = len(gmsh.model.mesh.getNodes()[0])")
    assert out["ok"] is True
    assert out["result"] == 4

    status = session_status()
    assert status["running"] is True
    assert status["gmsh_version"]
    assert status["n_calls"] >= 1

    stopped = session_shutdown()
    assert stopped["running"] is False
    assert session_status()["running"] is False


def test_exec_python_exception_returns_traceback():
    out = session_exec("1/0")
    assert out["ok"] is False
    assert "ZeroDivisionError" in out["error"]

    # the worker survived the exception
    out2 = session_exec("result = 'still alive'")
    assert out2["ok"] is True
    assert out2["result"] == "still alive"


def test_worker_crash_fails_loud_then_next_call_starts_fresh():
    pid1 = session_exec("result = 1")["session_pid"]
    with pytest.raises(GmshSessionError, match="died"):
        session_exec("import os; os._exit(3)", timeout_s=30.0)

    out = session_exec("result = 'fresh'")
    assert out["ok"] is True
    assert out["result"] == "fresh"
    assert out["session_pid"] != pid1


def test_timeout_kills_session_and_raises():
    with pytest.raises(GmshSessionError, match="did not answer"):
        session_exec("import time; time.sleep(30)", timeout_s=1.5)
    session_obj = GmshSession.peek()
    assert session_obj is None or not session_obj.alive()

    out = session_exec("result = 2 + 2", timeout_s=60.0)
    assert out["result"] == 4


def test_nonserializable_result_falls_back_to_repr():
    out = session_exec("result = object()")
    assert out["ok"] is True
    assert out.get("result_repr") is True
    assert "object object" in out["result"]
