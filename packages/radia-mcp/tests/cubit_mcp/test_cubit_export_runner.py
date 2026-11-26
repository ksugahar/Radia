"""Guards: the MCP export route starts the curving worker isolated and fails loudly.

Assumes: no Cubit and no worker run; ``subprocess.run`` is replaced by a
recorder, so these tests check the command line, the runtime probe verdicts
and the input refusals only.
"""
from __future__ import annotations

import json
import sys
from types import SimpleNamespace

import pytest

from radia_mcp.cubit import export_runner as runner
from radia_mcp.cubit import server


def _cubit_bin(tmp_path):
    (tmp_path / "bin" / "cubit").mkdir(parents=True, exist_ok=True)
    return tmp_path / "bin"


def _fake_run(probe, calls, status=None):
    def run(command, **kwargs):
        calls.append((command, kwargs))
        if "-c" in command:
            return SimpleNamespace(returncode=0, stdout=(json.dumps(probe) + "\n").encode(), stderr=b"")
        out = command[command.index("--result") + 1]
        if status is not None:
            with open(out, "w", encoding="utf-8") as handle:
                handle.write(json.dumps(status))
        return SimpleNamespace(returncode=0 if status and status["status"] == "ok" else 1,
                               stdout=b"", stderr=b"")
    return run


GOOD = {"worker": "vol_export.py", "netgen": "6.2", "geometric_refit": True}


def test_worker_runs_isolated_with_explicit_interpreter(tmp_path, monkeypatch):
    journal = tmp_path / "m.jou"
    journal.write_text("reset\n", encoding="utf-8")
    calls = []
    monkeypatch.setattr(runner.subprocess, "run", _fake_run(GOOD, calls, {"status": "ok"}))
    result = runner.run_export(journal, tmp_path / "m.vol", 3, cubit_bin=_cubit_bin(tmp_path),
                               python=sys.executable)
    assert result["status"] == "ok", result
    probe_cmd, export_cmd = calls[0][0], calls[1][0]
    assert probe_cmd[:3] == [sys.executable, "-I", "-B"], probe_cmd
    assert export_cmd[:5] == [sys.executable, "-I", "-B", "-m", "radia_mcp.cubit.vol_export"], export_cmd
    assert export_cmd[export_cmd.index("--order") + 1] == "3"
    assert calls[1][1]["stdin"] is runner.subprocess.DEVNULL
    assert not any("export netgen" in str(part) for part in export_cmd)


@pytest.mark.parametrize("probe, message", [
    ({"worker": None, "geometric_refit": True}, "not importable"),
    ({"worker": "x", "geometric_refit": False, "netgen": "6.2.2608"}, "lacks Mesh.GeometricRefit"),
    ({"worker": "x", "netgen_error": "ModuleNotFoundError: netgen"}, "lacks Mesh.GeometricRefit"),
])
def test_missing_worker_runtime_fails_loudly(tmp_path, monkeypatch, probe, message):
    journal = tmp_path / "m.jou"
    journal.write_text("reset\n", encoding="utf-8")
    calls = []
    monkeypatch.setattr(runner.subprocess, "run", _fake_run(probe, calls))
    with pytest.raises(runner.WorkerRuntimeError, match=message):
        runner.run_export(journal, tmp_path / "m.vol", 2, cubit_bin=_cubit_bin(tmp_path),
                          python=sys.executable)
    assert len(calls) == 1, "no export may start after a failed probe"


def test_missing_cubit_and_bad_inputs_are_refused(tmp_path, monkeypatch):
    journal = tmp_path / "m.jou"
    journal.write_text("reset\n", encoding="utf-8")
    monkeypatch.setattr(runner.subprocess, "run", _fake_run(GOOD, []))
    with pytest.raises(runner.WorkerRuntimeError, match="Cubit Python package not found"):
        runner.run_export(journal, tmp_path / "m.vol", 2, cubit_bin=tmp_path, python=sys.executable)
    with pytest.raises(runner.WorkerRuntimeError, match="absolute path"):
        runner.resolve_python("python")
    (tmp_path / "old.vol").write_text("x", encoding="utf-8")
    with pytest.raises(FileExistsError):
        runner.run_export(journal, tmp_path / "old.vol", 2, cubit_bin=_cubit_bin(tmp_path))
    with pytest.raises(ValueError, match="order"):
        runner.run_export(journal, tmp_path / "m.vol", 0, cubit_bin=_cubit_bin(tmp_path))


def test_worker_crash_without_status_is_an_error(tmp_path, monkeypatch):
    journal = tmp_path / "m.jou"
    journal.write_text("reset\n", encoding="utf-8")
    monkeypatch.setattr(runner.subprocess, "run", _fake_run(GOOD, [], None))
    result = runner.run_export(journal, tmp_path / "m.vol", 2, cubit_bin=_cubit_bin(tmp_path),
                               python=sys.executable)
    assert result["status"] == "error" and result["error_type"] == "WorkerCrash", result


def test_tool_reports_environment_errors_as_json(tmp_path, monkeypatch):
    journal = tmp_path / "m.jou"
    journal.write_text("reset\n", encoding="utf-8")
    monkeypatch.setattr(runner.subprocess, "run", _fake_run({"worker": None}, []))
    out = json.loads(server.cubit_export_vol(str(journal), str(tmp_path / "m.vol"), 2,
                                             cubit_bin=str(_cubit_bin(tmp_path)),
                                             worker_python=sys.executable))
    assert out["status"] == "error" and out["kind"] == "environment", out
