"""Guards: the MCP export route starts the native-plugin worker isolated and fails loudly.

Assumes: no Cubit, no Netgen and no worker run; ``subprocess.run`` is replaced
by a recorder, so these tests check the command line, the runtime and plugin
verdicts, the periodic/Kelvin argument contract and the input refusals only.
"""
from __future__ import annotations

import json
import sys
from types import SimpleNamespace

import pytest

from radia_mcp.cubit import export_runner as runner
from radia_mcp.cubit import kelvin
from radia_mcp.cubit import server
from radia_mcp.cubit import vol_export

PERIODIC = {"periodic_name": "kt", "periodic_master": "source", "periodic_slave": "target",
            "periodic_unit": "m", "periodic_patches": "6:7", "periodic_tx": 1, "periodic_ty": 0,
            "periodic_tz": 0, "periodic_match": "1e-6"}
GOOD = {"worker": "vol_export.py", "netgen": "6.2", "ngsolve": "6.2"}


def _cubit_bin(tmp_path):
    (tmp_path / "bin" / "cubit").mkdir(parents=True, exist_ok=True)
    return tmp_path / "bin"


def _plugin(tmp_path):
    (tmp_path / "plugin").mkdir(exist_ok=True)
    (tmp_path / "plugin" / vol_export.PLUGIN_DLL).write_bytes(b"")
    return tmp_path / "plugin"


def _journal(tmp_path):
    journal = tmp_path / "m.jou"
    journal.write_text("reset\n", encoding="utf-8")
    return journal


def _fake_run(probe, calls, status=None):
    def run(command, **kwargs):
        calls.append((command, kwargs))
        if "--result" not in command:
            return SimpleNamespace(returncode=0, stdout=(json.dumps(probe) + "\n").encode(), stderr=b"")
        out = command[command.index("--result") + 1]
        if status is not None:
            with open(out, "w", encoding="utf-8") as handle:
                handle.write(json.dumps(status))
        return SimpleNamespace(returncode=0 if status and status["status"] == "ok" else 1,
                               stdout=b"", stderr=b"")
    return run


def test_worker_runs_isolated_with_plugin_periodic_and_kelvin(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(runner.subprocess, "run", _fake_run(GOOD, calls, {"status": "ok"}))
    result = runner.run_export(_journal(tmp_path), tmp_path / "m.vol", 3, cubit_bin=_cubit_bin(tmp_path),
                               python=sys.executable, plugin_dir=_plugin(tmp_path), periodic=PERIODIC,
                               kelvin={"add_kelvin": True, "kelvin_sym_z": "bn"})
    assert result["status"] == "ok", result
    probe_cmd, export_cmd = calls[0][0], calls[1][0]
    assert probe_cmd[:3] == [sys.executable, "-I", "-B"], probe_cmd
    assert export_cmd[:4] == [sys.executable, "-I", "-B", "-c"], export_cmd
    assert "radia_mcp.cubit.vol_export" in export_cmd[4]
    assert export_cmd[export_cmd.index("--order") + 1] == "3"
    assert export_cmd[export_cmd.index("--plugin-dir") + 1] == str(tmp_path / "plugin")
    assert json.loads(export_cmd[export_cmd.index("--periodic") + 1]) == PERIODIC
    assert json.loads(export_cmd[export_cmd.index("--kelvin") + 1])["kelvin_sym_z"] == "bn"
    assert calls[1][1]["stdin"] is runner.subprocess.DEVNULL


def test_export_command_matches_the_plugin_syntax(tmp_path):
    command = vol_export.export_command(tmp_path / "a.vol", 3, PERIODIC)
    assert command.startswith(f'export netgen "{(tmp_path / "a.vol").as_posix()}" order 3 ')
    assert ('periodic_name "kt" periodic_master "source" periodic_slave "target" periodic_unit "m" '
            'periodic_patches "6:7" periodic_tx 1 periodic_ty 0 periodic_tz 0 periodic_match 1e-6 overwrite') in command
    assert vol_export.export_command(tmp_path / "a.vol", 2, None).endswith(" order 2 overwrite")


@pytest.mark.parametrize("periodic, message", [
    ({"periodic_name": "kt"}, "all nine"),
    ({**PERIODIC, "periodic_extra": 1}, "unknown periodic"),
    ({**PERIODIC, "periodic_name": 'k"t'}, "quotes"),
])
def test_periodic_contract_rejects(periodic, message):
    with pytest.raises(ValueError, match=message):
        vol_export.periodic_arguments(periodic)


@pytest.mark.parametrize("options, message", [
    ({"add_kelvin": True, "kelvin_sym_x": "on"}, "off\\|bn\\|ht"),
    ({"add_kelvin": True, "kelvin_mesh": -1}, "positive length"),
    ({"add_kelvin": True, "kelvin_radius": 1}, "unknown Kelvin"),
])
def test_kelvin_contract_rejects(options, message):
    with pytest.raises(ValueError, match=message):
        kelvin.validate_options(options)


def test_kelvin_symmetry_maps_to_reduction():
    assert kelvin.reduction_from_sym() is None
    assert kelvin.reduction_from_sym("ht", "off", "bn") == {"x": "ht=0", "z": "bn=0"}
    assert kelvin.validate_options({})["add_kelvin"] is False


def test_missing_plugin_runtime_and_cubit_fail_loudly(tmp_path, monkeypatch):
    journal = _journal(tmp_path)
    calls = []
    monkeypatch.setattr(runner.subprocess, "run", _fake_run(GOOD, calls))
    monkeypatch.delenv("RADIA_CUBIT_PLUGIN_DIR", raising=False)
    with pytest.raises(runner.WorkerRuntimeError, match="plugin directory not set"):
        runner.run_export(journal, tmp_path / "m.vol", 2, cubit_bin=_cubit_bin(tmp_path), python=sys.executable)
    with pytest.raises(runner.WorkerRuntimeError, match="cubit_mesh_export.dll not found"):
        runner.run_export(journal, tmp_path / "m.vol", 2, cubit_bin=_cubit_bin(tmp_path), python=sys.executable,
                          plugin_dir=tmp_path)
    with pytest.raises(runner.WorkerRuntimeError, match="Cubit Python package not found"):
        runner.run_export(journal, tmp_path / "m.vol", 2, cubit_bin=tmp_path, python=sys.executable,
                          plugin_dir=_plugin(tmp_path))
    assert calls == [], "no subprocess may start before the runtime checks pass"
    monkeypatch.setattr(runner.subprocess, "run", _fake_run({"worker": "x", "netgen_error": "ImportError"}, calls))
    with pytest.raises(runner.WorkerRuntimeError, match="ngsolve/netgen"):
        runner.run_export(journal, tmp_path / "m.vol", 2, cubit_bin=_cubit_bin(tmp_path), python=sys.executable,
                          plugin_dir=_plugin(tmp_path))
    assert len(calls) == 1, "no export may start after a failed probe"


def test_bad_inputs_are_refused(tmp_path, monkeypatch):
    journal = _journal(tmp_path)
    monkeypatch.setattr(runner.subprocess, "run", _fake_run(GOOD, []))
    with pytest.raises(runner.WorkerRuntimeError, match="absolute path"):
        runner.resolve_python("python")
    (tmp_path / "old.vol").write_text("x", encoding="utf-8")
    with pytest.raises(FileExistsError):
        runner.run_export(journal, tmp_path / "old.vol", 2, cubit_bin=_cubit_bin(tmp_path))
    with pytest.raises(ValueError, match="order"):
        runner.run_export(journal, tmp_path / "m.vol", 0, cubit_bin=_cubit_bin(tmp_path))


def test_worker_crash_without_status_is_an_error(tmp_path, monkeypatch):
    monkeypatch.setattr(runner.subprocess, "run", _fake_run(GOOD, [], None))
    result = runner.run_export(_journal(tmp_path), tmp_path / "m.vol", 2, cubit_bin=_cubit_bin(tmp_path),
                               python=sys.executable, plugin_dir=_plugin(tmp_path))
    assert result["status"] == "error" and result["error_type"] == "WorkerCrash", result


def test_tool_reports_environment_and_input_errors_as_json(tmp_path, monkeypatch):
    journal = _journal(tmp_path)
    monkeypatch.setattr(runner.subprocess, "run", _fake_run({"worker": None}, []))
    out = json.loads(server.cubit_export_vol(str(journal), str(tmp_path / "m.vol"), 2,
                                             cubit_bin=str(_cubit_bin(tmp_path)), worker_python=sys.executable,
                                             plugin_dir=str(_plugin(tmp_path))))
    assert out["status"] == "error" and out["kind"] == "environment", out
    out = json.loads(server.cubit_export_vol(str(journal), str(tmp_path / "m.vol"), 2, add_kelvin=True,
                                             kelvin_sym_y="sideways"))
    assert out["status"] == "error" and out["kind"] == "input", out
