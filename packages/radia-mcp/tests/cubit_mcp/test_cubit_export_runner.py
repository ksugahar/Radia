"""Guards: the MCP export route starts the native-plugin worker isolated and fails loudly.

Assumes: no Cubit, no Netgen and no worker run; ``subprocess.run`` is replaced
by a recorder, so these tests check the command line, the runtime and plugin
verdicts, the periodic/Kelvin argument contract and the input refusals only.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path
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


def _fake_run(probe, calls, status=None, returncode=None, before=None, companion=b"{}"):
    """Probe answers, then a mock worker; an ok status also writes the .vol pair.

    ``before`` runs inside the worker call (to synchronize concurrent runs);
    ``companion`` is the companion file written, while the status always
    reports the record ``{}`` and its sha256.
    """
    def run(command, **kwargs):
        calls.append((command, kwargs))
        if before is not None and "--result" in command:
            before()
        if "--result" not in command:
            return SimpleNamespace(returncode=0, stdout=(json.dumps(probe) + "\n").encode(), stderr=b"")
        result = command[command.index("--result") + 1]
        if status is not None:
            record = dict(status, token=command[command.index("--token") + 1])
            if status["status"] == "ok":
                vol = command[command.index("--out") + 1]
                with open(vol, "wb") as handle:
                    handle.write(b"mesh")
                with open(vol + ".json", "wb") as handle:
                    handle.write(companion)
                record.update(sha256=hashlib.sha256(b"mesh").hexdigest(), companion={},
                              companion_sha256=hashlib.sha256(b"{}").hexdigest())
            with open(result, "w", encoding="utf-8") as handle:
                handle.write(json.dumps(record))
        code = returncode if returncode is not None else (0 if status and status["status"] == "ok" else 1)
        return SimpleNamespace(returncode=code, stdout=b"", stderr=b"")
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
    with pytest.raises(runner.WorkerRuntimeError, match="radia_cubit_export.dll not found"):
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


def test_stale_status_and_silent_worker_are_not_success(tmp_path, monkeypatch):
    out = tmp_path / "m.vol"
    (tmp_path / "m.vol.status.json").write_text(json.dumps({"status": "ok"}), encoding="utf-8")
    monkeypatch.setattr(runner.subprocess, "run", _fake_run(GOOD, [], None, returncode=0))
    result = runner.run_export(_journal(tmp_path), out, 2, cubit_bin=_cubit_bin(tmp_path),
                               python=sys.executable, plugin_dir=_plugin(tmp_path))
    assert result["status"] == "error" and result["error_type"] == "WorkerCrash", result


def test_ok_status_without_outputs_is_an_error(tmp_path, monkeypatch):
    def run(command, **kwargs):
        if "--result" not in command:
            return SimpleNamespace(returncode=0, stdout=(json.dumps(GOOD) + "\n").encode(), stderr=b"")
        with open(command[command.index("--result") + 1], "w", encoding="utf-8") as handle:
            handle.write(json.dumps({"status": "ok", "sha256": "0" * 64,
                                     "token": command[command.index("--token") + 1]}))
        return SimpleNamespace(returncode=0, stdout=b"", stderr=b"")
    monkeypatch.setattr(runner.subprocess, "run", run)
    result = runner.run_export(_journal(tmp_path), tmp_path / "m.vol", 2, cubit_bin=_cubit_bin(tmp_path),
                               python=sys.executable, plugin_dir=_plugin(tmp_path))
    assert result["status"] == "error" and result["error_type"] == "WorkerOutputMismatch", result


def test_existing_companion_is_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(runner.subprocess, "run", _fake_run(GOOD, []))
    (tmp_path / "m.vol.json").write_text("{}", encoding="utf-8")
    with pytest.raises(FileExistsError):
        runner.run_export(_journal(tmp_path), tmp_path / "m.vol", 2, cubit_bin=_cubit_bin(tmp_path),
                          python=sys.executable, plugin_dir=_plugin(tmp_path))


def test_altered_companion_is_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(runner.subprocess, "run",
                        _fake_run(GOOD, [], {"status": "ok"}, companion=b'{"x": 1}'))
    out = tmp_path / "m.vol"
    result = runner.run_export(_journal(tmp_path), out, 2, cubit_bin=_cubit_bin(tmp_path),
                               python=sys.executable, plugin_dir=_plugin(tmp_path))
    assert result["status"] == "error" and "companion_sha256" in result["error"], result
    assert not out.exists() and not (tmp_path / "m.vol.json").exists()


def test_concurrent_exports_to_one_path_publish_once(tmp_path, monkeypatch):
    import threading
    inside, release = threading.Event(), threading.Event()

    def hold():
        if not inside.is_set():
            inside.set()
            assert release.wait(30)

    monkeypatch.setattr(runner.subprocess, "run", _fake_run(GOOD, [], {"status": "ok"}, before=hold))
    out, journal = tmp_path / "m.vol", _journal(tmp_path)
    kwargs = dict(cubit_bin=_cubit_bin(tmp_path), python=sys.executable, plugin_dir=_plugin(tmp_path))
    results = {}
    first = threading.Thread(target=lambda: results.setdefault(
        "first", runner.run_export(journal, out, 2, **kwargs)))
    first.start()
    assert inside.wait(30)
    try:
        with pytest.raises(FileExistsError, match="another export"):
            runner.run_export(journal, out, 2, **kwargs)
    finally:
        release.set()
        first.join(30)
    assert results["first"]["status"] == "ok", results
    assert out.read_bytes() == b"mesh"
    assert not list(tmp_path.glob(".m.vol.*.stage")) and not (tmp_path / "m.vol.lock").exists()


def test_worker_cwd_is_destination_parent_with_absolute_staged_outputs(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(runner.subprocess, "run", _fake_run(GOOD, calls, {"status": "ok"}))
    out = tmp_path / "sub" / "m.vol"
    out.parent.mkdir()
    result = runner.run_export(_journal(tmp_path), out, 2, cubit_bin=_cubit_bin(tmp_path),
                               python=sys.executable, plugin_dir=_plugin(tmp_path))
    assert result["status"] == "ok", result
    command, kwargs = calls[1]
    assert kwargs["cwd"] == str(out.parent), (kwargs["cwd"], str(out.parent))
    for flag in ("--out", "--result", "--journal"):
        value = command[command.index(flag) + 1]
        assert Path(value).is_absolute(), (flag, value)


def test_companion_planted_before_publish_rolls_back_and_retry_succeeds(tmp_path, monkeypatch):
    out = tmp_path / "m.vol"
    planted = tmp_path / "m.vol.json"

    def plant():
        planted.write_bytes(b'{"other": 1}')

    monkeypatch.setattr(runner.subprocess, "run", _fake_run(GOOD, [], {"status": "ok"}, before=plant))
    kwargs = dict(cubit_bin=_cubit_bin(tmp_path), python=sys.executable, plugin_dir=_plugin(tmp_path))
    with pytest.raises(FileExistsError):
        runner.run_export(_journal(tmp_path), out, 2, **kwargs)
    assert not out.exists(), "the .vol this run linked must be rolled back"
    assert planted.read_bytes() == b'{"other": 1}', "a file this run did not publish must stay"
    assert not list(tmp_path.glob(".m.vol.*.stage")) and not (tmp_path / "m.vol.lock").exists()
    planted.unlink()
    monkeypatch.setattr(runner.subprocess, "run", _fake_run(GOOD, [], {"status": "ok"}))
    result = runner.run_export(_journal(tmp_path), out, 2, **kwargs)
    assert result["status"] == "ok", result
    assert out.read_bytes() == b"mesh" and planted.read_bytes() == b"{}"


def test_companion_replaced_during_verification_keeps_it_and_removes_own_vol(tmp_path, monkeypatch):
    out = tmp_path / "m.vol"
    companion = tmp_path / "m.vol.json"
    original = runner._verify_outputs

    def swap_then_verify(status, token, vol, comp, returncode):
        if vol == out:  # post-publication verification of the real destination
            comp.unlink()
            comp.write_bytes(b'{"foreign": 1}')
        return original(status, token, vol, comp, returncode)

    monkeypatch.setattr(runner, "_verify_outputs", swap_then_verify)
    monkeypatch.setattr(runner.subprocess, "run", _fake_run(GOOD, [], {"status": "ok"}))
    result = runner.run_export(_journal(tmp_path), out, 2, cubit_bin=_cubit_bin(tmp_path),
                               python=sys.executable, plugin_dir=_plugin(tmp_path))
    assert result["status"] == "error", result
    assert "no longer this run's output" in result["error"], result
    assert not out.exists(), "this run's .vol must still be rolled back"
    assert companion.read_bytes() == b'{"foreign": 1}', "the foreign companion must stay"
    assert not list(tmp_path.glob(".m.vol.*.stage")) and not (tmp_path / "m.vol.lock").exists()


def test_worker_error_survives_a_failing_stage_cleanup(tmp_path, monkeypatch):
    class WorkerBoom(RuntimeError):
        pass

    def run(command, **kwargs):
        if "--result" not in command:
            return SimpleNamespace(returncode=0, stdout=(json.dumps(GOOD) + "\n").encode(), stderr=b"")
        stage = Path(command[command.index("--out") + 1]).parent
        (stage / "nested").mkdir()
        (stage / "nested" / "f").write_bytes(b"x")
        raise WorkerBoom("worker failed")

    def broken_rmtree(path, *a, **k):
        raise PermissionError(f"cannot remove {path}")

    monkeypatch.setattr(runner.subprocess, "run", run)
    monkeypatch.setattr(runner.shutil, "rmtree", broken_rmtree)
    with pytest.raises(WorkerBoom, match="worker failed") as info:
        runner.run_export(_journal(tmp_path), tmp_path / "m.vol", 2, cubit_bin=_cubit_bin(tmp_path),
                          python=sys.executable, plugin_dir=_plugin(tmp_path))
    notes = getattr(info.value, "__notes__", [])
    assert any("staging cleanup failed" in n for n in notes), notes
    assert not (tmp_path / "m.vol.lock").exists()


def test_worker_error_with_stage_subdirectory_cleans_up(tmp_path, monkeypatch):
    def run(command, **kwargs):
        if "--result" not in command:
            return SimpleNamespace(returncode=0, stdout=(json.dumps(GOOD) + "\n").encode(), stderr=b"")
        stage = Path(command[command.index("--out") + 1]).parent
        (stage / "nested").mkdir()
        raise OSError("worker launch failed")

    monkeypatch.setattr(runner.subprocess, "run", run)
    with pytest.raises(OSError, match="worker launch failed"):
        runner.run_export(_journal(tmp_path), tmp_path / "m.vol", 2, cubit_bin=_cubit_bin(tmp_path),
                          python=sys.executable, plugin_dir=_plugin(tmp_path))
    assert not list(tmp_path.glob(".m.vol.*.stage")) and not (tmp_path / "m.vol.lock").exists()


def test_rollback_continues_when_one_unlink_fails(tmp_path):
    staged_vol, staged_json = tmp_path / "s.vol", tmp_path / "s.json"
    staged_vol.write_bytes(b"mesh")
    staged_json.write_bytes(b"{}")
    dest_vol, dest_json = tmp_path / "m.vol", tmp_path / "m.vol.json"
    os.link(staged_vol, dest_vol)
    os.link(staged_json, dest_json)
    real_unlink = Path.unlink

    def failing_unlink(self, *args, **kwargs):
        if self == dest_json:
            raise PermissionError("locked by another process")
        return real_unlink(self, *args, **kwargs)

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(Path, "unlink", failing_unlink)
        with pytest.raises(OSError, match="rollback incomplete") as info:
            runner._rollback([(staged_vol, dest_vol), (staged_json, dest_json)])
    assert "locked by another process" in str(info.value), info.value
    assert not dest_vol.exists(), "the .vol must still be removed after the companion unlink fails"
    assert dest_json.exists()


def test_lock_release_failure_after_publish_keeps_the_result(tmp_path, monkeypatch):
    out = tmp_path / "m.vol"
    lock = tmp_path / "m.vol.lock"
    real_unlink = Path.unlink

    def failing_unlink(self, *args, **kwargs):
        if self == lock:
            raise PermissionError("lock held by a scanner")
        return real_unlink(self, *args, **kwargs)

    monkeypatch.setattr(runner.subprocess, "run", _fake_run(GOOD, [], {"status": "ok"}))
    monkeypatch.setattr(Path, "unlink", failing_unlink)
    result = runner.run_export(_journal(tmp_path), out, 2, cubit_bin=_cubit_bin(tmp_path),
                               python=sys.executable, plugin_dir=_plugin(tmp_path))
    assert result["status"] == "ok", result
    assert "lock held by a scanner" in result["lock_release_error"], result
    assert out.read_bytes() == b"mesh" and (tmp_path / "m.vol.json").exists()
