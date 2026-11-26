"""Guards: export claim, staging, worker status and MCP error mapping under injected faults.

Assumes: no Cubit, no Netgen and no worker run; ``subprocess.run`` is replaced
by a mock worker and filesystem calls are monkeypatched to fail at one step.
Each test checks that the primary error survives, this run's lock and stage
are removed, and the MCP response names the right error kind.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from radia_mcp.cubit import export_runner as runner
from radia_mcp.cubit import server
from radia_mcp.cubit import vol_export

GOOD = {"worker": "vol_export.py", "netgen": "6.2", "ngsolve": "6.2"}


def _setup(tmp_path):
    (tmp_path / "bin" / "cubit").mkdir(parents=True)
    (tmp_path / "plugin").mkdir()
    (tmp_path / "plugin" / vol_export.PLUGIN_DLL).write_bytes(b"")
    journal = tmp_path / "m.jou"
    journal.write_text("reset\n", encoding="utf-8")
    return journal, {"cubit_bin": tmp_path / "bin", "python": sys.executable, "plugin_dir": tmp_path / "plugin"}


def _worker(monkeypatch, write_status=None, raise_exc=None):
    """Probe answers; the export call writes ``write_status`` text or raises ``raise_exc``."""
    def run(command, **kwargs):
        if "--result" not in command:
            return SimpleNamespace(returncode=0, stdout=(json.dumps(GOOD) + "\n").encode(), stderr=b"")
        if raise_exc is not None:
            raise raise_exc
        if write_status is not None:
            target = Path(command[command.index("--result") + 1])
            if isinstance(write_status, bytes):
                target.write_bytes(write_status)
            else:
                target.write_text(write_status, encoding="utf-8")
        return SimpleNamespace(returncode=1, stdout=b"", stderr=b"")
    monkeypatch.setattr(runner.subprocess, "run", run)


def _leftovers(tmp_path):
    return sorted(p.name for p in tmp_path.iterdir() if p.name.endswith((".lock", ".stage")))


def _tool(journal, out, kw):
    return json.loads(server.cubit_export_vol(str(journal), str(out), 2, cubit_bin=str(kw["cubit_bin"]),
                                              worker_python=kw["python"], plugin_dir=str(kw["plugin_dir"])))


def test_claim_close_failure_releases_own_lock(tmp_path, monkeypatch):
    out = tmp_path / "m.vol"
    real_close = os.close

    def bad_close(fd):
        real_close(fd)
        raise OSError("injected close failure")
    monkeypatch.setattr(runner.os, "close", bad_close)
    with pytest.raises(OSError, match="injected close failure"):
        runner._claim(out)
    assert not Path(str(out) + ".lock").exists(), "own lock left after a failed close"


def test_claim_close_failure_reports_lock_when_cleanup_fails(tmp_path, monkeypatch):
    out = tmp_path / "m.vol"
    real_close = os.close
    monkeypatch.setattr(runner.os, "close", lambda fd: (real_close(fd), (_ for _ in ()).throw(OSError("close"))))
    monkeypatch.setattr(Path, "unlink", lambda self, missing_ok=False: (_ for _ in ()).throw(PermissionError("deny")))
    with pytest.raises(OSError, match="close") as info:
        runner._claim(out)
    notes = getattr(info.value, "__notes__", [])
    assert any(str(out) + ".lock" in n for n in notes), notes


def test_claim_recheck_error_releases_own_lock(tmp_path, monkeypatch):
    out = tmp_path / "m.vol"
    real_exists = Path.exists

    def exists(self):
        if self == out:
            raise PermissionError("injected exists failure")
        return real_exists(self)
    monkeypatch.setattr(Path, "exists", exists)
    with pytest.raises(PermissionError, match="injected exists failure"):
        runner._claim(out)
    assert not real_exists(Path(str(out) + ".lock")), "own lock left after a failed recheck"


def test_claim_collision_keeps_file_exists_when_unlink_denied(tmp_path, monkeypatch):
    out = tmp_path / "m.vol"
    out.write_bytes(b"other")
    monkeypatch.setattr(Path, "unlink", lambda self, missing_ok=False: (_ for _ in ()).throw(PermissionError("deny")))
    with pytest.raises(FileExistsError, match="refusing to overwrite") as info:
        runner._claim(out)
    notes = getattr(info.value, "__notes__", [])
    assert any(str(out) + ".lock" in n and "PermissionError" in n for n in notes), notes


def test_preparation_failure_after_stage_mkdir_removes_stage(tmp_path, monkeypatch):
    journal, kw = _setup(tmp_path)
    _worker(monkeypatch, write_status="{}")
    real_dumps = json.dumps

    def dumps(obj, *a, **k):
        if isinstance(obj, dict) and "periodic_name" in obj:
            raise TypeError("injected serialization failure")
        return real_dumps(obj, *a, **k)
    monkeypatch.setattr(runner.json, "dumps", dumps)
    periodic = {"periodic_name": "kt", "periodic_master": "a", "periodic_slave": "b", "periodic_unit": "m",
                "periodic_patches": "6:7", "periodic_tx": 1, "periodic_ty": 0, "periodic_tz": 0,
                "periodic_match": "1e-6"}
    with pytest.raises(TypeError, match="injected serialization"):
        runner.run_export(journal, tmp_path / "m.vol", 2, periodic=periodic, **kw)
    assert _leftovers(tmp_path) == [], _leftovers(tmp_path)


def test_worker_timeout_cleans_up_and_is_reported_as_timeout(tmp_path, monkeypatch):
    journal, kw = _setup(tmp_path)
    _worker(monkeypatch, raise_exc=subprocess.TimeoutExpired(["worker"], 5))
    with pytest.raises(subprocess.TimeoutExpired):
        runner.run_export(journal, tmp_path / "m.vol", 2, **kw)
    assert _leftovers(tmp_path) == []
    out = _tool(journal, tmp_path / "m.vol", kw)
    assert out["status"] == "error" and out["kind"] == "worker_timeout", out
    assert _leftovers(tmp_path) == []


@pytest.mark.parametrize("text, fragment", [
    ("{not json", "malformed"),
    ("[1, 2]", "not a JSON object"),
    (b'{"status": "ok' + bytes([0xFF]) + b'"}', "malformed"),
    ('{"token": "x"}', "expected 'ok' or 'error'"),
    ('{"status": "maybe"}', "expected 'ok' or 'error'"),
])
def test_bad_worker_status_is_a_protocol_error(tmp_path, monkeypatch, text, fragment):
    journal, kw = _setup(tmp_path)
    _worker(monkeypatch, write_status=text)
    with pytest.raises(runner.WorkerProtocolError, match=fragment):
        runner.run_export(journal, tmp_path / "m.vol", 2, **kw)
    assert _leftovers(tmp_path) == []
    out = _tool(journal, tmp_path / "m.vol", kw)
    assert out["kind"] == "worker_protocol" and fragment in out["error"], out
    assert _leftovers(tmp_path) == []


def test_lock_already_absent_at_release_is_not_a_warning(tmp_path):
    lock = tmp_path / "m.vol.lock"
    assert runner._release(lock) is None


@pytest.mark.parametrize("make, kind", [
    (lambda out: out.write_bytes(b"x"), "output_exists"),
    (lambda out: Path(str(out) + ".lock").write_bytes(b""), "output_exists"),
])
def test_tool_labels_existing_output_and_held_lock(tmp_path, make, kind):
    journal, kw = _setup(tmp_path)
    make(tmp_path / "m.vol")
    out = _tool(journal, tmp_path / "m.vol", kw)
    assert out["kind"] == kind and out["error_type"] == "FileExistsError", out


def test_tool_labels_filesystem_errors_and_keeps_notes(tmp_path, monkeypatch):
    journal, kw = _setup(tmp_path)
    _worker(monkeypatch, write_status="{}")
    real_mkdir = Path.mkdir
    real_unlink = Path.unlink

    def mkdir(self, *a, **k):
        if self.name.endswith(".stage"):
            raise PermissionError(13, "injected stage denial", str(self))
        return real_mkdir(self, *a, **k)

    def unlink(self, missing_ok=False):
        if self.name.endswith(".lock"):
            raise PermissionError("injected lock denial")
        return real_unlink(self, missing_ok=missing_ok)
    monkeypatch.setattr(Path, "mkdir", mkdir)
    monkeypatch.setattr(Path, "unlink", unlink)
    out = _tool(journal, tmp_path / "m.vol", kw)
    assert out["kind"] == "filesystem" and out["kind"] != "input", out
    assert "injected stage denial" in out["error"], out
    assert any("m.vol.lock" in n for n in out.get("notes", [])), out


def test_tool_labels_worker_output_mismatch_as_protocol(tmp_path, monkeypatch):
    journal, kw = _setup(tmp_path)
    _worker(monkeypatch, write_status=json.dumps({"status": "ok", "token": "stale"}))
    out = _tool(journal, tmp_path / "m.vol", kw)
    assert out["error_type"] == "WorkerOutputMismatch" and out["kind"] == "worker_protocol", out


def test_corrupt_companion_with_matching_hash_is_a_protocol_error(tmp_path, monkeypatch):
    import hashlib
    journal, kw = _setup(tmp_path)
    bad = b"{broken"

    def run(command, **kwargs):
        if "--result" not in command:
            return SimpleNamespace(returncode=0, stdout=(json.dumps(GOOD) + "\n").encode(), stderr=b"")
        vol =command[command.index("--out") + 1]
        Path(vol).write_bytes(b"mesh")
        Path(vol + ".json").write_bytes(bad)
        record = {"status": "ok", "token": command[command.index("--token") + 1],
                  "sha256": hashlib.sha256(b"mesh").hexdigest(), "companion": {},
                  "companion_sha256": hashlib.sha256(bad).hexdigest()}
        Path(command[command.index("--result") + 1]).write_text(json.dumps(record), encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout=b"", stderr=b"")

    monkeypatch.setattr(runner.subprocess, "run", run)
    out = _tool(journal, tmp_path / "m.vol", kw)
    assert out["status"] == "error" and out["kind"] == "worker_protocol", out
    assert "not valid JSON" in out["error"], out
    assert not (tmp_path / "m.vol").exists() and not (tmp_path / "m.vol.json").exists()
    assert _leftovers(tmp_path) == []
