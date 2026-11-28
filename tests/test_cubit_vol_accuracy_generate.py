"""The volume-accuracy generator collects only exports of the current run.

Guards: ``generate._run_cubit`` removes the artifacts it is about to write
before Cubit starts, and refuses a nonzero Cubit exit, so a reused
``--work-dir`` cannot hand stale ``.vol`` files to the dataset.
Assumes: no Cubit; the child process is replaced by a stub.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import types

import pytest

GENERATE = (Path(__file__).resolve().parents[1]
            / "validation_test" / "cubit" / "vol_accuracy_dataset" / "generate.py")


def _load():
    spec = importlib.util.spec_from_file_location("vol_accuracy_generate", GENERATE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_reused_work_dir_with_a_failed_cubit_run_is_refused(tmp_path, monkeypatch):
    generate = _load()
    case = {"case_id": "stale_case", "commands": ["reset"]}
    stale = generate._artifact_path(tmp_path, case["case_id"], 1)
    stale.parent.mkdir(parents=True, exist_ok=True)
    for path in (stale, Path(str(stale) + ".json")):
        path.write_text("from an earlier run", encoding="utf-8")

    seen_before_launch = {}

    def fake_run(command, **kwargs):
        seen_before_launch["stale_present"] = stale.exists()
        stale.write_text("written by a failing run", encoding="utf-8")
        Path(str(stale) + ".json").write_text("{}", encoding="utf-8")
        return types.SimpleNamespace(returncode=1, stdout="", stderr="startup failed")

    monkeypatch.setattr(generate, "_running_cubit_processes", lambda: [])
    monkeypatch.setattr(generate, "_find_cubit_exe", lambda: tmp_path / "bin" / "cubitx.exe")
    monkeypatch.setattr(generate, "export_plugin_directory", lambda root: tmp_path)
    monkeypatch.setattr(generate.subprocess, "run", fake_run)

    with pytest.raises(RuntimeError, match="Cubit exited 1"):
        generate._run_cubit(tmp_path, [case], [1])
    assert seen_before_launch == {"stale_present": False}, seen_before_launch
