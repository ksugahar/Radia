import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

PATH = Path(__file__).resolve().parents[1] / "tools/release_temp_shadows.py"
SPEC = importlib.util.spec_from_file_location("release_temp_shadows_test", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_exact_hosts_and_read_only_default(monkeypatch):
    calls = []
    def run(command, **kwargs):
        calls.append(command)
        assert "$apply = $false" in kwargs["input"]
        return SimpleNamespace(stdout=json.dumps({"passed": True}))
    monkeypatch.setattr(MODULE.subprocess, "run", run)
    assert all(x["passed"] for x in MODULE.inspect_shadows().values())
    assert [cmd[5] for cmd in calls] == ["mdx1", "mdx2", "hibino"]


def test_unreachable_host_is_not_reported_clean(monkeypatch):
    def run(*args, **kwargs):
        raise OSError("unreachable")
    monkeypatch.setattr(MODULE.subprocess, "run", run)
    assert all(not x["passed"] for x in MODULE.inspect_shadows(True).values())


def test_cleanup_has_literal_target_and_process_link_guards():
    assert "Remove-Item -LiteralPath 'C:\\temp\\radia-omega-test'" in MODULE.SCRIPT
    assert "ReparsePoint" in MODULE.SCRIPT
    assert "Get-CimInstance Win32_Process" in MODULE.SCRIPT
    assert "PYTHONPATH" in MODULE.SCRIPT
    assert "Stop-Process" not in MODULE.SCRIPT
    assert "Cleanup parent is a reparse point" in MODULE.SCRIPT


def test_missing_remote_report_fails_loudly(monkeypatch):
    monkeypatch.setattr(MODULE.subprocess, "run", lambda *a, **k: SimpleNamespace(stdout=""))
    results = MODULE.inspect_shadows()
    assert all(not record["passed"] for record in results.values())
    assert all("no JSON report" in record["error"] for record in results.values())


def test_active_host_remains_a_blocker(monkeypatch):
    report = {"passed": False, "exists": True, "blockers": ["active process"]}
    monkeypatch.setattr(MODULE.subprocess, "run", lambda *a, **k: SimpleNamespace(stdout=json.dumps(report)))
    assert all(record == report for record in MODULE.inspect_shadows(True).values())
