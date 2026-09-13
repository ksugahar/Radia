from __future__ import annotations

import importlib.util
import subprocess
import zipfile
from pathlib import Path
from types import SimpleNamespace

import pytest


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PACKAGE_ROOT / "tools" / "verify_wheel_contents.py"
SPEC = importlib.util.spec_from_file_location("verify_wheel_contents", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _probe_module():
    spec = importlib.util.spec_from_file_location("wheel_stdio_probe", PACKAGE_ROOT / "tools" / "smoke_mcp_stdio.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_wheel_probe_runs_every_catalog_server_and_selftest(monkeypatch):
    probe = _probe_module()
    monkeypatch.setattr(probe, "CATALOG", {"second": {}, "first": {}})
    calls = []

    def fake_probe(name, **kwargs):
        calls.append(("stdio", name))
        return {"server": name, "n_tools": 1, "protocol_version": "fixture"}

    monkeypatch.setattr(probe, "probe_server", fake_probe)
    monkeypatch.setattr(probe, "run_selftest", lambda name, *a, **kw: calls.append(("selftest", name)))
    assert probe.main(["--all", "--selftest"]) == 0
    assert calls == [("stdio", "first"), ("selftest", "first"),
                     ("stdio", "second"), ("selftest", "second")]


@pytest.mark.parametrize("defect", [None, "editable", "unknown", "reload", "source"])
def test_wheel_probe_requires_install_provenance_and_no_reload(tmp_path, defect):
    probe = _probe_module()
    root = tmp_path / "site-packages"
    result = {"server": "demo", "module_file": str(root / "radia_mcp/server.py"),
              "distribution": {"editable": False}, "reload_tools": []}
    if defect == "editable":
        result["distribution"]["editable"] = True
    elif defect == "unknown":
        result["distribution"] = {}
    elif defect == "reload":
        result["reload_tools"] = ["demo_reload_code"]
    elif defect == "source":
        result["module_file"] = str(tmp_path / "checkout/server.py")
    if defect:
        with pytest.raises(AssertionError):
            probe.check_wheel_result(result, str(root))
    else:
        probe.check_wheel_result(result, str(root))


@pytest.mark.parametrize("failure", [False, True, "timeout"])
def test_wheel_selftest_uses_isolated_launcher_and_propagates_failure(monkeypatch, failure):
    probe = _probe_module()

    def fake_run(command, **kwargs):
        assert command[1:] == ["-I", "-m", "radia_mcp.maintenance", "serve", "meta", "--selftest"]
        assert kwargs["timeout"] == 10
        if failure == "timeout":
            raise subprocess.TimeoutExpired(command, 10)
        return SimpleNamespace(returncode=int(failure), stdout="diagnostic", stderr="")

    monkeypatch.setattr(probe.subprocess, "run", fake_run)
    if failure:
        with pytest.raises((AssertionError, subprocess.TimeoutExpired)):
            probe.run_selftest("meta", 10, installed_wheel=True)
    else:
        probe.run_selftest("meta", 10, installed_wheel=True)


def _write_wheel(path: Path, names: set[str]) -> None:
    with zipfile.ZipFile(path, "w") as archive:
        for name in names:
            archive.writestr(name, "fixture")


def test_wheel_runtime_assets_are_complete(tmp_path):
    wheel = tmp_path / "radia_mcp-test.whl"
    _write_wheel(wheel, set(MODULE.REQUIRED_ASSETS))

    result = MODULE.verify_wheel_contents(wheel)

    assert result["ok"]
    assert result["missing"] == []


def test_wheel_runtime_asset_omission_fails(tmp_path):
    wheel = tmp_path / "radia_mcp-test.whl"
    required = set(MODULE.REQUIRED_ASSETS)
    required.remove("radia_mcp/paper_writing/skill.md")
    _write_wheel(wheel, required)

    result = MODULE.verify_wheel_contents(wheel)

    assert not result["ok"]
    assert result["missing"] == ["radia_mcp/paper_writing/skill.md"]
