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
    monkeypatch.setattr(probe, "CATALOG", {
        "second": {"entry_point": "radia-mcp serve second"},
        "first": {"entry_point": "radia-mcp serve first"},
        "domain": {"kind": "domain"},
    })
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


def _write_wheel(path: Path, names: set[str], scripts=None) -> None:
    scripts = MODULE.CONSOLE_SCRIPTS if scripts is None else scripts
    with zipfile.ZipFile(path, "w") as archive:
        for name in names:
            archive.writestr(name, "fixture")
        archive.writestr(
            "radia_mcp-test.dist-info/entry_points.txt",
            "[console_scripts]\n" + "".join(f"{s} = x:main\n" for s in sorted(scripts)),
        )


def _complete() -> set[str]:
    return set(MODULE.REQUIRED_ASSETS) | set(MODULE.served_server_modules())


def test_wheel_runtime_assets_are_complete(tmp_path):
    wheel = tmp_path / "radia_mcp-test.whl"
    _write_wheel(wheel, _complete())

    result = MODULE.verify_wheel_contents(wheel)

    assert result["ok"]
    assert result["missing"] == []


def test_wheel_runtime_asset_omission_fails(tmp_path):
    wheel = tmp_path / "radia_mcp-test.whl"
    required = _complete()
    required.remove("radia_mcp/matlab/skill.md")
    _write_wheel(wheel, required)

    result = MODULE.verify_wheel_contents(wheel)

    assert not result["ok"]
    assert result["missing"] == ["radia_mcp/matlab/skill.md"]


@pytest.mark.parametrize("member", [
    "radia_mcp/cubit/retired_native_plugin.py",
    "cubit_mesh_export/mcp/server.py", "cae_mcp_core/common/status.py",
    "radia_mcp/radia_ngsolve/comsol_v46_identity.py",
    "radia_mcp/obsolete_deleted_module.py",
])
def test_wheel_rejects_stale_build_output_from_other_owners(tmp_path, member):
    wheel = tmp_path / "candidate.whl"
    _write_wheel(wheel, _complete() | {member})
    result = MODULE.verify_wheel_contents(wheel)
    assert not result["ok"]
    assert result["unwanted"] == [member]


@pytest.mark.parametrize("extra", ["mcp-server-cubit", "mcp-server-radia-meta"])
def test_wheel_rejects_any_retired_server_script(tmp_path, extra):
    wheel = tmp_path / "candidate.whl"
    _write_wheel(wheel, _complete(), scripts=MODULE.CONSOLE_SCRIPTS | {extra})
    result = MODULE.verify_wheel_contents(wheel)
    assert not result["ok"]
    assert result["retired_entries"] == [extra]


def test_wheel_requires_the_single_launcher(tmp_path):
    wheel = tmp_path / "candidate.whl"
    _write_wheel(wheel, _complete(), scripts=MODULE.CONSOLE_SCRIPTS - {"radia-mcp"})
    result = MODULE.verify_wheel_contents(wheel)
    assert not result["ok"]
    assert result["missing_scripts"] == ["radia-mcp"]


def test_wheel_requires_every_served_server_module(tmp_path):
    wheel = tmp_path / "candidate.whl"
    names = _complete()
    names.remove("radia_mcp/radia_ih/server.py")
    _write_wheel(wheel, names)
    result = MODULE.verify_wheel_contents(wheel)
    assert not result["ok"]
    assert result["missing"] == ["radia_mcp/radia_ih/server.py"]


def test_every_tracked_runtime_asset_is_covered_by_package_data():
    # A tracked template, script or data file outside package-data is left out
    # of the wheel without error.
    assert MODULE.unpackaged_source_assets() == []
    assert "radia_mcp/acoustic_fembem/skill.md" in MODULE.REQUIRED_ASSETS
