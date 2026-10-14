import asyncio
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from radia_mcp.matlab import (
    matlab_extension_contract,
    matlab_extension_path,
    matlab_official_server_config,
    matlab_official_server_status,
)


def test_generic_extension():
    contract = matlab_extension_contract()

    assert contract["ok"]
    assert contract["tool_count"] == 43
    assert contract["matlab_function_count"] == 86


def test_official_server_composition():
    config = matlab_official_server_config(
        "new_nodesktop",
        include_generic_extension=True,
    )

    assert len(config["extension_files"]) == 1
    assert config["matlab_setup_code"].startswith("addpath(")


def test_server_tools():
    from radia_mcp.matlab.server import mcp, matlab_extension_contract as tool

    tool_names = {item.name for item in asyncio.run(mcp.list_tools())}

    assert "matlab_extension_contract" in tool_names
    assert "matlab_official_server_status" in tool_names
    assert json.loads(tool())["tool_count"] == 43


def test_extension_schema_matches_matlab_function_contract():
    extension_path = matlab_extension_path()
    payload = json.loads(extension_path.read_text(encoding="utf-8"))
    tools = payload["tools"]
    signatures = payload["signatures"]
    tool_names = [tool["name"] for tool in tools]
    matlab_root = extension_path.parent.parent / "matlab" / "+radia_mcp_matlab"
    matlab_functions = {path.stem for path in matlab_root.glob("*.m")}

    assert len(tool_names) == len(set(tool_names))
    assert set(tool_names) == set(signatures)

    for tool in tools:
        name = tool["name"]
        schema = tool["inputSchema"]
        properties = schema["properties"]
        required = schema["required"]
        argument_order = signatures[name]["input"]["order"]
        function_name = signatures[name]["function"].rsplit(".", maxsplit=1)[-1]

        assert argument_order == required
        assert set(argument_order) == set(properties)
        assert function_name in matlab_functions

        seed_schema = properties["seed"]
        assert seed_schema["type"] == "integer"
        assert seed_schema["minimum"] == 0


@pytest.fixture
def official_foundation(monkeypatch, tmp_path):
    from radia_mcp.matlab import runtime

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    monkeypatch.delenv("RADIA_MATLAB_MCP_SERVER", raising=False)
    monkeypatch.setattr(runtime.shutil, "which", lambda _: None)
    root = tmp_path / ".matlab" / "agentic-toolkits"
    (root / "bin").mkdir(parents=True)
    return runtime, root


def test_official_server_discovers_current_executable_name(official_foundation):
    _, root = official_foundation
    for name in ("matlab-mcp-server.exe", "matlab-mcp-server-windows-x64.exe"):
        (root / "bin" / name).touch()
    assert Path(matlab_official_server_config()["command"]).name == "matlab-mcp-server.exe"


def test_official_server_explicit_override_is_preserved(official_foundation, monkeypatch):
    _, root = official_foundation
    (root / "bin" / "matlab-mcp-server.exe").touch()
    monkeypatch.setenv("RADIA_MATLAB_MCP_SERVER", "custom-server")
    assert matlab_official_server_config()["command"] == "custom-server"


def test_missing_server_does_not_launch_a_process(official_foundation, monkeypatch):
    runtime, _ = official_foundation
    monkeypatch.setattr(runtime.subprocess, "run", lambda *a, **k: pytest.fail("unexpected launch"))
    result = matlab_official_server_status()
    assert result["status"] == "needs_attention"
    assert result["available"] is False


@pytest.mark.parametrize("probe", [
    SimpleNamespace(returncode=1, stdout="", stderr="failed"),
    SimpleNamespace(returncode=0, stdout="", stderr=""),
    subprocess.TimeoutExpired("server", 15),
    OSError("cannot execute"),
])
def test_official_version_probe_failure_is_diagnostic(official_foundation, monkeypatch, probe):
    runtime, root = official_foundation
    (root / "bin" / "matlab-mcp-server.exe").touch()

    def run(*args, **kwargs):
        assert args[0][-1] == "--version"
        assert kwargs["timeout"] == 15
        if isinstance(probe, Exception):
            raise probe
        return probe

    monkeypatch.setattr(runtime.subprocess, "run", run)
    result = matlab_official_server_status()
    assert result["available"] is True
    assert result["status"] == "needs_attention"
    assert result["version_probe_error"]


@pytest.mark.parametrize("contents", ["{broken", "[]", '{"matlab": null}'])
def test_invalid_toolkit_configuration_is_diagnostic(official_foundation, contents):
    runtime, root = official_foundation
    (root / "config.json").write_text(contents, encoding="utf-8")
    result = runtime._agentic_toolkit_status()
    assert result["available"] is False
    assert result["error"]


def test_official_status_reports_versions_through_server(official_foundation, monkeypatch):
    runtime, root = official_foundation
    (root / "bin" / "matlab-mcp-server.exe").touch()
    (root / "config.json").write_text(json.dumps({
        "mcpServerVersion": "0.11.2", "matlab": {"version": "R2026a"},
        "toolkits": {"simulink": {"version": "2026.07.22"}},
    }), encoding="utf-8-sig")
    monkeypatch.setattr(runtime.subprocess, "run", lambda *a, **k: SimpleNamespace(
        returncode=0, stdout="MATLAB MCP Server 0.11.2\n", stderr=""))
    from radia_mcp.matlab.server import matlab_official_server_status as tool

    result = json.loads(tool())
    assert result["status"] == "ready"
    assert result["version"] == "MATLAB MCP Server 0.11.2"
    assert result["agentic_toolkits"]["matlab_version"] == "R2026a"
    assert "does not claim upstream currency" in result["update_policy"]
