from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest


_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from radia_mcp.acoustic_fembem import (
    acoustic_fembem_extension_contract,
    acoustic_fembem_extension_path,
    acoustic_fembem_server_config,
)
from radia_mcp.acoustic_fembem.integration import main
from radia_mcp.matlab import (
    matlab_extension_contract,
    matlab_official_server_config,
    matlab_radia_acoustic_interface_contract,
)


def test_generic_matlab_extension_has_43_tools_and_86_self_contained_functions():
    contract = matlab_extension_contract()

    assert contract["ok"] is True
    assert contract["tool_count"] == 43
    assert contract["signature_count"] == 43
    assert contract["matlab_function_count"] == 86
    assert all(name.startswith("matlab_") for name in contract["tool_names"])
    assert not any("acoustic" in name for name in contract["tool_names"])

    config = matlab_official_server_config(
        profile="new_nodesktop", include_generic_extension=True
    )
    assert len(config["extension_files"]) == 1
    assert any(arg.startswith("--extension-file=") for arg in config["args"])
    assert "+radia_mcp_matlab" not in config["matlab_setup_code"]


def test_packaged_extension_is_curated_and_self_consistent():
    contract = acoustic_fembem_extension_contract()

    assert contract["ok"] is True
    assert contract["status"] == "ok"
    assert contract["tool_count"] == 10
    assert contract["signature_count"] == contract["tool_count"]
    assert len(contract["sha256"]) == 64
    assert Path(contract["extension_file"]).is_file()
    assert "acoustic_fembem_knowledge" in contract["tool_names"]
    assert "acoustic_fembem_cq_time_grid" in contract["tool_names"]
    assert not any("actor_critic" in name for name in contract["tool_names"])
    assert not any("quantum" in name for name in contract["tool_names"])

    manifest = json.loads(
        acoustic_fembem_extension_path().read_text(encoding="utf-8")
    )
    assert set(contract["tool_names"]) == set(manifest["signatures"])
    for tool in manifest["tools"]:
        name = tool["name"]
        order = manifest["signatures"][name]["input"]["order"]
        assert order == tool["inputSchema"]["required"]
        assert manifest["signatures"][name]["function"].startswith(
            "acoustic_fembem."
        )


def test_server_config_without_companion_root_is_still_actionable(monkeypatch):
    monkeypatch.delenv("ACOUSTIC_FEMBEM_ROOT", raising=False)

    config = acoustic_fembem_server_config(profile="existing")

    assert config["command_id"] == "matlab-mcp-server"
    assert config["command"]
    assert "--matlab-session-mode=existing" in config["args"]
    assert any(arg.startswith("--extension-file=") for arg in config["args"])
    assert not any(arg.startswith("--initial-working-folder=") for arg in config["args"])
    assert config["matlab_setup_code"] == ""
    assert config["matlab_entrypoints"] == []
    assert any("shareMATLABSession" in step for step in config["preflight"])


def test_server_config_validates_companion_matlab_entrypoints(tmp_path):
    root = tmp_path / "matlab-project"
    package = root / "+acoustic_fembem"
    package.mkdir(parents=True)
    (root / "matlab_api").mkdir()

    manifest = json.loads(
        acoustic_fembem_extension_path().read_text(encoding="utf-8")
    )
    for signature in manifest["signatures"].values():
        function_name = signature["function"].split(".", 1)[1]
        (package / f"{function_name}.m").write_text(
            f"function {function_name}()\nend\n", encoding="ascii"
        )

    config = acoustic_fembem_server_config(root, "auto_nodesktop")

    assert config["project_root"] == str(root.resolve())
    assert "--matlab-session-mode=auto" in config["args"]
    assert "--matlab-display-mode=nodesktop" in config["args"]
    assert any(arg.startswith("--initial-working-folder=") for arg in config["args"])
    assert "addpath(genpath(fullfile(root, 'matlab_api')))" in config["matlab_setup_code"]
    assert len(config["matlab_entrypoints"]) == config["tool_count"]


def test_server_config_rejects_incomplete_companion_project(tmp_path):
    root = tmp_path / "incomplete"
    (root / "+acoustic_fembem").mkdir(parents=True)
    (root / "matlab_api").mkdir()

    with pytest.raises(ValueError, match="does not implement"):
        acoustic_fembem_server_config(root)


def test_cli_prints_contract_json(capsys):
    assert main(["--contract-only"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["ok"] is True
    assert payload["tool_count"] == 10


def test_radia_ngsolve_mcp_exposes_integration_tools():
    from radia_mcp.radia_ngsolve.server import (
        matlab_acoustic_fembem_extension_contract as mcp_contract,
        matlab_acoustic_fembem_server_config as mcp_config,
    )

    contract = json.loads(mcp_contract())
    config = json.loads(mcp_config("", "new_nodesktop"))

    assert contract["ok"] is True
    assert config["status"] == "ok"
    assert "--matlab-session-mode=new" in config["args"]
    assert "--matlab-display-mode=nodesktop" in config["args"]


def test_environment_root_is_supported(tmp_path, monkeypatch):
    root = tmp_path / "matlab-project"
    package = root / "+acoustic_fembem"
    package.mkdir(parents=True)
    (root / "matlab_api").mkdir()
    manifest = json.loads(
        acoustic_fembem_extension_path().read_text(encoding="utf-8")
    )
    for signature in manifest["signatures"].values():
        function_name = signature["function"].split(".", 1)[1]
        (package / f"{function_name}.m").touch()
    monkeypatch.setenv("ACOUSTIC_FEMBEM_ROOT", os.fspath(root))

    config = acoustic_fembem_server_config()

    assert config["project_root"] == str(root.resolve())


def test_server_command_can_be_explicitly_configured(monkeypatch):
    monkeypatch.setenv("RADIA_MATLAB_MCP_SERVER", r"C:\tools\matlab-mcp.exe")
    monkeypatch.delenv("ACOUSTIC_FEMBEM_ROOT", raising=False)

    config = acoustic_fembem_server_config()

    assert config["command"] == r"C:\tools\matlab-mcp.exe"


def test_generic_matlab_namespace_has_no_fembem_domain_by_default():
    config = matlab_official_server_config(profile="new_nodesktop")
    contract = matlab_radia_acoustic_interface_contract()

    assert config["integration_owner"] == "radia-mcp.matlab"
    assert config["extension_files"] == []
    assert not any("acoustic-fembem" in arg for arg in config["args"])
    assert contract["production_owner"] == "radia.acoustics"
    assert contract["education_solver_owner"] == "radia_mcp.acoustic_fembem"


def test_old_combined_namespace_is_compatibility_only():
    from radia_mcp.matlab_acoustic_fembem import (
        matlab_acoustic_fembem_extension_contract as legacy_contract,
    )

    assert legacy_contract() == acoustic_fembem_extension_contract()


def test_separate_matlab_and_acoustic_fembem_servers_register_expected_tools():
    import asyncio

    from radia_mcp.acoustic_fembem.server import mcp as acoustic_mcp
    from radia_mcp.matlab.server import mcp as matlab_mcp

    matlab_tools = {tool.name for tool in asyncio.run(matlab_mcp.list_tools())}
    acoustic_tools = {tool.name for tool in asyncio.run(acoustic_mcp.list_tools())}

    assert "matlab_official_server_config" in matlab_tools
    assert "matlab_radia_acoustic_interface_contract" in matlab_tools
    assert "acoustic_fembem_extension_contract" not in matlab_tools
    assert "acoustic_fembem_extension_contract" in acoustic_tools
    assert "matlab_official_server_config" not in acoustic_tools
