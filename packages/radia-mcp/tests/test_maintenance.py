"""Maintenance must preserve user intent, fail closed, and stay idempotent."""
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from radia_mcp import maintenance as m
from radia_mcp import _maintenance_guard as guard


AUTH = dict(owner="test-task", reason="test change", clients_idle=True, targets="test client")


@pytest.fixture(autouse=True)
def isolated_guard(monkeypatch, tmp_path):
    monkeypatch.setattr(guard, "STATE_ROOT", tmp_path / "maintenance-state")


def test_json_plan_apply_is_explicit_preserves_policies_and_exact_repeat(tmp_path):
    path = tmp_path / "client.json"
    original = {"other": {"secret": "never-log-this"}, "mcpServers": {
        "external": {"url": "https://example.invalid"},
        "radia-motion": {"command": "python", "args": ["-m", "radia_mcp.radia_motion.server"],
                         "disabled": True, "tools": {"x": {"approval_policy": "always"}},
                         "env": {"CUSTOM": "kept"}},
    }}
    path.write_text(json.dumps(original), encoding="utf-8")
    before, after, report = m.plan_config(path, sys.executable)
    assert path.read_bytes() == before
    assert "never-log-this" not in json.dumps(report)
    backup = m.apply_config(path, before, after, **AUTH)
    assert Path(backup).read_bytes() == before
    updated = json.loads(path.read_bytes())
    assert updated["other"] == original["other"]
    assert updated["mcpServers"]["external"] == original["mcpServers"]["external"]
    motion = updated["mcpServers"]["radia-motion"]
    assert motion["disabled"] and motion["env"] == {"CUSTOM": "kept"}
    assert motion["tools"] == original["mcpServers"]["radia-motion"]["tools"]
    before2, after2, second = m.plan_config(path, sys.executable)
    assert before2 == after2 == after
    assert not second["changed"]
    assert m.apply_config(path, before2, after2) is None


@pytest.mark.parametrize("args", [["-m", "radia_mcp.radia_motion.server", "--profile", "ih"], ["-c", "custom()"]])
def test_custom_launch_arguments_block_the_entire_plan(tmp_path, args):
    path = tmp_path / "client.json"
    path.write_text(json.dumps({"mcpServers": {"radia-motion": {"command": "python", "args": args}}}))
    before, after, report = m.plan_config(path, sys.executable)
    assert report["conflicts"] == ["radia-motion"]
    assert not report["changed"] and before == after


def test_toml_preserves_comments_permissions_bom_and_newlines(tmp_path):
    tomlkit = pytest.importorskip("tomlkit")
    path = tmp_path / "config.toml"
    path.write_bytes(b'\xef\xbb\xbf# secret comment\r\nmodel = "kept"\r\n[mcp_servers.radia-meta]\r\ncommand = "python" # launcher\r\nargs = ["-m", "radia_mcp.meta.server"]\r\nenabled = false\r\n[mcp_servers.radia-meta.tools.x]\r\napproval_policy = "always"\r\n')
    before, after, report = m.plan_config(path, sys.executable)
    assert report["changed"] and not report["conflicts"]
    assert after.startswith(b'\xef\xbb\xbf# secret comment\r\n')
    assert b"# launcher" in after
    parsed = tomlkit.parse(after.decode("utf-8-sig"))
    assert parsed["model"] == "kept"
    assert parsed["mcp_servers"]["radia-meta"]["enabled"] is False
    assert parsed["mcp_servers"]["radia-meta"]["tools"]["x"]["approval_policy"] == "always"
    m.apply_config(path, before, after, **AUTH)
    assert m.plan_config(path, sys.executable)[1] == after


def test_apply_refuses_concurrent_change_and_does_not_overwrite(tmp_path):
    path = tmp_path / "client.json"
    before, after, _ = m.plan_config(path, sys.executable)
    path.write_bytes(b"another session")
    with pytest.raises(ValueError, match="changed"):
        m.apply_config(path, before, after, **AUTH)
    assert path.read_bytes() == b"another session"


def test_serve_calls_main_even_without_module_guard_and_restores_argv(monkeypatch, capsys):
    captured = []
    original = sys.argv
    monkeypatch.setattr(m.importlib, "import_module", lambda module: SimpleNamespace(main=lambda: captured.append(list(sys.argv))))
    assert m.main(["serve", "radia-matlab", "--selftest"]) == 0
    assert captured == [["radia_mcp.matlab.server", "--selftest"]]
    assert sys.argv is original
    assert capsys.readouterr().out == ""


def test_unknown_server_never_imports_domain(monkeypatch):
    monkeypatch.setattr(m.importlib, "import_module", lambda _: pytest.fail("must not import"))
    assert m.main(["serve", "unknown"]) == 2


def test_doctor_rejects_wrong_source_and_version(monkeypatch, tmp_path):
    monkeypatch.setattr(m.metadata, "distribution", lambda _: SimpleNamespace(
        version="old", read_text=lambda _: json.dumps({"dir_info": {"editable": True}, "url": "file:///old"})))
    report = m.doctor(tmp_path, "new")
    assert report["status"] == "fail"
    assert "loaded-source-mismatch" in report["issues"]
    assert "installed-version-mismatch" in report["issues"]


def test_every_baseline_launcher_is_a_catalog_module():
    for key in m.BASELINE.values():
        assert m.server_module(key).startswith("radia_mcp.")
        assert m.standard_entry(key, sys.executable)["args"][-1] == key


def test_versioned_python_executable_is_migratable():
    expected = {"args": ["-s", "-m", "radia_mcp.maintenance", "serve", "meta"]}
    assert m._migratable({"command": "/usr/bin/python3.12", **expected}, "meta", expected)


def test_duplicate_json_keys_are_rejected(tmp_path):
    path = tmp_path / "client.json"
    path.write_text('{"mcpServers": {}, "mcpServers": {}}')
    with pytest.raises(ValueError, match="Duplicate"):
        m.plan_config(path, sys.executable)
