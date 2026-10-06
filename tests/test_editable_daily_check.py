"""The old LAB editable command fails closed under the wheel-only policy."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def intent_file(tmp_path, monkeypatch):
    path = tmp_path / "editable-intent.json"
    monkeypatch.setenv("RADIA_EDITABLE_INTENT_FILE", str(path))
    return path


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_legacy_lab_editable_verifier_only_explains_new_command(capsys):
    tool = load(ROOT / "tools" / "verify_lab_editable.py", "retired_lab_editable")
    assert tool.main([]) == 2
    output = capsys.readouterr().err
    assert "LAB uses fixed wheels" in output
    assert "release_quad.py verify-editable" in output


def test_verify_editable_checks_lab_wheel_and_100_editable(monkeypatch):
    tool = load(ROOT / "tools" / "release_quad.py", "verify_runtime_roles")
    assert not hasattr(tool, "EDITABLE_REPO_LAB_ENV")
    assert not hasattr(tool, "_verify_lab_editable")
    calls = []
    monkeypatch.setattr(tool, "_verify_lab_wheel", lambda: calls.append("lab-wheel") or 0)
    monkeypatch.setattr(tool, "_verify_100_editable", lambda: calls.append("100-editable") or 0)
    assert tool.cmd_verify_editable(None) == 0
    assert calls == ["lab-wheel", "100-editable"]


def test_release_quad_repoint_refuses_lab_before_any_remote_or_local_action(monkeypatch, capsys):
    tool = load(ROOT / "tools" / "release_quad.py", "reject_lab_repoint")
    actions = []
    monkeypatch.setattr(tool.editable_intent, "main", lambda argv: actions.append(argv) or 0)
    monkeypatch.setattr(
        tool, "_remote_editable_intent", lambda *a, **k: actions.append(a) or (0, {}, "")
    )
    args = SimpleNamespace(
        host="lab",
        package=["radia"],
        source=["S:/tree"],
        reason="test",
        via=None,
        record_current=False,
        rollback=False,
        dry_run=False,
        require_pushed=False,
    )
    assert tool.cmd_repoint(args) == 2
    assert actions == []
    assert "only on 100号機" in capsys.readouterr().out


def test_release_quad_repoint_defaults_to_100():
    text = (ROOT / "tools" / "release_quad.py").read_text(encoding="utf-8")
    assert 'rp.add_argument("--host", default="100"' in text


@pytest.mark.parametrize(
    "args",
    [
        ["repoint", "--package", "radia", "--source", "{source}", "--reason", "policy test"],
        ["repoint", "--record-current", "--package", "radia", "--reason", "policy test"],
        ["repoint", "--rollback", "--package", "radia"],
    ],
)
def test_editable_intent_cli_refuses_mutation_from_lab(
    intent_file, monkeypatch, tmp_path, capsys, args
):
    tool = load(ROOT / "tools" / "editable_intent.py", "editable_intent_lab_refusal")
    monkeypatch.setattr(tool.platform, "node", lambda: "LAB")
    source = tmp_path / "source"
    source.mkdir()
    argv = [part.format(source=source) for part in args]
    assert tool.main(argv) == tool.EXIT_PRECONDITION
    assert "only on 100号機" in capsys.readouterr().out
    assert not Path(tool.intent_file_path()).exists()
