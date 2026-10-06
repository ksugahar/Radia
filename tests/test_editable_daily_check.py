"""The old LAB editable command fails closed under the wheel-only policy."""

from __future__ import annotations

import importlib.util
from pathlib import Path
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
