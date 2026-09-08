"""Daily checks share release verification without running deployment."""
import importlib.util
from pathlib import Path
import sys


def load_checker():
    tools = Path(__file__).resolve().parents[1] / "tools"
    sys.path.insert(0, str(tools))
    try:
        spec = importlib.util.spec_from_file_location(
            "daily_editable_check", tools / "verify_lab_editable.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.pop(0)


def test_explicit_mcp_source_does_not_repoint_physics():
    module = load_checker()
    canonical = dict(module.expected_packages())
    actual = dict(module.expected_packages("S:/Radia/mcp-runtime/packages/radia-mcp"))
    assert actual.pop("radia-mcp") == "S:/Radia/mcp-runtime/packages/radia-mcp"
    canonical.pop("radia-mcp")
    assert actual == canonical


def test_failure_exit_code(monkeypatch):
    module = load_checker()
    monkeypatch.setattr(module.release_quad, "_verify_lab_editable", lambda packages: 1)
    assert module.main([]) == 4


def test_noneditable_install_cannot_pass(monkeypatch, capsys):
    module = load_checker().release_quad
    monkeypatch.setattr(module, "_pip_show", lambda name: {
        "version": "1", "location": "S:/Radia/01_GitHub"})
    assert module._verify_lab_editable([("radia", "S:/Radia/01_GitHub")]) == 1
    output = capsys.readouterr().out
    assert "NOT EDITABLE" in output
    assert "Stop-Process" not in output
