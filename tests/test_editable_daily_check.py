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


def test_source_behind_origin_main_is_drift(monkeypatch, capsys):
    """A correct pointer at a stale tree must not read as clean.

    Measured 2026-09-10: every LAB editable satisfied the path check while
    running older code than origin/main (radia-mcp 1.4.39 vs 1.4.53).
    """
    module = load_checker()
    monkeypatch.setattr(module, "_running_version",
                        lambda name: ("1.4.39", "S:/Radia/01_GitHub/.../__init__.py"))
    monkeypatch.setattr(module, "_origin_main_version", lambda path: "1.4.53")
    assert module.verify_against_origin_main([("radia-mcp", "S:/Radia/01_GitHub")]) == 1
    assert "origin/main carries 1.4.53" in capsys.readouterr().out


def test_source_ahead_of_origin_main_is_not_drift(monkeypatch):
    """Mid-development the checkout leads origin/main; that is not drift."""
    module = load_checker()
    monkeypatch.setattr(module, "_running_version", lambda name: ("1.4.60", "x"))
    monkeypatch.setattr(module, "_origin_main_version", lambda path: "1.4.53")
    assert module.verify_against_origin_main([("radia-mcp", "S:/Radia/01_GitHub")]) == 0


def test_release_order_compares_numerically_not_lexically():
    """4.95.9 precedes 4.95.81; string order would invert that."""
    module = load_checker()
    assert module._release_order("4.95.9") < module._release_order("4.95.81")


def test_unknown_version_is_skipped_rather_than_passed(monkeypatch, capsys):
    """A missing answer must be visible, never counted as agreement."""
    module = load_checker()
    monkeypatch.setattr(module, "_running_version", lambda name: (None, None))
    monkeypatch.setattr(module, "_origin_main_version", lambda path: "1.4.53")
    assert module.verify_against_origin_main([("radia-mcp", "S:/Radia/01_GitHub")]) == 0
    assert "skipped" in capsys.readouterr().out


def test_skip_origin_check_bypasses_the_comparison(monkeypatch):
    module = load_checker()
    monkeypatch.setattr(module.release_quad, "_verify_lab_editable", lambda packages: 0)

    def must_not_run(packages):
        raise AssertionError("origin/main comparison should have been skipped")

    monkeypatch.setattr(module, "verify_against_origin_main", must_not_run)
    assert module.main(["--skip-origin-check"]) == 0


def test_noneditable_install_cannot_pass(monkeypatch, capsys):
    module = load_checker().release_quad
    monkeypatch.setattr(module, "_pip_show", lambda name: {
        "version": "1", "location": "S:/Radia/01_GitHub"})
    assert module._verify_lab_editable([("radia", "S:/Radia/01_GitHub")]) == 1
    output = capsys.readouterr().out
    assert "NOT EDITABLE" in output
    assert "Stop-Process" not in output
