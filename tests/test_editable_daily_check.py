"""Daily checks share release verification without running deployment."""
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest


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


def test_daily_checker_is_connected_to_impact_ci():
    root = Path(__file__).resolve().parents[1]
    rules = json.loads((root / "tests/test_tier_manifest.json").read_text())["impact_rules"]
    assert rules["tools/verify_lab_editable.py"] == ["tests/test_editable_daily_check.py"]


def test_explicit_mcp_source_does_not_repoint_physics():
    module = load_checker()
    canonical = dict(module.expected_packages())
    actual = dict(module.expected_packages("S:/Radia/mcp-runtime/packages/radia-mcp"))
    assert actual.pop("radia-mcp") == "S:/Radia/mcp-runtime/packages/radia-mcp"
    canonical.pop("radia-mcp")
    assert actual == canonical


def test_explicit_source_root_repoints_only_monorepo_packages():
    module = load_checker()
    canonical = dict(module.expected_packages())
    root = Path("S:/Radia/release-quad/main-current")
    actual = dict(module.expected_packages(source_root=root))

    assert actual["radia"] == str(root)
    assert actual["cubit-mesh-export"] == str(root / "packages" / "cubit-mesh-export")
    assert actual["radia-mcp"] == str(root / "packages" / "radia-mcp")
    assert actual["mcp-server-document"] == canonical["mcp-server-document"]


def test_explicit_mcp_source_can_override_source_root():
    module = load_checker()
    root = Path("S:/Radia/release-quad/main-current")
    mcp = "S:/Radia/mcp-runtime/packages/radia-mcp"
    actual = dict(module.expected_packages(mcp_source=mcp, source_root=root))

    assert actual["radia"] == str(root)
    assert actual["cubit-mesh-export"] == str(root / "packages" / "cubit-mesh-export")
    assert actual["radia-mcp"] == mcp


def test_failure_exit_code(monkeypatch):
    module = load_checker()
    monkeypatch.setattr(module.release_quad, "_verify_lab_editable", lambda packages: 1)
    monkeypatch.setattr(module, "verify_against_origin_main", lambda packages: 0)
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


@pytest.mark.parametrize("running,reference", [(None, "1.4.53"), ("1.4.53", None), (None, None)])
def test_unknown_version_fails_even_when_editable_path_matches(monkeypatch, capsys, running, reference):
    """A missing answer must be visible, never counted as agreement."""
    module = load_checker()
    monkeypatch.setattr(module, "_running_version", lambda name: (running, None))
    monkeypatch.setattr(module, "_origin_main_version", lambda path: reference)
    monkeypatch.setattr(module, "expected_packages", lambda **kwargs: [("radia-mcp", "source")])
    monkeypatch.setattr(module.release_quad, "_verify_lab_editable", lambda packages: 0)
    assert module.main([]) == 4
    assert "UNVERIFIED" in capsys.readouterr().out


def test_matching_version_and_path_pass(monkeypatch):
    module = load_checker()
    monkeypatch.setattr(module, "_running_version", lambda name: ("1.4.53", "source"))
    monkeypatch.setattr(module, "_origin_main_version", lambda path: "1.4.53")
    monkeypatch.setattr(module.release_quad, "_verify_lab_editable", lambda packages: 0)
    assert module.main([]) == 0


@pytest.mark.parametrize("returncode,stdout,expected", [(0, '__version__ = "1.4.53"', "1.4.53"),
                                                     (1, "", None), (0, "no version", None)])
def test_reference_git_access_is_scoped_and_missing_data_is_unknown(monkeypatch, returncode, stdout, expected):
    module = load_checker()
    def run(command, **kwargs):
        assert command == ["git", "-c", f"safe.directory={module._REPO.as_posix()}",
                           "-C", str(module._REPO), "show", "origin/main:file.py"]
        assert kwargs["encoding"] == "utf-8"
        return SimpleNamespace(returncode=returncode, stdout=stdout)
    monkeypatch.setattr(module.subprocess, "run", run)
    assert module._origin_main_version("file.py") == expected


def test_missing_git_is_unknown_not_an_uncaught_exception(monkeypatch):
    module = load_checker()
    def missing(*args, **kwargs):
        raise FileNotFoundError("git")
    monkeypatch.setattr(module.subprocess, "run", missing)
    assert module._origin_main_version("file.py") is None


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
