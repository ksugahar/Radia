"""Regression tests for impact-scoped preflight path discovery."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[1]


def _load_module():
    path = ROOT / "tools" / "ci_preflight.py"
    spec = importlib.util.spec_from_file_location("ci_preflight", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _git(repo: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=repo, text=True, encoding="utf-8"
    ).strip()


def test_changed_files_never_include_git_stderr(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-b", "main")
    _git(repo, "config", "user.email", "ci@example.invalid")
    _git(repo, "config", "user.name", "CI Test")
    _git(repo, "config", "core.autocrlf", "true")
    path = repo / "sample.txt"
    path.write_text("base\n", encoding="utf-8")
    _git(repo, "add", "sample.txt")
    _git(repo, "commit", "-m", "base")
    path.write_text("changed\n", encoding="utf-8")

    module = _load_module()
    module.REPO = str(repo)

    assert module._changed_since("HEAD") == ["sample.txt"]


def test_mcp_impact_lane_disables_external_pytest_plugins():
    source = (ROOT / "tools" / "ci_preflight.py").read_text(encoding="utf-8")

    assert '"PYTEST_DISABLE_PLUGIN_AUTOLOAD": "1"' in source
    assert '*pytest_targets' in source


def test_meta_catalog_loads_with_package_relative_imports():
    module = _load_module()

    catalog = module._load_meta_catalog()

    assert "meta" in catalog


def test_lab_unc_remap_supports_synthetic_current_and_historical_routes(monkeypatch):
    module = _load_module()
    monkeypatch.setattr(module, "_LAB_UNC_HOSTS", ("192.0.2.10", "192.0.2.11"))
    assert module._remap_lab_unc(r"\\192.0.2.10\share\repo\Radia\main-checkout") == r"S:\Radia\main-checkout"
    assert module._remap_lab_unc(r"\\192.0.2.11\share\repo\Radia\candidate-worktree") == r"S:\Radia\candidate-worktree"
    assert module._remap_lab_unc(r"\\192.0.2.10\share\repo\RADIA\candidate") == r"S:\RADIA\candidate"
    unrelated = r"\\unrelated-server\share\repo\Radia\candidate"
    assert module._remap_lab_unc(unrelated) == unrelated
