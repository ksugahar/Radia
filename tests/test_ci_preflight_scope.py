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


def test_lab_unc_remap_supports_current_and_historical_nas_addresses():
    module = _load_module()

    assert module._remap_lab_unc(
        r"\\192.168.121.100\work\00_CAE\Radia\01_GitHub"
    ) == r"S:\Radia\01_GitHub"
    assert module._remap_lab_unc(
        r"\\192.168.11.100\work\00_CAE\Radia\02_Worktrees\candidate"
    ) == r"S:\Radia\02_Worktrees\candidate"
    unrelated = r"\\server\share\Radia\01_GitHub"
    assert module._remap_lab_unc(unrelated) == unrelated
