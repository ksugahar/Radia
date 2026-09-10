"""Wheel evidence must not come from an editable or shadowing checkout."""

import importlib.util
import subprocess
from pathlib import Path

import pytest


def _make_junction(link: Path, target: Path) -> bool:
    """Create a Windows directory junction, or report that we cannot."""
    link.parent.mkdir(parents=True, exist_ok=True)
    target.mkdir(parents=True, exist_ok=True)
    try:
        completed = subprocess.run(
            ["cmd", "/c", "mklink", "/J", str(link), str(target)],
            capture_output=True, text=True)
    except OSError:
        return False
    return completed.returncode == 0 and link.exists()


@pytest.fixture
def source_guard():
    repo = Path(__file__).resolve().parents[2]
    path = repo / "validation_test/quadrupole_cefc2020/qmag_source.py"
    spec = importlib.util.spec_from_file_location("_qmag_source_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Distribution:
    version = "4.95.81"

    def __init__(self, root, *, editable=False, files=True, wheel=True):
        self.root = root
        self.editable = editable
        self.files = [Path("radia/__init__.py")] if files else None
        self.wheel = wheel

    def read_text(self, name):
        if name == "direct_url.json" and self.editable:
            return '{"dir_info": {"editable": true}, "url": "file:///other-tree"}'
        if name == "WHEEL" and self.wheel:
            return "Wheel-Version: 1.0\n"
        return None

    def locate_file(self, path):
        return self.root / path


def test_repo_accepts_only_requested_checkout(source_guard, tmp_path):
    checkout = tmp_path / "repo/src"
    result = source_guard.require_radia_source(
        "repo", checkout / "radia/__init__.py", checkout)
    assert result["resolved_from_checkout"] is True
    with pytest.raises(RuntimeError, match="not the requested checkout"):
        source_guard.require_radia_source(
            "repo", tmp_path / "another/src/radia/__init__.py", checkout)


def test_repo_accepts_a_checkout_staged_as_a_junction(source_guard, tmp_path):
    """mdx1 and hibino stage <checkout>/src/radia as a junction into the venv.

    Path.resolve() follows it, so a plain equality test rejected the very tree
    the runner asked for and the default --radia-source repo aborted on both
    heavy validation hosts.
    """
    site_packages = tmp_path / "venv/Lib/site-packages/radia"
    checkout = tmp_path / "repo/src"
    link = checkout / "radia"
    if not _make_junction(link, site_packages):
        pytest.skip("cannot create a directory junction on this filesystem")
    (site_packages / "__init__.py").write_text("", encoding="utf-8")

    result = source_guard.require_radia_source(
        "repo", link / "__init__.py", checkout)
    assert result["imported_as"] == str((link / "__init__.py").absolute())
    # The junction is what the interpreter went through, so an 'installed'
    # claim from the same name must still be refused.
    with pytest.raises(RuntimeError, match="which is the checkout"):
        source_guard.require_radia_source(
            "installed", link / "__init__.py", checkout)


def test_installed_rejects_current_checkout(source_guard, tmp_path):
    checkout = tmp_path / "repo/src"
    with pytest.raises(RuntimeError, match="which is the checkout"):
        source_guard.require_radia_source(
            "installed", checkout / "radia/__init__.py", checkout)


@pytest.mark.parametrize("editable", [False, True])
def test_installed_rejects_another_checkout(source_guard, tmp_path, monkeypatch, editable):
    dist = Distribution(tmp_path / "site-packages", editable=editable)
    monkeypatch.setattr(source_guard.metadata, "distribution", lambda name: dist)
    with pytest.raises(RuntimeError, match="editable|does not match"):
        source_guard.require_radia_source(
            "installed", tmp_path / "another/src/radia/__init__.py", tmp_path / "repo/src")


def test_installed_accepts_recorded_wheel(source_guard, tmp_path, monkeypatch):
    dist = Distribution(tmp_path / "site-packages")
    monkeypatch.setattr(source_guard.metadata, "distribution", lambda name: dist)
    loaded = dist.root / "radia/__init__.py"
    result = source_guard.require_radia_source("installed", loaded, tmp_path / "repo/src")
    assert result["distribution_module"] == str(loaded.resolve())
    assert result["distribution_version"] == dist.version
    assert result["editable"] is False


@pytest.mark.parametrize("missing", ["files", "wheel", "distribution"])
def test_installed_requires_distribution_evidence(source_guard, tmp_path, monkeypatch, missing):
    dist = Distribution(tmp_path / "site-packages", files=missing != "files",
                        wheel=missing != "wheel")

    def lookup(name):
        if missing == "distribution":
            raise source_guard.metadata.PackageNotFoundError(name)
        return dist

    monkeypatch.setattr(source_guard.metadata, "distribution", lookup)
    with pytest.raises(RuntimeError):
        source_guard.require_radia_source(
            "installed", dist.root / "radia/__init__.py", tmp_path / "repo/src")
