"""Guards: radia.cubit_locator picks the newest Cubit by version number and
fails loudly on an environment override without Cubit's Python package.

Assumes: Cubit installs are named "<prefix> <major>.<minor>[...]/bin".
"""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _locator():
    spec = importlib.util.spec_from_file_location("cubit_locator", ROOT / "src/radia/cubit_locator.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _install(root, name):
    bin_dir = root / name / "bin"
    bin_dir.mkdir(parents=True)
    (bin_dir / "cubit.py").write_text("")
    return bin_dir


def test_newest_version_wins_numerically(tmp_path, monkeypatch):
    locator = _locator()
    for name in ("Coreform Cubit 2026.9", "Coreform Cubit 2026.10", "Coreform Cubit 2025.12"):
        _install(tmp_path, name)
    monkeypatch.delenv("CUBIT_BIN_DIR", raising=False)
    monkeypatch.delenv("CUBIT_INSTALL_DIR", raising=False)
    monkeypatch.setattr(locator, "_DEFAULT_PATTERNS", (str(tmp_path / "Coreform Cubit */bin"),))
    found = locator.find_cubit_bin()
    expected = str((tmp_path / "Coreform Cubit 2026.10" / "bin").resolve())
    assert found == expected, (found, expected)


def test_override_without_cubit_python_raises(tmp_path, monkeypatch):
    locator = _locator()
    monkeypatch.setenv("CUBIT_BIN_DIR", str(tmp_path))
    with pytest.raises(FileNotFoundError, match="CUBIT_BIN_DIR"):
        locator.find_cubit_bin()
