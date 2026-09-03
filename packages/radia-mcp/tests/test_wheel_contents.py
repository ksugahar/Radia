from __future__ import annotations

import importlib.util
import zipfile
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PACKAGE_ROOT / "tools" / "verify_wheel_contents.py"
SPEC = importlib.util.spec_from_file_location("verify_wheel_contents", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _write_wheel(path: Path, names: set[str]) -> None:
    with zipfile.ZipFile(path, "w") as archive:
        for name in names:
            archive.writestr(name, "fixture")


def test_wheel_runtime_assets_are_complete(tmp_path):
    wheel = tmp_path / "radia_mcp-test.whl"
    _write_wheel(wheel, set(MODULE.REQUIRED_ASSETS))

    result = MODULE.verify_wheel_contents(wheel)

    assert result["ok"]
    assert result["missing"] == []


def test_wheel_runtime_asset_omission_fails(tmp_path):
    wheel = tmp_path / "radia_mcp-test.whl"
    required = set(MODULE.REQUIRED_ASSETS)
    required.remove("radia_mcp/paper_writing/skill.md")
    _write_wheel(wheel, required)

    result = MODULE.verify_wheel_contents(wheel)

    assert not result["ok"]
    assert result["missing"] == ["radia_mcp/paper_writing/skill.md"]
