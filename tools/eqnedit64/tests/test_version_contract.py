"""The PE numeric version must not lag the visible/package version."""
from pathlib import Path
import re
import tomllib

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT.parents[1] / "packages/eqnedit64"


def test_version_contract():
    version = tomllib.loads((PACKAGE / "pyproject.toml").read_text("utf-8"))["project"]["version"]
    major, minor, patch = map(int, version.split("."))
    header = (ROOT / "src/eqnedit64_version.h").read_text("utf-8")
    for key, value in (("MAJOR", major), ("MINOR", minor), ("PATCH", patch)):
        assert re.search(rf"#define EQNEDIT64_VERSION_{key}\s+{value}\b", header)
    assert f'#define EQNEDIT64_VERSION_TUPLE {major},{minor},{patch},0' in header
    assert f'#define EQNEDIT64_VERSION_TEXT "{version}"' in header
    assert f'#define EQNEDIT64_VERSION_TEXT_W L"{version}"' in header
    assert f'__version__ = "{version}"' in (PACKAGE / "src/eqnedit64/__init__.py").read_text("utf-8")
    assert f'project(Eqnedit64 VERSION {version} ' in (ROOT / "CMakeLists.txt").read_text("utf-8")
    assert f'var BUILD = "{version} (' in (ROOT / "web/equation-editor.js").read_text("utf-8")
    for changelog in [ROOT / "CHANGELOG.md", PACKAGE / "CHANGELOG.md"]:
        assert re.search(rf"^## {re.escape(version)}(?: |$)", changelog.read_text("utf-8"), re.M)
    for readme in [ROOT / "README.md", PACKAGE / "README.md"]:
        assert f"Source version: **{version}**" in readme.read_text("utf-8")
    assert f"- 対象: Eqnedit64 {version}" in (ROOT / "docs/GUI_SPEC.md").read_text("utf-8")
    # The installed-wheel check pins the version too; a stale pin stops the
    # hosted job before the Office bit gate runs (3.1.3 candidate, 2026-10-08).
    assert f'assert eqnedit64.__version__ == "{version}"' in (
        PACKAGE / "tests/verify_installed.py").read_text("utf-8")
