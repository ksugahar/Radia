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
