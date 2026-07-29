"""Keep package metadata and minimal CI on the supported FastMCP SDK line."""

from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10
    from importlib import import_module

    tomllib = import_module("tomli")


ROOT = Path(__file__).resolve().parents[3]
PACKAGE_ROOT = ROOT / "packages" / "radia-mcp"
SDK_REQUIREMENT = "mcp>=1.0,<2"


def test_mcp_sdk_dependency_and_minimal_ci_are_synchronized():
    project = tomllib.loads(
        (PACKAGE_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )
    assert SDK_REQUIREMENT in project["project"]["dependencies"]

    workflow = (ROOT / ".github" / "workflows" / "radia-mcp-matrix.yml").read_text(
        encoding="utf-8"
    )
    assert f'"{SDK_REQUIREMENT}"' in workflow


def test_supported_sdk_still_exposes_fastmcp():
    from mcp.server.fastmcp import FastMCP

    assert FastMCP is not None
