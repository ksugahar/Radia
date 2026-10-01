"""Verify non-Python runtime assets in a built ``radia-mcp`` wheel."""
from __future__ import annotations

import argparse
import configparser
import pathlib
import re
import subprocess
import zipfile

PACKAGE_ROOT = pathlib.Path(__file__).resolve().parents[1]
PACKAGE_SOURCE = PACKAGE_ROOT / "src" / "radia_mcp"
# Tracked non-Python files that document the source tree and are not runtime
# assets. Everything else under src/radia_mcp must be shipped by package-data.
SOURCE_ONLY_ASSETS = frozenset({"radia_mcp/radia_ngsolve/README.md"})


def package_data_patterns(pyproject: pathlib.Path = PACKAGE_ROOT / "pyproject.toml") -> dict:
    """``[tool.setuptools.package-data]`` as {package: [glob, ...]} (stdlib only)."""
    text = pyproject.read_text(encoding="utf-8")
    table = re.search(r"^\[tool\.setuptools\.package-data\]\s*$(.*?)(?=^\[|\Z)",
                      text, re.MULTILINE | re.DOTALL)
    if table is None:
        raise ValueError(f"{pyproject}: no [tool.setuptools.package-data] table")
    patterns = {}
    for package, body in re.findall(r'^"([\w.]+)"\s*=\s*\[(.*?)\]', table.group(1),
                                    re.MULTILINE | re.DOTALL):
        patterns[package] = re.findall(r'"([^"]+)"', body)
    return patterns


def required_assets() -> frozenset:
    """Every source file that package-data promises to ship in the wheel."""
    found = set()
    for package, globs in package_data_patterns().items():
        base = PACKAGE_ROOT / "src" / pathlib.Path(*package.split("."))
        for pattern in globs:
            matches = [p for p in base.glob(pattern) if p.is_file()]
            if not matches:
                raise ValueError(f"package-data {package}:{pattern} matches no source file")
            found.update(p.relative_to(PACKAGE_ROOT / "src").as_posix() for p in matches)
    return frozenset(found)


def _source_files() -> list[str]:
    try:
        listed = subprocess.run(
            ["git", "-C", str(PACKAGE_SOURCE), "ls-files", "--", "."],
            capture_output=True, text=True, check=True, timeout=60,
        ).stdout.splitlines()
        return [f"radia_mcp/{line}" for line in listed if line]
    except (OSError, subprocess.SubprocessError):
        return [p.relative_to(PACKAGE_ROOT / "src").as_posix()
                for p in PACKAGE_SOURCE.rglob("*")
                if p.is_file() and "__pycache__" not in p.parts]


def unpackaged_source_assets() -> list[str]:
    """Non-Python source files that the wheel would silently leave out."""
    shipped = required_assets()
    return sorted(
        name for name in _source_files()
        if not name.endswith((".py", ".pyc"))
        and name not in shipped and name not in SOURCE_ONLY_ASSETS
    )


REQUIRED_ASSETS = required_assets()


def verify_wheel_contents(wheel_path: str | pathlib.Path) -> dict:
    """Return required asset coverage for one wheel."""
    wheel = pathlib.Path(wheel_path)
    if not wheel.is_file():
        raise FileNotFoundError(f"wheel not found: {wheel}")
    with zipfile.ZipFile(wheel) as archive:
        names = set(archive.namelist())
        unwanted = {name for name in names if name.startswith((
            "radia_mcp/cubit/",
            "cubit_mesh_export/", "cae_mcp_core/",
        ))}
        # setuptools' reusable build/lib tree can retain deleted Python modules
        # and silently put them back into a later wheel. Every packaged Python
        # module must therefore still exist in the selected source tree.
        for name in names:
            if name.startswith("radia_mcp/") and name.endswith(".py"):
                relative = pathlib.PurePosixPath(name).relative_to("radia_mcp")
                if not PACKAGE_SOURCE.joinpath(*relative.parts).is_file():
                    unwanted.add(name)
        unwanted = sorted(unwanted)
        retired_entries = []
        for name in names:
            if name.endswith(".dist-info/entry_points.txt"):
                entries = configparser.ConfigParser(interpolation=None)
                entries.read_string(archive.read(name).decode("utf-8"))
                if entries.has_option("console_scripts", "mcp-server-cubit"):
                    retired_entries.append("mcp-server-cubit")
    missing = sorted(REQUIRED_ASSETS - names)
    return {
        "wheel": str(wheel),
        "required": sorted(REQUIRED_ASSETS),
        "missing": missing,
        "unwanted": unwanted,
        "retired_entries": retired_entries,
        "ok": not (missing or unwanted or retired_entries),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("wheel")
    args = parser.parse_args()
    result = verify_wheel_contents(args.wheel)
    if not result["ok"]:
        print("radia-mcp wheel violates its asset/ownership contract:")
        for category in ("missing", "unwanted", "retired_entries"):
            for name in result[category]:
                print(f"  {category}: {name}")
        return 1
    print(f"OK: {len(result['required'])} required runtime assets are present")
    return 0


def check_source_assets() -> int:
    """Fail when a tracked runtime asset is not covered by package-data."""
    missing = unpackaged_source_assets()
    for name in missing:
        print(f"  not in package-data: {name}")
    return 1 if missing else 0


if __name__ == "__main__":
    raise SystemExit(main())
