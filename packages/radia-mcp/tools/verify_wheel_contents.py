"""Verify non-Python runtime assets in a built ``radia-mcp`` wheel."""
from __future__ import annotations

import argparse
import pathlib
import zipfile


REQUIRED_ASSETS = frozenset(
    {
        "radia_mcp/bibliography/data/references.bib",
        "radia_mcp/grant_writing/skill.md",
        "radia_mcp/paper_writing/skill.md",
        "radia_mcp/poster/skill.md",
        "radia_mcp/presentation/skill.md",
    }
)


def verify_wheel_contents(wheel_path: str | pathlib.Path) -> dict:
    """Return required asset coverage for one wheel."""
    wheel = pathlib.Path(wheel_path)
    if not wheel.is_file():
        raise FileNotFoundError(f"wheel not found: {wheel}")
    with zipfile.ZipFile(wheel) as archive:
        names = set(archive.namelist())
    missing = sorted(REQUIRED_ASSETS - names)
    return {
        "wheel": str(wheel),
        "required": sorted(REQUIRED_ASSETS),
        "missing": missing,
        "ok": not missing,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("wheel")
    args = parser.parse_args()
    result = verify_wheel_contents(args.wheel)
    if result["missing"]:
        print("radia-mcp wheel is missing required runtime assets:")
        for name in result["missing"]:
            print(f"  {name}")
        return 1
    print(f"OK: {len(result['required'])} required runtime assets are present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
