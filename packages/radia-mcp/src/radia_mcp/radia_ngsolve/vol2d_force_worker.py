"""Owned process entry point for thread-sensitive dimension-2 force solves."""

from __future__ import annotations

import json
import sys
from pathlib import Path


def _prefer_sibling_radia_source() -> None:
    package_file = Path(__file__).resolve()
    if len(package_file.parents) >= 6:
        source = package_file.parents[5] / "src"
        if source.is_dir() and str(source) not in sys.path:
            sys.path.insert(0, str(source))


def main() -> None:
    _prefer_sibling_radia_source()
    from .vol2d_force import FORCE_SCHEMA, analyze_vol2d_force

    try:
        request = json.loads(sys.stdin.read())
        result = analyze_vol2d_force(request)
    except (json.JSONDecodeError, TypeError, ValueError, RuntimeError) as exc:
        result = {"schema": FORCE_SCHEMA, "status": "invalid_input", "error": str(exc)}
    sys.stdout.write(json.dumps(result, separators=(",", ":")))


if __name__ == "__main__":
    main()
