"""Owned process entry point for thread-sensitive NGSolve 2-D assembly."""

from __future__ import annotations

import json
import sys
from pathlib import Path


def _prefer_sibling_radia_source() -> None:
    package_file = Path(__file__).resolve()
    if len(package_file.parents) < 6:
        return
    source = package_file.parents[5] / "src"
    if source.is_dir() and str(source) not in sys.path:
        sys.path.insert(0, str(source))


def main() -> None:
    _prefer_sibling_radia_source()
    from .vol2d_circuit import RESULT_SCHEMA, analyze_vol2d_circuit

    try:
        request = json.loads(sys.stdin.read())
        result = analyze_vol2d_circuit(request)
    except (json.JSONDecodeError, TypeError, ValueError, RuntimeError) as exc:
        result = {
            "schema": RESULT_SCHEMA,
            "status": "invalid_input",
            "error": str(exc),
        }
    sys.stdout.write(json.dumps(result, separators=(",", ":")))


if __name__ == "__main__":
    main()
