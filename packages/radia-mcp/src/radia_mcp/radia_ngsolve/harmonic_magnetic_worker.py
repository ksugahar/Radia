"""Owned worker for harmonic magnetic validation."""

from __future__ import annotations

import json
import sys

from .harmonic_magnetic_validation import SCHEMA, run_harmonic_magnetic_validation


def main() -> None:
    try:
        result = run_harmonic_magnetic_validation(json.loads(sys.stdin.read()))
    except (json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        result = {
            "schema": SCHEMA,
            "status": "invalid_input",
            "pass": False,
            "solver_launched": False,
            "error": str(exc),
        }
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
