"""Owned worker for native Kelvin-transform validation."""

from __future__ import annotations

import json
import sys

from .kelvin_open_boundary_validation import SCHEMA, run_kelvin_open_boundary_validation


def main() -> None:
    try:
        request = json.loads(sys.stdin.read())
        result = run_kelvin_open_boundary_validation(request)
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
