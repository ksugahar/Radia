"""Owned process entry point for AGE periodic-motion analysis."""

from __future__ import annotations

import json
import sys


def main() -> None:
    from .age_periodic_motion import SCHEMA, analyze_age_periodic_motion

    try:
        request = json.loads(sys.stdin.read())
        result = analyze_age_periodic_motion(request)
    except (json.JSONDecodeError, TypeError, ValueError, RuntimeError) as exc:
        result = {"schema": SCHEMA, "status": "invalid_input", "error": str(exc)}
    sys.stdout.write(json.dumps(result, separators=(",", ":")))


if __name__ == "__main__":
    main()
