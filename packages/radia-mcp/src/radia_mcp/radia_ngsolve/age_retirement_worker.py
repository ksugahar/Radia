from __future__ import annotations

import json
import sys

from .age_retirement_validation import SCHEMA, run_age_retirement_validation


def main() -> None:
    try:
        result = run_age_retirement_validation(json.loads(sys.stdin.read()))
    except (json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
        result = {"schema": SCHEMA, "status": "invalid_input", "pass": False, "solver_launched": False, "error": str(exc)}
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
