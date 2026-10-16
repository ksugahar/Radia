#!/usr/bin/env python3
"""Check every branch candidate from Git pre-push stdin on the wheel test host."""
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).absolute().parents[1]
ZERO = "0" * 40


def candidates(lines):
    """Validate the entire input before running any remote preflight."""
    result = []
    for number, line in enumerate(lines, 1):
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(f"invalid pre-push record at line {number}")
        local_ref, head, remote_ref, base = fields
        if not all(re.fullmatch(r"[0-9a-f]{40}", sha) for sha in (head, base)):
            raise ValueError(f"invalid commit ID at line {number}")
        if not remote_ref.startswith("refs/heads/") or head == ZERO:
            continue  # Release tags have separate release gates; deletes have no candidate.
        if not local_ref:
            raise ValueError(f"missing local ref at line {number}")
        pair = (base, head)
        if pair not in result:
            result.append(pair)
    return result


def main(lines=None):
    try:
        updates = candidates(sys.stdin if lines is None else lines)
        helper = ROOT / "tools" / "ci_preflight_lab.py"
        if updates and not helper.is_file():
            raise ValueError("candidate preflight helper is missing")
        for base, head in updates:
            print(f"pre-push: impact-scoped candidate preflight {head}", flush=True)
            completed = subprocess.run(
                [sys.executable, str(helper), "--base", base, "--head", head],
                cwd=ROOT,
            )
            if completed.returncode:
                return completed.returncode
        return 0
    except (ValueError, OSError) as exc:
        print(f"pre-push: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
