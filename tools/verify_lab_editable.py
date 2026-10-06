"""Compatibility notice for the retired LAB editable verifier.

LAB is a fixed-wheel test host. Run ``python tools/release_quad.py
verify-editable`` from 100号機 to check LAB's wheel and 100号機's development
editable together. This legacy entry point is read-only and never inspects or
changes an editable installation on LAB.
"""

from __future__ import annotations

import sys


def main(argv: list[str] | None = None) -> int:
    print(
        "This verifier is retired: LAB uses fixed wheels, not editable installs. "
        "From 100号機 run `python tools/release_quad.py verify-editable` to "
        "check LAB's wheel and 100号機's development editable.",
        file=sys.stderr,
    )
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
