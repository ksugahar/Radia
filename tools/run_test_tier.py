#!/usr/bin/env python3
"""Run one bounded Radia pytest tier from the checked manifest."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "tests" / "test_tier_manifest.json"


def load_profile(name: str) -> tuple[list[str], float | None]:
    data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    profiles = data.get("profiles", {})
    if name not in profiles:
        choices = ", ".join(sorted(profiles))
        raise ValueError(f"unknown test tier {name!r}; choose one of {choices}")

    def collect(profile_name: str, stack: tuple[str, ...] = ()) -> list[str]:
        if profile_name in stack:
            chain = " -> ".join((*stack, profile_name))
            raise ValueError(f"cyclic test-tier inheritance: {chain}")
        try:
            profile = profiles[profile_name]
        except KeyError as exc:
            raise ValueError(f"unknown parent test tier {profile_name!r}") from exc

        paths: list[str] = []
        parent = profile.get("extends")
        if parent:
            paths.extend(collect(parent, (*stack, profile_name)))
        paths.extend(profile.get("paths", []))
        return paths

    paths = collect(name)
    duplicates = sorted({path for path in paths if paths.count(path) > 1})
    if duplicates:
        raise ValueError(f"test-tier manifest repeats files: {', '.join(duplicates)}")
    missing = [path for path in paths if not (ROOT / path).is_file()]
    if missing:
        raise ValueError(f"test-tier manifest names missing files: {', '.join(missing)}")

    budget = profiles[name].get("max_elapsed_seconds")
    if budget is not None and (not isinstance(budget, (int, float)) or budget <= 0):
        raise ValueError(f"profile {name!r} has an invalid max_elapsed_seconds")
    return paths, float(budget) if budget is not None else None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", default="fast-contracts")
    parser.add_argument("--junitxml", type=Path)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--collect-only", action="store_true")
    args = parser.parse_args(argv)

    try:
        paths, budget = load_profile(args.profile)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"test-tier configuration error: {exc}", file=sys.stderr)
        return 2

    command = [sys.executable, "-m", "pytest", *paths]
    if args.collect_only:
        command.append("--collect-only")
    if not args.verbose:
        command.append("-q")
    if args.junitxml:
        args.junitxml.parent.mkdir(parents=True, exist_ok=True)
        command.append(f"--junitxml={args.junitxml}")
    if not args.collect_only:
        command.extend(("--durations=10", "--durations-min=0.25"))

    print(f"Radia test tier {args.profile}: {len(paths)} files", flush=True)
    started = time.monotonic()
    result = subprocess.run(command, cwd=ROOT)
    elapsed = time.monotonic() - started
    print(f"Radia test tier {args.profile}: elapsed={elapsed:.2f}s", flush=True)
    if result.returncode == 0 and not args.collect_only and budget is not None:
        if elapsed > budget:
            print(
                f"test-tier budget exceeded: {args.profile} took {elapsed:.2f}s "
                f"(budget {budget:.2f}s); move heavy evidence to validation_test/",
                file=sys.stderr,
            )
            return 1
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
