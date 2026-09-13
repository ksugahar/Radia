#!/usr/bin/env python3
"""Run one bounded Radia pytest tier from the checked manifest."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

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
    parser.add_argument("--since", help="Add affected lightweight contracts since this commit; empty means all registered impacts")
    args = parser.parse_args(argv)

    try:
        paths, budget = load_profile(args.profile)
        if args.since is not None:
            changed = None
            previous_manifest = None
            if args.since:
                diff = subprocess.run(
                    ['git', '-c', f'safe.directory={ROOT}', 'diff', '--name-only', '-z',
                     args.since, 'HEAD', '--'], cwd=ROOT, capture_output=True)
                if diff.returncode == 0:
                    changed = diff.stdout.decode('utf-8').split('\0')
                    if 'tests/test_tier_manifest.json' in changed:
                        previous_manifest = read_previous_manifest(args.since)
            paths = select_impact_tests(paths, changed, previous_manifest=previous_manifest)
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


def read_previous_manifest(ref: str) -> dict | None:
    """Read the exact comparison tree; missing or invalid evidence stays unknown."""
    try:
        result = subprocess.run(
            ['git', '-c', f'safe.directory={ROOT}', 'show',
             f'{ref}:tests/test_tier_manifest.json'],
            cwd=ROOT, capture_output=True, timeout=30, check=False,
        )
        if result.returncode != 0:
            return None
        data = json.loads(result.stdout.decode('utf-8'))
        return data if isinstance(data, dict) else None
    except (OSError, ValueError, subprocess.TimeoutExpired):
        return None


def changed_impact_tests(current: dict, previous: dict | None) -> set[str] | None:
    """Return tests owned by changed rules; structural changes remain broad."""
    if not isinstance(previous, dict):
        return None
    if ({key: value for key, value in current.items() if key != 'impact_rules'}
            != {key: value for key, value in previous.items() if key != 'impact_rules'}):
        return None
    before = previous.get('impact_rules')
    after = current.get('impact_rules')
    if not isinstance(before, dict) or not isinstance(after, dict):
        return None
    selected: set[str] = set()
    for source in before.keys() | after.keys():
        old_tests = before.get(source, [])
        new_tests = after.get(source, [])
        if old_tests == new_tests:
            continue
        if not isinstance(old_tests, list) or not isinstance(new_tests, list):
            return None
        selected.update(old_tests)
        selected.update(new_tests)
    return selected


def select_impact_tests(
    paths: list[str], changed: list[str] | None, *, previous_manifest: dict | None = None,
) -> list[str]:
    """Match files exactly and trailing-slash directories recursively; unknown bases fail broad."""
    manifest = json.loads(MANIFEST.read_text(encoding='utf-8'))
    rules = manifest.get('impact_rules', {})
    selected = list(paths)
    manifest_tests = set()
    if changed is not None and 'tests/test_tier_manifest.json' in changed:
        manifest_tests = changed_impact_tests(manifest, previous_manifest)
        if manifest_tests is None:
            changed = None
            manifest_tests = set()
        else:
            selected.extend(sorted(manifest_tests))
    for source, tests in rules.items():
        if (changed is None or source in changed
                or (source.endswith('/') and any(path.startswith(source) for path in changed))
                or any(test in changed for test in tests)):
            selected.extend(tests)
    selected = list(dict.fromkeys(selected))
    for path in selected:
        if not (ROOT / path).is_file():
            raise ValueError(f'impact rule names missing test: {path}')
    return selected


if __name__ == "__main__":
    raise SystemExit(main())
