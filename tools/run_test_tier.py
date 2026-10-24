#!/usr/bin/env python3
"""Run one bounded Radia pytest tier from the checked manifest."""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "tests" / "test_tier_manifest.json"


def read_tier_data(path):
    """DATA is checked separately from executable sealing; duplicate keys fail."""
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"{path}: duplicate key {key!r}")
            result[key] = value
        return result
    return json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=pairs)


def _safe_path(root, path):
    if (not isinstance(path, str) or not path or "\\" in path or ":" in path
            or "\0" in path or PurePosixPath(path).is_absolute()
            or ".." in PurePosixPath(path).parts or any(c in path for c in "*?[]!")):
        raise ValueError(f"unsafe manifest path: actual={path!r}; expected relative literal path")
    if not (root / path).resolve().is_relative_to(root.resolve()):
        raise ValueError(f"manifest path escapes repository: {path}")


def profile_paths(data, name, stack=()):
    if name in stack:
        raise ValueError(f"cyclic test-tier inheritance: {(*stack, name)}")
    if name not in data["profiles"]:
        raise ValueError(f"unknown test tier or parent {name!r}")
    profile = data["profiles"][name]
    paths = profile_paths(data, profile["extends"], (*stack, name)) if profile.get("extends") else []
    paths += profile["paths"]
    if len(set(paths)) != len(paths):
        raise ValueError(f"test-tier manifest repeats files: {name}")
    return paths


def validate_tier_data(root, data):
    """Validate inert selector DATA using the same resolver as the CLI."""
    if (not isinstance(data, dict) or set(data) != {"schema", "profiles", "impact_rules"}
            or data["schema"] != "radia.test-tier-manifest.v1"):
        raise ValueError("tier DATA: expected schema/profiles/impact_rules contract")
    if not isinstance(data["profiles"], dict) or not data["profiles"] or not isinstance(data["impact_rules"], dict):
        raise ValueError("tier DATA profiles/impact_rules: expected objects")
    for name, profile in data["profiles"].items():
        if (not name or not isinstance(profile, dict) or "paths" not in profile
                or set(profile) - {"paths", "description", "extends", "max_elapsed_seconds"}):
            raise ValueError(f"profile {name!r}: expected paths and supported metadata only")
        if "description" in profile and not isinstance(profile["description"], str):
            raise ValueError(f"profile {name}: description must be text")
        if "extends" in profile and (not isinstance(profile["extends"], str) or not profile["extends"]):
            raise ValueError(f"profile {name}: extends must name a parent")
        budget = profile.get("max_elapsed_seconds")
        if "max_elapsed_seconds" in profile:
            try:
                valid = type(budget) in (int, float) and math.isfinite(budget) and budget > 0
            except OverflowError:
                valid = False
            if not valid:
                raise ValueError(f"profile {name!r}: invalid max_elapsed_seconds budget {budget!r}")
        if not isinstance(profile["paths"], list):
            raise ValueError(f"profile {name}: expected path array")
        for path in profile["paths"]:
            _safe_path(root, path)
            if not (root / path).is_file():
                raise ValueError(f"test-tier manifest names missing file: {path}")
    for name in data["profiles"]:
        profile_paths(data, name)
    for source, paths in data["impact_rules"].items():
        _safe_path(root, source)
        if not isinstance(paths, list) or any(not isinstance(path, str) for path in paths) or len(set(paths)) != len(paths):
            raise ValueError(f"impact rule {source}: expected unique path array")
        for path in paths:
            _safe_path(root, path)
            if not (root / path).is_file():
                raise ValueError(f"impact rule names missing test: {path}")
    return data


def load_profile(name: str) -> tuple[list[str], float | None]:
    data = validate_tier_data(ROOT, read_tier_data(MANIFEST))
    paths = profile_paths(data, name)
    budget = data["profiles"][name].get("max_elapsed_seconds")
    return paths, float(budget) if budget is not None else None


def check_light_selection(paths):
    """A changed/retired manifest cannot admit code outside the sealed union."""
    base, _ = load_profile("fast-contracts")
    permitted = set(select_impact_tests(base, None, profile_name="fast-contracts"))
    if not set(paths) <= permitted:
        raise ValueError(f"selected light tests outside audited alternatives: {sorted(set(paths) - permitted)}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", default="fast-contracts")
    parser.add_argument("--junitxml", type=Path)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--collect-only", action="store_true")
    parser.add_argument("--since", help="Add affected lightweight contracts since this commit; empty means all registered impacts")
    parser.add_argument("--numerical-plan", type=Path, help="Recompute and verify exact-source per-case selection")
    args = parser.parse_args(argv)

    try:
        paths, budget = load_profile(args.profile)
        if args.numerical_plan:
            if args.since is not None:
                raise ValueError('numerical plan cannot be mixed with legacy --since selection')
            from ci_numerical_runtime import verified_plan, selected_paths
            paths = selected_paths(verified_plan(ROOT, args.numerical_plan), args.profile)
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
            paths = select_impact_tests(paths, changed, previous_manifest=previous_manifest,
                                        profile_name=args.profile)
        if args.profile == "fast-contracts":
            check_light_selection(paths)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"test-tier configuration error: {exc}", file=sys.stderr)
        return 2

    command = [sys.executable, "-m", "pytest", *paths]
    if args.profile == "repo-contracts":
        command.extend(("-p", "repository_contract_scratch"))
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
    environment = dict(os.environ)
    if args.profile == "repo-contracts":
        environment["PYTHONPATH"] = str(ROOT / "tools") + os.pathsep + environment.get("PYTHONPATH", "")
    result = subprocess.run(command, cwd=ROOT, env=environment)
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


def changed_impact_tests(
    current: dict, previous: dict | None, *, profile_name: str | None = None,
) -> set[str] | None:
    """Select rule/membership deltas; budgets, inheritance and unknowns stay broad."""
    if not isinstance(previous, dict):
        return None
    if profile_name is not None:
        def scoped(data):
            profiles = data.get('profiles')
            if not isinstance(profiles, dict):
                return None
            selected = {}
            name = profile_name
            while name:
                if not isinstance(name, str) or name in selected or not isinstance(profiles.get(name), dict):
                    return None
                selected[name] = profiles[name]
                name = profiles[name].get('extends')
            return {**data, 'profiles': selected}
        current, previous = scoped(current), scoped(previous)
        if current is None or previous is None:
            return None
    excluded = {'impact_rules', 'profiles'}
    if ({key: value for key, value in current.items() if key not in excluded}
            != {key: value for key, value in previous.items() if key not in excluded}):
        return None
    before_profiles, after_profiles = previous.get('profiles'), current.get('profiles')
    if (not isinstance(before_profiles, dict) or not isinstance(after_profiles, dict)
            or before_profiles.keys() != after_profiles.keys()):
        return None
    selected: set[str] = set()
    for name, after_profile in after_profiles.items():
        before_profile = before_profiles[name]
        if not isinstance(before_profile, dict) or not isinstance(after_profile, dict):
            return None
        if ({key: value for key, value in before_profile.items() if key != 'paths'}
                != {key: value for key, value in after_profile.items() if key != 'paths'}):
            return None
        old_paths, new_paths = before_profile.get('paths', []), after_profile.get('paths', [])
        for paths in (old_paths, new_paths):
            if (not isinstance(paths, list) or not all(isinstance(path, str) for path in paths)
                    or len(paths) != len(set(paths))):
                return None
        selected.update(set(old_paths) ^ set(new_paths))
    before = previous.get('impact_rules')
    after = current.get('impact_rules')
    if not isinstance(before, dict) or not isinstance(after, dict):
        return None
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
    profile_name: str | None = None,
) -> list[str]:
    """Match files exactly and trailing-slash directories recursively; unknown bases fail broad."""
    manifest = json.loads(MANIFEST.read_text(encoding='utf-8'))
    rules = manifest.get('impact_rules', {})
    selected = list(paths)
    manifest_tests = set()
    if changed is not None and 'tests/test_tier_manifest.json' in changed:
        manifest_tests = changed_impact_tests(
            manifest, previous_manifest, profile_name=profile_name)
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
    current_references = {test for tests in rules.values() for test in tests}
    current_references.update(test for profile in manifest.get('profiles', {}).values()
                              for test in profile.get('paths', []))
    repository_paths = set(manifest.get("profiles", {}).get("repo-contracts", {}).get("paths", []))
    policy_paths = {"tests/test_sparsesolv_ci_contract.py"}
    runnable = []
    numerical_paths = set(manifest.get('profiles', {}).get('solver-numerics', {}).get('paths', []))
    fast_candidates = (set(profile_paths(manifest, "fast-contracts"))
                       if profile_name == "fast-contracts" else set())
    fast_candidates.update(test for tests in rules.values() for test in tests)
    fast_candidates.difference_update(repository_paths | policy_paths | numerical_paths)
    for path in selected:
        if profile_name == "fast-contracts" and path in repository_paths | policy_paths:
            continue  # These exact tests run in the always-on repository lane.
        if not (ROOT / path).is_file():
            if path in manifest_tests and path not in current_references:
                print(f'Retired test removed from manifest and checkout: {path}')
                continue
            raise ValueError(f'impact rule names missing test: {path}')
        if profile_name == 'fast-contracts' and path in numerical_paths:
            if path in paths:
                raise ValueError(f'fast-contracts explicitly contains numerical test: {path}')
            print(f'Numerical contract belongs to solver-numerics: {path}')
            continue
        if (profile_name == "fast-contracts" and path in manifest_tests
                and path not in fast_candidates and path in current_references):
            print(f'Test no longer a current fast alternative; retained in another profile: {path}')
            continue
        runnable.append(path)
    return runnable


if __name__ == "__main__":
    raise SystemExit(main())
