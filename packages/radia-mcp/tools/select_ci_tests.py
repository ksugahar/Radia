"""Select impact-scoped radia-mcp tests and server selftests for CI.

The complete package suite remains available for explicit full audits. Normal
pull-request and main-push CI always runs a compact contract set, then adds
tests that reference the changed package family and selftests only the affected
servers. The emitted JSON is retained as CI evidence.
"""

from __future__ import annotations

import argparse
import ast
from functools import lru_cache
import json
from pathlib import Path
import re
import runpy
import subprocess
import sys
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[3]
PACKAGE_ROOT = REPO_ROOT / "packages" / "radia-mcp"
TEST_ROOT = PACKAGE_ROOT / "tests"
SOURCE_PREFIX = "packages/radia-mcp/src/radia_mcp/"
TEST_PREFIX = "packages/radia-mcp/tests/"

ALWAYS_TESTS = (
    "tests/test_ci_selection.py",
    "tests/test_collection_dependency_gate.py",
    "tests/test_mcp_naming_contract.py",
    "tests/test_mcp_sdk_dependency_contract.py",
    "tests/test_meta_health.py::test_meta_catalog_has_at_least_30_servers",
    "tests/test_meta_health.py::test_every_cataloged_server_has_register_status_tool",
    "tests/test_meta_health.py::test_meta_overview_returns_expected_shape",
)


def _test_file(selector: str) -> str:
    return selector.split("::", 1)[0]


def _normalize(path: str) -> str:
    normalized = path.strip().replace("\\", "/")
    while normalized.startswith("./"):
        normalized = normalized[2:]
    return normalized


@lru_cache(maxsize=1)
def _catalog() -> dict[str, dict]:
    namespace = runpy.run_path(
        str(PACKAGE_ROOT / "src" / "radia_mcp" / "meta" / "catalog.py")
    )
    return namespace["CATALOG"]


def _family_to_server(catalog: dict[str, dict]) -> dict[str, str]:
    return {
        str(info["subpackage"]).removeprefix("radia_mcp."): short
        for short, info in catalog.items()
    }


@lru_cache(maxsize=1)
def _test_sources() -> frozenset[str]:
    return frozenset(
        path.relative_to(PACKAGE_ROOT).as_posix()
        for path in TEST_ROOT.rglob("test_*.py")
    )


@lru_cache(maxsize=None)
def _tests_containing(token: str) -> frozenset[str]:
    matches = set()
    for relative in _test_sources():
        path = PACKAGE_ROOT / relative
        try:
            source = path.read_text(encoding="utf-8")
        except (OSError, UnicodeError):
            continue
        if token in source:
            matches.add(relative)
    return frozenset(matches)


def _related_tests(
    sources: frozenset[str],
    *,
    family: str,
    module_parts: tuple[str, ...],
    changed_symbols: set[str] | None = None,
) -> set[str]:
    package_token = f"radia_mcp.{family}"
    stem = module_parts[-1] if module_parts else ""
    broad_family_change = not module_parts or stem == "__init__"
    if broad_family_change:
        search_token = package_token
        hints = {family, family.replace("_", "-")}
    else:
        search_token = ".".join((package_token, *module_parts))
        hints = {stem, stem.replace("_", "-")}
        if stem == "server":
            family_hint = family.replace("-", "_")
            hints = {f"{family_hint}_mcp", f"{family_hint}_server"}

    if changed_symbols and "*" not in changed_symbols and stem != "server":
        selected: set[str] = set()
        for symbol in changed_symbols:
            selected.update(_tests_containing(symbol))
        return selected

    # Server modules are registration/composition boundaries. Searching every
    # direct server import selects nearly the whole historical physics-gate
    # corpus (radia_ngsolve alone exceeded 100 files). Dedicated MCP contract
    # tests plus the affected server selftest own that boundary; numerical
    # behavior is selected when its implementation module changes.
    selected = set() if stem == "server" else set(_tests_containing(search_token))
    for relative in sources:
        filename = Path(relative).name
        if any(hint and hint in filename for hint in hints):
            selected.add(relative)
    return selected


def _symbols_for_line_ranges(
    text: str, line_ranges: list[tuple[int, int]]
) -> set[str]:
    tree = ast.parse(text)
    symbol_ranges = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            start = min(
                [node.lineno, *(item.lineno for item in node.decorator_list)]
            )
            symbol_ranges.append((start, node.end_lineno, node.name))
    symbols: set[str] = set()
    for changed_start, changed_end in line_ranges:
        matches = {
            name
            for symbol_start, symbol_end, name in symbol_ranges
            if symbol_start <= changed_end and changed_start <= symbol_end
        }
        symbols.update(matches or {"*"})
    return symbols


def _changed_symbols_by_file(
    changed_files: Iterable[str], base: str
) -> dict[str, set[str]]:
    paths = sorted(
        {
            _normalize(path)
            for path in changed_files
            if _normalize(path).startswith(SOURCE_PREFIX)
            and _normalize(path).endswith(".py")
        }
    )
    if not paths:
        return {}
    completed = subprocess.run(
        [
            "git",
            "-C",
            str(REPO_ROOT),
            "diff",
            "--unified=0",
            "--no-color",
            base,
            "--",
            *paths,
        ],
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=60,
    )
    if completed.returncode != 0:
        return {}

    hunks: dict[str, list[tuple[int, int]]] = {}
    current_path = ""
    hunk_pattern = re.compile(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@")
    for line in completed.stdout.splitlines():
        if line.startswith("+++ "):
            value = line[4:].strip()
            current_path = "" if value == "/dev/null" else value.removeprefix("b/")
            continue
        match = hunk_pattern.match(line)
        if not current_path or match is None:
            continue
        start = int(match.group(1))
        count = int(match.group(2) or "1")
        if count == 0:
            hunks.setdefault(current_path, []).append((0, 0))
        else:
            hunks.setdefault(current_path, []).append((start, start + count - 1))

    result: dict[str, set[str]] = {}
    for path, line_ranges in hunks.items():
        if any(start == 0 for start, _ in line_ranges):
            result[path] = {"*"}
            continue
        current = REPO_ROOT / path
        try:
            result[path] = _symbols_for_line_ranges(
                current.read_text(encoding="utf-8"), line_ranges
            )
        except (OSError, SyntaxError, UnicodeError):
            result[path] = {"*"}
    return result


def build_plan(
    changed_files: Iterable[str],
    *,
    full: bool = False,
    changed_symbols_by_file: dict[str, set[str]] | None = None,
) -> dict:
    """Return the deterministic package-test and server-selftest selection."""

    changed = sorted({_normalize(path) for path in changed_files if path.strip()})
    catalog = _catalog()
    if full:
        return {
            "schema": "radia-mcp.ci-selection.v1",
            "mode": "full",
            "changed_files": changed,
            "package_tests": ["tests"],
            "server_selftests": sorted(catalog),
            "run_mcp_response_tests": True,
        }

    sources = _test_sources()
    selected = {path for path in ALWAYS_TESTS if _test_file(path) in sources}
    servers: set[str] = set()
    family_servers = _family_to_server(catalog)

    for path in changed:
        if path == "packages/radia-mcp/tests/conftest.py":
            selected.add("tests/test_ci_collection_scope.py")
            continue

        if path.startswith(TEST_PREFIX) and path.endswith(".py"):
            relative = path.removeprefix("packages/radia-mcp/")
            if relative in sources:
                selected.add(relative)
                selected.add("tests/test_validation_lane_separation.py")
            continue

        if path == "packages/radia-mcp/pyproject.toml":
            selected.update(
                {
                    "tests/test_mcp_sdk_dependency_contract.py",
                    "tests/test_optional_dependency_imports.py",
                    "tests/test_meta_health.py::test_readme_mcp_entrypoints_are_cataloged_or_external",
                }
            )
            continue

        if path.endswith("packages/radia-mcp/tools/policy_lint.py"):
            selected.add("tests/test_policy_lint.py")
            continue

        if not path.startswith(SOURCE_PREFIX):
            continue

        relative = path.removeprefix(SOURCE_PREFIX)
        parts = tuple(part for part in relative.split("/") if part)
        if len(parts) < 2:
            servers.update(catalog)
            continue

        family = parts[0]
        module_parts = parts[1:]
        if module_parts and module_parts[-1].endswith(".py"):
            module_parts = (*module_parts[:-1], module_parts[-1][:-3])

        stem = module_parts[-1] if module_parts else ""
        if family == "common":
            servers.update(catalog)
        elif family in family_servers:
            servers.add(family_servers[family])

        if family == "common" and stem == "__init__":
            selected.update(
                {
                    "tests/test_coarse_tool_registry.py",
                    "tests/test_optional_dependency_imports.py",
                }
            )
            continue

        selected.update(
            _related_tests(
                sources,
                family=family,
                module_parts=module_parts,
                changed_symbols=(changed_symbols_by_file or {}).get(path),
            )
        )

    whole_files = {selector for selector in selected if "::" not in selector}
    deduplicated = {
        selector
        for selector in selected
        if "::" not in selector or _test_file(selector) not in whole_files
    }

    return {
        "schema": "radia-mcp.ci-selection.v1",
        "mode": "targeted",
        "changed_files": changed,
        "package_tests": sorted(
            selector for selector in deduplicated if _test_file(selector) in sources
        ),
        "server_selftests": sorted(servers),
        "run_mcp_response_tests": any(
            path.startswith("tests/mcp_server/") for path in changed
        ),
    }


def _git_changed_files(base: str, head: str) -> list[str]:
    completed = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "diff", "--name-only", base, head],
        check=True,
        capture_output=True,
        text=True,
        timeout=60,
    )
    return completed.stdout.splitlines()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--changed-files-json")
    parser.add_argument("--base")
    parser.add_argument("--head", default="HEAD")
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.changed_files_json is not None:
        changed_files = json.loads(args.changed_files_json)
        if not isinstance(changed_files, list) or not all(
            isinstance(path, str) for path in changed_files
        ):
            parser.error("--changed-files-json must encode a list of strings")
    elif args.base:
        changed_files = _git_changed_files(args.base, args.head)
    elif args.full:
        changed_files = []
    else:
        parser.error("provide --changed-files-json, --base, or --full")

    changed_symbols = (
        _changed_symbols_by_file(changed_files, args.base)
        if args.base and not args.full
        else None
    )
    plan = build_plan(
        changed_files,
        full=args.full,
        changed_symbols_by_file=changed_symbols,
    )
    text = json.dumps(plan, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
