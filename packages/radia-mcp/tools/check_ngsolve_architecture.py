"""Keep new gate responsibilities out of the legacy radia_ngsolve monolith.

Compare definitions, not line counts: bug fixes to existing gates remain
possible, while a new named responsibility must live in a domain module.
"""

from __future__ import annotations

import argparse
import ast
from pathlib import Path
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[3]
LEGACY_PATH = Path(
    "packages/radia-mcp/src/radia_mcp/radia_ngsolve/slot_gates.py"
)
DEFINITION_NODES = (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)


def definitions(source: str) -> set[str]:
    """Collect qualified function/class names, including nested helpers."""
    found: set[str] = set()

    def visit(node: ast.AST, scope: str = "") -> None:
        for child in ast.iter_child_nodes(node):
            if isinstance(child, DEFINITION_NODES):
                name = f"{scope}.{child.name}" if scope else child.name
                found.add(name)
                visit(child, name)
            else:
                visit(child, scope)

    visit(ast.parse(source))
    return found


def new_definitions(base_source: str, current_source: str) -> set[str]:
    return definitions(current_source) - definitions(base_source)


def _git_source(ref: str) -> str:
    result = subprocess.run(
        [
            "git", "-c", f"safe.directory={REPO_ROOT}", "-C", str(REPO_ROOT),
            "show", f"{ref}:{LEGACY_PATH.as_posix()}",
        ],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return result.stdout


def _selftest() -> None:
    base = "def existing():\n    return 1\n"
    assert not new_definitions(base, "def existing():\n    return 2\n")
    assert new_definitions(base, base + "\ndef added():\n    return 1\n") == {
        "added"
    }
    nested = "def existing():\n    def added():\n        return 1\n    return added()\n"
    assert new_definitions(base, nested) == {"existing.added"}
    assert not new_definitions(base, "")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", help="Git commit to compare with the worktree")
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        _selftest()
        print("ngsolve architecture checker selftest: PASS")
        if not args.base:
            return 0
    if not args.base:
        parser.error("--base is required unless running --selftest")

    before = _git_source(args.base)
    after = (REPO_ROOT / LEGACY_PATH).read_text(encoding="utf-8")
    added = sorted(new_definitions(before, after))
    if added:
        print(f"New definitions in legacy {LEGACY_PATH}:")
        for name in added:
            print(f"  {name}")
        print(
            "Put new gate logic in a cohesive radia_ngsolve domain module. "
            "Keep only a compatibility import or thin adapter in slot_gates.py."
        )
        return 1
    print("ngsolve legacy gate boundary: PASS (no new definitions)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
