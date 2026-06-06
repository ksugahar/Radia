#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""policy_lint -- enforce the Sugahara-lab PUBLISH BOUNDARY on radia-mcp.

Per CLAUDE.md "公開境界 (OSS / PyPI / GitHub)":

  * MCP servers that TARGET / WRAP a commercial tool (COMSOL / FEMM / JMAG)
    must NOT be public -- regardless of who wrote the code. "It's my own code"
    is not an exemption. This explicitly includes ``comsol_converter`` (the
    COMSOL-conversion MCP server). The FEMM (S:\\FEMM) and JMAG (S:\\JMAG)
    servers are likewise lab-internal.
  * Publishable = OPEN-system servers only (radia-ngsolve / cubit / gmsh /
    build123d / ...), plus their engineering models / helpers / knowledge /
    own lint.
  * Commercial tools are used INTERNALLY as a verification *benchmark*; their
    content / models / docs / bench numbers must not be mixed into public
    artifacts (public showcase is analytic-solution-led).

This is the executable form of that policy ("失敗からも lint に学ぶ"): it turns
the near-miss of shipping a commercial wrapper into a permanent, CI-enforceable
guard. Run it before any public commit / push / publish.

  ERROR  (exit 1): a commercial-wrapper server is wired for publication
                   (public ``[project.scripts]`` entry, or shipped in the wheel
                   because it is not excluded from ``packages.find``).
  WARN   (exit 0): the wrapper source still lives under the public ``src/`` tree
                   (relocate to a lab-private package for full GitHub
                   compliance), or a public subpackage carries a module named
                   for a commercial tool (review for content / bench leakage).

Usage:  python tools/policy_lint.py [--root <radia-mcp dir>] [--strict]
        (--strict promotes WARN to ERROR.)
"""
from __future__ import annotations

import argparse
import re
import sys
import tomllib
from pathlib import Path

# --- policy data -----------------------------------------------------------
# Subpackages that WRAP a commercial tool -> must never be published.
HARD_DENY = {
    "comsol_converter",   # COMSOL .mph -> NGSolve/Java converter (named in policy)
}
# Subpackages that CONTAIN commercial-tool content and need a human call on
# whether they belong in the public package (warn, don't hard-fail).
REVIEW = {
    "interop",            # radia<->COMSOL/MATLAB LiveLink interop tips
}
# Public modules whose FILENAME is dedicated to a commercial tool warrant a
# human review for content / bench-number leakage. Filename-based on purpose:
# a free-text token scan flags the legitimate benchmark *mentions* (".mph",
# "JMAG-Designer", "COMSOL multilingual RAG") that pervade an open, FEMM-parity
# tool -- those are allowed; only their content / models / bench numbers are not.
COMMERCIAL_NAME_RE = re.compile(r"(comsol|jmag)", re.IGNORECASE)


def _subpackage_of(target: str) -> str | None:
    """'radia_mcp.comsol_converter.server:main' -> 'comsol_converter'."""
    mod = target.split(":", 1)[0]
    parts = mod.split(".")
    if len(parts) >= 2 and parts[0] == "radia_mcp":
        return parts[1]
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description="radia-mcp publish-boundary lint")
    ap.add_argument("--root", default=str(Path(__file__).resolve().parents[1]),
                    help="radia-mcp package dir (contains pyproject.toml)")
    ap.add_argument("--strict", action="store_true",
                    help="promote warnings to errors")
    args = ap.parse_args()

    root = Path(args.root)
    pyproject = root / "pyproject.toml"
    src = root / "src" / "radia_mcp"
    errors: list[str] = []
    warns: list[str] = []

    with open(pyproject, "rb") as f:
        cfg = tomllib.load(f)

    scripts = cfg.get("project", {}).get("scripts", {})
    find = cfg.get("tool", {}).get("setuptools", {}).get("packages", {}).get("find", {})
    excludes = find.get("exclude", [])

    def _is_excluded(sub: str) -> bool:
        return any(sub in pat for pat in excludes)

    # 1) public console-script entry points must not point at a wrapper
    for name, target in scripts.items():
        sub = _subpackage_of(target)
        if sub in HARD_DENY:
            errors.append(
                f"[scripts] '{name} = {target}' publishes commercial-wrapper "
                f"subpackage '{sub}'. Remove from [project.scripts].")
        elif sub in REVIEW:
            warns.append(
                f"[scripts] '{name}' exposes REVIEW subpackage '{sub}' "
                f"(contains commercial content -- confirm it may be public).")

    # 2) wrapper subpackages must be excluded from the wheel (packages.find)
    for sub in sorted(HARD_DENY):
        if (src / sub).is_dir():
            if not _is_excluded(sub):
                errors.append(
                    f"[packaging] subpackage '{sub}' is shipped in the wheel "
                    f"(not in [tool.setuptools.packages.find] exclude). Add "
                    f"'radia_mcp.{sub}*' to exclude, or relocate it out of src/.")
            # 3) even if excluded from the wheel, source in the public repo leaks on GitHub
            warns.append(
                f"[source] commercial-wrapper '{sub}' source still lives under "
                f"the public src/ tree -> relocate to a lab-private package for "
                f"full public-GitHub compliance.")

    # 4) dedicated commercial-tool modules in PUBLIC subpackages -> review for
    #    content / bench-number leakage (benchmark *mentions* are fine, so this
    #    is filename-based, not a free-text token scan).
    public_subs = set()
    for target in scripts.values():
        s = _subpackage_of(target)
        if s and s not in HARD_DENY and s not in REVIEW:
            public_subs.add(s)
    for sub in sorted(public_subs):
        d = src / sub
        if not d.is_dir():
            continue
        for py in sorted(d.rglob("*.py")):
            if COMMERCIAL_NAME_RE.search(py.stem):
                warns.append(
                    f"[content] {py.relative_to(root)} is a commercial-tool "
                    f"module in public subpackage '{sub}' -- review for content "
                    f"/ bench-number leakage (benchmark mentions are OK).")

    # --- report ---
    print("=" * 70)
    print("radia-mcp publish-boundary lint (CLAUDE.md 公開境界)")
    print("=" * 70)
    if not errors and not warns:
        print("OK -- no publish-boundary issues.")
        return 0
    if errors:
        print(f"\nERRORS ({len(errors)}) -- block public commit/push/publish:")
        for e in errors:
            print(f"  [X] {e}")
    if warns:
        print(f"\nWARNINGS ({len(warns)}):")
        for w in warns:
            print(f"  [!] {w}")
    fail = bool(errors) or (args.strict and bool(warns))
    print("\nRESULT:", "FAIL" if fail else "PASS (warnings only)")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
