"""Bind repository integration checks to one checkout, never an installed wheel.

This test-runner helper does not change normal MCP runtime package selection.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys


def source_environment(root: Path, inherited: dict[str, str]) -> dict[str, str]:
    root = root.resolve()
    env = dict(inherited)
    env.pop("RADIA_MCP_CI_SELECTION_JSON", None)
    env["PYTHONPATH"] = os.pathsep.join(str(root / "packages" / package / "src")
                                       for package in ("radia-mcp", "radia-optuna"))
    env["PYTHONNOUSERSITE"] = "1"
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
    env["RADIA_REPO_ROOT"] = str(root)
    return env


def identity_command(root: Path, report: Path) -> list[str]:
    return [sys.executable, str(Path(__file__).resolve()), "--root", str(root.resolve()),
            "--report", str(report.resolve())]


def checked_identity(root: Path) -> dict:
    root = root.resolve()
    imports = {}
    for module_name, distribution in (("radia_mcp", "radia-mcp"),
                                      ("radia_optuna", "radia-optuna")):
        module = importlib.import_module(module_name)
        actual = Path(module.__file__).resolve()
        expected = root / "packages" / distribution / "src" / module_name / "__init__.py"
        if actual != expected.resolve() or not actual.is_relative_to(root):
            raise ValueError(f"Wrong integration source for {module_name}: actual={actual}, expected={expected}")
        imports[module_name] = str(actual)

    # Use the existing evidence verifier: no producer execution or source rescue.
    from radia_mcp.matlab.optuna_evidence import checked_evidence
    package = root / "packages/radia-optuna"
    evidence = checked_evidence(package)
    evidence["distribution_manifest"] = package / "src/radia_optuna/manifest.json"
    records = {key: {"path": str(path), "sha256": hashlib.sha256(
        path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()}
        for key, path in evidence.items()}
    sha = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"],
                                  text=True).strip()
    return {"schema": "radia.integration-source.v1", "source_sha": sha,
            "interpreter": sys.executable, "imports": imports,
            "hash_basis": "utf8-lf-bytes", "evidence": records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    # An exception exits nonzero before any integration test is launched.
    result = checked_identity(args.root)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
