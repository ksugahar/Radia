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


OPTUNA_SOURCE_SHA = "274e312cf084e1975f32e4d9dbf35b63d3f4fcd7"


def optuna_source(environment: dict[str, str]) -> Path:
    configured = environment.get("RADIA_OPTUNA_SOURCE", "").strip()
    if not configured:
        raise ValueError("Set RADIA_OPTUNA_SOURCE to the pinned standalone radia-optuna checkout")
    source = Path(configured).resolve()
    sha = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"],
                                  text=True).strip()
    if sha != OPTUNA_SOURCE_SHA:
        raise ValueError(f"Wrong standalone Optuna source: actual={sha}, expected={OPTUNA_SOURCE_SHA}")
    dirty = subprocess.check_output(["git", "-C", str(source), "status",
                                    "--porcelain", "--untracked-files=all"], text=True)
    if dirty:
        raise ValueError("Pinned standalone Optuna checkout is not clean")
    return source


def source_environment(root: Path, inherited: dict[str, str]) -> dict[str, str]:
    root = root.resolve()
    env = dict(inherited)
    env.pop("RADIA_MCP_CI_SELECTION_JSON", None)
    source = optuna_source(env)
    env["RADIA_OPTUNA_SOURCE"] = str(source)
    env["PYTHONPATH"] = os.pathsep.join((str(root / "packages/radia-mcp/src"), str(source / "src")))
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
    package = optuna_source(os.environ)
    for module_name, expected in (
        ("radia_mcp", root / "packages/radia-mcp/src/radia_mcp/__init__.py"),
        ("radia_optuna", package / "src/radia_optuna/__init__.py"),
    ):
        module = importlib.import_module(module_name)
        actual = Path(module.__file__).resolve()
        if actual != expected.resolve():
            raise ValueError(f"Wrong integration source for {module_name}: actual={actual}, expected={expected}")
        imports[module_name] = str(actual)

    # Use the existing evidence verifier: no producer execution or source rescue.
    from radia_mcp.matlab.optuna_evidence import checked_evidence
    evidence = checked_evidence(package)
    evidence["distribution_manifest"] = package / "src/radia_optuna/manifest.json"
    records = {key: {"path": str(path), "sha256": hashlib.sha256(
        path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()}
        for key, path in evidence.items()}
    sha = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"],
                                  text=True).strip()
    return {"schema": "radia.integration-source.v1", "source_sha": sha, "optuna_source_sha": OPTUNA_SOURCE_SHA,
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
