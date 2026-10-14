#!/usr/bin/env python
"""Verify that required native modules load from the requested source tree.

The check runs each import in a fresh process.  This isolates DLL loader state
between modules and turns an NGSolve/Netgen ABI mismatch into a concise build
failure instead of letting a stale ``.pyd`` reach tests or a release package.
"""

from __future__ import annotations

import argparse
import importlib
import importlib.machinery
import importlib.metadata
import json
import os
import subprocess
import sys
from pathlib import Path


class ProbeError(RuntimeError):
    """A native module failed its isolated import or provenance check."""


def _is_native_extension(path: Path) -> bool:
    name = path.name.lower()
    return any(
        name.endswith(suffix.lower())
        for suffix in importlib.machinery.EXTENSION_SUFFIXES
    )


def _probe_in_process(module_name: str, source_root: Path) -> dict[str, str]:
    source_root = source_root.resolve()
    sys.path.insert(0, str(source_root))
    module = importlib.import_module(module_name)
    module_file = Path(module.__file__).resolve()
    expected_parent = source_root.joinpath(*module_name.split(".")[:-1]).resolve()

    if module_file.parent != expected_parent:
        raise ProbeError(
            f"{module_name} resolved outside the requested source tree: "
            f"{module_file} (expected under {expected_parent})"
        )
    if not _is_native_extension(module_file):
        raise ProbeError(
            f"{module_name} resolved to {module_file}, which is not a native extension"
        )

    try:
        ngsolve_version = importlib.metadata.version("ngsolve")
    except importlib.metadata.PackageNotFoundError:
        ngsolve_version = "not-installed"
    return {
        "module": module_name,
        "path": str(module_file),
        "python": sys.version.split()[0],
        "ngsolve": ngsolve_version,
    }


def probe_module(
    module_name: str,
    source_root: Path,
    *,
    python_executable: str = sys.executable,
    timeout: float = 30.0,
) -> dict[str, str]:
    """Import one module in an isolated interpreter and return its inventory."""
    command = [
        python_executable,
        str(Path(__file__).resolve()),
        "--probe",
        "--source-root",
        str(source_root),
        "--module",
        module_name,
    ]
    environment = os.environ.copy()
    environment["PYTHONIOENCODING"] = "utf-8"
    existing = environment.get("PYTHONPATH")
    environment["PYTHONPATH"] = (
        str(source_root) if not existing else str(source_root) + os.pathsep + existing
    )
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=environment,
        timeout=timeout,
    )
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        if "DLL load failed" in detail:
            prefix = detail.split("DLL load failed", 1)[0]
            detail = prefix + "DLL load failed (ABI or native dependency mismatch)"
        raise ProbeError(
            f"{module_name} failed to load with {python_executable}: "
            f"{detail or f'exit code {completed.returncode}'}"
        )
    try:
        return json.loads(completed.stdout.strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError) as error:
        raise ProbeError(
            f"{module_name} produced an invalid ABI probe response: "
            f"{completed.stdout.strip()!r}"
        ) from error


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--source-root",
        required=True,
        type=Path,
        help="directory placed first on sys.path, normally <repo>/src",
    )
    parser.add_argument(
        "--module",
        action="append",
        required=True,
        help="fully qualified native module name; repeat for multiple modules",
    )
    parser.add_argument("--probe", action="store_true", help=argparse.SUPPRESS)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.probe:
        if len(args.module) != 1:
            raise SystemExit("--probe accepts exactly one --module")
        try:
            result = _probe_in_process(args.module[0], args.source_root)
        except Exception as error:  # noqa: BLE001 - probe boundary reports any import failure
            print(f"{type(error).__name__}: {error}", file=sys.stderr)
            return 1
        print(json.dumps(result, sort_keys=True))
        return 0

    failed = False
    for module_name in args.module:
        try:
            result = probe_module(module_name, args.source_root)
        except (ProbeError, subprocess.TimeoutExpired) as error:
            print(f"FAIL {error}", file=sys.stderr)
            failed = True
            continue
        print(
            f"PASS {result['module']}: {result['path']} "
            f"(Python {result['python']}, NGSolve {result['ngsolve']})"
        )
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
