"""Run the CI-equivalent source checks locally before submitting or merging.

The fast-contracts steps run in a fresh isolated virtual environment that holds
only the lightweight dependencies declared for that job, exactly as the light
CI runner does. Repository contracts and policy checks then run with the
caller's interpreter. Every step is recorded in a JSON report; the command exits
non-zero if any step fails. Nothing is installed into the caller's environment.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import tempfile
import time
import venv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIGHT_REQUIREMENTS = ["pytest", "pyyaml", "numpy", "scipy", "setuptools>=45", "wheel"]
# Caller settings that would change what runs or where packages land.
_DROPPED_PREFIXES = ("PYTEST_", "PIP_", "PYTHON")
_DROPPED_NAMES = {"VIRTUAL_ENV", "CONDA_PREFIX"}


def _clean_env(**extra):
    env = {key: value for key, value in os.environ.items()
           if not key.upper().startswith(_DROPPED_PREFIXES) and key.upper() not in _DROPPED_NAMES}
    env.update(PYTHONNOUSERSITE="1", **extra)
    return env


def _run(name, command, env, steps):
    started = time.perf_counter()
    try:
        result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True)
        returncode = result.returncode
        tail = (result.stdout + result.stderr).strip().splitlines()[-3:]
    except OSError as error:  # the step could not even start
        returncode, tail = 1, [repr(error)]
    steps.append({
        "name": name,
        "command": [str(part) for part in command],
        "returncode": returncode,
        "seconds": round(time.perf_counter() - started, 3),
        "tail": tail,
    })
    status = "PASS" if returncode == 0 else "FAIL"
    print(f"[{status}] {name}: {tail[-1] if tail else ''}", flush=True)
    return returncode == 0


def _light_python(directory, steps):
    started = time.perf_counter()
    try:
        venv.EnvBuilder(with_pip=True, clear=True).create(directory)
    except Exception as error:  # recorded, then reported as a failed step
        steps.append({"name": "create isolated venv", "command": [], "returncode": 1,
                      "seconds": round(time.perf_counter() - started, 3), "tail": [repr(error)]})
        print(f"[FAIL] create isolated venv: {error!r}", flush=True)
        return None
    scripts = "Scripts" if os.name == "nt" else "bin"
    python = Path(directory) / scripts / ("python.exe" if os.name == "nt" else "python")
    # --isolated ignores pip configuration files and PIP_* variables.
    ok = _run("install light requirements",
              [python, "-m", "pip", "--isolated", "install", "--quiet",
               "--disable-pip-version-check", *LIGHT_REQUIREMENTS], _clean_env(), steps)
    return python if ok else None


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--since", default="origin/main",
                        help="base revision for impact-selected fast contracts")
    parser.add_argument("--report", type=Path, default=ROOT / "ci-precheck-report.json",
                        help="JSON report path (default: ci-precheck-report.json at the repository root)")
    parser.add_argument("--skip-repo-contracts", action="store_true",
                        help="run only the fast-contracts steps")
    args = parser.parse_args(argv)

    steps = []
    with tempfile.TemporaryDirectory(prefix="radia-ci-precheck-") as scratch:
        light = _light_python(Path(scratch) / "venv", steps)
        env = _clean_env(RADIA_SOURCE_ONLY_CONTRACTS="1")
        fast = [] if light is None else [
            ("runner routing", [light, "tools/ci_runner_routing.py"]),
            ("routing contract", [light, "-m", "pytest", "tests/ci_routing/test_runner_routing.py",
                                  "--confcutdir=tests/ci_routing", "-q", "-o", "addopts="]),
            ("compile radia sources", [light, "-m", "compileall", "-q", "src/radia"]),
            ("TaskManager audit", [light, "tools/audit_taskmanager.py", "--distribution", "radia"]),
            ("CPU version", [light, "-m", "pytest",
                             "tests/test_version_consistency.py::test_pyproject_matches_init[radia]", "-q"]),
            ("fast contracts", [light, "tools/run_test_tier.py", "--profile", "fast-contracts",
                                "--since", args.since]),
            ("showcase readers", [light, "-m", "pytest", "tests/showcase/test_result_readers.py",
                                  "--confcutdir=tests/showcase", "-q", "-o", "addopts="]),
        ]
        if light is not None:
            for name, command in fast:
                _run(name, command, env, steps)

    caller = [
        ("policy lint", [sys.executable, "tools/policy_lint.py"]),
        ("policy tests", [sys.executable, "-m", "pytest", "tests/policy", "-q"]),
    ]
    if not args.skip_repo_contracts:
        caller.append(("repository contracts",
                       [sys.executable, "tools/run_test_tier.py", "--profile", "repo-contracts"]))
    for name, command in caller:
        _run(name, command, _clean_env(), steps)

    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    report = {
        "schema": "radia.ci-precheck.v1",
        "head": head,
        "since": args.since,
        "python": platform.python_version(),
        "passed": bool(steps) and all(step["returncode"] == 0 for step in steps),
        "steps": steps,
    }
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("PRECHECK", "PASS" if report["passed"] else "FAIL")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
