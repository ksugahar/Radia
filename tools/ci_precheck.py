"""Run impact-scoped Python source checks before specialist submission.

CI policy/repository contracts always run. Selected Radia checks use a fresh
light environment; selected MCP pytest paths use the isolated caller and this
checkout's source. --full executes every source step and full MCP pytest plan.
Native/platform/release and server stdio qualification remain separate CI gates.
Nothing is installed into the caller's environment.
"""

from __future__ import annotations

import argparse
import importlib.util
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
_integration_spec = importlib.util.spec_from_file_location(
    "mcp_integration_source", ROOT / "packages/radia-mcp/tools/integration_source.py")
integration_source = importlib.util.module_from_spec(_integration_spec)
_integration_spec.loader.exec_module(integration_source)
LIGHT_REQUIREMENTS = ["pytest", "pyyaml", "numpy", "scipy", "setuptools>=45", "wheel"]
# Caller settings that would change what runs or where packages land.
_DROPPED_PREFIXES = ("PYTEST_", "PIP_", "PYTHON")
_DROPPED_NAMES = {"VIRTUAL_ENV", "CONDA_PREFIX", "RADIA_MCP_CI_SELECTION_JSON",
                  "RADIA_SOURCE_ONLY_CONTRACTS", "RADIA_TESTS_ALLOW_PARTIAL",
                  "RADIA_MCP_FORCE_MINIMAL", "RADIA_REPO_ROOT"}


def _clean_env(**extra):
    env = {key: value for key, value in os.environ.items()
           if not key.upper().startswith(_DROPPED_PREFIXES) and key.upper() not in _DROPPED_NAMES}
    env.update(PYTHONNOUSERSITE="1", PYTHONDONTWRITEBYTECODE="1",
               RADIA_TESTS_ALLOW_PARTIAL="0", **extra)
    return env


LOG_DIR = None  # set by main(); full step output is written here


def _run(name, command, env, steps, cwd=None):
    started = time.perf_counter()
    try:
        child_env = dict(env, PYTHONUTF8="1", PYTHONIOENCODING="utf-8")
        result = subprocess.run(command, cwd=cwd or ROOT, env=child_env,
                                capture_output=True, text=True,
                                encoding="utf-8", errors="replace")
        returncode = result.returncode
        output = (result.stdout or "") + (result.stderr or "")
    except OSError as error:  # the step could not even start
        returncode, output = 1, repr(error)
    tail = output.strip().splitlines()[-3:]
    log = None
    if LOG_DIR is not None:
        LOG_DIR.mkdir(parents=True, exist_ok=True)
        slug = "".join(ch if ch.isalnum() else "-" for ch in name).strip("-")
        log = LOG_DIR / f"{len(steps) + 1:02d}-{slug}.log"
        log.write_text(output, encoding="utf-8", errors="replace")
    steps.append({
        "name": name,
        "cwd": str(cwd or ROOT),
        "command": [str(part) for part in command],
        "returncode": returncode,
        "seconds": round(time.perf_counter() - started, 3),
        "tail": tail,
        "log": str(log) if log else None,
    })
    status = "PASS" if returncode == 0 else "FAIL"
    message = f"[{status}] {name}: {tail[-1] if tail else ''}"
    # Preserve the UTF-8 log even when the parent's console cannot display it.
    console_encoding = sys.stdout.encoding or "utf-8"
    message = message.encode(console_encoding, errors="backslashreplace").decode(console_encoding)
    print(message, flush=True)
    return returncode == 0


def _light_python(directory, steps):
    started = time.perf_counter()
    try:
        venv.EnvBuilder(with_pip=True, clear=True).create(directory)
    except Exception as error:  # recorded, then reported as a failed step
        log = None
        if LOG_DIR is not None:
            LOG_DIR.mkdir(parents=True, exist_ok=True)
            log = LOG_DIR / f"{len(steps) + 1:02d}-create-isolated-venv.log"
            log.write_text(repr(error), encoding="utf-8", errors="replace")
        steps.append({"name": "create isolated venv", "command": [], "returncode": 1,
                      "seconds": round(time.perf_counter() - started, 3), "tail": [repr(error)],
                      "log": str(log) if log else None})
        print(f"[FAIL] create isolated venv: {error!r}", flush=True)
        return None
    scripts = "Scripts" if os.name == "nt" else "bin"
    python = Path(directory) / scripts / ("python.exe" if os.name == "nt" else "python")
    # --isolated ignores pip configuration files and PIP_* variables.
    ok = _run("install light requirements",
              [python, "-m", "pip", "--isolated", "install", "--quiet",
               "--disable-pip-version-check", *LIGHT_REQUIREMENTS], _clean_env(), steps)
    return python if ok else None


POLICY_TESTS = ("tests/policy/", "tests/test_sparsesolv_ci_contract.py")
POLICY_COMMANDS = (
    ("policy lint", ("tools/policy_lint.py",)),
    ("policy tests", ("-m", "pytest", *POLICY_TESTS, "-q")),
)
REPOSITORY_COMMANDS = (
    ("repository contracts", ("tools/run_test_tier.py", "--profile", "repo-contracts")),
    ("repository TaskManager audit", ("tools/audit_taskmanager.py",)),
    ("CI install audit", ("tools/audit_ci_no_system_install.py",)),
    ("repository versions", ("tools/ci_preflight.py", "--only", "version")),
)
MCP_STATIC_COMMANDS = (
    ("MCP publish lint selftest", ("packages/radia-mcp/tools/policy_lint.py", "--selftest")),
    ("MCP publish lint", ("packages/radia-mcp/tools/policy_lint.py",)),
    ("MCP architecture selftest", ("packages/radia-mcp/tools/check_ngsolve_architecture.py", "--selftest")),
    ("MCP compile", ("-m", "compileall", "-q", "packages/radia-mcp/src")),
)


def dirty_paths():
    return subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=all"], cwd=ROOT, text=True
    ).splitlines()


_CONTROL_INPUTS = {
    "tools/ci_precheck.py", "tools/ci_distribution_contract.py",
    "tools/ci_distribution_manifest.json", "tools/run_test_tier.py",
    "tests/test_tier_manifest.json", "tests/conftest.py",
    "tests/test_ci_precheck.py",
}


def _ci_only_lanes(changed, data=None):
    if data is None:
        spec = importlib.util.spec_from_file_location("precheck_lanes", ROOT / "tools/ci_distribution_contract.py")
        contract = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(contract)
        data = contract.load(ROOT)
    else:
        spec = importlib.util.spec_from_file_location("precheck_lanes", ROOT / "tools/ci_distribution_contract.py")
        contract = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(contract)
    return sorted(lane for lane, rule in data["lanes"].items()
                  if lane not in {"radia", "radia-mcp-matrix"}
                  and "pull_request" in rule["events"]
                  and (changed is None or any(contract.selected(
                      path, rule["events"]["pull_request"]) for path in changed)))


def build_source_plan(since, *, full=False, git=None):
    """Use the checked workflow projections, never a second path allowlist."""
    def run_git(*args):
        return subprocess.check_output(["git", *args], cwd=ROOT, text=True, stderr=subprocess.PIPE)
    git = git or run_git
    if full:
        base = git("merge-base", since, "HEAD").strip()
        if not base:
            raise ValueError("comparison base could not be resolved")
        return {"mode": "full", "base": base, "changed_paths": [],
                "radia": True, "mcp": True, "ci_only_lanes": _ci_only_lanes(None),
                "reason": "explicit full source audit"}
    base = git("merge-base", since, "HEAD").strip()
    if not base:
        raise ValueError("comparison base could not be resolved")
    changed = sorted(set(filter(None, git("diff", "--name-only", "-z",
                                           "--no-renames", base, "HEAD", "--").split("\0"))))
    spec = importlib.util.spec_from_file_location("precheck_contract", ROOT / "tools/ci_distribution_contract.py")
    contract = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(contract)
    data = contract.load(ROOT)
    retired_deleted = []
    deleted, previous = None, None
    for path in changed:
        try:
            contract.owner_rule(data, path)
        except ValueError:
            # A removed distribution may retire its ownership rules in this
            # same change. Only a Git-proven deletion may use checked-base
            # ownership; added/modified unknown inputs still fail closed.
            if deleted is None:
                deleted = set(filter(None, git("diff", "--name-only", "-z",
                    "--no-renames", "--diff-filter=D", base, "HEAD", "--").split("\0")))
            if path not in deleted:
                raise
            try:
                if previous is None:
                    previous = json.loads(git("show", f"{base}:tools/ci_distribution_manifest.json"))
                contract.owner_rule(previous, path)
            except (ValueError, KeyError, subprocess.CalledProcessError) as error:
                raise ValueError(f"deleted input has no checked-base ownership: {path}") from error
            retired_deleted.append(path)
    control = bool(retired_deleted) or any(path in _CONTROL_INPUTS or path.startswith(".github/workflows/")
                  for path in changed)
    def selected(lane):
        return control or any(contract.selected(path, data["lanes"][lane]["events"]["pull_request"])
                              for path in changed)
    return {"mode": "affected", "base": base, "changed_paths": changed,
            "retired_deleted_paths": retired_deleted,
            "radia": not changed or selected("radia"),
            "mcp": selected("radia-mcp-matrix"),
            "ci_only_lanes": _ci_only_lanes(changed, data),
            "reason": "shared execution control" if control else
                      "verified empty diff: prior source baseline" if not changed else
                      "checked distribution workflow projections"}


def _mcp_checks(plan, steps, output):
    """Execute the selector's Python test fields, using CI collection semantics."""
    artifacts = output.with_name(output.stem + "-logs")
    artifacts.mkdir(parents=True, exist_ok=True)
    selection_file = artifacts / "mcp-selection.json"
    command = [sys.executable, "packages/radia-mcp/tools/select_ci_tests.py",
               "--output", str(selection_file)]
    command += ["--full"] if plan["mode"] == "full" else ["--base", plan["base"]]
    if not _run("MCP test selection", command, _clean_env(), steps):
        return None
    selection = json.loads(selection_file.read_text(encoding="utf-8"))
    if (selection.get("schema") != "radia-mcp.ci-selection.v1"
            or not isinstance(selection.get("package_tests"), list)
            or not selection["package_tests"]
            or not isinstance(selection.get("integration_tests"), list)
            or type(selection.get("run_mcp_response_tests")) is not bool
            or not isinstance(selection.get("server_selftests"), list)):
        raise ValueError("invalid MCP selection: refusing an empty or malformed test plan")
    for field, prefix in (("package_tests", "tests"), ("integration_tests", "tests/mcp_integration")):
        for value in selection[field]:
            if (not isinstance(value, str) or "\\" in value or ":" in value.split("::", 1)[0]
                    or value.startswith("/") or ".." in value.split("/")
                    or not (value == prefix or value.startswith(prefix + "/"))):
                raise ValueError(f"unsafe MCP test selector: {value!r}")
    package = ROOT / "packages/radia-mcp"
    env = _clean_env(PYTHONPATH=str(package / "src"),
                     PYTEST_DISABLE_PLUGIN_AUTOLOAD="1")
    code = ("from pathlib import Path; import radia_mcp; "
            f"expected=Path({str(package / 'src/radia_mcp')!r}).resolve(); "
            "actual=Path(radia_mcp.__file__).resolve(); "
            "assert actual.parent == expected, (actual, expected); print(actual)")
    _run("MCP source identity", [sys.executable, "-c", code], env, steps)
    package_env = dict(env)
    if selection["package_tests"] != ["tests"]:
        package_env["RADIA_MCP_CI_SELECTION_JSON"] = json.dumps(selection["package_tests"])
    else:
        package_env.pop("RADIA_MCP_CI_SELECTION_JSON", None)
    # Probe the actual collection policy in a fresh child before pytest can
    # collect_ignore selected files. Missing dependencies never become PASS.
    dependency_report = artifacts / "mcp-not-collected.json"
    probe = (
        "import importlib.util,json,runpy; from pathlib import Path; "
        "required=['pytest','numpy','scipy','h5py']; "
        + ("required+=['markdown','pymdownx']; " if selection["integration_tests"] else "")
        + "missing=[name for name in required if importlib.util.find_spec(name) is None]; "
        "state=runpy.run_path('tests/conftest.py'); "
        "ignored=state['_DEPENDENCY_IGNORED']; "
        f"Path({str(dependency_report)!r}).write_text(json.dumps(dict(missing=missing,not_collected=ignored),indent=2)); "
        "print(dict(missing=missing,not_collected=ignored)); "
        "raise SystemExit(bool(missing or ignored))"
    )
    if not _run("MCP collection dependencies", [sys.executable, "-c", probe],
                package_env, steps, cwd=package):
        return selection
    for name, arguments in MCP_STATIC_COMMANDS:
        _run(name, [sys.executable, *arguments], env, steps)
    # --full still validates the checked base for this changed-definition gate.
    _run("MCP architecture diff", [sys.executable,
         "packages/radia-mcp/tools/check_ngsolve_architecture.py", "--base", plan["base"]], env, steps)
    _run("MCP selected package tests",
         [sys.executable, "-m", "pytest", "-v", "-r", "aR",
          "-m", "not slow", "tests"], package_env, steps, cwd=package)
    if selection["run_mcp_response_tests"]:
        _run("MCP response tests", [sys.executable, "-m", "pytest", "tests/mcp_server", "-q"],
             env, steps)
    if selection["integration_tests"]:
        integration_env = integration_source.source_environment(ROOT, env)
        if not _run("MCP integration source identity",
                    integration_source.identity_command(ROOT, artifacts / "integration-source.json"),
                    integration_env, steps):
            return selection
        _run("MCP selected integration tests",
             [sys.executable, "-m", "pytest", "-v", "-ra",
              "--confcutdir=tests/mcp_integration", *selection["integration_tests"]], integration_env, steps)
        _run("MCP package test boundary",
             [sys.executable, "packages/radia-mcp/tools/check_package_test_boundary.py",
              "--scratch-root", str(artifacts / "boundary-scratch")], env, steps)
    return selection


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--since", default="origin/main", help="Git base for impact selection")
    parser.add_argument("--full", action="store_true", help="all source steps and full MCP pytest plan")
    parser.add_argument("--plan-only", action="store_true", help="record selection without executing checks")
    parser.add_argument("--report", type=Path, default=ROOT / "ci-precheck-report.json")
    parser.add_argument("--skip-repo-contracts", action="store_true",
                        help="explicit partial check; report cannot qualify the full source gate")
    args = parser.parse_args(argv)
    # Child checks use package-specific working directories. Bind output paths
    # to the caller's directory before handing them to any child process.
    args.report = args.report.resolve()
    if args.full and args.skip_repo_contracts:
        parser.error("--full cannot omit repository contracts")
    global LOG_DIR
    LOG_DIR = args.report.with_name(args.report.stem + "-logs")
    steps, errors = [], []
    plan, selection = None, None
    dirty = []
    try:
        dirty = dirty_paths()
        if dirty:
            raise ValueError("dirty working tree: commit all tracked/untracked inputs before precheck")
        plan = build_source_plan(args.since, full=args.full)
    except Exception as error:
        detail = str(error)
        if isinstance(error, subprocess.CalledProcessError):
            diagnostic = error.stderr or error.stdout
            if diagnostic:
                if isinstance(diagnostic, bytes):
                    diagnostic = diagnostic.decode("utf-8", errors="replace")
                detail += f": {diagnostic.strip()}"
        errors.append(f"selection failed: {type(error).__name__}: {detail}")
    if plan is not None and not args.plan_only:
        if plan["radia"]:
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
                                        *([] if args.full else ["--since", plan["base"]])]),
                    ("showcase readers", [light, "-m", "pytest", "tests/showcase/test_result_readers.py",
                                          "--confcutdir=tests/showcase", "-q", "-o", "addopts="]),
                ]
                for name, command in fast:
                    _run(name, command, env, steps)
        caller = [(name, [sys.executable, *arguments])
                  for name, arguments in POLICY_COMMANDS]
        if not args.skip_repo_contracts:
            caller += [(name, [sys.executable, *arguments])
                       for name, arguments in REPOSITORY_COMMANDS]
        for name, command in caller:
            _run(name, command, _clean_env(), steps)
        if plan["mcp"]:
            try:
                selection = _mcp_checks(plan, steps, args.report)
            except Exception as error:
                errors.append(f"MCP checks failed: {type(error).__name__}: {error}")
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                          capture_output=True, text=True).stdout.strip()
    passed = bool(steps) and not errors and all(step["returncode"] == 0 for step in steps)
    report = {
        "schema": "radia.ci-precheck.v2", "head": head, "since": args.since,
        "python": platform.python_version(), "interpreter": sys.executable,
        "passed": None if (args.plan_only and not errors) or (args.skip_repo_contracts and passed) else passed,
        "qualified": passed and not args.skip_repo_contracts,
        "status": "failed" if errors or (not args.plan_only and not passed) else
                  "planned" if args.plan_only else "partial" if args.skip_repo_contracts else "passed",
        "dirty_paths": dirty,
        "ci_only_lanes": plan.get("ci_only_lanes", []) if plan else [],
        "ci_only_lane_note": "Distribution/platform qualification remains CI-only; policy-lint source commands are also executed locally.",
        "plan": plan, "steps": steps, "errors": errors, "mcp_selection": selection,
        "mcp_ci_only_server_selftests": selection["server_selftests"] if selection else [],
        "scope": "Python 3.12 source precheck; native/platform/release and server stdio qualification remain CI",
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for error in errors:
        print(f"[FAIL] {error}", flush=True)
    for step in steps:
        if step["returncode"]:
            print(f"[FAIL] step={step['name']} rc={step['returncode']} log={step.get('log')}", flush=True)
    print("PRECHECK", report["status"].upper())
    return 0 if report["status"] in {"passed", "partial", "planned"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
