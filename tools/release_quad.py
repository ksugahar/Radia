#!/usr/bin/env python
"""release_quad.py — CPU numerical Radia solver two-machine gate.

Walks the numerical solver release phases in order, gating each
phase on the success of the previous one. Refuses to skip steps that
have caused real outages (2026-04-14 incident series).

Usage:
    python tools/release_quad.py preflight
        Read-only: report current state and consistency. Use anytime.

    python tools/release_quad.py phase8 [--target validation|100|all]
        Install the bound Radia source editable on either operating role,
        with metadata-pinned GPU dependencies and same-source native payloads.
        Exact-wheel QA runs in separate explicitly configured environments.
        Never change radia-mcp or cubit-mesh-export.

    python tools/release_quad.py phase8e
        Reject the retired compute refresh route without contacting a host.

    python tools/release_quad.py phase9
        Cross-machine consistency probe. Final gate.

    python tools/release_quad.py optuna-candidate --ci-run-id <id> --target all
        Download the exact radia-optuna wheel from one successful main CI run
        and execute its installed-wheel MATLAB/Simulink test on both machines.

    python tools/release_quad.py optuna-done --wheel <path>
        Require the retained wheel bytes, source commit, CI run, version, and
        two-machine candidate state to agree before tagging/publication.

    python tools/release_quad.py all
        phase8 -> phase9 with all preconditions enforced.

        Both normal roles require explicit host, Python, source and release
        commit bindings; isolated QA interpreters are configured separately.

    python tools/release_quad.py done [--release-source <path>]
        Require the CPU solver release gate. The release source is the controller (which must then
        sit at the tag) or a separate exact-tag checkout; with the latter the
        controller must be a tracked-clean descendant declaring the same
        version. No runtime is changed.

    python tools/release_quad.py verify-editable
        Read-only: verify both bound normal editable installations.

    python tools/release_quad.py repoint --host <role> --package radia --source <path> --reason "..."
        Explicitly update one bound normal Radia source and verify it before
        recording intent. --record-current adopts an already verified pointer.
        Updates are forward-only; no process is stopped or package uninstalled.

    python tools/release_quad.py restore-editable
        Removed. The development tier has no default target; use repoint
        with a reason (policy P02, 2026-09-11).

The independent radia-optuna lane does not install or deploy Radia, Cubit,
radia-mcp, or NGSolve. It gates one exact CI wheel before publication.

Exit codes:
    0   success
    2   precondition failure (skip detected)
    3   action failure (external command)
    4   verification mismatch
    5   unverified: no recorded editable intent and no explicit expectation
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import importlib.util
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
import tomllib
import zipfile
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen
from urllib.request import url2pathname


def _load_editable_intent():
    """Load tools/editable_intent.py beside this file; its text is shipped to remotes."""
    source_path = Path(__file__).resolve().parent / "editable_intent.py"
    spec = importlib.util.spec_from_file_location("radia_editable_intent", source_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, source_path.read_text(encoding="utf-8")


editable_intent, _EDITABLE_INTENT_SOURCE = _load_editable_intent()


def gh_get(path):
    """Return one GitHub REST JSON response, using a token when available."""
    url = path if path.startswith("https://") else (
        "https://api.github.com/" + path.lstrip("/")
    )
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "radia-release-quad",
        "X-GitHub-Api-Version": "2022-11-28",
    }
    token = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"
    request = Request(url, headers=headers)
    try:
        with urlopen(request, timeout=30) as response:
            return json.loads(response.read()), dict(response.headers.items())
    except HTTPError as exc:
        raise RuntimeError(f"GitHub API HTTP {exc.code}: {exc.reason}") from exc
    except URLError as exc:
        raise RuntimeError(f"GitHub API network error: {exc.reason}") from exc

# Force UTF-8 stdout/stderr so en/em dashes and CJK in messages do not
# crash the script on ja-JP cp932 consoles.
if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

GIT_EXE = shutil.which("git")
if GIT_EXE is None:
    raise RuntimeError("Git is required by tools/release_quad.py")

REPO = Path(__file__).resolve().parent.parent
EDITABLE_REPO_100_ENV = "RADIA_RELEASE_EDITABLE_REPO_100"
# The three repository packages installed editable on the development tier.
EDITABLE_PACKAGES = ("radia", "cubit-mesh-export", "radia-mcp")
# Resolve 100号機 through the machine SSH configuration.  Its LAN address may
# change between lab network segments, while the supported host alias remains
# stable and carries the correct user/key settings.
SSH_100 = "100"
# 100号機's editable lives only in this dedicated development venv, separate
# from the machine's release runtime (host policy, 2026-10-04).
MATLAB_EXE = r"C:\Program Files\MATLAB\R2026a\bin\matlab.exe"
SIMULINK_GATE_ROOT = Path(r"C:\temp\radia-release-quad")
OPTUNA_GATE_ROOT = SIMULINK_GATE_ROOT / "radia-optuna"
OPTUNA_SUCCESS_MARKER = "RADIA_OPTUNA_WHEEL_SIMULINK_OK"
OPTUNA_TARGETS = {"100": ("build", SSH_100, "python")}
OPTUNA_WHEEL_RUNNER = REPO / "packages/radia-optuna/tests/run_installed_wheel_simulink.ps1"
_KNOWN_SOLVER_TARGETS = {
    "validation": ("validation", None, None),
    "100": ("build", SSH_100, "python"),
}
# The promotion gate imports the same tuple, so the machines that produce
# release evidence and the machines the gate demands evidence from cannot
# drift apart again.  release_acceptance.py is the single owner of the list.
try:
    from release_acceptance import RELEASE_ACCEPTANCE_HOSTS
    from release_acceptance import RELEASE_ACCEPTANCE_HOST_LABELS
except ModuleNotFoundError:  # Loaded by path in unit tests rather than run as a script.
    _acceptance_spec = importlib.util.spec_from_file_location(
        "radia_release_acceptance", Path(__file__).resolve().with_name("release_acceptance.py"))
    _acceptance_module = importlib.util.module_from_spec(_acceptance_spec)
    _acceptance_spec.loader.exec_module(_acceptance_module)
    RELEASE_ACCEPTANCE_HOSTS = _acceptance_module.RELEASE_ACCEPTANCE_HOSTS
    RELEASE_ACCEPTANCE_HOST_LABELS = _acceptance_module.RELEASE_ACCEPTANCE_HOST_LABELS

SOLVER_TARGETS = {host: _KNOWN_SOLVER_TARGETS[host]
                    for host in RELEASE_ACCEPTANCE_HOSTS}

assert tuple(SOLVER_TARGETS) == RELEASE_ACCEPTANCE_HOSTS, (
    "release_quad targets and the release acceptance hosts have drifted")
# The promotion gate checks each stored row label against the same table.
assert {key: row[0] for key, row in SOLVER_TARGETS.items()} == RELEASE_ACCEPTANCE_HOST_LABELS, (
    "Simulink target labels and the release acceptance host labels have drifted")


VALIDATION_SSH_ENV = "RADIA_RELEASE_SSH_VALIDATION"
VALIDATION_PYTHON_ENV = "RADIA_RELEASE_PYTHON_VALIDATION"


def _solver_target(role):
    """Resolve private normal-runtime routing without a PATH/host fallback.

    The selected absolute Python may be standard Python or an owned runtime venv.
    Reading configuration never provisions or changes that runtime.
    """
    label = SOLVER_TARGETS[role][0]
    suffix = role.upper()
    ssh_env = "RADIA_RELEASE_SSH_" + suffix
    python_env = "RADIA_RELEASE_PYTHON_" + suffix
    host = os.environ.get(ssh_env, "")
    python = os.environ.get(python_env, "")
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", host):
        raise ValueError(ssh_env + " must name one configured SSH alias")
    _absolute_windows_path(python, python_env, python=True)
    return label, host, python


def _absolute_windows_path(value, setting, *, python=False):
    from pathlib import PureWindowsPath
    path = PureWindowsPath(value)
    if (not value or not path.is_absolute() or (python and path.name.lower() != "python.exe")
            or any(ord(c) < 32 or c == '"' for c in value)
            or "windowsapps" in [part.lower() for part in path.parts]):
        raise ValueError(setting + " must name an absolute non-WindowsApps " + ("python.exe" if python else "source path"))
    return value


def _editable_target(role):
    label, ssh, python = _solver_target(role)
    suffix = role.upper()
    hostname = os.environ.get("RADIA_RELEASE_HOSTNAME_" + suffix, "")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", hostname):
        raise ValueError("RADIA_RELEASE_HOSTNAME_" + suffix + " must pin the actual target hostname")
    source_env = "RADIA_RELEASE_EDITABLE_REPO_" + suffix
    source = _absolute_windows_path(os.environ.get(source_env, ""), source_env)
    commit = os.environ.get("RADIA_RELEASE_EDITABLE_COMMIT_" + suffix, "")
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("RADIA_RELEASE_EDITABLE_COMMIT_" + suffix + " must pin the source commit")
    aliases = json.loads(os.environ.get("RADIA_RELEASE_EDITABLE_PATH_ALIASES_" + suffix, "[]"))
    if not isinstance(aliases, list) or any(not isinstance(pair, list) or len(pair) != 2
            or not all(isinstance(x, str) and x for x in pair) for pair in aliases):
        raise ValueError("Private editable path aliases must be JSON prefix pairs")
    binding = dict(schema=editable_intent.BINDING_SCHEMA, role=label, hostname=hostname,
                   python=python, sources={"radia": dict(source=source, commit=commit)})
    return dict(role=role, label=label, ssh=ssh, python=python, source=source,
                commit=commit, binding=binding, aliases=aliases)


def _solver_qa_target(role):
    target = _editable_target(role)
    setting = "RADIA_RELEASE_QA_PYTHON_" + role.upper()
    python = _absolute_windows_path(os.environ.get(setting, ""), setting, python=True)
    import ntpath
    if ntpath.normcase(ntpath.normpath(python)) == ntpath.normcase(ntpath.normpath(target["python"])):
        raise ValueError(setting + " must be separate from the normal editable runtime")
    return target["label"], target["ssh"], python


def _target_check_code(target, *, python=None, source=True, qa=False):
    """Read-only checks run in the actual target interpreter before mutation."""
    values = dict(hostname=target["binding"]["hostname"], python=python or target["python"], normal_python=target["python"],
                  source=target["source"], commit=target["commit"], check_source=source, qa=qa)
    token = base64.b64encode(json.dumps(values).encode()).decode("ascii")
    return '''import base64,json,os,pathlib,platform,subprocess,sys
t=json.loads(base64.b64decode("TOKEN"))
if platform.node().casefold()!=t["hostname"].casefold(): raise SystemExit("Selected hostname mismatch")
expected=pathlib.Path(t["python"]).resolve()
if "windowsapps" in [p.casefold() for p in expected.parts]: raise SystemExit("WindowsApps runtime refused")
if os.path.normcase(os.path.realpath(sys.executable))!=os.path.normcase(str(expected)): raise SystemExit("Selected Python mismatch")
if t["qa"]:
 if sys.prefix==sys.base_prefix: raise SystemExit("Exact-wheel QA requires an isolated venv")
 normal_prefix=subprocess.check_output([t["normal_python"],"-I","-c","import sys;print(sys.prefix)"],text=True).strip()
 if pathlib.Path(normal_prefix).resolve()==pathlib.Path(sys.prefix).resolve(): raise SystemExit("QA and normal runtimes share one environment")
if t["check_source"]:
 root=pathlib.Path(t["source"]).resolve(strict=True)
 cmd=["git","-c","safe.directory="+root.as_posix(),"-C",str(root)]
 head=subprocess.check_output(cmd+["rev-parse","HEAD"],text=True).strip().lower()
 if head!=t["commit"]: raise SystemExit("Selected source commit mismatch")
 dirty=subprocess.check_output(cmd+["status","--porcelain","--untracked-files=no"],text=True)
 if dirty.strip(): raise SystemExit("Selected source has tracked changes")
 import tomllib
 if tomllib.loads((root/"pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]!="radia": raise SystemExit("Selected source is not Radia")
'''.replace("TOKEN", token)


def _ps_python(command):
    return "& " + " ".join(_powershell_literal(str(word)) for word in command)


def _target_check_powershell(target, *, python=None, source=True, qa=False):
    python = python or target["python"]
    code = _target_check_code(target, python=python, source=source, qa=qa)
    token = base64.b64encode(code.encode()).decode("ascii")
    return (_runtime_guard_powershell(python) +
            _ps_python([python, "-I", "-c", f"import base64;exec(base64.b64decode('{token}'))"]) +
            "\nif ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }\n")


def _check_remote_target(target, *, python=None, source=True, qa=False):
    script = _target_check_powershell(target, python=python, source=source, qa=qa) + "exit $LASTEXITCODE"
    encoded = base64.b64encode(script.encode("utf-16le")).decode("ascii")
    return run(["ssh", target["ssh"], "pwsh", "-NoProfile", "-EncodedCommand", encoded], check=False).returncode


def _powershell_literal(value):
    return "'" + value.replace("'", "''") + "'"


def _runtime_guard_powershell(python_path):
    """Refuse missing/redirected runtimes before install or verification."""
    literal = _powershell_literal(python_path)
    code = ("import ntpath,sys; assert ntpath.normcase(ntpath.abspath(sys.executable)) == "
            + "ntpath.normcase(ntpath.abspath(" + repr(python_path) + ")), 'Selected Python identity mismatch'")
    token = base64.b64encode(code.encode('utf-8')).decode('ascii')
    return (f"if (-not (Test-Path -LiteralPath {literal} -PathType Leaf)) {{ throw 'Selected Python is missing' }}\n"
            f"& {literal} -c \"import base64; exec(base64.b64decode('{token}'))\"\n"
            "if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }\n")


SIMULINK_TARGETS = {"100": SOLVER_TARGETS["100"]}


def _editable_repo_100():
    """Return the explicitly selected build-role source, never a historic path."""
    return _editable_target("100")["source"]


def _radia_version() -> str:
    return _read_repo_versions({"radia"})["radia"]


def _release_commit() -> str:
    """Return the peeled commit for the version being deployed."""
    version = _radia_version()
    result = _git("rev-list", "-n", "1", f"v{version}", check=False)
    commit = result.stdout.strip().lower()
    if result.returncode or len(commit) != 40:
        raise RuntimeError(f"annotated or lightweight release tag v{version} is missing")
    return commit



def _release_head():
    return _release_commit()


def cmd_deployment_plan(args):
    """Read-only normal editable routing; QA runtimes remain separate."""
    plans = []
    for role in SOLVER_TARGETS:
        target = _editable_target(role)
        plans.append(dict(host=role, solver_runtime="editable", python=target["python"],
                          source=target["source"], source_commit=target["commit"],
                          solver_action="verify-and-install-radia-gpu", isolated_wheel_qa_required=True))
    print(json.dumps(plans, indent=2))
    return 0


def _color(c, s):
    return f"\033[{c}m{s}\033[0m"


def info(msg):  print("  " + msg)
def ok(msg):    print("  " + _color("32;1", "[OK]    ") + msg)
def warn(msg):  print("  " + _color("33;1", "[WARN]  ") + msg)
def fail(msg):  print("  " + _color("31;1", "[FAIL]  ") + msg)
def step(msg):  print("\n" + _color("36;1", f"=== {msg} ==="))


def run(cmd, *, check=True, capture=False, shell=False, **kw):
    """Subprocess wrapper that prints what it runs."""
    print("  $ " + (cmd if shell else " ".join(str(c) for c in cmd)))
    p = subprocess.run(cmd, shell=shell, capture_output=capture,
                        text=True, **kw)
    if check and p.returncode != 0:
        fail(f"command failed (exit {p.returncode})")
        if capture and p.stderr:
            print(p.stderr.strip())
        sys.exit(3)
    return p


def copy_file(src: Path, dst: Path):
    """Copy a file without relying on POSIX cp being on PATH."""
    print(f"  $ copy-file {src} {dst}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


# ============================================================
# state inspectors
# ============================================================

def _read_repo_versions(labels=None):
    """Parse versions out of pyproject.toml + __init__.py (no toml dep)."""
    out = {}
    import re
    for label, path in [
        ("radia",             REPO / "pyproject.toml"),
        ("radia.__version__", REPO / "src/radia/__init__.py"),
        ("cubit-mesh-export", REPO / "packages/cubit-mesh-export/pyproject.toml"),
        ("cme.__version__",   REPO / "packages/cubit-mesh-export/src/cubit_mesh_export/__init__.py"),
        ("radia-mcp",         REPO / "packages/radia-mcp/pyproject.toml"),
        ("radia-optuna",      REPO / "packages/radia-optuna/pyproject.toml"),
        ("optuna.__version__", REPO / "packages/radia-optuna/src/radia_optuna/__init__.py"),
    ]:
        if labels is not None and label not in labels:
            continue
        text = path.read_text(encoding="utf-8")
        m = re.search(r'(?:^version|^__version__)\s*=\s*"([^"]+)"', text, re.M)
        out[label] = m.group(1) if m else None
    return out


def _solver_dependency_pins() -> dict[str, str]:
    """Return the exact native solver ABI pins owned by pyproject.toml."""
    with (REPO / "pyproject.toml").open("rb") as stream:
        dependencies = tomllib.load(stream)["project"]["dependencies"]
    pins = {}
    for dependency in dependencies:
        for package in ("ngsolve", "netgen-mesher"):
            prefix = f"{package}=="
            if dependency.startswith(prefix):
                pins[package] = dependency.removeprefix(prefix)
    if set(pins) != {"ngsolve", "netgen-mesher"}:
        raise ValueError("pyproject.toml must exactly pin ngsolve and netgen-mesher")
    return pins


def _solver_hdf5_requirement() -> str:
    """Read the required HDF5 reader dependency from package metadata."""
    with (REPO / "pyproject.toml").open("rb") as stream:
        dependencies = tomllib.load(stream)["project"]["dependencies"]
    requirements = [item for item in dependencies if item.startswith("h5py")]
    if len(requirements) != 1 or not re.fullmatch(r"h5py[<>=!0-9.,]+", requirements[0]):
        raise ValueError("pyproject.toml must declare one bounded h5py requirement")
    return requirements[0]


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _simulink_manifest(package: Path) -> dict:
    verifier = REPO / "tools/verify_simulink_release.py"
    result = subprocess.run(
        [sys.executable, str(verifier), str(package), "--manifest-only"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or result.stdout.strip())
    with zipfile.ZipFile(package) as bundle:
        return json.loads(bundle.read("manifest.json"))


def _simulink_state_path(package_sha256: str) -> Path:
    return SIMULINK_GATE_ROOT / f"simulink-{package_sha256}.json"


def _simulink_mex_sha256(manifest: dict) -> str:
    matches = [row.get("sha256") for row in manifest.get("files", [])
               if row.get("path") == "packages/radia-simulink/matlab/radia_mex.mexw64"]
    if len(matches) != 1 or not isinstance(matches[0], str) or len(matches[0]) != 64:
        raise ValueError("Simulink manifest must identify exactly one radia_mex hash")
    return matches[0].lower()


def _release_tag_commit(version: object) -> str | None:
    """Return the peeled commit for the public ``v<version>`` tag, if any."""
    if not isinstance(version, str) or not version.strip():
        return None
    result = _git(
        "rev-parse", "--verify", f"refs/tags/v{version}^{{}}", check=False
    )
    if result.returncode != 0:
        return None
    return result.stdout.strip().lower()


def _verify_head_release_tag() -> int:
    """Require the declared version's ``v<version>`` tag to name the release commit."""
    version = _read_repo_versions()["radia"]
    head = _release_head()
    tag_commit = _release_tag_commit(version)
    if tag_commit is None:
        fail(f"Radia {version} has no local v{version} release tag")
        return 4
    if tag_commit != head:
        fail(
            f"release commit {head} is not the v{version} tag commit {tag_commit}; "
            "release_quad done must verify the exact tagged source"
        )
        return 4
    ok(f"release commit is anchored at v{version} ({head[:12]})")
    return 0


def _simulink_candidate_commit_is_release_anchored(
        manifest: dict, head: str) -> tuple[bool, str]:
    """Accept an exact tagged candidate when release-tooling later advances main.

    The archive bytes are built and tested from the tagged release source.  The
    definition-of-done controller can acquire a later release-only repair, so
    requiring the archive commit to equal the controller's ``HEAD`` would reject
    the already-tested public artifact.  A candidate is therefore accepted only
    when its commit is the peeled ``v<version>`` tag and remains an ancestor of
    the controller commit.
    """
    commit = manifest.get("commit")
    version = manifest.get("radia_version") or manifest.get("version")
    if not isinstance(commit, str) or not commit.strip():
        return False, "Simulink manifest has no source commit"
    commit = commit.strip().lower()
    tag_commit = _release_tag_commit(version)
    if tag_commit is None:
        return False, f"Simulink manifest version {version!r} has no v<version> release tag"
    if commit != tag_commit:
        return False, "Simulink manifest commit does not match its peeled release tag"
    ancestry = _git("merge-base", "--is-ancestor", commit, head, check=False)
    if ancestry.returncode != 0:
        return False, "Simulink tagged commit is not an ancestor of the release controller HEAD"
    return True, f"Simulink candidate is anchored at v{version} ({commit[:12]})"


def _write_simulink_state(path: Path, state: dict, target: str) -> None:
    """Merge only the completed target under a process lock, then publish atomically."""
    path.parent.mkdir(parents=True, exist_ok=True)
    # Keep this lock file: unlinking it could split concurrent writers' locks.
    with path.with_suffix(path.suffix + '.lock').open('a+b') as lock:
        deadline = time.monotonic() + 30
        while True:
            try:
                if os.name == 'nt':
                    import msvcrt
                    lock.seek(0)
                    msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
                else:
                    import fcntl
                    fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except OSError:
                if time.monotonic() >= deadline:
                    raise TimeoutError(f'Timed out locking candidate state: {path}')
                time.sleep(0.05)
        temporary = None
        try:
            previous = json.loads(path.read_text(encoding='utf-8')) if path.exists() else {}
            if previous:
                for field in ('schema', 'package_sha256', 'wheel_sha256', 'commit',
                              'ci_run_id', 'mex_sha256'):
                    if previous.get(field) != state.get(field):
                        raise ValueError(f'Candidate state identity mismatch: {field}')
            merged = {**state, 'targets': {**previous.get('targets', {})}}
            merged['targets'][target] = state['targets'][target]
            with tempfile.NamedTemporaryFile(
                    mode='w', encoding='utf-8', dir=path.parent,
                    prefix=path.name + '.', suffix='.tmp', delete=False) as output:
                temporary = Path(output.name)
                json.dump(merged, output, indent=2)
                output.write('\n')
                output.flush()
                os.fsync(output.fileno())
            os.replace(temporary, path)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
            if os.name == 'nt':
                lock.seek(0)
                msvcrt.locking(lock.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def _run_ssh_powershell(host: str, script: str, timeout: int):
    encoded = base64.b64encode(script.encode("utf-16le")).decode("ascii")
    command = ["ssh", "-v", "-o", "BatchMode=yes", "-o", "ConnectTimeout=15",
               host, "pwsh", "-NoProfile",
               "-EncodedCommand", encoded]
    # Use a file rather than inherited output pipes: Windows SSH/Engine
    # descendants can keep a pipe open after the command has exited.
    with tempfile.TemporaryDirectory(prefix="radia-ssh-", dir=tempfile.gettempdir()) as scratch:
        with (Path(scratch) / "output.log").open("w+b") as output:
            result = subprocess.run(command, stdin=subprocess.DEVNULL,
                                    stdout=output, stderr=subprocess.STDOUT,
                                    timeout=timeout)
            output.seek(0)
            return subprocess.CompletedProcess(
                command, result.returncode,
                output.read().decode("utf-8", errors="replace"), "")


def _engine_sessions(args, requested: list[str]) -> dict[str, str] | None:
    """Per-target Engine names or com:PID selectors; None after a bad pair."""
    sessions = {}
    for entry in getattr(args, "engine_session", None) or []:
        key, separator, name = entry.partition("=")
        if not separator or key not in requested or not name.strip() or key in sessions:
            fail("engine sessions must be unique selected-target=shared-name pairs")
            return None
        sessions[key] = name.strip()
    return sessions


def _local_com_target_error(key: str, session: str | None) -> str | None:
    """COM borrows the chosen host's desktop, never another host's session."""
    if not session or not session.startswith("com:"):
        return None
    if not re.fullmatch(r"com:[1-9][0-9]*", session):
        return "COM session must be com:<positive PID>"
    expected = os.environ.get("RADIA_ACCEPTANCE_HOSTNAME_" + key.upper(), "").strip().lower()
    if not expected or platform.node().strip().lower() != expected:
        return f"Run COM acceptance locally in the selected {key} desktop logon"
    return None


def _run_simulink_candidate_target(
        key: str, package: Path, package_sha256: str,
        success_marker: str, engine_session: str | None = None,
        python_executable: str | None = None) -> tuple[bool, str]:
    if key not in SIMULINK_TARGETS:
        return False, "unsupported Simulink target"
    label, host, python_command = SIMULINK_TARGETS[key]
    if engine_session and engine_session.startswith("com:"):
        error = _local_com_target_error(key, engine_session)
        if error:
            return False, error
        host = None
    if python_executable:
        # The verifier hands its own interpreter to MATLAB as RADIA_PYTHON_EXECUTABLE.
        python_command = "'" + python_executable.replace("'", "''") + "'"
    verifier = REPO / "tools/verify_simulink_release.py"
    if host is None:
        command = [
            python_executable or sys.executable, str(verifier), str(package),
            "--matlab", MATLAB_EXE,
        ]
        if engine_session:
            command.extend(["--engine-session", engine_session])
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
    else:
        remote_root_posix = f"C:/temp/radia-release-quad/{package_sha256[:16]}"
        remote_root_windows = remote_root_posix.replace("/", "\\")
        prepare = f"New-Item -ItemType Directory -Force -Path '{remote_root_windows}' | Out-Null\n"
        created = _run_ssh_powershell(host, prepare, timeout=60)
        if created.returncode != 0:
            return False, created.stderr.strip() or created.stdout.strip()
        remote_package = f"{remote_root_posix}/{package.name}"
        remote_verifier = f"{remote_root_posix}/{verifier.name}"
        for source, destination in (
                (package, remote_package), (verifier, remote_verifier)):
            copied = subprocess.run(
                ["scp", str(source), f"{host}:{destination}"],
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=180,
            )
            if copied.returncode != 0:
                return False, copied.stderr.strip() or copied.stdout.strip()
        session_option = (
            " --engine-session '" + engine_session.replace("'", "''") + "'"
            if engine_session else "")
        invocation = (
            f"& {python_command} '{remote_verifier}' '{remote_package}' "
            f"--matlab '{MATLAB_EXE}'{session_option}\nexit $LASTEXITCODE\n"
        )
        result = _run_ssh_powershell(host, invocation, timeout=420)
    output = ((result.stdout or "") + (result.stderr or "")).strip()
    if result.returncode != 0:
        return False, output
    if success_marker not in output or '"status": "passed"' not in output:
        return False, f"{label} did not emit both release success markers:\n{output}"
    return True, output


def cmd_simulink_candidate(args):
    """Reject the retired CPU MEX release route."""
    fail("CPU MEX releases are unsupported; use the independent radia-simulink release workflow.")
    return 4


def _verify_simulink_candidate_state(package_arg: str) -> int:
    package = Path(package_arg).resolve()
    try:
        manifest = _simulink_manifest(package)
        package_sha256 = _sha256_file(package)
        mex_sha256 = _simulink_mex_sha256(manifest)
    except (OSError, RuntimeError, ValueError, zipfile.BadZipFile) as error:
        fail(f"invalid Simulink candidate: {error}")
        return 2
    state_path = _simulink_state_path(package_sha256)
    if not state_path.is_file():
        fail("Simulink candidate has no release-quad state. Run "
             "`release_quad simulink-candidate --package <zip> --target all`.")
        return 4
    state = json.loads(state_path.read_text(encoding="utf-8"))
    if state.get("package_sha256") != package_sha256 or \
            state.get("commit") != manifest.get("commit") or \
            state.get("mex_sha256") != mex_sha256:
        fail("Simulink candidate state does not match the supplied archive")
        return 4
    missing = [key for key in SIMULINK_TARGETS
               if state.get("targets", {}).get(key, {}).get("status") != "passed"
               or state.get("targets", {}).get(key, {}).get("mex_sha256") != mex_sha256]
    if missing:
        fail(f"Simulink candidate has not passed: {', '.join(missing)}")
        return 4
    head = _release_head()
    anchored, message = _simulink_candidate_commit_is_release_anchored(
        manifest, head)
    if not anchored:
        fail(message)
        return 4
    info(message)
    ok("supplied Simulink candidate passed on the build target")
    return 0


# ============================================================
# Standalone radia-optuna exact-wheel candidate gate
# ============================================================

def _optuna_state_path(wheel_sha256: str) -> Path:
    return OPTUNA_GATE_ROOT / f"candidate-{wheel_sha256}.json"


def _verify_optuna_wheel(wheel: Path) -> tuple[dict | None, str]:
    verifier = REPO / "packages/radia-optuna/verify_wheel.py"
    result = subprocess.run(
        [
            sys.executable,
            str(verifier),
            str(wheel),
            "--json",
            "--release-candidate",
        ],
        capture_output=True, text=True,
    )
    output = ((result.stdout or "") + (result.stderr or "")).strip()
    if result.returncode != 0:
        return None, output
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError as error:
        return None, f"wheel verifier returned invalid JSON: {error}\n{output}"
    if payload.get("ok") is not True:
        return None, f"wheel verifier did not report ok=true:\n{output}"
    return payload, output


def _optuna_release_source_ready() -> tuple[bool, str]:
    """Return the checked-out radia-optuna candidate commit when it is main.

    radia-optuna is released independently of the solver, so its source
    identity is the checkout HEAD, never the solver ``v<radia version>`` tag.
    """
    current = _git("rev-parse", "HEAD", check=False)
    if current.returncode != 0:
        return False, current.stderr.strip() or "git rev-parse HEAD failed"
    head = current.stdout.strip().lower()
    status = _git("status", "--porcelain", "--untracked-files=no", check=False)
    if status.returncode != 0:
        return False, status.stderr.strip() or "git status failed"
    if status.stdout.strip():
        return False, "tracked release source is dirty"
    fetched = _git("fetch", "origin", "main", "--quiet", check=False)
    if fetched.returncode != 0:
        return False, fetched.stderr.strip() or "git fetch origin main failed"
    origin_main = _git("rev-parse", "origin/main", check=False)
    if origin_main.returncode != 0:
        return False, origin_main.stderr.strip() or "origin/main is unavailable"
    remote_head = origin_main.stdout.strip().lower()
    if head != remote_head:
        return False, f"HEAD {head} differs from origin/main {remote_head}"
    return True, head


def _download_verified_optuna_ci_wheel(ci_run_id: str) -> tuple[dict | None, str]:
    ready, source = _optuna_release_source_ready()
    if not ready:
        return None, source
    head = source
    result = subprocess.run(
        [
            "gh", "run", "view", str(ci_run_id), "--json",
            "databaseId,headSha,headBranch,event,status,conclusion,workflowName,url,jobs",
        ],
        cwd=str(REPO), capture_output=True, text=True,
    )
    if result.returncode != 0:
        return None, result.stderr.strip() or result.stdout.strip()
    try:
        run_info = json.loads(result.stdout)
    except json.JSONDecodeError as error:
        return None, f"gh run view returned invalid JSON: {error}"
    expected = {
        "headSha": head,
        "headBranch": "main",
        "event": "push",
        "status": "completed",
        "conclusion": "success",
        "workflowName": "radia-optuna",
    }
    if run_info.get("event") == "workflow_dispatch":
        expected["event"] = "workflow_dispatch"
    mismatches = [
        f"{key}={run_info.get(key)!r} (expected {value!r})"
        for key, value in expected.items()
        if str(run_info.get(key, "")).lower() != str(value).lower()
    ]
    if mismatches:
        return None, "CI run is not the exact successful main run: " + "; ".join(mismatches)
    jobs = {
        row.get("name"): (row.get("status"), row.get("conclusion"))
        for row in run_info.get("jobs", [])
    }
    required_jobs = (
        "build-test",
        "radia-optuna installed-wheel MATLAB/Simulink E2E",
    )
    missing_jobs = [name for name in required_jobs
                    if jobs.get(name) != ("completed", "success")]
    if missing_jobs:
        return None, "CI run lacks successful required jobs: " + ", ".join(missing_jobs)

    Path(r"C:\temp").mkdir(parents=True, exist_ok=True)
    download_dir = Path(tempfile.mkdtemp(
        prefix=f"radia-optuna-ci-{ci_run_id}-", dir=tempfile.gettempdir()
    ))
    downloaded = subprocess.run(
        [
            "gh", "run", "download", str(ci_run_id),
            "--name", "radia-optuna-wheel", "--dir", str(download_dir),
        ],
        cwd=str(REPO), capture_output=True, text=True,
    )
    if downloaded.returncode != 0:
        return None, downloaded.stderr.strip() or downloaded.stdout.strip()
    wheels = sorted(download_dir.rglob("*.whl"))
    if len(wheels) != 1:
        return None, f"expected one radia-optuna wheel artifact, found {len(wheels)}"
    verification, output = _verify_optuna_wheel(wheels[0])
    if verification is None:
        return None, output
    wheel_sha256 = _sha256_file(wheels[0])
    retained_dir = OPTUNA_GATE_ROOT / f"{head[:12]}-ci-{ci_run_id}"
    retained_dir.mkdir(parents=True, exist_ok=True)
    retained_wheel = retained_dir / wheels[0].name
    if retained_wheel.exists() and _sha256_file(retained_wheel) != wheel_sha256:
        return None, f"retained candidate path contains different bytes: {retained_wheel}"
    if not retained_wheel.exists():
        copy_file(wheels[0], retained_wheel)
    return {
        "ci_run_id": str(ci_run_id),
        "ci_url": run_info.get("url"),
        "commit": head,
        "wheel": str(retained_wheel),
        "wheel_sha256": wheel_sha256,
        "version": verification["version"],
    }, output


def _optuna_remote_invocation(
        runner: str, wheel: str, worker: str, python_command: str,
        wheel_sha256: str, engine_session: str | None = None) -> str:
    """PowerShell script that exits nonzero on any failure.

    A terminating error inside the runner leaves $LASTEXITCODE at the last
    native command's value, often 0, so the exit code cannot come from it
    alone: the runner's errors exit 1 and a native failure keeps its code.
    """
    session_option = (
        " -EngineSession '" + engine_session.replace("'", "''") + "'"
        if engine_session else "")
    return (
        "$ErrorActionPreference = 'Stop'\n"
        "try {\n"
        f"    $pythonExe = (& {python_command} -c \"import sys; print(sys.executable)\").Trim()\n"
        "    if ($LASTEXITCODE -ne 0) { throw \"Python discovery failed with exit code $LASTEXITCODE\" }\n"
        f"    & '{runner}' -Wheel '{wheel}' "
        f"-MatlabExecutable '{MATLAB_EXE}' -PythonExecutable $pythonExe "
        f"-PreverifiedWheelSha256 '{wheel_sha256}' "
        f"-EngineWorker '{worker}'{session_option}\n"
        "    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }\n"
        "} catch {\n"
        "    [Console]::Error.WriteLine(($_ | Out-String))\n"
        "    exit 1\n"
        "}\n"
        "exit 0\n"
    )


def _run_optuna_candidate_target(
        key: str, wheel: Path, wheel_sha256: str,
        engine_session: str | None = None) -> tuple[bool, str]:
    if key not in OPTUNA_TARGETS:
        return False, f"Unsupported Optuna CUDA acceptance target: {key}"
    label, host, python_command = OPTUNA_TARGETS[key]
    if key == "100" and platform.node().strip().lower() == "t7920":
        host = None
    if engine_session and engine_session.startswith("com:"):
        error = _local_com_target_error(key, engine_session)
        if error:
            return False, error
        host = None
    runner = OPTUNA_WHEEL_RUNNER
    matlab_tests = REPO / "packages/radia-optuna/tests/matlab"
    # The solver release Engine worker owns MATLAB session reuse and launch.
    worker = REPO / "tools/verify_simulink_release.py"
    if host is None:
        command = [
            "pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File",
            str(runner), "-Wheel", str(wheel),
            "-MatlabExecutable", MATLAB_EXE,
            "-PythonExecutable", sys.executable,
            "-PreverifiedWheelSha256", wheel_sha256,
            "-EngineWorker", str(worker),
        ]
        if engine_session:
            command.extend(["-EngineSession", engine_session])
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
    else:
        remote_root_posix = f"C:/temp/radia-release-quad/optuna-{wheel_sha256[:16]}"
        remote_root_windows = remote_root_posix.replace("/", "\\")
        remote_matlab_windows = remote_root_windows + r"\matlab"
        prepare = (
            "$ErrorActionPreference = 'Stop'\n"
            f"New-Item -ItemType Directory -Force -Path '{remote_matlab_windows}' | Out-Null\n"
        )
        created = _run_ssh_powershell(host, prepare, timeout=60)
        if created.returncode != 0:
            return False, created.stderr.strip() or created.stdout.strip()
        sources = [runner, worker, wheel, *sorted(matlab_tests.glob("*.m"))]
        for source_file in sources:
            subdir = "matlab/" if source_file.suffix == ".m" else ""
            destination = f"{remote_root_posix}/{subdir}{source_file.name}"
            copied = subprocess.run(
                ["scp", str(source_file), f"{host}:{destination}"],
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
            if copied.returncode != 0:
                return False, copied.stderr.strip() or copied.stdout.strip()
        invocation = _optuna_remote_invocation(
            remote_root_windows + "\\" + runner.name,
            remote_root_windows + "\\" + wheel.name,
            remote_root_windows + "\\" + worker.name,
            python_command, wheel_sha256, engine_session,
        )
        result = _run_ssh_powershell(host, invocation, timeout=1800)
    output = ((result.stdout or "") + (result.stderr or "")).strip()
    if result.returncode != 0:
        return False, output
    if OPTUNA_SUCCESS_MARKER not in output:
        return False, f"{label} did not emit {OPTUNA_SUCCESS_MARKER}:\n{output}"
    return True, output


def cmd_optuna_candidate(args):
    """Download one main-CI wheel and run it on both MATLAB machines."""
    step("radia-optuna exact-wheel candidate gate (build)")
    candidate, output = _download_verified_optuna_ci_wheel(args.ci_run_id)
    if candidate is None:
        fail(f"invalid radia-optuna CI candidate: {output}")
        return 2
    wheel = Path(candidate["wheel"])
    wheel_sha256 = candidate["wheel_sha256"]
    state_path = _optuna_state_path(wheel_sha256)
    state = {
        "schema": "radia.release-quad.optuna-candidate.v1",
        **candidate,
        "targets": {},
    }
    if state_path.is_file():
        previous = json.loads(state_path.read_text(encoding="utf-8"))
        if (previous.get("wheel_sha256") == wheel_sha256
                and previous.get("commit") == candidate["commit"]
                and str(previous.get("ci_run_id")) == str(args.ci_run_id)):
            state["targets"] = previous.get("targets", {})

    requested = [part.strip().lower() for part in args.target.split(",")]
    if "all" in requested:
        requested = list(OPTUNA_TARGETS)
    unknown = sorted(set(requested) - set(OPTUNA_TARGETS))
    if unknown:
        fail(f"unknown radia-optuna target(s): {', '.join(unknown)}")
        return 2
    sessions = _engine_sessions(args, requested)
    if sessions is None:
        return 2
    failed = 0
    for key in requested:
        label = OPTUNA_TARGETS[key][0]
        info(f"verifying exact installed wheel on {label}")
        session_args = [sessions[key]] if key in sessions else []
        passed, target_output = _run_optuna_candidate_target(
            key, wheel, wheel_sha256, *session_args
        )
        state["targets"][key] = {
            "label": label,
            "status": "passed" if passed else "failed",
            "verified_at_utc": datetime.now(timezone.utc).isoformat(),
            "engine_session": sessions.get(key),
            "output_tail": target_output[-4000:],
        }
        _write_simulink_state(state_path, state, key)
        if passed:
            ok(f"radia-optuna wheel passed on {label}")
        else:
            failed += 1
            fail(f"radia-optuna wheel failed on {label}")
            if target_output:
                print(target_output[-4000:])
    if failed:
        return 4
    ok(f"retained wheel: {wheel}")
    ok(f"wheel SHA256: {wheel_sha256}")
    ok(f"candidate state: {state_path}")
    return 0


def _verify_optuna_candidate_state(wheel_arg: str) -> tuple[int, dict | None]:
    wheel = Path(wheel_arg).resolve()
    verification, output = _verify_optuna_wheel(wheel)
    if verification is None:
        fail(f"invalid radia-optuna wheel: {output}")
        return 2, None
    wheel_sha256 = _sha256_file(wheel)
    state_path = _optuna_state_path(wheel_sha256)
    if not state_path.is_file():
        fail("radia-optuna wheel has no release-quad state. Run "
             "`release_quad optuna-candidate --ci-run-id <id> --target all`.")
        return 4, None
    state = json.loads(state_path.read_text(encoding="utf-8"))
    ready, source = _optuna_release_source_ready()
    if not ready:
        fail(source)
        return 4, None
    if (state.get("wheel_sha256") != wheel_sha256
            or state.get("commit") != source
            or state.get("version") != verification.get("version")):
        fail("radia-optuna candidate state differs from wheel, version, or HEAD")
        return 4, None
    missing = [key for key in OPTUNA_TARGETS
               if state.get("targets", {}).get(key, {}).get("status") != "passed"]
    if missing:
        fail(f"radia-optuna candidate has not passed: {', '.join(missing)}")
        return 4, None
    ok("exact radia-optuna wheel passed 100号機")
    return 0, state


def cmd_optuna_done(args):
    """Definition-of-done gate for the independent radia-optuna wheel."""
    step("radia-optuna definition of done")
    rc, state = _verify_optuna_candidate_state(args.wheel)
    if rc != 0 or state is None:
        return rc
    print("")
    ok("RADIA-OPTUNA DEFINITION OF DONE met")
    info(f"CI run ID: {state['ci_run_id']}")
    info(f"commit: {state['commit']}")
    info(f"wheel: {state['wheel']}")
    info(f"SHA256: {state['wheel_sha256']}")
    info("Next: create the matching radia-optuna-v<version> tag, then dispatch "
         "release-radia-optuna.yml with this CI run ID and SHA256.")
    return 0


# ============================================================
# Phases
# ============================================================

def cmd_preflight(args):
    """Read-only numerical-solver state report. Always safe to run."""
    step("Solver preflight: state report")

    # Versions in repo
    v = _read_repo_versions({"radia", "radia.__version__"})
    info("Repo versions:")
    info(f"  radia              pyproject={v['radia']}  __version__={v['radia.__version__']}")
    pp_radia = (v["radia"] == v["radia.__version__"])
    if pp_radia: ok("radia pyproject == __init__")
    else:        fail("radia pyproject != __init__ — fix before any release")
    if not pp_radia:
        return 2

    _check_main_synced(hard=False)

    return 0


SOLVER_INSTALL_GUARD = r'''
import importlib.metadata as metadata
import importlib.util
import os
from pathlib import Path
import sys
import sysconfig
if os.name == "nt":
    candidates = set()
    try:
        distribution = metadata.distribution("radia")
    except metadata.PackageNotFoundError:
        distribution = None
    if distribution is not None:
        scripts = Path(sysconfig.get_path("scripts"))
        candidates.update(scripts / (ep.name + ".exe") for ep in distribution.entry_points
                          if ep.group == "console_scripts")
    spec = importlib.util.find_spec("radia")
    for location in (spec.submodule_search_locations or []) if spec else []:
        root = Path(location)
        candidates.update(root.glob("*.pyd"))
        candidates.update(root.glob("*.dll"))
    blocked = []
    for path in sorted(candidates):
        if path.is_file():
            try:
                with path.open("r+b"):
                    pass
            except OSError:
                blocked.append(str(path))
    if blocked:
        sys.exit("Radia install blocked before pip: binaries are locked or not writable: "
                 + "; ".join(blocked))
'''


def _solver_install_guard_powershell(python_command="python"):
    encoded = base64.b64encode(SOLVER_INSTALL_GUARD.encode("utf-8")).decode("ascii")
    return (f'{python_command} -c "import base64; exec(base64.b64decode(\'{encoded}\'))"\n'
            'if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }\n')


EDITABLE_RELEASE_VERIFY = r'''import hashlib, importlib.metadata as md, importlib.util, json, pathlib, sys, sysconfig
repo = pathlib.Path(sys.argv[1]).resolve()
expected_commit, expected_ngsolve, expected_netgen = sys.argv[2:5]
preinstall = len(sys.argv) > 5 and sys.argv[5] == "preinstall"
package = repo / "src" / "radia"
manifest_path = package / "release_native_payloads.json"
if not manifest_path.is_file():
    raise SystemExit(f"missing exact-wheel native manifest: {manifest_path}")
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
if manifest.get("schema") != "radia.release-native-payloads.v1":
    raise SystemExit("invalid exact-wheel native manifest schema")
if manifest.get("source_commit", "").lower() != expected_commit:
    raise SystemExit("native manifest source commit differs from the release tag")
files = manifest.get("files")
if not isinstance(files, dict):
    raise SystemExit("native manifest files must be an object")
actual = {p.name: p for pattern in ("*.pyd", "*.dll") for p in package.glob(pattern)}
required = {"_radia_pybind.pyd", "axifem.pyd", "peec_matrices.pyd", "sparsesolv_ngsolve.pyd", "radia_motor_rom.dll"}
if not required <= files.keys():
    raise SystemExit("native manifest omits required NGSolve ABI payloads")
if set(actual) != set(files):
    raise SystemExit(f"native payload set differs: actual={sorted(actual)} expected={sorted(files)}")
for name, path in actual.items():
    lowered = name.lower()
    if name == "c" + "ln_core.pyd" or "locked-old" in lowered or ".pre_" in lowered:
        raise SystemExit(f"retired native payload remains import-visible: {name}")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if files[name] != digest:
        raise SystemExit(f"native payload hash differs for {name}")
    if name in required and path.suffix == ".pyd":
        provenance = json.loads(pathlib.Path(str(path) + ".build.json").read_text(encoding="utf-8"))
        if (provenance.get("schema") != "radia.native-build-provenance.v1"
                or provenance.get("source_commit") != expected_commit
                or provenance.get("source_dirty") is not False
                or provenance.get("binary_name") != name
                or provenance.get("binary_sha256") != digest
                or provenance.get("binary_bytes") != path.stat().st_size):
            raise SystemExit(f"native source provenance differs for {name}")
spec = importlib.util.spec_from_file_location("radia_gpu_payload_check", repo / "tools/python_gpu_payload.py")
payload = importlib.util.module_from_spec(spec)
spec.loader.exec_module(payload)
for module in payload.MODULES:
    name = module + sysconfig.get_config_var("EXT_SUFFIX")
    binary = package / name
    receipt = json.loads(pathlib.Path(str(binary) + ".build.json").read_text(encoding="utf-8"))
    payload.verify_payload(name, binary.read_bytes(), receipt, expected_commit)
if not preinstall:
    if md.version("ngsolve") != expected_ngsolve or md.version("netgen-mesher") != expected_netgen:
        raise SystemExit("installed NGSolve/Netgen differs from package metadata pins")
    import ngsolve, radia
    if ngsolve.__version__ != expected_ngsolve:
        raise SystemExit("loaded NGSolve differs from package metadata pin")
    distribution = md.distribution("radia")
    import tomllib
    expected_version = tomllib.loads((repo / "pyproject.toml").read_text(encoding="utf-8"))["project"]["version"]
    if distribution.version != expected_version or radia.__version__ != expected_version:
        raise SystemExit("editable metadata/import version differs from bound source")
    registration = json.loads(distribution.read_text("direct_url.json") or "{}")
    if registration.get("dir_info", {}).get("editable") is not True:
        raise SystemExit("normal runtime is not registered editable")
    if repo not in pathlib.Path(radia.__file__).resolve().parents:
        raise SystemExit(f"radia imported outside fixed release checkout: {radia.__file__}")
    from radia._cuda_runtime import _pins
    for name, version in _pins().items():
        if md.version(name) != version:
            raise SystemExit("GPU runtime metadata pin mismatch: " + name)
print(json.dumps({"source_commit": expected_commit, "native_sha256": files}, sort_keys=True))
'''


def _editable_release_verify_command(repo, expected_sha, python_command="python", *, preinstall=False):
    pins = _solver_dependency_pins()
    encoded = base64.b64encode(EDITABLE_RELEASE_VERIFY.encode("utf-8")).decode("ascii")
    command = [python_command, "-I", "-c", f"import base64; exec(base64.b64decode('{encoded}'))",
               repo, expected_sha, pins["ngsolve"], pins["netgen-mesher"]]
    return command + (["preinstall"] if preinstall else [])


def _solver_abi_probe_command(python_command="python"):
    pins = _solver_dependency_pins()
    code = ("import importlib.metadata as m,ngsolve;"
            f"assert ngsolve.__version__ == '{pins['ngsolve']}';"
            f"assert m.version('netgen-mesher') == '{pins['netgen-mesher']}'")
    return [python_command, "-c", code]


def _verify_final_pip_checks() -> int:
    """Check both normal editable and separate exact-wheel QA runtimes."""
    targets = [_editable_target(role) for role in SOLVER_TARGETS]
    checks = [(t["ssh"], t["python"]) for t in targets]
    checks += [(host, python) for _, host, python in (_solver_qa_target(role) for role in SOLVER_TARGETS)]
    for host, python in checks:
        script = _runtime_guard_powershell(python) + f"& {_powershell_literal(python)} -I -m pip check; exit $LASTEXITCODE"
        encoded = base64.b64encode(script.encode("utf-16le")).decode("ascii")
        if run(["ssh", host, "pwsh", "-NoProfile", "-EncodedCommand", encoded], check=False).returncode:
            fail("pip check failed in a normal or isolated QA runtime")
            return 4
    return 0


def _verify_qa_wheel(role):
    """Exact-wheel qualification is isolated from normal editable operation."""
    target = _editable_target(role)
    label, host, python = _solver_qa_target(role)
    if _check_remote_target(target, python=python, source=False, qa=True):
        return 4
    return _verify_wheel_runtime(host, label + " QA", python_path=python, isolated=True)


def _wheel_runtime_probe_code(version, label):
    """Build the target-side check for wheel metadata, import origin and RECORD."""
    return """import base64, hashlib, importlib.metadata as m, json, pathlib, sysconfig
import radia
version = VERSION_PLACEHOLDER
d = m.distribution('radia')
assert d.version == version and radia.__version__ == version, 'Radia version mismatch'
assert not json.loads(d.read_text('direct_url.json') or '{}').get('dir_info', {}).get('editable', False), 'LABEL_PLACEHOLDER must use a wheel'
origin = pathlib.Path(radia.__file__).resolve()
assert origin == pathlib.Path(d.locate_file('radia/__init__.py')).resolve(), 'Shadowed Radia import'
install_paths = sysconfig.get_paths()
wheel_roots = [pathlib.Path(install_paths[key]).resolve() for key in ('purelib', 'platlib')]
assert any(origin.is_relative_to(root) for root in wheel_roots), 'Radia wheel must be under site-packages'
assert d.files, 'Missing wheel RECORD'
checked = 0
for f in d.files:
    if f.hash is not None:
        digest = hashlib.new(f.hash.mode, pathlib.Path(d.locate_file(f)).read_bytes()).digest()
        actual = base64.urlsafe_b64encode(digest).rstrip(b'=').decode()
        assert actual == f.hash.value, 'Wheel payload mismatch: ' + str(f)
        checked += 1
assert checked, 'Empty wheel hash manifest'
print('LABEL_PLACEHOLDER wheel verified:', version, checked)
""".replace("VERSION_PLACEHOLDER", repr(version)).replace("LABEL_PLACEHOLDER", label)


def _verify_wheel_runtime(ssh_host, label, *, python_path, isolated=True):
    """Verify the explicitly selected QA Python; never select a PATH runtime."""
    _absolute_windows_path(python_path, "QA Python", python=True)
    code = _wheel_runtime_probe_code(_radia_version(), label)
    if isolated:
        code = "import sys\nif sys.prefix == sys.base_prefix: raise SystemExit('QA runtime must be isolated')\n" + code
    encoded_code = base64.b64encode(code.encode('utf-8')).decode('ascii')
    invocation = '& ' + _powershell_literal(python_path)
    guard = _runtime_guard_powershell(python_path)
    script = guard + f'{invocation} -c "import base64; exec(base64.b64decode(\'{encoded_code}\'))"; exit $LASTEXITCODE'
    encoded = base64.b64encode(script.encode('utf-16le')).decode('ascii')
    result = run(['ssh', ssh_host, 'pwsh', '-NoProfile', '-EncodedCommand', encoded], check=False)
    if result.returncode != 0:
        fail(f"{label} does not hold the verified Radia {_radia_version()} wheel")
        return 4
    return 0


def _deploy_validation():
    return _deploy_editable_role("validation")


# Remote editable-deploy refusals, all before any install.  The guards write to
# stderr and exit: Write-Error under ErrorActionPreference Stop would end the
# script with exit 1 before the documented code.
SOURCE_SHA_MISMATCH_EXIT = 41
SOURCE_DIRTY_EXIT = 42
DEV_VENV_MISSING_EXIT = 43


def _dev_venv_guard_powershell(python_path):
    """Exit DEV_VENV_MISSING_EXIT before any install when the venv is absent."""
    literal = python_path.replace("'", "''")
    return (f"if (-not (Test-Path -LiteralPath '{literal}' -PathType Leaf)) {{\n"
            f"  [Console]::Error.WriteLine('Development venv interpreter is missing: {literal}')\n"
            f"  exit {DEV_VENV_MISSING_EXIT}\n"
            "}\n")


def _verify_live_editable(target):
    script = _target_check_powershell(target)
    script += _ps_python(_editable_release_verify_command(target["source"], target["commit"], target["python"]))
    script += "\nexit $LASTEXITCODE"
    encoded = base64.b64encode(script.encode("utf-16le")).decode("ascii")
    return run(["ssh", target["ssh"], "pwsh", "-NoProfile", "-EncodedCommand", encoded], check=False).returncode


def _deploy_editable_role(role, reason=None):
    """One guarded normal Radia installation; exact-wheel QA is independent."""
    target = _editable_target(role)
    if _check_remote_target(target):
        return 3
    python = target["python"]
    # Compression is transport only: avoid the Windows command-line limit while
    # executing the same inspectable checks, not a different install route.
    import zlib
    program = _editable_install_program(target)
    token = base64.b64encode(zlib.compress(program.encode("utf-8"))).decode("ascii")
    invocation = [python, "-I", "-c", f"import base64,zlib;exec(zlib.decompress(base64.b64decode('{token}')))"]
    script = _runtime_guard_powershell(python) + _ps_python(invocation) + "\nexit $LASTEXITCODE"
    encoded = base64.b64encode(script.encode("utf-16le")).decode("ascii")
    if len(encoded) > 28000:
        raise ValueError("Editable installer exceeds the Windows command-line transport budget")
    if run(["ssh", target["ssh"], "pwsh", "-NoProfile", "-EncodedCommand", encoded], check=False).returncode:
        fail(target["label"] + " editable installation/ABI/native/pip verification failed")
        return 3
    rc = _record_release_intent_remote(role, reason)
    return rc or _verify_role_editable(role)


def _editable_install_program(target):
    """All target checks precede pip; all payload checks precede intent writes."""
    pins = _solver_dependency_pins()
    args = ["native-check", target["source"], target["commit"], pins["ngsolve"], pins["netgen-mesher"]]
    install = [target["python"], "-I", "-m", "pip", "install", "--no-cache-dir", "--only-binary=:all:", "-e", target["source"] + "[gpu]"]
    check = [target["python"], "-I", "-m", "pip", "check"]
    return ("import subprocess,sys\n" +
            "target_check=" + repr(_target_check_code(target)) + "\n" +
            "native_check=" + repr(EDITABLE_RELEASE_VERIFY) + "\n" +
            "native_args=" + repr(args) + "\n" +
            "def verify_native(preinstall):\n"
            " old=sys.argv\n"
            " sys.argv=native_args+(['preinstall'] if preinstall else [])\n"
            " try: exec(native_check, {'__name__':'__native_payload_check__'})\n"
            " finally: sys.argv=old\n"
            "exec(target_check)\nverify_native(True)\n" +
            "exec(" + repr(SOLVER_INSTALL_GUARD) + ")\n" +
            "subprocess.run(" + repr(install) + ",check=True)\n" +
            "exec(target_check)\nverify_native(False)\n" +
            "subprocess.run(" + repr(check) + ",check=True)\n")


def _deploy_100():
    return _deploy_editable_role("100")


def _deploy_mdx(host):
    """Reject every retired refresh route before runtime or network changes."""
    fail("Retired compute deployment routes are disabled")
    return 2


def cmd_phase8(args):
    """Deploy and verify the Radia solver on selected QUAD targets."""
    targets = [t.strip().lower() for t in (args.target or "validation,100").split(",")]
    if any(t not in (*SOLVER_TARGETS, "all") for t in targets):
        fail("unsupported solver target")
        return 2
    selected = list(SOLVER_TARGETS) if "all" in targets else list(dict.fromkeys(targets))
    resolved = [_editable_target(role) for role in selected]
    expected = _release_head()
    if any(target["commit"] != expected for target in resolved):
        raise ValueError("Normal editable source bindings must match the exact release tag")
    rc = cmd_preflight(args)
    if rc != 0:
        return rc
    for t in targets:
        if t == "validation":
            rc = _deploy_validation()
            if rc != 0:
                return rc
        elif t == "100":
            rc = _deploy_100()
            if rc != 0:
                return rc
        elif t == "all":
            rc = _deploy_validation()
            if rc != 0:
                return rc
            rc = _deploy_100()
            if rc != 0:
                return rc
        else:
            fail(f"unknown target: {t!r}")
            return 2
    return 0


def cmd_phase8e(args):
    """Reject the retired optional compute refresh command."""
    fail("Retired compute refresh is disabled; release targets are validation and 100")
    return 2


CROSS_MACHINE_PROBE = '''import hashlib, importlib.metadata as md, os

# This is an identity probe, not a throughput benchmark.  Keep its import
# footprint bounded so it can coexist with an authorized MATLAB computation.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

def hsh_text(p):
    h = hashlib.sha256()
    try:
        d = open(p, "rb").read().replace(b"\\r\\n", b"\\n").replace(b"\\r", b"\\n")
        h.update(d); return h.hexdigest()[:12]
    except Exception: return "MISSING"

import radia
rad = os.path.dirname(radia.__file__)
print(f"VER radia              = {radia.__version__}")
print(f"ABI ngsolve            = {md.version('ngsolve')}")
print(f"ABI netgen-mesher      = {md.version('netgen-mesher')}")
native = sorted(p for pattern in ("*.pyd", "*.dll") for p in __import__('pathlib').Path(rad).glob(pattern))
h = hashlib.sha256()
for p in native:
    h.update(p.name.encode()); h.update(b"\\0"); h.update(hashlib.sha256(p.read_bytes()).digest())
print(f"NATIVE radia           = {h.hexdigest()[:12]}")
for r in ["panels/calc_inductance.py",
          "panels/calc_fem_kelvin.py",
          "panels/calc_fem_coilmesh.py"]:
    print(f"SHA radia/{r:35s} = {hsh_text(os.path.join(rad,r))}")
'''


def _normal_probe_prefix(target):
    token = base64.b64encode(target["source"].encode()).decode("ascii")
    return _target_check_code(target) + '''import base64,importlib.metadata as md,json,pathlib
d=md.distribution("radia")
if json.loads(d.read_text("direct_url.json") or "{}").get("dir_info",{}).get("editable") is not True: raise SystemExit("Normal runtime must be editable")
import radia
expected=pathlib.Path(base64.b64decode("TOKEN").decode()).resolve()/"src"/"radia"
if pathlib.Path(radia.__file__).resolve().parent!=expected: raise SystemExit("Normal runtime import differs from source binding")
'''.replace("TOKEN", token)


def _qa_probe_prefix(target, python):
    return _target_check_code(target, python=python, source=False, qa=True) + '''import importlib.metadata as md,json,pathlib
d=md.distribution("radia")
if json.loads(d.read_text("direct_url.json") or "{}").get("dir_info",{}).get("editable",False): raise SystemExit("QA runtime must contain the exact wheel")
import radia
if pathlib.Path(radia.__file__).resolve()!=pathlib.Path(d.locate_file("radia/__init__.py")).resolve(): raise SystemExit("QA wheel import is shadowed")
'''


def _run_script_with_file_output(command, source, timeout):
    """Bound remote Python runs without waiting for inherited output pipes."""
    if command[0] == "ssh" and source is not None:
        # Python's stdin route can remain open under Windows OpenSSH. Stage
        # the script, then execute with stdin closed and clean it remotely.
        host, python, dash, *arguments = command[1:]
        if dash != "-":
            raise ValueError("Expected a Python stdin command")
        with tempfile.TemporaryDirectory(prefix="radia-probe-", dir=tempfile.gettempdir()) as scratch:
            source_path = Path(scratch) / "probe.py"
            source_path.write_text(source, encoding="utf-8")
            remote = "C:/temp/" + Path(scratch).name + ".py"
            subprocess.run(["scp", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
                            str(source_path), f"{host}:{remote}"],
                           stdin=subprocess.DEVNULL, timeout=60, check=True)
            quoted_args = " ".join("'" + arg.replace("'", "''") + "'" for arg in arguments)
            runtime_guard = ""
            if python.startswith("'") and python.endswith("'"):
                runtime_guard = _runtime_guard_powershell(python[1:-1].replace("''", "'"))
            script = (f"try {{ {runtime_guard}& {python} -X utf8 '{remote}' {quoted_args}; "
                      "$result = $LASTEXITCODE } finally { "
                      f"Remove-Item -LiteralPath '{remote}' -Force }}; exit $result")
            encoded = base64.b64encode(script.encode("utf-16-le")).decode("ascii")
            return _run_script_with_file_output(
                ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", host,
                 "pwsh", "-NoProfile", "-EncodedCommand", encoded], None, timeout)
    with tempfile.TemporaryFile(mode="w+b", dir=tempfile.gettempdir()) as output, \
            tempfile.TemporaryFile(mode="w+b", dir=tempfile.gettempdir()) as errors:
        result = subprocess.run(command, input=source.encode("utf-8") if source is not None else None,
                                **({"stdin": subprocess.DEVNULL} if source is None else {}),
                                stdout=output, stderr=errors, timeout=timeout,
                                check=False)
        output.seek(0)
        errors.seek(0)
        return subprocess.CompletedProcess(
            command, result.returncode,
            output.read().decode("utf-8", errors="replace"),
            errors.read().decode("utf-8", errors="replace"))


def _probe(host_label, cmd_prefix, probe_src=CROSS_MACHINE_PROBE):
    """Run the probe on a target (cmd_prefix is the python invocation).

    probe_src defaults to the consumer probe (hashes the installed wheel
    files).  Normal editable probes check actual working files after exact-source guards.
    Isolated QA probes check installed wheel files with explicit runtime guards.
    """
    try:
        p = _run_script_with_file_output(cmd_prefix, probe_src, timeout=120)
    except (OSError, subprocess.TimeoutExpired) as exc:
        fail(f"probe could not run on {host_label}: {exc}")
        return None
    if p.returncode != 0:
        fail(f"probe failed on {host_label}: {p.stderr.strip()}")
        return None
    return p.stdout


_PHASE9_FIELDS = (
    "VER radia",
    "ABI ngsolve",
    "ABI netgen-mesher",
    "NATIVE radia",
    "SHA radia/panels/calc_inductance.py",
    "SHA radia/panels/calc_fem_kelvin.py",
    "SHA radia/panels/calc_fem_coilmesh.py",
)


def _parse_phase9_probe(label, output):
    """Require complete unique fields, PEP 440 versions and ordered inclusive bounds.

    The packaging parser is mandatory; missing support fails the gate. Bounds
    certify a well-formed declared interval, not installed-version membership.
    Cross-host comparison still requires identical reported declarations.
    """
    try:
        from packaging.version import InvalidVersion, Version
    except ImportError as exc:
        raise ValueError("phase9 requires packaging; install it with python -m pip install packaging") from exc
    result = {}
    for line in output.splitlines():
        if not line.strip():
            continue
        key, separator, value = line.partition("=")
        key, value = " ".join(key.split()), value.strip()
        if not separator or key not in _PHASE9_FIELDS or not value:
            raise ValueError(f"{label}: invalid probe row {line!r}")
        if key in result:
            raise ValueError(f"{label}: duplicate field {key}")
        if key.startswith(("SHA ", "NATIVE ")):
            valid = len(value) == 12 and all(c in "0123456789abcdef" for c in value)
        else:
            try:
                Version(value)
                valid = True
            except InvalidVersion:
                valid = False
        if not valid:
            raise ValueError(f"{label}: invalid value for {key}: {value!r}")
        result[key] = value
    missing = set(_PHASE9_FIELDS) - result.keys()
    if missing:
        raise ValueError(f"{label}: missing fields: {', '.join(sorted(missing))}")
    return result


def cmd_phase9(args):
    """Cross-machine consistency probe."""
    # Resolve every private binding before probing any target.
    normal = [_editable_target(role) for role in SOLVER_TARGETS]
    qa = [_solver_qa_target(role) for role in SOLVER_TARGETS]
    targets = []
    for target, (label, host, python) in zip(normal, qa):
        targets.append((label + " QA", ["ssh", host, _powershell_literal(python), "-"],
                        _qa_probe_prefix(target, python) + CROSS_MACHINE_PROBE))
        targets.append((label, ["ssh", target["ssh"], _powershell_literal(target["python"]), "-"],
                        _normal_probe_prefix(target) + CROSS_MACHINE_PROBE))
    step("Phase 9: normal editable and isolated exact-wheel QA consistency")
    outputs = []
    for label, cmd_prefix, probe_src in targets:
        out = _probe(label, cmd_prefix, probe_src)
        if not out:
            fail(f"could not collect probe data from {label}")
            return 4
        try:
            outputs.append((label, _parse_phase9_probe(label, out)))
        except ValueError as exc:
            fail(str(exc))
            return 4

    if not outputs:
        fail("could not collect probe data from any machine")
        return 4

    drift = 0

    labels = [label for label, _ in outputs]
    header = f"\n  {'field':<44} | " + " | ".join(f"{label:<18}" for label in labels)
    print(header)
    print("  " + "-" * max(110, len(header) - 2))
    for key in _PHASE9_FIELDS:
        values = [rows[key] for _, rows in outputs]
        comparable = [value for value in values if value != "N/A"]
        match = bool(comparable) and all(value == comparable[0] for value in comparable)
        marker = _color("32;1", "OK") if match else _color("31;1", "DRIFT")
        if not match:
            drift += 1
        row = f"  {key:<44} | " + " | ".join(f"{value:<18}" for value in values)
        print(f"{row}  [{marker}]")

    if drift:
        print("")
        fail(f"{drift} field(s) drift across machines — release NOT done.")
        return 4
    print("")
    ok("all fields match across validation / 100号機 — release verified.")
    return 0


def cmd_all(args):
    """Run the full deploy + verify chain (phase8 validation+100, phase9)."""
    rc = cmd_phase8(argparse.Namespace(target="validation,100"))
    if rc != 0: return rc
    return cmd_phase9(args)


def _remote_100_editable_packages():
    root = _editable_repo_100()
    return [("radia", root)]


def _verify_release_controller(controller, release_sha, version):
    """Accept a later tooling controller only on top of the exact release.

    The controller may carry release-tool repairs made after the tag, so its
    HEAD may differ from the release commit.  It must then be tracked-clean,
    descend from the release commit, and declare the same Radia version that
    the release commit declares; the release source itself is checked
    separately by the strict exact-SHA gate.
    """
    safe = str(controller).replace("\\", "/")
    git = [GIT_EXE, "-c", f"safe.directory={safe}", "-C", str(controller)]
    head = run([*git, "rev-parse", "HEAD"], capture=True, check=False)
    got = (head.stdout or "").strip().lower()
    if head.returncode != 0 or len(got) != 40:
        fail(f"release controller HEAD is unavailable: {controller}")
        return 4
    dirty = run([*git, "status", "--porcelain", "--untracked-files=no"],
                capture=True, check=False)
    if dirty.returncode != 0 or (dirty.stdout or "").strip():
        fail("release controller has tracked changes: "
             f"{(dirty.stdout or '').strip() or '<status unavailable>'}")
        return 4
    ancestry = run([*git, "merge-base", "--is-ancestor", release_sha, got],
                   capture=True, check=False)
    if ancestry.returncode != 0:
        fail(f"release controller {got[:12]} does not descend from the release "
             f"commit {release_sha[:12]}")
        return 4
    released = run([*git, "show", f"{release_sha}:pyproject.toml"],
                   capture=True, check=False, encoding="utf-8")
    try:
        declared = tomllib.loads(released.stdout or "")["project"]["version"]
    except (tomllib.TOMLDecodeError, KeyError, TypeError):
        declared = None
    if released.returncode != 0 or declared != version:
        fail(f"release controller declares radia {version}, but the release commit "
             f"declares {declared or '<unreadable>'}")
        return 4
    ok(f"release controller {got[:10]} is clean and descends from v{version} "
       f"({release_sha[:10]})")
    return 0


def _verify_local_release_source(repo, expected_sha):
    """Fail unless the release source is the exact tracked-clean release SHA."""
    path = Path(repo)
    if not path.is_dir():
        fail(f"release source does not exist: {repo}")
        return 2
    safe_repo = repo.replace("\\", "/")
    head = run(
        [GIT_EXE, "-c", f"safe.directory={safe_repo}",
         "-C", repo, "rev-parse", "HEAD"],
        capture=True, check=False,
    )
    got_sha = (head.stdout or "").strip().lower()
    if head.returncode != 0 or got_sha != expected_sha.lower():
        fail(
            "release source SHA mismatch: "
            f"expected {expected_sha}, got {got_sha or '<unavailable>'}"
        )
        return 4
    dirty = run(
        [GIT_EXE, "-c", f"safe.directory={safe_repo}",
         "-C", repo, "status", "--porcelain", "--untracked-files=no"],
        capture=True, check=False,
    )
    tracked_changes = (dirty.stdout or "").strip()
    if dirty.returncode != 0 or tracked_changes:
        fail(
            "release source has tracked changes: "
            f"{tracked_changes or '<status unavailable>'}"
        )
        return 4
    ok(f"release source is exact and clean ({got_sha[:10]})")
    return 0


_PATH_ALIASES = ()  # Private namespace mappings are operator configuration.

def _norm_path(p):
    """Lower-case + forward-slashes + strip trailing slash so a UNC and
    a drive-letter form of the same NAS path compare equal."""
    p = (p or "").replace("\\", "/").rstrip("/").lower()
    # Normalize equivalent configured path spellings.
    aliases = json.loads(os.environ.get("RADIA_RELEASE_CONTROLLER_PATH_ALIASES", "[]"))
    for alias, canonical in (*_PATH_ALIASES, *aliases):
        p = p.replace(alias, canonical)
    return p


def _pip_show(pkg):
    """Return parsed dict from `pip show <pkg>`, or None if not
    installed.  Keys are lower-cased; values are stripped strings."""
    try:
        dist = importlib_metadata.distribution(pkg)
    except importlib_metadata.PackageNotFoundError:
        return None

    d = {
        "version": dist.version,
        "location": str(dist.locate_file("")),
    }
    direct_url = dist.read_text("direct_url.json")
    if direct_url:
        try:
            data = json.loads(direct_url)
        except json.JSONDecodeError:
            data = {}
        if data.get("dir_info", {}).get("editable") and data.get("url"):
            parsed = urlparse(data["url"])
            if parsed.scheme == "file":
                path = url2pathname(parsed.path)
                if parsed.netloc:
                    path = f"//{parsed.netloc}{path}"
                elif os.name == "nt" and path.startswith("/") and len(path) > 2 and path[2] == ":":
                    path = path[1:]
                d["editable project location"] = path
            else:
                d["editable project location"] = data["url"]
    return d


_EDITABLE_IMPORT_MODULES = {
    "radia": "radia",
}


def _fresh_import_origin(pkg):
    """Return a package's imported __file__ from a fresh Python process."""
    module = _EDITABLE_IMPORT_MODULES.get(pkg)
    if not module:
        return None
    probe = (
        "import importlib; "
        f"m=importlib.import_module({module!r}); "
        "print(m.__file__ or '')"
    )
    result = run([sys.executable, "-c", probe], capture=True, check=False)
    if result.returncode != 0:
        return ""
    return (result.stdout or "").strip()


# tools/editable_intent.py is piped to the remote interpreter as a script, so
# 100号機 runs the same verification and repoint code. Arguments travel
# as one base64 JSON token so the remote shell cannot re-split them.
REMOTE_EDITABLE_VERIFY = _EDITABLE_INTENT_SOURCE


def _remote_editable_intent(target, argv, timeout=900):
    """Ship one bound Radia intent command with explicit Python and namespace."""
    rc = _check_remote_target(target)
    if rc:
        return rc, None, "Target identity/source check failed before staging the intent helper"
    options = ["--json", "--path-aliases", json.dumps(target["aliases"]),
               "--target-binding", json.dumps(target["binding"]), *argv]
    token = base64.b64encode(json.dumps(options).encode("utf-8")).decode("ascii")
    proc = _run_script_with_file_output(
        ["ssh", target["ssh"], _powershell_literal(target["python"]), "-", "--argv-b64", token],
        REMOTE_EDITABLE_VERIFY, timeout=timeout)
    text = (proc.stdout or "").strip()
    report = None
    if "{" in text:
        try:
            report = json.loads(text[text.index("{"):])
        except ValueError:
            pass
    if proc.stderr and proc.stderr.strip():
        text += ("\n" if text else "") + proc.stderr.strip()
    return proc.returncode, report, text


def _verify_role_editable(role, report=None):
    target = _editable_target(role)
    if _verify_live_editable(target):
        return 1
    argv = ["verify", "--package", "radia", "--expect", "radia=" + target["source"]]
    rc, result, text = _remote_editable_intent(target, argv)
    if result is None:
        fail(target["label"] + " editable verification returned no report")
        if report is not None:report.update(drift=1, unverified=0)
        return 1
    print(editable_intent.format_verify(result))
    counts = result.get("counts", {})
    drift, unverified = int(counts.get("drift", 0)), int(counts.get("unverified", 0))
    if report is not None:report.update(drift=drift, unverified=unverified, ok=int(counts.get("match", 0)))
    return drift + unverified + int(rc != 0 and not (drift or unverified))


def cmd_verify_editable(args):
    """Verify both normal editable roles; isolated QA has its own wheel checks."""
    for role in SOLVER_TARGETS:_editable_target(role)
    return 0 if sum(_verify_role_editable(role) for role in SOLVER_TARGETS) == 0 else 1


def _record_release_intent_remote(role, reason=None):
    target = _editable_target(role)
    argv = ["repoint", "--record-current", "--package", "radia", "--reason",
            reason or "verified normal editable deployment " + target["commit"], "--via", "release controller"]
    rc, result, text = _remote_editable_intent(target, argv)
    if result is not None:print(editable_intent.format_action(result))
    elif text:print(text)
    return rc


def _done_release_source(args):
    """The checkout `done` verifies as the exact release source.

    By default that is the controller itself, which must then sit at the tag.
    ``--release-source`` names a separate exact-tag checkout so a later,
    repaired tooling controller can run `done`; validation intent never selects it.
    """
    selected = getattr(args, "release_source", None)
    return str(Path(selected).resolve()) if selected else str(REPO)


def cmd_repoint(args):
    """Explicit, bound Radia update on either supported normal operating role."""
    role = (getattr(args, "host", None) or "100").strip().lower()
    if role not in SOLVER_TARGETS:
        fail("Unsupported editable target role")
        return 2
    if getattr(args, "rollback", False):
        fail("Editable updates are forward-only")
        return 2
    if not getattr(args, "reason", ""):
        fail("--reason is required for an explicit source move")
        return 2
    packages = list(getattr(args, "package", None) or ["radia"])
    if packages != ["radia"]:
        fail("This solver controller updates radia only")
        return 2
    target = _editable_target(role)
    sources = list(getattr(args, "source", None) or [])
    adopt = getattr(args, "record_current", False)
    if (adopt and sources) or (not adopt and sources != [target["source"]]):
        fail("Source arguments must match the explicit target binding")
        return 2
    if getattr(args, "dry_run", False):
        print(json.dumps(dict(role=role, binding=target["binding"], mode="editable",
                             install=[target["python"], "-m", "pip", "install", "-e", target["source"] + "[gpu]"]), indent=2))
        return 0
    if adopt:
        if _verify_live_editable(target):return 4
        return _record_release_intent_remote(role, args.reason) or _verify_role_editable(role)
    return _deploy_editable_role(role, args.reason)


def cmd_restore_editable(args):
    """Removed: the development tier has no default source to restore."""
    fail("restore-editable was removed (policy P02, 2026-09-11): no canonical target "
         "is assumed, nothing is uninstalled and no process is stopped.")
    print("  If the installed pointer is the intended one, record it:")
    print('    python tools/release_quad.py repoint --record-current --reason "<why>"')
    print("  To move one package explicitly:")
    print('    python tools/release_quad.py repoint --package <name> --source <path> --reason "<why>"')
    return 2


# ============================================================
# Phase 5.5 gate: CI-green BEFORE tagging
# ============================================================
# Routine CI now runs on development/build host. Verify GitHub check-runs by immutable commit SHA;
# never infer remote CI state from a local Runner.Worker process or C:\temp.


def _git_repo_owner_name():
    """Extract 'owner/name' from `git config remote.origin.url`."""
    import re
    url = _git("config", "--get", "remote.origin.url").stdout.strip()
    m = re.search(r"github\.com[:/]([^/\s]+)/([^/\s.]+?)(?:\.git)?$", url)
    if not m:
        raise RuntimeError(f"cannot parse GitHub repo from origin {url!r}")
    return f"{m.group(1)}/{m.group(2)}"


def _check_github_hosted_workflows(
    sha, *, required_names=None, require_present=None,
    excluded_workflows=None,
    timeout_sec=1800, poll_sec=20, registration_grace_sec=90
):
    """Wait for SHA-bound check-runs and require their latest attempts green.

    When ``required_names`` is omitted, every check registered for the commit
    is included. This handles path-scoped monorepo CI: only workflows affected
    by the commit need to exist, but every workflow that does exist must pass.

    ``require_present`` names checks that must EXIST as well as pass, while
    still holding every other registered check to green. Path-scoped CI makes
    absence normal, so absence cannot be an error by default -- but a workflow
    that never runs on push has no registration to inspect, and treating that
    silence as success is what let four release tags burn (4.95.86..4.95.89).

    ``excluded_workflows`` scopes a distribution release by the owning Actions
    workflow, rather than ambiguous job names. Unknown/external checks remain
    required. Workflow lookup errors or mismatched commit identities fail closed.

    Returns (ok: bool, message: str).
    """
    import time as _time
    import re

    repo = _git_repo_owner_name()
    path = f"repos/{repo}/commits/{sha}/check-runs?per_page=100"

    started = _time.time()
    deadline = started + timeout_sec
    last_pending = []
    workflow_by_run_id = {}
    while True:
        try:
            data, _headers = gh_get(path)
        except Exception as e:
            return False, f"GitHub API error: {type(e).__name__}: {e}"

        all_runs = data.get("check_runs", [])
        latest_by_name = {}
        for run in all_runs:
            if excluded_workflows:
                match = re.fullmatch(
                    rf"https://github\.com/{re.escape(repo)}/actions/runs/(\d+)(?:/job/\d+)?/?",
                    run.get("details_url") or "",
                )
                if match:
                    run_id = match.group(1)
                    if run_id not in workflow_by_run_id:
                        try:
                            workflow, _ = gh_get(f"repos/{repo}/actions/runs/{run_id}")
                        except Exception as e:
                            return False, f"Workflow ownership lookup failed: {type(e).__name__}: {e}"
                        if workflow.get("head_sha") != sha or not workflow.get("path"):
                            return False, f"Unverified workflow identity for check-run {run.get('name')}"
                        workflow_by_run_id[run_id] = workflow["path"].split("@", 1)[0]
                    if workflow_by_run_id[run_id] in excluded_workflows:
                        continue
            name = run.get("name")
            if not name:
                continue
            if name not in latest_by_name or run.get("id", 0) > latest_by_name[name].get("id", 0):
                latest_by_name[name] = run

        if required_names is None:
            runs = list(latest_by_name.values())
            missing_names = set() if runs else {"any check-run"}
        else:
            selected_names = set(required_names) | set(require_present or ())
            runs = [latest_by_name[name] for name in selected_names if name in latest_by_name]
            missing_names = set(required_names) - set(latest_by_name)
        if require_present:
            missing_names = missing_names | (
                set(require_present) - set(latest_by_name))

        # Push-triggered workflows take time to register. A commit with no
        # applicable/registered CI is not release evidence.
        if missing_names:
            if _time.time() - started >= registration_grace_sec:
                return False, ("required check-runs not registered for "
                               f"{sha[:8]} after {registration_grace_sec} s: "
                               + ", ".join(sorted(missing_names)))
            _time.sleep(poll_sec)
            continue

        pending = [r for r in runs if r["status"] != "completed"]
        if not pending:
            break
        last_pending = pending
        if _time.time() > deadline:
            return False, ("timeout: github-hosted runs still pending: "
                           + ", ".join(r["name"] for r in last_pending))
        _time.sleep(poll_sec)

    # Every run completed. Green conclusions: success / skipped / neutral.
    must_succeed = set(require_present or ())
    failures = [r for r in runs
                if (r["conclusion"] != "success" if r["name"] in must_succeed
                    else r["conclusion"] not in ("success", "skipped", "neutral"))]
    if failures:
        msg = "; ".join(f"{r['name']}: {r['conclusion']}" for r in failures)
        return False, "github-hosted CI RED -- " + msg
    names = sorted(set(r["name"] for r in runs))
    return True, (f"all {len(runs)} latest SHA-bound check-runs GREEN "
                  f"({', '.join(names)})")


# Job name of the release build in .github/workflows/build-test.yml
# ("Radia Native Release"). That workflow triggers only on a v* tag or a
# manual dispatch -- never on a push to main -- so a green main proves nothing
# about it. ci-verify requires it by name: a SHA that never ran it carries no
# evidence, which is not the same as passing evidence.
RELEASE_CHECK_RUN = "radia-native-release-build"

# These independently released owners cannot block the CPU Radia distribution.
# Keep shared/unknown checks required and the exact native release check mandatory.
INDEPENDENT_DISTRIBUTION_WORKFLOWS = frozenset({
    ".github/workflows/gpu-validation.yml",
    ".github/workflows/radia-simulink.yml",
    ".github/workflows/radia-simulink-release.yml",
    ".github/workflows/ltspice-data-safety.yml",
    ".github/workflows/radia-optuna.yml",
    ".github/workflows/release-radia-optuna.yml",
    ".github/workflows/radia-mcp-matrix.yml",
    ".github/workflows/cubit-mesh-export.yml",
})


def cmd_ci_verify(args):
    """Require CPU Radia/shared checks for the release-tag commit to be green.

    The commit is the one the existing ``v<version>`` tag names, so this runs
    after tagging; it fails while that tag is missing. It also requires the
    release build itself to exist for this SHA, which the dispatched
    'Radia Native Release' run provides before the tag is pushed -- measured
    4.95.86..4.95.89, four tags spent on four missing names without it.
    """

    step("Phase 5.5: CI verify -- SHA-bound GitHub check-runs")
    head_sha = _release_head()
    print(f"  release tag commit = {head_sha[:8]}")
    gh_ok, gh_msg = _check_github_hosted_workflows(
        head_sha, required_names=None, require_present=[RELEASE_CHECK_RUN],
        excluded_workflows=INDEPENDENT_DISTRIBUTION_WORKFLOWS,
        timeout_sec=3600
    )
    print("  " + gh_msg)
    if not gh_ok:
        fail("CI is RED or unverified -- inspect at "
             f"github.com/{_git_repo_owner_name()}/actions, fix-forward.\n"
             f"        If '{RELEASE_CHECK_RUN}' is simply absent, dispatch\n"
             "        'Radia Native Release' on this exact SHA and re-run:\n"
             "        that workflow does not run on a push to main.")
        return 4

    print("")
    ok("CI is GREEN for the exact release-tag commit, including the native release check.")
    return 0


def _run_retired_standalone_pyside_guard():
    """Guard that retired non-Cubit PySide panels have not been reintroduced."""
    step("retired standalone PySide panel guard")
    retired = [
        REPO / "src/radia/radia_gui_base.py",
        REPO / "src/radia/radia_ih.py",
        REPO / "src/radia/radia_em.py",
        REPO / "src/radia/radia_pcb.py",
        REPO / "src/radia/radia_motor.py",
        REPO / "src/radia/radia_streamfunction.py",
        REPO / "src/radia/_heat_panel.py",
    ]
    present = [str(p.relative_to(REPO)) for p in retired if p.exists()]
    pyproject = (REPO / "pyproject.toml").read_text(encoding="utf-8")
    forbidden_tokens = [
        "PySide6>=6.5",
        "radia-ih = \"radia.radia_ih:main\"",
        "radia-em = \"radia.radia_em:main\"",
        "radia-pcb = \"radia.radia_pcb:main\"",
        "radia-streamfunction = \"radia.radia_streamfunction:main\"",
    ]
    present += [f"pyproject token: {t}" for t in forbidden_tokens if t in pyproject]
    if present:
        for item in present:
            fail(f"retired standalone PySide panel still present: {item}")
        return 1
    ok("retired standalone PySide panel files and entry points are absent")
    return 0


def cmd_done(args):
    """Verify both normal editable roles and their isolated exact-wheel QA runtimes.

    Exact-tag, native/source and cross-role checks remain mandatory. Independent
    distributions retain their own release gates.
    """
    step("Definition-of-done check "
         "(preflight + editable tier + phase9 + retired standalone panel guard)")
    rc = cmd_preflight(args)
    if rc != 0:
        fail("preflight failed — repo state not release-ready.")
        return rc

    release_sha = _release_head()
    release_source = _done_release_source(args)
    if _norm_path(release_source) != _norm_path(str(REPO)):
        rc = _verify_release_controller(REPO, release_sha, _radia_version())
        if rc != 0:
            fail("release controller is not a clean descendant of the release commit.")
            return rc
    rc = _verify_local_release_source(release_source, release_sha)
    if rc != 0:
        fail("release source is not the exact clean release SHA.")
        return rc


    rc = _verify_head_release_tag()
    if rc != 0:
        fail("release HEAD is not anchored by its declared Radia version tag.")
        return rc

    if any(_editable_target(role)["commit"] != release_sha for role in SOLVER_TARGETS):
        fail("Normal editable source bindings do not match the exact release tag")
        return 4
    drift = sum(_verify_qa_wheel(role) + _verify_role_editable(role) for role in SOLVER_TARGETS)
    if drift > 0:
        fail(f"{drift} runtime-tier check(s) drifted.  "
             "Run the printed recovery commands, then re-run "
             "`release_quad done`.")
        return 1

    rc = cmd_phase9(args)
    if rc != 0:
        fail("phase9 drift detected — at least one machine is out of sync.")
        return rc

    rc = _run_retired_standalone_pyside_guard()
    if rc != 0:
        fail("retired standalone PySide panel surface reappeared. Remove it, then re-run "
             "`release_quad done`.")
        return rc

    rc = _check_main_synced(hard=True)
    if rc != 0:
        fail("NAS main diverged from origin/main — sync before calling the "
             "release done (this is the root cause of the recurring "
             "rebase-conflict sessions).")
        return rc

    if getattr(args, "simulink_package", None):
        fail("MATLAB/Simulink artifacts belong to radia-simulink, not the CPU Radia release.")
        return 4

    rc = _verify_final_pip_checks()
    if rc != 0:
        return rc

    print("")
    ok("Windows role checks passed. The four-case wrapper still requires its CPU/GPU compute-build receipts. "
       "both normal operating roles use their bound editable source, isolated QA runtimes hold the verified wheel, "
       "and the retired standalone PySide panel surface is absent. "
       "The recorded editable intent names this release source; move it later "
       "only with an explicit `release_quad repoint --reason`.")
    return 0


# ============================================================
# main-sync gate + sync-main + evidence-motor (2026-08-07)
# ============================================================

def _git(*argv, capture=True, check=True):
    safe_repo = REPO.resolve().as_posix()
    return subprocess.run([GIT_EXE, "-c", f"safe.directory={safe_repo}",
                           "-C", str(REPO), *argv],
                          capture_output=capture, text=True, check=check)


def _main_ahead_behind(fetch=True):
    """Return (ahead, behind) of local main vs origin/main."""
    if fetch:
        _git("fetch", "origin", "--quiet", check=False)
    ahead = int(_git("rev-list", "--count", "origin/main..main").stdout.strip())
    behind = int(_git("rev-list", "--count", "main..origin/main").stdout.strip())
    return ahead, behind


def _check_main_synced(*, hard: bool) -> int:
    """NAS main must equal origin/main.

    History divergence is the root cause of the recurring 40-minute
    merge/rebase archaeology (2026-08-07 analysis): releases cut from
    clones push rebased twins while NAS main is left behind.  `done`
    now refuses until main is fast-forward-identical to origin/main.
    """
    ahead, behind = _main_ahead_behind()
    if ahead == 0 and behind == 0:
        ok("NAS main == origin/main (no divergence)")
        return 0
    msg = f"NAS main is ahead {ahead} / behind {behind} of origin/main"
    if not hard:
        warn(msg + " (informational at preflight; `done` enforces it)")
        return 0
    fail(msg + " — the next session will pay for this in rebase conflicts.")
    if ahead == 0:
        info("recovery: git -C " + str(REPO) + " merge --ff-only origin/main")
    else:
        info("recovery: python tools/release_quad.py sync-main")
    return 4


def cmd_sync_main(args):
    """Deterministic fetch -> twin-aware rebase -> preflight -> push.

    Automates the 2026-08-07 manual dance: git rebase drops commits whose
    patches already landed on origin as rebased twins; --empty=drop removes
    ones that become empty.  Genuine conflicts stop the rebase IN PLACE
    with instructions (this tool never resolves content on its own).
    """
    step("sync-main: fetch -> rebase (twin-aware) -> preflight -> push")

    if (REPO / ".git/rebase-merge").exists() or (REPO / ".git/rebase-apply").exists():
        fail("a rebase is already in progress — finish it first "
             "(git rebase --continue / --abort), then re-run sync-main.")
        return 2
    dirty = _git("status", "--porcelain").stdout.strip()
    if dirty:
        fail("working tree not clean — commit this-session files by name "
             "(or stash) before sync-main:")
        for line in dirty.splitlines()[:12]:
            info("  " + line)
        return 2

    ahead, behind = _main_ahead_behind()
    info(f"main vs origin/main: ahead {ahead} / behind {behind}")

    if behind > 0:
        p = _git("rebase", "origin/main", "--empty=drop", check=False)
        print((p.stdout or "").strip())
        if p.returncode != 0:
            print((p.stderr or "").strip()[-1500:])
            conflicts = _git("diff", "--name-only", "--diff-filter=U",
                             check=False).stdout.strip()
            fail("rebase stopped on conflicts — resolve, `git rebase "
                 "--continue`, then re-run sync-main for the push leg:")
            for c in conflicts.splitlines():
                info("  UU " + c)
            return 3
        ok("rebased onto origin/main (patch-identical twins dropped)")

    if getattr(args, "no_push", False):
        ok("sync-main --no-push: stopping before preflight/push as requested")
        return 0

    ahead, _behind = _main_ahead_behind(fetch=False)
    if ahead == 0:
        ok("nothing to push — main already equals origin/main")
        return 0

    p = run([sys.executable, str(REPO / "tools/ci_preflight.py")],
            check=False)
    if p.returncode != 0:
        fail("ci_preflight is RED — fix the printed gate, "
             "then re-run sync-main.")
        return 3

    p = _git("push", "origin", "main", capture=False, check=False)
    if p.returncode != 0:
        fail("push failed (pre-push hook output above names the gate).")
        return 3
    ok("main pushed — origin/main == NAS main")
    return 0


# ---- evidence-motor -----------------------------------------------------

_MOTOR_ARTIFACT = ("validation_test/radia_mcp/artifacts/"
                   "annular_motor_dual_lane_v1/native_motor_angle_family.json")
_MOTOR_PROOF = ("validation_test/radia_mcp/artifacts/"
                "annular_motor_dual_lane_v1/production_replacement_proof.json")
_MOTOR_PROOF_GENERATOR = ("validation_test/radia_mcp/"
                          "generate_production_replacement_proof.py")
# Existing acceptance checks for a regenerated artifact. They are run, never
# edited: a failure means the evidence or the source is wrong.
_MOTOR_ACCEPTANCE_TESTS = (
    "validation_test/radia_mcp/test_motor_angle_source_freshness.py",
    "tests/mcp_integration/test_annular_motor_dual_lane_artifact.py",
    "tests/mcp_integration/test_matlab_source_contract.py",
)
# Complete root-relative closure the HIBINO run needs (learned 2026-08-07:
# missing pyproject.toml / maglev data / team28 docs each cost one round trip).
_MOTOR_SNAPSHOT_ROOTS = (
    "matlab",
    "src/matlab",
    "tests/matlab",
    "validation_test/radia_mcp",
    "validation_test/maglev",
    "docs/maglev/demos/team28",
)
_MOTOR_SNAPSHOT_FILES = ("pyproject.toml",)
_HIBINO_DEST = r"C:\temp\radia_motor_evidence_quad"


def _motor_text_sha(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _motor_pin_status():
    art = json.loads((REPO / _MOTOR_ARTIFACT).read_text(encoding="utf-8"))
    drifted = []
    for path_key, sha_key in (("source_relative_path", "source_sha256"),
                              ("setup_relative_path", "setup_sha256"),
                              ("generator_relative_path", "generator_sha256")):
        rel = art[path_key]
        if _motor_text_sha(REPO / rel) != art[sha_key]:
            drifted.append(rel)
    return art, drifted


def cmd_evidence_motor(args):
    """Check historical motor evidence; CPU MEX regeneration is retired."""
    if getattr(args, "check", False):
        _, drifted = _motor_pin_status()
        for rel in drifted:
            fail("pin drift: " + rel)
        return 4 if drifted else 0
    fail("CPU MEX evidence regeneration is retired. CPU uses Python/pybind11; "
         "GPU capabilities require independent mexcuda validation.")
    return 2


def _motor_post_regeneration(run_fn=None):
    """Regenerate the production proof, then accept; stop at the first failure."""
    run_fn = run if run_fn is None else run_fn
    p = run_fn([sys.executable, _MOTOR_PROOF_GENERATOR], check=False, cwd=str(REPO))
    if p.returncode != 0:
        fail(f"production replacement proof failed (exit {p.returncode})")
        return 4
    ok("production replacement proof regenerated: " + _MOTOR_PROOF)
    return _motor_acceptance(run_fn)


def _motor_acceptance(run_fn=None):
    """Run the motor acceptance tests unchanged; return 0 or a failure code."""
    run_fn = run if run_fn is None else run_fn
    missing = [rel for rel in _MOTOR_ACCEPTANCE_TESTS if not (REPO / rel).is_file()]
    if missing:
        for rel in missing:
            fail("motor acceptance test missing: " + rel)
        return 4
    p = run_fn([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
                *_MOTOR_ACCEPTANCE_TESTS], check=False, cwd=str(REPO))
    if p.returncode != 0:
        fail(f"motor acceptance tests failed (exit {p.returncode})")
        return 4
    ok("motor acceptance tests passed")
    return 0


# ============================================================
# CLI
# ============================================================

def cmd_temp_shadows(args):
    # cmd_done reaches this from callers that never put tools/ on sys.path.
    _tools_dir = str(Path(__file__).resolve().parent)
    if _tools_dir not in sys.path:
        sys.path.insert(0, _tools_dir)
    from release_temp_shadows import inspect_shadows
    if getattr(args, "apply", False):
        merged = run(["git", "-c", f"safe.directory={REPO.as_posix()}",
                      "merge-base", "--is-ancestor", "7ffda79d3", "origin/main"],
                     check=False, capture=True)
        if merged.returncode:
            fail("Kelvin source commits are not confirmed in origin/main; refusing cleanup")
            return 2
    results = inspect_shadows(apply=getattr(args, "apply", False))
    print(json.dumps(results, indent=2))
    return 0 if all(record.get("passed") for record in results.values()) else 2


def main():
    p = argparse.ArgumentParser(prog="release_quad",
                                 description="Verify normal editable deployment and isolated exact-wheel release QA.")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("preflight",
                    help="read-only state report (always safe)")
    s8 = sub.add_parser("phase8",
                         help="deploy and verify the Radia solver on validation / 100号機")
    s8.add_argument("--target", default="validation,100",
                     help="comma list: validation, 100, all (default validation,100)")
    sub.add_parser("phase8e",
                    help="retired compute refresh; always rejected")
    sub.add_parser("phase9",
                    help="cross-machine consistency probe")
    ss = sub.add_parser(
        "simulink-candidate",
        help="verify one extracted Simulink package on the build target")
    ss.add_argument("--package", required=True,
                    help="path to an IH preview or full Radia Simulink ZIP")
    ss.add_argument("--target", default="all",
                    help="comma list: 100, all")
    ss.add_argument("--engine-session", action="append", metavar="HOST=NAME",
                    help="shared Engine name or com:<PID> in the selected host's desktop logon (repeatable)")
    ss.add_argument("--python", action="append", metavar="HOST=PATH",
                    help="candidate interpreter on this target, which runs the verifier and "
                         "becomes MATLAB's RADIA_PYTHON_EXECUTABLE (repeatable); it must "
                         "provide the package's radia release")
    optuna_candidate = sub.add_parser(
        "optuna-candidate",
        help="download one main-CI radia-optuna wheel and verify it on the build target")
    optuna_candidate.add_argument(
        "--ci-run-id", required=True,
        help="successful main-push CI run containing radia-optuna-wheel")
    optuna_candidate.add_argument(
        "--target", default="all",
        help="comma list: 100, all")
    optuna_candidate.add_argument(
        "--engine-session", action="append", metavar="HOST=NAME",
        help="shared Engine name or com:<PID> in the selected host's desktop logon (repeatable); "
             "without it a target runs only when no MATLAB exists there")
    optuna_done = sub.add_parser(
        "optuna-done",
        help="require the exact radia-optuna wheel to have passed both machines")
    optuna_done.add_argument(
        "--wheel", required=True,
        help="retained wheel path emitted by optuna-candidate")
    sub.add_parser("all",
                    help="phase8 -> phase9 in one shot")
    sub.add_parser("verify-editable",
                    help="verify validation's fixed Radia wheel and 100号機's development editable")
    rp = sub.add_parser(
        "repoint",
        help=("explicit editable move: record previous pointer, pip install -e, "
              "verify, record intent (no process stop, no uninstall)"))
    rp.add_argument("--host", default="100",
                    help="100 only; editable installs are not allowed on validation")
    rp.add_argument("--package", action="append", default=[],
                    help="package name (repeatable)")
    rp.add_argument("--source", action="append", default=[],
                    help=("source tree for the package in the same position, "
                          "spelled in the target host's own path namespace"))
    rp.add_argument("--reason", default="", help="recorded reason for the move")
    rp.add_argument("--via", default="release_quad repoint", help="recorded_via label")
    rp.add_argument("--record-current", action="store_true",
                    help="adopt the installed pointers as the recorded intent (no pip)")
    rp.add_argument("--rollback", action="store_true",
                    help="reinstall the previous pointer kept in the record")
    rp.add_argument("--dry-run", action="store_true", help="plan only")
    rp.add_argument("--require-pushed", action="store_true",
                    help=("refuse a source commit no remote ref contains "
                          "(formal handoff / completion evidence only)"))
    sub.add_parser(
        "restore-editable",
        help="removed: no default development source; use repoint with a reason")
    sub.add_parser("ci-verify",
                    help="Phase 5.5: CI-green gate for the commit the existing "
                         "v<version> tag names (requires the tag)")
    sm = sub.add_parser("sync-main",
                        help="fetch -> twin-aware rebase -> preflight -> push")
    sm.add_argument("--no-push", action="store_true",
                    help="stop after the rebase (prep-only)")
    em = sub.add_parser("evidence-motor",
                        help="regenerate the motor MEX evidence artifact on HIBINO")
    em.add_argument("--check", action="store_true",
                    help="verify the SHA pins only (fast, no MATLAB)")
    em.add_argument("--force", action="store_true",
                    help="regenerate even when the pins are current")
    done = sub.add_parser(
        "done",
        help="non-mutating definition-of-done: exact source + phase9 + guards")
    done.add_argument(
        "--simulink-package",
        help="also require a matching two-machine Simulink candidate pass")
    done.add_argument(
        "--release-source",
        help="separate exact-tag checkout to verify as the release source; the "
             "controller running `done` must then be a tracked-clean descendant of "
             "the release commit declaring the same version (default: the "
             "controller itself, which must sit at the tag)")

    shadows = sub.add_parser("temp-shadows", help="verify retired Omega overrides on mdx1/mdx2/hibino")
    shadows.add_argument("--apply", action="store_true", help="remove only unused, non-linked known override trees")
    args = p.parse_args()
    handler = {
        "preflight":        cmd_preflight,
        "temp-shadows":     cmd_temp_shadows,
        "phase8":           cmd_phase8,
        "phase8e":          cmd_phase8e,
        "phase9":           cmd_phase9,
        "simulink-candidate": cmd_simulink_candidate,
        "optuna-candidate":  cmd_optuna_candidate,
        "optuna-done":       cmd_optuna_done,
        "all":              cmd_all,
        "verify-editable":  cmd_verify_editable,
        "repoint":          cmd_repoint,
        "restore-editable": cmd_restore_editable,
        "ci-verify":        cmd_ci_verify,
        "sync-main":        cmd_sync_main,
        "evidence-motor":   cmd_evidence_motor,
        "done":             cmd_done,
    }[args.cmd]
    raise SystemExit(handler(args))


if __name__ == "__main__":
    main()
