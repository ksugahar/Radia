#!/usr/bin/env python
"""release_quad.py — Radia release candidate verification and explicit runtime operations.

Walks the numerical solver release phases in order, gating each
phase on the success of the previous one. Refuses to skip steps that
have caused real outages (2026-04-14 incident series).

Usage:
    python tools/release_quad.py preflight
        Read-only: report current state and consistency. Use anytime.

    python tools/release_quad.py deployment-plan
        Read-only: print the four release lanes (build and validation roles)
        and the solver runtime role. Never contacts a host.

    python tools/release_quad.py phase8 [--target validation]
        Install the bound Radia source editable on the validation role, with
        metadata-pinned GPU dependencies and same-source native payloads.
        Exact-wheel QA runs in a separate explicitly configured environment.
        The CPU-only build host is never a target: it prepares CPU candidates
        and installs no GPU dependencies. Never change radia-mcp or
        cubit-mesh-export.

    python tools/release_quad.py phase9
        Editable and isolated exact-wheel QA consistency probe.


    python tools/release_quad.py done --source-sha SHA --manifest FILE [--release-source PATH] [--public-summary FILE]
    python tools/release_quad.py all --source-sha SHA --manifest FILE [--release-source PATH]
        Verify four candidate lanes and their exact-source, exact-artifact
        acceptance evidence for release-owner handoff. Does not deploy.

    python tools/release_quad.py verify-editable
        Read-only: verify the bound normal editable installation.

    python tools/release_quad.py repoint --host <role> --package radia --source <path> --reason "..."
        Explicitly update one bound normal Radia source and verify it before
        recording intent. --record-current adopts an already verified pointer.
        Updates are forward-only; no process is stopped or package uninstalled.

    python tools/release_quad.py restore-editable
        Removed. The development tier has no default target; use repoint
        with a reason (policy P02, 2026-09-11).


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
import re
import shutil
import subprocess
import sys
import tempfile
import tomllib
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
# The three repository packages installed editable on the development tier.
EDITABLE_PACKAGES = ("radia", "cubit-mesh-export", "radia-mcp")
# The solver runtime installs the "[gpu]" extra and verifies GPU payloads, so
# it runs only on the validation role. The CPU-only build host prepares CPU
# candidates; it is neither a deployment nor an acceptance target.
_KNOWN_SOLVER_TARGETS = {
    "validation": ("validation", None, None),
}
# Exact-wheel release acceptance runs only on the validation role
# (release_acceptance.py owns that list). SOLVER_TARGETS are the explicit
# runtime-maintenance targets.
try:
    from release_acceptance import RELEASE_ACCEPTANCE_HOSTS
    import release_acceptance as _release_acceptance
except ModuleNotFoundError:  # Loaded by path in unit tests rather than run as a script.
    _acceptance_spec = importlib.util.spec_from_file_location(
        "radia_release_acceptance", Path(__file__).resolve().with_name("release_acceptance.py"))
    _acceptance_module = importlib.util.module_from_spec(_acceptance_spec)
    _acceptance_spec.loader.exec_module(_acceptance_module)
    RELEASE_ACCEPTANCE_HOSTS = _acceptance_module.RELEASE_ACCEPTANCE_HOSTS
    _release_acceptance = _acceptance_module

SOLVER_TARGETS = dict(_KNOWN_SOLVER_TARGETS)

assert set(RELEASE_ACCEPTANCE_HOSTS) <= set(SOLVER_TARGETS), (
    "release acceptance names a role release_quad cannot reach")


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
    """Read-only: the four release lanes and the solver runtime role; no host access."""
    plan = dict(
        lanes=[dict(lane=name, build_role=build, validation_role=validation)
               for name, (build, validation) in _release_acceptance.RELEASE_LANES.items()],
        published_wheel_lane=_release_acceptance.PUBLISHED_WHEEL_LANE,
        solver_runtime_roles=list(SOLVER_TARGETS),
        acceptance_roles=list(RELEASE_ACCEPTANCE_HOSTS),
    )
    print(json.dumps(plan, indent=2))
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
        ("radia-mcp",         REPO / "packages/radia-mcp/pyproject.toml"),
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
    """Check the normal editable and separate exact-wheel QA runtimes."""
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
    """Check before pip, then verify new editable registration at fresh startup."""
    pins = _solver_dependency_pins()
    args = ["native-check", target["source"], target["commit"], pins["ngsolve"], pins["netgen-mesher"]]
    install = [target["python"], "-I", "-m", "pip", "install", "--no-cache-dir", "--only-binary=:all:", "-e", target["source"] + "[gpu]"]
    check = [target["python"], "-I", "-m", "pip", "check"]
    # This process predates pip's new .pth registration and may have imported
    # the previous source. Only a fresh target interpreter may certify it.
    postinstall = ("import sys\nsys.argv=" + repr(args) + "\n" +
                   _target_check_code(target) + "\n" + EDITABLE_RELEASE_VERIFY)
    verify = [target["python"], "-I", "-c", postinstall]
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
            "subprocess.run(" + repr(verify) + ",check=True)\n" +
            "subprocess.run(" + repr(check) + ",check=True)\n")


def cmd_phase8(args):
    """Deploy and verify the Radia solver on the validation role."""
    target = (args.target or "validation").strip().lower()
    if target not in SOLVER_TARGETS:
        fail("unsupported solver target; the CPU-only build host is never a deployment target")
        return 2
    if _editable_target(target)["commit"] != _release_head():
        raise ValueError("Normal editable source bindings must match the exact release tag")
    rc = cmd_preflight(args)
    if rc != 0:
        return rc
    return _deploy_validation()


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
        fail(f"{drift} field(s) drift across runtimes — release NOT done.")
        return 4
    print("")
    ok("all fields match across the normal and isolated QA runtimes — release verified.")
    return 0


def cmd_all(args):
    """Verify release candidates without deploying them."""
    return cmd_done(args)


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
# the target runs the same verification and repoint code. Arguments travel
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
    """Verify the normal editable role; isolated QA has its own wheel checks."""
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
    """Explicit, bound Radia update on the supported normal operating role."""
    role = (getattr(args, "host", None) or "validation").strip().lower()
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
    ".github/workflows/radia-simulink-a100.yml",
    ".github/workflows/radia-simulink-release.yml",
    ".github/workflows/ltspice-data-safety.yml",
    ".github/workflows/radia-mcp-matrix.yml",
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


RELEASE_LANES = _release_acceptance.RELEASE_LANES
REQUIRED_VALIDATION_CHECKS = _release_acceptance.REQUIRED_VALIDATION_CHECKS


def _verify_candidate_receipts(manifest, source_sha):
    base = manifest.resolve().parent

    def read(rel):
        path = (base / rel).resolve()
        if not path.is_relative_to(base) or not path.is_file():
            raise ValueError(f"Missing file or path outside receipt directory: {rel}")
        return path.read_bytes()

    data = json.loads(manifest.read_text(encoding="utf-8"))
    return _release_acceptance.verify_lane_receipts(data, source_sha, read)


def cmd_done(args):
    """Verify the fixed source and artifacts for release-owner handoff."""
    source_sha = args.source_sha
    if not re.fullmatch(r"[0-9a-f]{40}", source_sha):
        fail("A full lowercase source SHA is required")
        return 4
    rc = _verify_local_release_source(str(Path(args.release_source or REPO).resolve()), source_sha)
    if rc:
        return rc
    try:
        validated = _verify_candidate_receipts(Path(args.manifest), source_sha)
        summary = getattr(args, "public_summary", None)
        if summary:
            data = json.loads(Path(args.manifest).read_text(encoding="utf-8"))
            Path(summary).write_text(json.dumps(
                _release_acceptance.lane_summary(data, validated), indent=2) + "\n", encoding="utf-8")
    except (OSError, KeyError, TypeError, ValueError, AttributeError) as exc:
        fail(f"Release candidate evidence failed: {exc}")
        return 4
    ok(f"All four artifact lanes verified for {source_sha}; ready for release-owner review.")
    if summary:
        info(f"Public lane summary for the acceptance commit: {summary}")
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

def main():
    p = argparse.ArgumentParser(prog="release_quad",
                                 description="Verify normal editable deployment and isolated exact-wheel release QA.")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("preflight",
                    help="read-only state report (always safe)")
    sub.add_parser("deployment-plan",
                    help="print the four release lanes and the solver runtime role (read-only)")
    s8 = sub.add_parser("phase8",
                         help="deploy and verify the Radia solver on the validation role")
    s8.add_argument("--target", choices=tuple(SOLVER_TARGETS), default="validation")
    sub.add_parser("phase9",
                    help="editable and isolated exact-wheel QA consistency probe")
    sub.add_parser("verify-editable",
                    help="verify the validation role's bound Radia editable installation")
    rp = sub.add_parser(
        "repoint",
        help=("explicit editable move: record previous pointer, pip install -e, "
              "verify, record intent (no process stop, no uninstall)"))
    rp.add_argument("--host", choices=tuple(SOLVER_TARGETS), default="validation",
                    help="solver runtime role (validation only)")
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
                        help="check the motor evidence SHA pins; CPU MEX regeneration is retired")
    em.add_argument("--check", action="store_true",
                    help="verify the SHA pins only (fast, no MATLAB)")
    em.add_argument("--force", action="store_true",
                    help="regenerate even when the pins are current")
    for name in ("done", "all"):
        handoff = sub.add_parser(name, help="verify exact candidate receipts for release handoff")
        handoff.add_argument("--source-sha", required=True, help="full fixed source commit SHA")
        handoff.add_argument("--manifest", required=True, help="candidate evidence manifest JSON")
        handoff.add_argument("--release-source", help="clean candidate checkout (default: controller checkout)")
        handoff.add_argument("--public-summary",
                             help="write the public four-lane summary that the publish gate reads")
    args = p.parse_args()
    handler = {
        "preflight":        cmd_preflight,
        "deployment-plan":  cmd_deployment_plan,
        "phase8":           cmd_phase8,
        "phase9":           cmd_phase9,
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
