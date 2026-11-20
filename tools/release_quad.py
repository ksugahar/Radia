#!/usr/bin/env python
"""release_quad.py — Radia release candidate verification and verified-wheel deployment.

Walks the numerical solver release phases in order, gating each
phase on the success of the previous one. Refuses to skip steps that
have caused real outages (2026-04-14 incident series).

The build host keeps the only editable development installs; the validation
role runs verified wheels only. No command here installs an editable source.

Usage:
    python tools/release_quad.py preflight
        Read-only: report current state and consistency. Use anytime.

    python tools/release_quad.py deployment-plan
        Read-only: print the four release lanes (build and validation roles)
        and the solver runtime role. Never contacts a host.

    python tools/release_quad.py phase8 --source-sha SHA --manifest FILE [--target validation]
        Install the exact release wheel of the published lane into the
        validation role's normal runtime. The wheel's SHA-256 must match the
        verified candidate receipts before and after staging; a fresh process
        then checks non-editable registration, import origin, RECORD hashes
        bound to the candidate wheel, metadata-pinned NGSolve/Netgen and GPU
        runtime pins, and pip check. The CPU-only build host is never a
        target. Never change radia-mcp or cubit-mesh-export.

    python tools/release_quad.py phase9
        Read-only verified-wheel runtime consistency probe.

    python tools/release_quad.py done --source-sha SHA --manifest FILE [--release-source PATH] [--public-summary FILE]
    python tools/release_quad.py all --source-sha SHA --manifest FILE [--release-source PATH]
        Verify four candidate lanes and their exact-source, exact-artifact
        acceptance evidence for release-owner handoff. Does not deploy.


Exit codes:
    0   success
    2   precondition failure (skip detected)
    3   action failure (external command)
    4   verification mismatch
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
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


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


def _wheel_target(role):
    """Resolve the role's normal verified-wheel runtime from private configuration."""
    label, ssh, python = _solver_target(role)
    setting = "RADIA_RELEASE_HOSTNAME_" + role.upper()
    hostname = os.environ.get(setting, "")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", hostname):
        raise ValueError(setting + " must pin the actual target hostname")
    return dict(role=role, label=label, ssh=ssh, python=python, hostname=hostname)


def _target_check_code(target):
    """Read-only identity checks run in the actual target interpreter before mutation."""
    values = dict(hostname=target["hostname"], python=target["python"])
    token = base64.b64encode(json.dumps(values).encode()).decode("ascii")
    return '''import base64,json,os,pathlib,platform,sys
t=json.loads(base64.b64decode("TOKEN"))
if platform.node().casefold()!=t["hostname"].casefold(): raise SystemExit("Selected hostname mismatch")
expected=pathlib.Path(t["python"]).resolve()
if "windowsapps" in [p.casefold() for p in expected.parts]: raise SystemExit("WindowsApps runtime refused")
if os.path.normcase(os.path.realpath(sys.executable))!=os.path.normcase(str(expected)): raise SystemExit("Selected Python mismatch")
'''.replace("TOKEN", token)


def _ps_python(command):
    return "& " + " ".join(_powershell_literal(str(word)) for word in command)


def _ssh_pwsh(ssh, script):
    encoded = base64.b64encode(script.encode("utf-16le")).decode("ascii")
    return run(["ssh", ssh, "pwsh", "-NoProfile", "-EncodedCommand", encoded], check=False).returncode


def _target_check_powershell(target):
    python = target["python"]
    token = base64.b64encode(_target_check_code(target).encode()).decode("ascii")
    return (_runtime_guard_powershell(python) +
            _ps_python([python, "-I", "-c", f"import base64;exec(base64.b64decode('{token}'))"]) +
            "\nif ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }\n")


def _check_remote_target(target):
    return _ssh_pwsh(target["ssh"], _target_check_powershell(target) + "exit $LASTEXITCODE")


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


def _wheel_runtime_probe_code(version, label, wheel=None):
    """Target-side check: non-editable wheel metadata, import origin and RECORD.

    With ``wheel`` (the staged, hash-verified candidate file), every hashed
    entry of the wheel's own RECORD must equal the installed RECORD hash, so
    the runtime is bound to the exact candidate rather than to its version.
    """
    return """import base64, csv, hashlib, importlib.metadata as m, io, json, pathlib, sysconfig, zipfile
version = VERSION_PLACEHOLDER
wheel = WHEEL_PLACEHOLDER
d = m.distribution('radia')
assert not json.loads(d.read_text('direct_url.json') or '{}').get('dir_info', {}).get('editable', False), 'LABEL_PLACEHOLDER must use a verified wheel, not an editable install'
import radia
assert d.version == version and radia.__version__ == version, 'Radia version mismatch'
origin = pathlib.Path(radia.__file__).resolve()
assert origin == pathlib.Path(d.locate_file('radia/__init__.py')).resolve(), 'Shadowed Radia import'
install_paths = sysconfig.get_paths()
wheel_roots = [pathlib.Path(install_paths[key]).resolve() for key in ('purelib', 'platlib')]
assert any(origin.is_relative_to(root) for root in wheel_roots), 'Radia wheel must be under site-packages'
assert d.files, 'Missing wheel RECORD'
installed = {}
for f in d.files:
    if f.hash is not None:
        digest = hashlib.new(f.hash.mode, pathlib.Path(d.locate_file(f)).read_bytes()).digest()
        actual = base64.urlsafe_b64encode(digest).rstrip(b'=').decode()
        assert actual == f.hash.value, 'Wheel payload mismatch: ' + str(f)
        installed[str(f).replace(chr(92), '/')] = f.hash.mode + '=' + f.hash.value
assert installed, 'Empty wheel hash manifest'
if wheel is not None:
    with zipfile.ZipFile(wheel) as archive:
        records = [n for n in archive.namelist() if n.endswith('.dist-info/RECORD')]
        assert len(records) == 1, 'Candidate wheel has no unique RECORD'
        rows = list(csv.reader(io.TextIOWrapper(archive.open(records[0]), encoding='utf-8')))
    expected = {r[0]: r[1] for r in rows if len(r) >= 2 and r[1] and '.data/' not in r[0]}
    differs = sorted(path for path, value in expected.items() if installed.get(path) != value)
    assert expected and not differs, 'Installed files differ from the candidate wheel: ' + ', '.join(differs[:5])
print('LABEL_PLACEHOLDER wheel verified:', version, len(installed))
""".replace("VERSION_PLACEHOLDER", repr(version)).replace(
        "WHEEL_PLACEHOLDER", repr(wheel)).replace("LABEL_PLACEHOLDER", label)


def _solver_pin_check_code():
    """Target-side check of the metadata-pinned ABI and GPU runtime packages."""
    pins = _solver_dependency_pins()
    return f"""import importlib.metadata as m
assert m.version('ngsolve') == {pins['ngsolve']!r} and m.version('netgen-mesher') == {pins['netgen-mesher']!r}, 'Installed NGSolve/Netgen differ from package metadata pins'
import ngsolve
assert ngsolve.__version__ == {pins['ngsolve']!r}, 'Loaded NGSolve differs from package metadata pin'
from radia._cuda_runtime import _pins
for name, pin in _pins().items():
    assert m.version(name) == pin, 'GPU runtime metadata pin mismatch: ' + name
"""


def _release_wheel(manifest, source_sha):
    """Return the one radia wheel of the published lane, re-hashed after receipt checks."""
    validated = _verify_candidate_receipts(manifest, source_sha)
    base = manifest.resolve().parent
    prefix = "radia-" + _radia_version() + "-"
    wheels = [(rel, digest) for rel, digest in validated[_release_acceptance.PUBLISHED_WHEEL_LANE].items()
              if rel.endswith(".whl") and rel.rsplit("/", 1)[-1].startswith(prefix)]
    if len(wheels) != 1:
        raise ValueError(f"Expected exactly one {prefix}*.whl in the published lane, found {len(wheels)}")
    rel, digest = wheels[0]
    path = (base / rel).resolve()
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError(f"SHA-256 mismatch for {rel}")
    return path, digest


def _wheel_install_program(target, wheel, sha256):
    """Hash-check the staged wheel, install it, then verify in a fresh process."""
    python = target["python"]
    dependencies = [python, "-I", "-m", "pip", "install", "--no-cache-dir", "--only-binary=:all:",
                    wheel + "[gpu]"]
    exact = [python, "-I", "-m", "pip", "install", "--no-cache-dir", "--only-binary=:all:",
             "--force-reinstall", "--no-deps", wheel]
    # This process may predate the install; only a fresh target interpreter certifies it.
    verify = [python, "-I", "-c", _target_check_code(target)
              + _wheel_runtime_probe_code(_radia_version(), target["label"], wheel)
              + _solver_pin_check_code()]
    check = [python, "-I", "-m", "pip", "check"]
    return ("import hashlib,subprocess\n" + _target_check_code(target) +
            f"with open({wheel!r},'rb') as staged: staged_sha256=hashlib.sha256(staged.read()).hexdigest()\n"
            f"if staged_sha256!={sha256!r}:"
            " raise SystemExit('Staged wheel SHA-256 differs from the release candidate')\n" +
            "exec(" + repr(SOLVER_INSTALL_GUARD) + ")\n" +
            "".join(f"subprocess.run({command!r},check=True)\n"
                    for command in (dependencies, exact, verify, check)))


def _deploy_wheel(role, wheel, sha256):
    """Stage, install and verify one exact wheel; remove and verify the staging copy."""
    target = _wheel_target(role)
    if _check_remote_target(target):
        return 3
    staging = "C:/temp/radia-release-wheel-" + sha256[:16]
    remote = staging + "/" + wheel.name
    literal = _powershell_literal(staging)
    if _ssh_pwsh(target["ssh"], f"New-Item -ItemType Directory -Force -Path {literal} | Out-Null; exit 0"):
        fail(target["label"] + ": could not create the wheel staging directory")
        return 3
    try:
        if run(["scp", "-o", "BatchMode=yes", str(wheel), f"{target['ssh']}:{remote}"], check=False).returncode:
            fail(target["label"] + ": wheel transfer failed")
            return 3
        # Compression is transport only: avoid the Windows command-line limit
        # while executing the same inspectable checks.
        import zlib
        program = _wheel_install_program(target, remote, sha256)
        token = base64.b64encode(zlib.compress(program.encode("utf-8"))).decode("ascii")
        invocation = [target["python"], "-I", "-c",
                      f"import base64,zlib;exec(zlib.decompress(base64.b64decode('{token}')))"]
        script = _runtime_guard_powershell(target["python"]) + _ps_python(invocation) + "\nexit $LASTEXITCODE"
        if len(base64.b64encode(script.encode("utf-16le"))) > 28000:
            raise ValueError("Wheel installer exceeds the Windows command-line transport budget")
        if _ssh_pwsh(target["ssh"], script):
            fail(target["label"] + " wheel hash/installation/ABI/pip verification failed")
            return 3
    finally:
        if _ssh_pwsh(target["ssh"], f"if (Test-Path -LiteralPath {literal}) {{ Remove-Item -LiteralPath {literal}"
                     f" -Recurse -Force }}; if (Test-Path -LiteralPath {literal}) {{ exit 1 }}; exit 0"):
            fail(target["label"] + ": staged wheel copy was not removed")
    ok(f"{target['label']} holds the verified Radia wheel {wheel.name} ({sha256[:12]})")
    return 0


def cmd_phase8(args):
    """Install and verify the exact release wheel on the validation role."""
    role = (args.target or "validation").strip().lower()
    if role not in SOLVER_TARGETS:
        fail("unsupported solver target; the CPU-only build host is never a deployment target")
        return 2
    _wheel_target(role)
    source_sha = args.source_sha
    if not re.fullmatch(r"[0-9a-f]{40}", source_sha or ""):
        fail("A full lowercase source SHA is required")
        return 4
    if source_sha != _release_head():
        fail("The candidate source SHA differs from the exact release tag")
        return 4
    rc = cmd_preflight(args)
    if rc != 0:
        return rc
    try:
        wheel, digest = _release_wheel(Path(args.manifest), source_sha)
    except (OSError, KeyError, TypeError, ValueError, AttributeError) as exc:
        fail(f"Release candidate wheel failed verification: {exc}")
        return 4
    return _deploy_wheel(role, wheel, digest)


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


def _wheel_probe_prefix(target):
    return _target_check_code(target) + '''import importlib.metadata as md,json,pathlib
d=md.distribution("radia")
if json.loads(d.read_text("direct_url.json") or "{}").get("dir_info",{}).get("editable",False): raise SystemExit("The validation runtime must hold the verified wheel, not an editable install")
import radia
if pathlib.Path(radia.__file__).resolve()!=pathlib.Path(d.locate_file("radia/__init__.py")).resolve(): raise SystemExit("Wheel import is shadowed")
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
    files) after the explicit runtime identity guards.
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
    """Verified-wheel runtime consistency probe."""
    # Resolve every private binding before probing any target.
    targets = [_wheel_target(role) for role in SOLVER_TARGETS]
    expected = _radia_version()
    step("Phase 9: verified-wheel runtime consistency")
    outputs = []
    for target in targets:
        label = target["label"]
        out = _probe(label, ["ssh", target["ssh"], _powershell_literal(target["python"]), "-"],
                     _wheel_probe_prefix(target) + CROSS_MACHINE_PROBE)
        if not out:
            fail(f"could not collect probe data from {label}")
            return 4
        try:
            rows = _parse_phase9_probe(label, out)
        except ValueError as exc:
            fail(str(exc))
            return 4
        if rows["VER radia"] != expected:
            fail(f"{label} runs Radia {rows['VER radia']}, not the release {expected}")
            return 4
        outputs.append((label, rows))

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
    ok("all fields match across the verified-wheel runtimes — release verified.")
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
                                 description="Verify release candidates and deploy verified wheels to the validation role.")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("preflight",
                    help="read-only state report (always safe)")
    sub.add_parser("deployment-plan",
                    help="print the four release lanes and the solver runtime role (read-only)")
    s8 = sub.add_parser("phase8",
                         help="install and verify the exact release wheel on the validation role")
    s8.add_argument("--target", choices=tuple(SOLVER_TARGETS), default="validation")
    s8.add_argument("--source-sha", required=True, help="full release source commit SHA")
    s8.add_argument("--manifest", required=True, help="candidate evidence manifest JSON")
    sub.add_parser("phase9",
                    help="read-only verified-wheel runtime consistency probe")
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
        "ci-verify":        cmd_ci_verify,
        "sync-main":        cmd_sync_main,
        "evidence-motor":   cmd_evidence_motor,
        "done":             cmd_done,
    }[args.cmd]
    raise SystemExit(handler(args))


if __name__ == "__main__":
    main()
