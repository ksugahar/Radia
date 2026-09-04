#!/usr/bin/env python
"""release_quad.py — orchestrator for the 4-distribution / 4-machine release flow.

Walks Phase 0 -> 9 of the release-quad skill in order, gating each
phase on the success of the previous one. Refuses to skip steps that
have caused real outages (2026-04-14 incident series).

Usage:
    python tools/release_quad.py preflight
        Read-only: report current state and consistency. Use anytime.

    python tools/release_quad.py phase0
        Mandatory clean rebuild of the Cubit plugin (~3-4 min).

    python tools/release_quad.py phase8 [--target lab|100|hibino|all]
        Run Phase 8a..8d on each target: kill Cubit, install by the
        target's tier (LAB/100 editable, hibino PyPI), cubit-plugin-install,
        --verify-only, cubit-smoke-test. Refuses
        to start if Phase 0 has not been done since the last source
        change in src/cubit_plugin/.

    python tools/release_quad.py phase8e
        Upgrade mdx from PyPI. Refuses to run if pip index versions
        radia / cubit-mesh-export don't match the local repo
        (i.e. PyPI hasn't propagated yet). radia-mcp is intentionally
        not installed on mdx -- and is actively uninstalled if a prior
        release left it behind (mdx is a compute consumer, no MCP).

    python tools/release_quad.py phase9
        Cross-machine consistency probe. Final gate.

    python tools/release_quad.py simulink-candidate --package <zip> --target all
        Extract and execute the exact Simulink package on all four MATLAB machines.

    python tools/release_quad.py optuna-candidate --ci-run-id <id> --target all
        Download the exact radia-optuna wheel from one successful main CI run
        and execute its installed-wheel MATLAB/Simulink test on all four machines.

    python tools/release_quad.py optuna-done --wheel <path>
        Require the retained wheel bytes, source commit, CI run, version, and
        four-machine candidate state to agree before tagging/publication.

    python tools/release_quad.py all
        phase8 -> phase8e -> phase9 with all preconditions enforced.

        When the canonical LAB worktree contains parallel WIP, set
        RADIA_RELEASE_EDITABLE_REPO_LAB and RADIA_RELEASE_EDITABLE_REPO_100
        to the LAB and 100-machine views of one clean NAS release worktree.
        Both editable deployments then fail before install unless that
        worktree is tracked-clean and matches the invoking release SHA.

    python tools/release_quad.py done --simulink-package <zip>
        Require both the normal release gate and the matching four-machine
        Simulink candidate state. The active LAB/100号機 editable sources stay
        unchanged so a verified release worktree cannot be replaced by an
        older canonical WIP tree after the final gate.

    python tools/release_quad.py restore-editable
        Recovery command for a failed or interrupted release. Stop active
        MCP transports, reinstall all three packages from canonical
        01_GitHub worktrees, and verify both machines.

The independent radia-optuna lane does not install or deploy Radia, Cubit,
radia-mcp, or NGSolve. It gates one exact CI wheel before publication.

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
import json
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen
from urllib.request import url2pathname


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
NAS_REPO_LAB = "S:/Radia/01_GitHub"
NAS_REPO_100 = r"W:\00_CAE\Radia\01_GitHub"
EDITABLE_REPO_LAB_ENV = "RADIA_RELEASE_EDITABLE_REPO_LAB"
EDITABLE_REPO_100_ENV = "RADIA_RELEASE_EDITABLE_REPO_100"
# Resolve 100号機 through the machine SSH configuration.  Its LAN address may
# change between lab network segments, while the supported host alias remains
# stable and carries the correct user/key settings.
SSH_100 = "100"
SSH_MDX = "mdx"
SSH_HIBINO = "hibino"
PY_HIBINO = "py -3.12"
MATLAB_EXE = r"C:\Program Files\MATLAB\R2026a\bin\matlab.exe"
SIMULINK_GATE_ROOT = Path(r"C:\temp\radia-release-quad")
OPTUNA_GATE_ROOT = SIMULINK_GATE_ROOT / "radia-optuna"
OPTUNA_SUCCESS_MARKER = "RADIA_OPTUNA_WHEEL_SIMULINK_OK"
SIMULINK_TARGETS = {
    "lab": ("LAB", None, "python"),
    "100": ("100号機", SSH_100, "python"),
    "mdx": ("mdx", SSH_MDX, "python"),
    "hibino": ("hibino", SSH_HIBINO, "py -3.12"),
}


def _editable_repo_lab():
    """Return the LAB editable source, allowing an exact release worktree."""
    return os.environ.get(EDITABLE_REPO_LAB_ENV, NAS_REPO_LAB).strip().rstrip("/\\")


def _editable_repo_100():
    """Return the 100-machine view of the same release worktree."""
    return os.environ.get(EDITABLE_REPO_100_ENV, NAS_REPO_100).strip().rstrip("/\\")


def _release_head():
    return _git("rev-parse", "HEAD").stdout.strip().lower()


# ============================================================
# tiny io helpers
# ============================================================

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


def remove_tree(path: Path):
    """Remove a build directory without relying on POSIX rm being on PATH."""
    print(f"  $ remove-tree {path}")
    shutil.rmtree(path)


def copy_file(src: Path, dst: Path):
    """Copy a file without relying on POSIX cp being on PATH."""
    print(f"  $ copy-file {src} {dst}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


# ============================================================
# state inspectors
# ============================================================

def _read_repo_versions():
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
        text = path.read_text(encoding="utf-8")
        m = re.search(r'(?:^version|^__version__)\s*=\s*"([^"]+)"', text, re.M)
        out[label] = m.group(1) if m else None
    return out


def _newest_mtime(root: Path, suffixes):
    latest = 0.0
    for r, dirs, files in os.walk(root):
        # skip build dirs
        dirs[:] = [d for d in dirs
                    if d not in ("build-pyd", "build-ccm", "build", "compact_netgen")]
        for f in files:
            if Path(f).suffix.lower() in suffixes:
                p = Path(r) / f
                try:
                    mt = p.stat().st_mtime
                    if mt > latest:
                        latest = mt
                except OSError:
                    pass
    return latest


def _bundled_plugin_mtime():
    """Newest mtime of bundled .ccm in cubit-mesh-export package.

    Note: the retired Qt5 .ccl target is gone. The Cubit-embedded PySide
    toolbar is Python package data and is checked by deploy probes, not by this
    compiled-plugin freshness gate.
    """
    pkg = REPO / "packages/cubit-mesh-export/src/cubit_mesh_export"
    times = []
    for name in ("cubit_mesh_export.ccm",):
        p = pkg / name
        if p.is_file():
            times.append(p.stat().st_mtime)
    return max(times) if times else 0.0


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
    """Require ``HEAD`` to be the peeled tag for the declared Radia version."""
    version = _read_repo_versions()["radia"]
    head = _release_head()
    tag_commit = _release_tag_commit(version)
    if tag_commit is None:
        fail(f"Radia {version} has no local v{version} release tag")
        return 4
    if tag_commit != head:
        fail(
            f"HEAD {head} is not the v{version} release commit {tag_commit}; "
            "release_quad done must run on the exact tagged source"
        )
        return 4
    ok(f"HEAD is anchored at v{version} ({head[:12]})")
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


def _write_simulink_state(path: Path, state: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")


def _run_simulink_candidate_target(
        key: str, package: Path, package_sha256: str,
        success_marker: str) -> tuple[bool, str]:
    label, host, python_command = SIMULINK_TARGETS[key]
    verifier = REPO / "tools/verify_simulink_release.py"
    if host is None:
        command = [
            sys.executable, str(verifier), str(package),
            "--matlab", MATLAB_EXE,
        ]
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
        created = subprocess.run(
            ["ssh", host, "pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass",
             "-Command", "-"],
            input=prepare,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
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
            )
            if copied.returncode != 0:
                return False, copied.stderr.strip() or copied.stdout.strip()
        invocation = (
            f"& {python_command} '{remote_verifier}' '{remote_package}' "
            f"--matlab '{MATLAB_EXE}'\nexit $LASTEXITCODE\n"
        )
        result = subprocess.run(
            ["ssh", host, "pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass",
             "-Command", "-"],
            input=invocation,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
    output = ((result.stdout or "") + (result.stderr or "")).strip()
    if result.returncode != 0:
        return False, output
    if success_marker not in output or '"status": "passed"' not in output:
        return False, f"{label} did not emit both release success markers:\n{output}"
    return True, output


def cmd_simulink_candidate(args):
    """Verify one extracted Simulink archive on the requested MATLAB machines."""
    package = Path(args.package).resolve()
    step("Simulink candidate gate (LAB / 100号機 / mdx / hibino)")
    try:
        manifest = _simulink_manifest(package)
    except (OSError, RuntimeError, ValueError, zipfile.BadZipFile) as error:
        fail(f"invalid Simulink candidate: {error}")
        return 2
    package_sha256 = _sha256_file(package)
    success_marker = (
        "RADIA_SIMULINK_RELEASE_OK"
        if manifest.get("schema") in {
            "radia.simulink.library-release-manifest.v1",
            "radia.simulink.library-release-manifest.v2",
            "radia.simulink.library-release-manifest.v3",
        }
        else "RADIA_IH_RELEASE_OK"
    )
    state_path = _simulink_state_path(package_sha256)
    state = {
        "schema": "radia.release-quad.simulink-candidate.v1",
        "package": str(package),
        "package_sha256": package_sha256,
        "version": manifest.get("version"),
        "commit": manifest.get("commit"),
        "targets": {},
    }
    if state_path.is_file():
        previous = json.loads(state_path.read_text(encoding="utf-8"))
        if previous.get("package_sha256") == package_sha256:
            state["targets"] = previous.get("targets", {})

    requested = [part.strip().lower() for part in args.target.split(",")]
    if "all" in requested:
        requested = list(SIMULINK_TARGETS)
    unknown = sorted(set(requested) - set(SIMULINK_TARGETS))
    if unknown:
        fail(f"unknown Simulink target(s): {', '.join(unknown)}")
        return 2

    failed = 0
    for key in requested:
        label = SIMULINK_TARGETS[key][0]
        info(f"verifying extracted package on {label}")
        passed, output = _run_simulink_candidate_target(
            key, package, package_sha256, success_marker)
        state["targets"][key] = {
            "label": label,
            "status": "passed" if passed else "failed",
            "verified_at_utc": datetime.now(timezone.utc).isoformat(),
            "output_tail": output[-4000:],
        }
        _write_simulink_state(state_path, state)
        if passed:
            ok(f"Simulink package passed on {label}")
        else:
            failed += 1
            fail(f"Simulink package failed on {label}")
            if output:
                print(output[-4000:])
    if failed:
        return 4
    ok(f"Simulink candidate state: {state_path}")
    return 0


def _verify_simulink_candidate_state(package_arg: str) -> int:
    package = Path(package_arg).resolve()
    try:
        manifest = _simulink_manifest(package)
        package_sha256 = _sha256_file(package)
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
            state.get("commit") != manifest.get("commit"):
        fail("Simulink candidate state does not match the supplied archive")
        return 4
    missing = [key for key in SIMULINK_TARGETS
               if state.get("targets", {}).get(key, {}).get("status") != "passed"]
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
    ok("supplied Simulink candidate passed LAB / 100号機 / mdx / hibino")
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
    head = _release_head()
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
        prefix=f"radia-optuna-ci-{ci_run_id}-", dir=r"C:\temp"
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


def _run_optuna_candidate_target(
        key: str, wheel: Path, wheel_sha256: str) -> tuple[bool, str]:
    label, host, python_command = SIMULINK_TARGETS[key]
    runner = REPO / "packages/radia-optuna/tests/run_installed_wheel_simulink.ps1"
    matlab_tests = REPO / "packages/radia-optuna/tests/matlab"
    if host is None:
        command = [
            "pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File",
            str(runner), "-Wheel", str(wheel),
            "-MatlabExecutable", MATLAB_EXE,
            "-PythonExecutable", sys.executable,
            "-PreverifiedWheelSha256", wheel_sha256,
        ]
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
        created = subprocess.run(
            ["ssh", host, "pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass",
             "-Command", "-"],
            input=prepare,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        if created.returncode != 0:
            return False, created.stderr.strip() or created.stdout.strip()
        sources = [runner, wheel, *sorted(matlab_tests.glob("*.m"))]
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
        remote_runner = remote_root_windows + "\\" + runner.name
        remote_wheel = remote_root_windows + "\\" + wheel.name
        python_parts = python_command.split()
        python_line = (
            "$pythonExe = (& " + " ".join(python_parts)
            + " -c \"import sys; print(sys.executable)\").Trim()\n"
        )
        invocation = (
            "$ErrorActionPreference = 'Stop'\n"
            + python_line
            + f"& '{remote_runner}' -Wheel '{remote_wheel}' "
              f"-MatlabExecutable '{MATLAB_EXE}' -PythonExecutable $pythonExe "
              f"-PreverifiedWheelSha256 '{wheel_sha256}'\n"
            + "exit $LASTEXITCODE\n"
        )
        result = subprocess.run(
            ["ssh", host, "pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass",
             "-Command", "-"],
            input=invocation,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
    output = ((result.stdout or "") + (result.stderr or "")).strip()
    if result.returncode != 0:
        return False, output
    if OPTUNA_SUCCESS_MARKER not in output:
        return False, f"{label} did not emit {OPTUNA_SUCCESS_MARKER}:\n{output}"
    return True, output


def cmd_optuna_candidate(args):
    """Download one main-CI wheel and run it on all four MATLAB machines."""
    step("radia-optuna exact-wheel candidate gate (LAB / 100号機 / mdx / hibino)")
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
        requested = list(SIMULINK_TARGETS)
    unknown = sorted(set(requested) - set(SIMULINK_TARGETS))
    if unknown:
        fail(f"unknown radia-optuna target(s): {', '.join(unknown)}")
        return 2
    failed = 0
    for key in requested:
        label = SIMULINK_TARGETS[key][0]
        info(f"verifying exact installed wheel on {label}")
        passed, target_output = _run_optuna_candidate_target(
            key, wheel, wheel_sha256
        )
        state["targets"][key] = {
            "label": label,
            "status": "passed" if passed else "failed",
            "verified_at_utc": datetime.now(timezone.utc).isoformat(),
            "output_tail": target_output[-4000:],
        }
        _write_simulink_state(state_path, state)
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
    missing = [key for key in SIMULINK_TARGETS
               if state.get("targets", {}).get(key, {}).get("status") != "passed"]
    if missing:
        fail(f"radia-optuna candidate has not passed: {', '.join(missing)}")
        return 4, None
    ok("exact radia-optuna wheel passed LAB / 100号機 / mdx / hibino")
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
    """Read-only state report. Always safe to run."""
    step("Phase preflight: state report")

    # Versions in repo
    v = _read_repo_versions()
    info(f"Repo versions:")
    info(f"  radia              pyproject={v['radia']}  __version__={v['radia.__version__']}")
    info(f"  cubit-mesh-export  pyproject={v['cubit-mesh-export']}  __version__={v['cme.__version__']}")
    info(f"  radia-mcp          pyproject={v['radia-mcp']}")
    info(f"  radia-optuna       pyproject={v['radia-optuna']}  __version__={v['optuna.__version__']}")

    pp_radia = (v["radia"] == v["radia.__version__"])
    pp_cme   = (v["cubit-mesh-export"] == v["cme.__version__"])
    pp_optuna = (v["radia-optuna"] == v["optuna.__version__"])
    if pp_radia: ok("radia pyproject == __init__")
    else:        fail("radia pyproject != __init__ — fix before any release")
    if pp_cme: ok("cubit-mesh-export pyproject == __init__")
    else:      fail("cubit-mesh-export pyproject != __init__ — fix before any release")
    if pp_optuna: ok("radia-optuna pyproject == __init__")
    else:         fail("radia-optuna pyproject != __init__ — fix before any release")
    if not (pp_radia and pp_cme and pp_optuna):
        return 2

    # Cubit plugin freshness
    src_dir = REPO / "src/cubit_plugin"
    src_mtime = _newest_mtime(src_dir, {".cpp", ".cc", ".cxx", ".c", ".h",
                                          ".hpp", ".hh", ".hxx", ".cmake", ".txt"})
    bin_mtime = _bundled_plugin_mtime()
    if src_mtime == 0:
        warn("could not measure src/cubit_plugin/ mtime")
    elif bin_mtime == 0:
        fail("bundled .ccm missing — Phase 0 not done")
    elif bin_mtime + 1 < src_mtime:
        from datetime import datetime
        fail(f"bundled .ccm ({datetime.fromtimestamp(bin_mtime)}) older than "
              f"src/cubit_plugin/ ({datetime.fromtimestamp(src_mtime)}). "
              "Run `python tools/release_quad.py phase0`.")
        return 2
    else:
        ok("bundled plugin .ccm >= src/cubit_plugin/ mtime")

    _check_main_synced(hard=False)

    return 0


def cmd_phase0(args):
    """Clean rebuild of Cubit plugin (.ccm + .pyd).

    Note: the retired Qt5 .ccl target is gone. Phase 0 builds only the
    C++/APREPRO plugin (.ccm + .pyd); the Cubit-embedded PySide toolbar is
    shipped as Python package data.
    """
    step("Phase 0: clean rebuild of Cubit plugin (~2-3 min)")
    build_pyd = REPO / "src/cubit_plugin/build-pyd"
    build_ccm = REPO / "src/cubit_plugin/build-ccm"
    if build_pyd.exists():
        remove_tree(build_pyd)
    if build_ccm.exists():
        remove_tree(build_ccm)

    # Build via the same ps1 we used in the 2026-04-14 manual run.
    ps1 = REPO / "tools/_build_cubit_plugin.ps1"
    if not ps1.is_file():
        fail(f"missing helper script: {ps1}")
        return 3
    run(["pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", str(ps1)])

    # Propagate to the cubit-mesh-export package ONLY (Tier-2, 2026-06-01):
    # cme is the sole shipper of the Cubit plugin binary; radia no longer
    # bundles cubit_mesh_export.ccm, so radia + cme release fully independently.
    # Only .ccm here (the .ccl target was removed in radia 4.80.0; the .pyd
    # is propagated by the full Build.ps1, not this fast .ccm-only phase0).
    for src_name, dst_dirs in [
        ("build-ccm/cubit_mesh_export.ccm", ["packages/cubit-mesh-export/src/cubit_mesh_export"]),
    ]:
        src = REPO / "src/cubit_plugin" / src_name
        if not src.is_file():
            fail(f"build did not produce {src}")
            return 3
        for d in dst_dirs:
            dst = REPO / d / src.name
            copy_file(src, dst)

    ok("Phase 0 complete; .ccm propagated to cubit-mesh-export pkg (radia no longer bundles it)")
    return 0


_CONSOLE_SCRIPT_WRAPPER_RE = (
    r'(?i)^\s*"?[^\"]*pythonw?\.exe"?\s+"?'
    r'[^\"]*(?:mcp-server-|radia[-_])[^\"\s]*\.exe(?:[\"\s]|$)'
)


def _kill_cubit_local():
    info("force-kill any local Cubit process")
    run(["pwsh", "-NoProfile", "-Command",
         "Get-Process -ErrorAction SilentlyContinue | Where-Object { "
         "$_.ProcessName -eq 'coreform_cubit' -or $_.ProcessName -eq 'cubit' "
         "} | ForEach-Object { Stop-Process -Id $_.Id -Force }; "
         "Start-Sleep -Seconds 2"], check=False)


def _kill_mcp_local():
    info("force-kill local MCP and Radia console-script processes")
    run(["pwsh", "-NoProfile", "-Command",
         "Get-CimInstance Win32_Process -ErrorAction SilentlyContinue | "
         "Where-Object { $_.ProcessId -ne $PID -and ("
         "$_.Name -like 'mcp-*' -or $_.Name -like 'radia-*' -or "
         "$_.Name -like 'radia_*' -or "
         "((($_.Name -eq 'python.exe') -or ($_.Name -eq 'pythonw.exe')) -and "
         f"$_.CommandLine -match '{_CONSOLE_SCRIPT_WRAPPER_RE}') "
         ") } | ForEach-Object { Stop-Process -Id $_.ProcessId -Force "
         "-ErrorAction SilentlyContinue }; "
         "Start-Sleep -Seconds 2"], check=False)


def _deploy_lab():
    step("Phase 8 (LAB): kill, install from NAS, plugin install, verify, smoke")
    repo = _editable_repo_lab()
    rc = _verify_local_release_source(repo, _release_head())
    if rc != 0:
        return rc
    _kill_cubit_local()
    _kill_mcp_local()
    run([sys.executable, "-m", "pip", "uninstall", "-y",
         "radia", "cubit-mesh-export", "radia-mcp"], check=False)
    run([sys.executable, "-m", "pip", "install", "--no-deps",
         "--no-cache-dir", "--no-build-isolation",
         "-e", repo,
         "-e", repo + "/packages/cubit-mesh-export",
         "-e", repo + "/packages/radia-mcp"])
    run(["cubit-plugin-install"])
    run(["cubit-plugin-install", "--verify-only"])
    run(["cubit-smoke-test"])
    run(["cubit-toolbar-smoke-test", "--restarts", "2"])
    ok("Phase 8 complete on LAB")
    return 0


def _deploy_editable_remote(ssh_host, label, repo):
    """Editable-install recipe for machines that should read NAS source.

    LAB and 100号機 are the editable tier.  PyPI propagation is not a
    precondition for this tier; hibino/mdx remain the wheel-consumer
    verification tier.
    """
    step(f"Phase 8 ({label}): kill + NAS editable install + plugin install + verify + smoke (over SSH)")
    expected_sha = _release_head()
    safe_repo = repo.replace("\\", "/")
    ps_block = f"""
$ErrorActionPreference = 'Continue'
$sourceHead = (& git -c "safe.directory={safe_repo}" -C "{repo}" rev-parse HEAD).Trim().ToLowerInvariant()
if ($LASTEXITCODE -ne 0 -or $sourceHead -ne "{expected_sha}") {{
  Write-Error "Release source SHA mismatch: expected {expected_sha}, got $sourceHead"
  exit 41
}}
$sourceDirty = (& git -c "safe.directory={safe_repo}" -C "{repo}" status --porcelain --untracked-files=no) -join "`n"
if ($LASTEXITCODE -ne 0 -or $sourceDirty) {{
  Write-Error "Release source has tracked changes: $sourceDirty"
  exit 42
}}
Get-CimInstance Win32_Process -ErrorAction SilentlyContinue | Where-Object {{
  $_.ProcessId -ne $PID -and (
    $_.Name -eq 'coreform_cubit.exe' -or $_.Name -eq 'cubit.exe' -or
    $_.Name -like 'mcp-*' -or
    $_.Name -like 'radia-*' -or $_.Name -like 'radia_*' -or
    ((($_.Name -eq 'python.exe') -or ($_.Name -eq 'pythonw.exe')) -and
      $_.CommandLine -match '{_CONSOLE_SCRIPT_WRAPPER_RE}')
  )
}} | ForEach-Object {{
  Write-Host "Stopping $($_.Name) pid=$($_.ProcessId)"
  Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue
}}
Start-Sleep -Seconds 2
python -m pip uninstall -y radia cubit-mesh-export radia-mcp
python -m pip install --no-deps --no-cache-dir --no-build-isolation -e "{repo}" -e "{repo}\\packages\\cubit-mesh-export" -e "{repo}\\packages\\radia-mcp"
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
cubit-plugin-install --all-users
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
cubit-plugin-install --verify-only
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
cubit-smoke-test
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
"""
    encoded = base64.b64encode(ps_block.encode("utf-16le")).decode("ascii")
    run(["ssh", ssh_host, "pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass",
         "-EncodedCommand", encoded])
    ok(f"Phase 8 complete on {label}")
    return 0


def _check_pypi_propagation(versions, *, include_mcp=True):
    """Refuse to deploy if PyPI hasn't propagated to repo's current versions.

    Returns 0 on success, 2 if any package is stale.  Used by every PyPI
    install target (hibino + mdx) to prevent installing the OLD version
    while CI is still publishing the new one.
    """
    info("checking PyPI propagation...")
    packages = [("radia", versions["radia"]),
                ("cubit-mesh-export", versions["cubit-mesh-export"])]
    if include_mcp:
        packages.append(("radia-mcp", versions["radia-mcp"]))
    for pkg, want in packages:
        p = run(["python", "-m", "pip", "index", "versions", pkg],
                capture=True, check=False)
        first = p.stdout.splitlines()[0] if p.stdout else ""
        if want and want in first:
            ok(f"PyPI {pkg} live at {want}")
        else:
            fail(f"PyPI {pkg} not yet at {want} (got {first!r}). "
                 "Wait for CI / PyPI propagation, then retry.")
            return 2
    return 0


def _deploy_pypi(ssh_host, label, *, include_mcp=True, python_cmd="python", cubit_optional=False):
    """PyPI-install recipe for downstream Cubit-equipped machines.

    Used for hibino and mdx.  hibino gets radia-mcp; mdx is a compute
    consumer and intentionally skips the MCP server package -- and
    actively uninstalls radia-mcp if a prior release left it behind.

    Recipe:
      1. PyPI propagation check (refuse if stale)
      2. force-kill Cubit + mcp-server-*.exe (otherwise pip install blocks
         on locked Scripts/mcp-server-*.exe)
      3. pip install --upgrade --no-cache-dir from PyPI, pinned to the
         repo's current versions
      4. cubit-plugin-install --all-users (regular-file deploy of the
         freshly-installed wheel's plugin; skipped on optional-Cubit
         targets when Cubit is not installed)
      5. cubit-plugin-install --verify-only (sha256 sanity)
      6. cubit-smoke-test (Cubit 2025.12+ -batch run on ih_bem_sample.jou)
    """
    step(f"Phase 8 ({label}): kill + PyPI install + plugin install + verify + smoke (over SSH)")
    v = _read_repo_versions()
    rc = _check_pypi_propagation(v, include_mcp=include_mcp)
    if rc != 0:
        return rc

    v_radia = v["radia"]
    v_cme   = v["cubit-mesh-export"]
    v_mcp   = v["radia-mcp"]
    mcp_pin = f' "radia-mcp=={v_mcp}"' if include_mcp else ""
    # mdx is a compute consumer: radia-mcp must NOT be present there.  Older
    # (pre-policy) releases left radia-mcp installed, so on any
    # include_mcp=False target actively uninstall it rather than merely
    # skipping the install -- "absent" is the enforced invariant, not a
    # passive side effect.  `pip uninstall -y` on an already-absent package
    # is a no-op (exit 0), so this is safe to run unconditionally there.
    mcp_uninstall_ps = (
        "" if include_mcp
        else f"{python_cmd} -m pip uninstall -y radia-mcp\n"
    )
    cubit_optional_ps = "$true" if cubit_optional else "$false"

    ps_block = f"""
$ErrorActionPreference = 'Continue'
Get-CimInstance Win32_Process -ErrorAction SilentlyContinue | Where-Object {{
  $_.ProcessId -ne $PID -and (
    $_.Name -eq 'coreform_cubit.exe' -or $_.Name -eq 'cubit.exe' -or
    $_.Name -like 'mcp-server*' -or
    $_.Name -like 'radia-*' -or $_.Name -like 'radia_*' -or
    ((($_.Name -eq 'python.exe') -or ($_.Name -eq 'pythonw.exe')) -and
      $_.CommandLine -match '{_CONSOLE_SCRIPT_WRAPPER_RE}')
  )
}} | ForEach-Object {{
  Write-Host "Stopping $($_.Name) pid=$($_.ProcessId)"
  Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue
}}
Start-Sleep -Seconds 2
{mcp_uninstall_ps}        # --force-reinstall is mandatory: `pip install --upgrade X==Y` on a
        # machine already at X==Y is a NO-OP and leaves the on-disk files
        # untouched. That bit us 2026-05-26: 100号機 was already at
        # cubit-mesh-export==0.10.1 from a prior NAS-source install, so
        # `--upgrade 0.10.1` did nothing and the worktree binary
        # (1fd45675...) stayed on 100号機 even though PyPI 0.10.1 had a
        # different binary (ef49da18...). Force-reinstall guarantees the
        # PyPI wheel's bytes overwrite whatever is on disk, which is the
        # whole point of "PyPI is the canonical channel" in the 2-tier policy.
{python_cmd} -m pip install --upgrade --force-reinstall --no-deps --no-cache-dir "radia[cubit]=={v_radia}" "cubit-mesh-export=={v_cme}"{mcp_pin}
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
$pyExe = (& {python_cmd} -c "import sys; print(sys.executable)").Trim()
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
$scriptDir = Join-Path (Split-Path -Parent $pyExe) "Scripts"
$env:PATH = $scriptDir + ";" + $env:PATH
$cubitOptional = {cubit_optional_ps}
$cubitExeCandidates = @()
if ($env:CUBIT_PATH) {{ $cubitExeCandidates += (Join-Path $env:CUBIT_PATH "coreform_cubit.exe") }}
$cubitExeCandidates += "C:\\Program Files\\Coreform Cubit 2025.12\\bin\\coreform_cubit.exe"
$cubitFound = $false
foreach ($candidate in $cubitExeCandidates) {{
  if (Test-Path $candidate) {{ $cubitFound = $true; break }}
}}
if (-not $cubitFound -and $cubitOptional) {{
  Write-Host "[WARN] Coreform Cubit 2025.12+ not found; skipping Cubit plugin install on optional-Cubit target."
  exit 0
}}
cubit-plugin-install --all-users
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
cubit-plugin-install --verify-only
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
cubit-smoke-test
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
"""
    encoded = base64.b64encode(ps_block.encode("utf-16le")).decode("ascii")
    run(["ssh", ssh_host, "pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass",
         "-EncodedCommand", encoded])
    ok(f"Phase 8 complete on {label}")
    return 0


def _deploy_100():
    return _deploy_editable_remote(
        SSH_100, "100号機", _editable_repo_100()
    )


def _deploy_hibino():
    return _deploy_pypi(SSH_HIBINO, "hibino", python_cmd=PY_HIBINO, cubit_optional=True)


def cmd_phase8(args):
    """Deploy + verify + smoke on LAB, 100号機, and/or hibino."""
    # precondition: Phase 0 freshness
    rc = cmd_preflight(args)
    if rc != 0:
        fail("preflight failed; refusing Phase 8")
        return rc

    targets = args.target.split(",") if args.target else ["lab", "100"]
    for t in targets:
        t = t.strip().lower()
        if t == "lab":
            rc = _deploy_lab()
            if rc != 0:
                return rc
        elif t in ("100", "100号機", "100goki"):
            rc = _deploy_100()
            if rc != 0:
                return rc
        elif t == "hibino":
            rc = _deploy_hibino()
            if rc != 0:
                return rc
        elif t == "all":
            rc = _deploy_lab()
            if rc != 0:
                return rc
            rc = _deploy_100()
            if rc != 0:
                return rc
            rc = _deploy_hibino()
            if rc != 0:
                return rc
        else:
            fail(f"unknown target: {t!r}")
            return 2
    return 0


def cmd_phase8e(args):
    """Upgrade mdx from PyPI (radia + cubit-mesh-export only)."""
    return _deploy_pypi(SSH_MDX, "mdx", include_mcp=False)


CROSS_MACHINE_PROBE = '''import hashlib, os
import importlib.metadata as md

def hsh_text(p):
    h = hashlib.sha256()
    try:
        d = open(p, "rb").read().replace(b"\\r\\n", b"\\n").replace(b"\\r", b"\\n")
        h.update(d); return h.hexdigest()[:12]
    except Exception: return "MISSING"

def ver(n):
    try: return md.version(n)
    except Exception: return "MISSING"

import radia, cubit_mesh_export
rad = os.path.dirname(radia.__file__)
print(f"VER radia              = {radia.__version__}")
print(f"VER cubit-mesh-export  = {cubit_mesh_export.__version__}")
print(f"VER radia-mcp          = {ver('radia-mcp')}")
print(f"COMPAT cme  -> radia   = [{cubit_mesh_export.COMPAT_RADIA_MIN}, {cubit_mesh_export.COMPAT_RADIA_MAX}]")
print(f"COMPAT rad  -> cme     = [{radia.COMPAT_CUBIT_MESH_EXPORT_MIN}, {radia.COMPAT_CUBIT_MESH_EXPORT_MAX}]")
for r in ["panels/register_toolbar.py",
          "panels/radia_export_menu.py",
          "simulink/application.py",
          "panels/calc_inductance.py",
          "panels/calc_fem_kelvin.py",
          "panels/calc_fem_coilmesh.py"]:
    print(f"SHA radia/{r:35s} = {hsh_text(os.path.join(rad,r))}")
'''


CROSS_MACHINE_PROBE_NO_MCP = CROSS_MACHINE_PROBE.replace(
    "print(f\"VER radia-mcp          = {ver('radia-mcp')}\")",
    'print("VER radia-mcp          = N/A")',
)


# Editable-tier probe (2026-05-28 fix).  LAB/100号機 are editable DEV checkouts, so
# os.path.dirname(radia.__file__) is the *working tree* -- full of
# uncommitted dev WIP that has nothing to do with the released wheel.
# Hashing that guaranteed perpetual false-positive drift (the whole reason
# this variant exists).  Instead, hash each tracked panel file as it exists
# at the RELEASE TAG (v{radia.__version__}) via `git show`: byte-identical to
# what the consumers' wheel was built from, immune to (a) uncommitted WIP and
# (b) post-release commits on main.  Versions/COMPAT come from installed
# metadata (re-synced to the release in Phase 8), identical to the consumer
# probe so the 10 rows line up 1:1 for the row-by-row comparison.
CROSS_MACHINE_PROBE_LAB = '''import hashlib, os, shutil, subprocess
import importlib.metadata as md

def ver(n):
    try: return md.version(n)
    except Exception: return "MISSING"

import radia, cubit_mesh_export
root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(radia.__file__))))
tag = "v" + radia.__version__
git_exe = shutil.which("git")
if git_exe is None:
    raise RuntimeError("Git is required by the editable release probe")

def hsh_git(relpath):
    result = subprocess.run(
        [git_exe, "-c", "safe.directory=" + root, "-C", root,
         "show", tag + ":" + relpath],
        capture_output=True,
    )
    if result.returncode != 0:
        detail = result.stderr.decode(errors="replace").strip()
        raise RuntimeError(f"git show failed for {relpath}: {detail}")
    NL = bytes([10]); CR = bytes([13])
    d = result.stdout.replace(CR + NL, NL).replace(CR, NL)
    h = hashlib.sha256(); h.update(d); return h.hexdigest()[:12]

print(f"VER radia              = {radia.__version__}")
print(f"VER cubit-mesh-export  = {cubit_mesh_export.__version__}")
print(f"VER radia-mcp          = {ver('radia-mcp')}")
print(f"COMPAT cme  -> radia   = [{cubit_mesh_export.COMPAT_RADIA_MIN}, {cubit_mesh_export.COMPAT_RADIA_MAX}]")
print(f"COMPAT rad  -> cme     = [{radia.COMPAT_CUBIT_MESH_EXPORT_MIN}, {radia.COMPAT_CUBIT_MESH_EXPORT_MAX}]")
for r in ["panels/register_toolbar.py",
          "panels/radia_export_menu.py",
          "simulink/application.py",
          "panels/calc_inductance.py",
          "panels/calc_fem_kelvin.py",
          "panels/calc_fem_coilmesh.py"]:
    print(f"SHA radia/{r:35s} = {hsh_git('src/radia/' + r)}")
'''


def _probe(host_label, cmd_prefix, probe_src=CROSS_MACHINE_PROBE):
    """Run the probe on a target (cmd_prefix is the python invocation).

    probe_src defaults to the consumer probe (hashes the installed wheel
    files).  LAB passes CROSS_MACHINE_PROBE_LAB (hashes tracked files at the
    release tag via git) -- see those probe strings for the rationale.
    """
    p = subprocess.run(cmd_prefix, input=probe_src,
                        capture_output=True, text=True, shell=False)
    if p.returncode != 0:
        fail(f"probe failed on {host_label}: {p.stderr.strip()}")
        return None
    return p.stdout


def cmd_phase9(args):
    """Cross-machine consistency probe."""
    targets = [
        ("LAB", ["python", "-"], CROSS_MACHINE_PROBE_LAB),
        ("100号機", ["ssh", SSH_100, "python", "-"], CROSS_MACHINE_PROBE_LAB),
        ("mdx", ["ssh", SSH_MDX, "python", "-"], CROSS_MACHINE_PROBE_NO_MCP),
        ("hibino", ["ssh", SSH_HIBINO, "py", "-3.12", "-"], CROSS_MACHINE_PROBE),
    ]
    step("Phase 9: cross-machine consistency (LAB / 100号機 / mdx / hibino)")
    outputs = []
    for label, cmd_prefix, probe_src in targets:
        out = _probe(label, cmd_prefix, probe_src)
        if not out:
            fail(f"could not collect probe data from {label}")
            return 4
        outputs.append((label, out.strip().splitlines()))

    if not outputs:
        fail("could not collect probe data from any machine")
        return 4

    n = min(len(rows) for _, rows in outputs)
    drift = 0

    def split(s):
        k, _, v = s.partition("=")
        k = k.strip()
        # strip leading "VER " / "SHA " / "COMPAT " from key
        for prefix in ("VER ", "SHA ", "COMPAT "):
            if k.startswith(prefix):
                k = k[len(prefix):]
        return k.strip(), v.strip()

    labels = [label for label, _ in outputs]
    header = f"\n  {'field':<44} | " + " | ".join(f"{label:<18}" for label in labels)
    print(header)
    print("  " + "-" * max(110, len(header) - 2))
    for i in range(n):
        parsed = [split(rows[i]) for _, rows in outputs]
        key = parsed[0][0]
        values = [value for _, value in parsed]
        comparable = [value for value in values if value != "N/A"]
        match = bool(comparable) and all(value == comparable[0] for value in comparable)
        marker = ok.__name__.upper() if match else "DRIFT"
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
    ok("all fields match across LAB / 100号機 / mdx / hibino — release verified.")
    return 0


def cmd_all(args):
    """Run the full deploy + verify chain (phase8 LAB+100+hibino, phase8e mdx, phase9)."""
    rc = cmd_phase8(argparse.Namespace(target="lab,100,hibino"))
    if rc != 0: return rc
    rc = cmd_phase8e(args)
    if rc != 0:
        warn("Phase 8e failed (PyPI not yet live or mdx unreachable). "
              "Phase 9 will be skipped — re-run later when PyPI propagates.")
        return rc
    return cmd_phase9(args)


# ============================================================
# LAB editable-install verifier (POLICY 2026-05-27, CLAUDE.md
# "release 後の LAB editable 再確認")
# ============================================================

# (package name, expected Editable project location prefix on LAB)
# Path-prefix match (case-insensitive, slash-normalised) so a UNC vs
# drive-letter representation of the same NAS path is accepted.
def _lab_editable_packages():
    root = _editable_repo_lab()
    return [
        ("radia", root),
        ("cubit-mesh-export", root + "/packages/cubit-mesh-export"),
        ("radia-mcp", root + "/packages/radia-mcp"),
        # LAB-private; tolerated as missing if pip show says not installed.
        ("mcp-server-document", "S:/mcp-server"),
    ]


def _canonical_lab_editable_packages():
    """Return LAB development pointers, ignoring release-worktree overrides."""
    root = NAS_REPO_LAB.rstrip("/\\")
    return [
        ("radia", root),
        ("cubit-mesh-export", root + "/packages/cubit-mesh-export"),
        ("radia-mcp", root + "/packages/radia-mcp"),
        ("mcp-server-document", "S:/mcp-server"),
    ]


def _remote_100_editable_packages():
    root = _editable_repo_100()
    return [
        ("radia", root),
        ("cubit-mesh-export", root + r"\packages\cubit-mesh-export"),
        ("radia-mcp", root + r"\packages\radia-mcp"),
    ]


def _canonical_remote_100_editable_packages():
    """Return 100号機 development pointers, ignoring release overrides."""
    root = NAS_REPO_100.rstrip("/\\")
    return [
        ("radia", root),
        ("cubit-mesh-export", root + r"\packages\cubit-mesh-export"),
        ("radia-mcp", root + r"\packages\radia-mcp"),
    ]


def _verify_local_release_source(repo, expected_sha):
    """Fail before install unless the editable source is the exact clean SHA."""
    path = Path(repo)
    if not path.is_dir():
        fail(f"LAB editable release source does not exist: {repo}")
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
            "LAB editable release source SHA mismatch: "
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
            "LAB editable release source has tracked changes: "
            f"{tracked_changes or '<status unavailable>'}"
        )
        return 4
    ok(f"LAB editable release source is exact and clean ({got_sha[:10]})")
    return 0


def _norm_path(p):
    """Lower-case + forward-slashes + strip trailing slash so a UNC and
    a drive-letter form of the same NAS path compare equal."""
    p = (p or "").replace("\\", "/").rstrip("/").lower()
    # Both UNC spellings resolve to the LAB S: drive.  Normalize the whole
    # Radia namespace, not just 01_GitHub: release worktrees live under the
    # same share and editable imports report their real UNC location.
    for unc_root in (
        "//192.168.11.100/work/00_cae/radia/",
        "//192.168.121.100/work/00_cae/radia/",
    ):
        p = p.replace(unc_root, "s:/radia/")
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
    "cubit-mesh-export": "cubit_mesh_export",
    "radia-mcp": "radia_mcp",
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


def _verify_lab_editable(packages=None):
    """Check the 4 LAB-editable packages still point at NAS source.

    Returns (n_ok, n_drift, n_missing, details) tuple.  Drift is the
    one we care about for the "release-done" gate -- it means LAB
    cannot dev-loop on that package because edits won't reflect.

    Missing for `mcp-server-document` is OK (LAB-private, may not be
    installed on every developer's box).  Missing for `radia` /
    `cubit-mesh-export` / `radia-mcp` is itself a drift state (LAB
    should always have these editable).
    """
    step("LAB editable verify (POLICY 2026-05-27)")
    n_ok = n_drift = n_missing = 0
    details = []

    packages = packages or _lab_editable_packages()
    for pkg, want_prefix in packages:
        d = _pip_show(pkg)
        if d is None:
            if pkg == "mcp-server-document":
                info(f"{pkg:25s}  not installed  (LAB-private; tolerated)")
                n_missing += 1
            else:
                fail(f"{pkg:25s}  NOT INSTALLED  -- expected editable @ "
                     f"{want_prefix}")
                n_drift += 1
                details.append((pkg, "not_installed", want_prefix))
            continue

        version = d.get("version", "?")
        editable = d.get("editable project location") \
                   or d.get("editable-project-location") \
                   or d.get("location")
        if not editable:
            fail(f"{pkg:25s}  v{version}  -- no Location field in "
                 f"pip show output")
            n_drift += 1
            details.append((pkg, "no_location", "?"))
            continue

        if d.get("editable project location"):
            kind = "editable"
        else:
            kind = "non-editable (Location only)"

        got_norm = _norm_path(editable)
        want_norm = _norm_path(want_prefix)
        if got_norm != want_norm:
            fail(f"{pkg:25s}  v{version}  DRIFT  {kind}\n"
                 f"        got:      {editable}\n"
                 f"        expected: {want_prefix}")
            n_drift += 1
            details.append((pkg, "drift", editable))
            continue

        origin = _fresh_import_origin(pkg)
        if origin is not None and not _norm_path(origin).startswith(want_norm + "/"):
            fail(f"{pkg:25s}  v{version}  IMPORT DRIFT\n"
                 f"        imported: {origin or '<import failed>'}\n"
                 f"        expected below: {want_prefix}")
            n_drift += 1
            details.append((pkg, "import_drift", origin or "<import failed>"))
            continue

        suffix = f"; import -> {origin}" if origin else ""
        ok(f"{pkg:25s}  v{version}  {kind}  -> {editable}{suffix}")
        n_ok += 1

    print("")
    if n_drift == 0:
        ok(f"All {n_ok} LAB-editable package(s) point at LAB source"
           + (f" ({n_missing} LAB-private skipped)" if n_missing else ""))
    else:
        fail(f"{n_drift} LAB-editable package(s) DRIFTED.  Fix:")
        for pkg, why, _got in details:
            print(f"        # {pkg} ({why})")
            print(f"        Get-Process | Where-Object {{ $_.Name -like "
                  f"'mcp-server*' }} | Stop-Process -Force")
            print(f"        pip uninstall -y {pkg}")
            print(f"        pip install -e {dict(packages)[pkg]} "
                  f"--no-deps --no-cache-dir")
        print("")
        print("        See CLAUDE.md \"POLICY (2026-05-27): release 後の "
              "LAB editable 再確認\" for the full recovery procedure.")
    return n_drift


REMOTE_EDITABLE_VERIFY = r'''
import importlib
import importlib.metadata as md
import json
import sys
from urllib.parse import urlparse
from urllib.request import url2pathname

EXPECT = __EXPECT__
MODULES = {
    "radia": "radia",
    "cubit-mesh-export": "cubit_mesh_export",
    "radia-mcp": "radia_mcp",
}

def norm(p):
    p = (p or "").replace("\\", "/").rstrip("/").lower()
    for unc_root in (
        "//192.168.11.100/work/00_cae/radia/",
        "//192.168.121.100/work/00_cae/radia/",
    ):
        p = p.replace(unc_root, "s:/radia/")
    return p

def editable_location(dist):
    direct_url = dist.read_text("direct_url.json")
    if not direct_url:
        return None
    try:
        data = json.loads(direct_url)
    except json.JSONDecodeError:
        return None
    if not data.get("dir_info", {}).get("editable"):
        return None
    url = data.get("url")
    if not url:
        return None
    parsed = urlparse(url)
    if parsed.scheme == "file":
        path = url2pathname(parsed.path)
        if parsed.netloc:
            path = f"//{parsed.netloc}{path}"
        elif sys.platform.startswith("win") and path.startswith("/") and len(path) > 2 and path[2] == ":":
            path = path[1:]
        return path
    return url

bad = 0
for pkg, want in EXPECT:
    try:
        dist = md.distribution(pkg)
    except md.PackageNotFoundError:
        print(f"FAIL {pkg}: not installed; expected editable @ {want}")
        bad += 1
        continue
    got = editable_location(dist)
    if not got:
        print(f"FAIL {pkg}: installed v{dist.version}, but not editable")
        bad += 1
        continue
    if norm(got) != norm(want):
        print(f"FAIL {pkg}: editable drift")
        print(f"  got:      {got}")
        print(f"  expected: {want}")
        bad += 1
        continue
    try:
        origin = importlib.import_module(MODULES[pkg]).__file__ or ""
    except Exception as exc:
        print(f"FAIL {pkg}: import failed: {exc!r}")
        bad += 1
        continue
    if not norm(origin).startswith(norm(want) + "/"):
        print(f"FAIL {pkg}: import drift")
        print(f"  imported: {origin}")
        print(f"  expected below: {want}")
        bad += 1
        continue
    print(f"OK   {pkg}: v{dist.version} editable -> {got}; import -> {origin}")

sys.exit(1 if bad else 0)
'''


def _verify_remote_editable(ssh_host, label, expected):
    """Check remote editable installs by reading package direct_url.json."""
    step(f"{label} editable verify")
    script = REMOTE_EDITABLE_VERIFY.replace("__EXPECT__", repr(expected))
    p = subprocess.run(["ssh", ssh_host, "python", "-"], input=script,
                       capture_output=True, text=True)
    if p.stdout:
        print(p.stdout.rstrip())
    if p.stderr:
        print(p.stderr.rstrip())
    if p.returncode == 0:
        ok(f"{label} editable packages point at NAS source")
        return 0
    fail(f"{label} editable package drift detected")
    return 1


def _verify_100_editable(expected=None):
    return _verify_remote_editable(
        SSH_100, "100号機", expected or _remote_100_editable_packages()
    )


def cmd_verify_editable(args):
    """Standalone editable verifier (no preflight, no phase9).

    Use this between deploys, or any time pip operations have run on
    LAB/100号機 and you want to confirm editable-loop integrity is intact.
    """
    drift = _verify_lab_editable()
    drift += _verify_100_editable()
    return 1 if drift else 0


def _restore_lab_canonical_editable():
    """Reinstall the three-package development tier from canonical LAB source."""
    step("Restore LAB canonical editable installs")
    _kill_mcp_local()
    packages = _canonical_lab_editable_packages()[:3]
    uninstall = run(
        [sys.executable, "-m", "pip", "uninstall", "-y",
         *(pkg for pkg, _path in packages)],
        check=False,
    )
    if uninstall.returncode != 0:
        warn("pip uninstall reported an error; attempting a clean editable install")
    install_cmd = [sys.executable, "-m", "pip", "install", "--no-deps",
                   "--no-cache-dir", "--no-build-isolation"]
    for _pkg, path in packages:
        install_cmd.extend(["-e", path])
    installed = run(install_cmd, check=False)
    if installed.returncode != 0:
        fail("LAB canonical editable reinstall failed")
        return 3
    smoke = run(["mcp-server-grant-writing", "--selftest"], check=False)
    if smoke.returncode != 0:
        fail("LAB grant-writing MCP self-test failed after restore")
        return 4
    ok("LAB canonical editable reinstall and grant-writing self-test passed")
    return 0


def _restore_100_canonical_editable():
    """Reinstall the three-package development tier from canonical 100 source."""
    step("Restore 100号機 canonical editable installs")
    repo = NAS_REPO_100.rstrip("/\\")
    ps_block = fr"""
$ErrorActionPreference = 'Stop'
$env:PYTHONUTF8 = '1'
Get-CimInstance Win32_Process -ErrorAction SilentlyContinue | Where-Object {{
  $_.ProcessId -ne $PID -and (
    $_.Name -like 'mcp-*' -or $_.Name -like 'radia-*' -or $_.Name -like 'radia_*' -or
    ((($_.Name -eq 'python.exe') -or ($_.Name -eq 'pythonw.exe')) -and
      $_.CommandLine -match '{_CONSOLE_SCRIPT_WRAPPER_RE}')
  )
}} | ForEach-Object {{
  Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue
}}
python -m pip uninstall -y radia cubit-mesh-export radia-mcp
python -m pip install --no-deps --no-cache-dir --no-build-isolation -e "{repo}" -e "{repo}\packages\cubit-mesh-export" -e "{repo}\packages\radia-mcp"
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
mcp-server-grant-writing --selftest
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
"""
    encoded = base64.b64encode(ps_block.encode("utf-16le")).decode("ascii")
    restored = run(
        ["ssh", SSH_100, "pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass",
         "-EncodedCommand", encoded],
        check=False,
    )
    if restored.returncode != 0:
        fail("100号機 canonical editable reinstall failed")
        return 3
    ok("100号機 canonical editable reinstall and grant-writing self-test passed")
    return 0


def cmd_restore_editable(args):
    """Restore canonical editable installs after a release or interrupted deploy."""
    failures = 0
    failures += int(_restore_lab_canonical_editable() != 0)
    failures += int(_restore_100_canonical_editable() != 0)
    if failures:
        fail(f"canonical editable restore failed on {failures} machine(s)")
        return 3

    drift = _verify_lab_editable(_canonical_lab_editable_packages())
    drift += _verify_100_editable(_canonical_remote_100_editable_packages())
    if drift:
        fail("canonical editable verification failed after reinstall")
        return 4
    ok("LAB and 100号機 are back on canonical 01_GitHub editable sources")
    return 0


# ============================================================
# Phase 5.5 gate: CI-green BEFORE tagging
# ============================================================
# Routine CI now runs on mdx. Verify GitHub check-runs by immutable commit SHA;
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
    sha, *, required_names=None, timeout_sec=1800, poll_sec=20
):
    """Wait for SHA-bound check-runs and require their latest attempts green.

    When ``required_names`` is omitted, every check registered for the commit
    is included. This handles path-scoped monorepo CI: only workflows affected
    by the commit need to exist, but every workflow that does exist must pass.

    Returns (ok: bool, message: str).
    """
    import time as _time

    repo = _git_repo_owner_name()
    path = f"repos/{repo}/commits/{sha}/check-runs?per_page=100"

    started = _time.time()
    deadline = started + timeout_sec
    last_pending = []
    while True:
        try:
            data, _headers = gh_get(path)
        except Exception as e:
            return False, f"GitHub API error: {type(e).__name__}: {e}"

        all_runs = data.get("check_runs", [])
        latest_by_name = {}
        for run in all_runs:
            name = run.get("name")
            if not name:
                continue
            if name not in latest_by_name or run.get("id", 0) > latest_by_name[name].get("id", 0):
                latest_by_name[name] = run

        if required_names is None:
            runs = list(latest_by_name.values())
            missing_names = set() if runs else {"any check-run"}
        else:
            runs = [latest_by_name[name] for name in required_names if name in latest_by_name]
            missing_names = set(required_names) - set(latest_by_name)

        # Push-triggered workflows take time to register. A commit with no
        # applicable/registered CI is not release evidence.
        if missing_names:
            if _time.time() - started > 90:
                return False, ("required check-runs not registered for "
                               f"{sha[:8]} after 90 s: "
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
    failures = [r for r in runs
                if r["conclusion"] not in ("success", "skipped", "neutral")]
    if failures:
        msg = "; ".join(f"{r['name']}: {r['conclusion']}" for r in failures)
        return False, "github-hosted CI RED -- " + msg
    names = sorted(set(r["name"] for r in runs))
    return True, (f"all {len(runs)} latest SHA-bound check-runs GREEN "
                  f"({', '.join(names)})")


def cmd_ci_verify(args):
    """Require every latest GitHub check-run for HEAD to be green."""

    step("Phase 5.5: CI verify -- SHA-bound GitHub check-runs")
    head_sha = _release_head()
    print(f"  HEAD = {head_sha[:8]}")
    gh_ok, gh_msg = _check_github_hosted_workflows(
        head_sha, required_names=None, timeout_sec=3600
    )
    print("  " + gh_msg)
    if not gh_ok:
        fail("CI is RED or unverified -- inspect at "
             f"github.com/{_git_repo_owner_name()}/actions, fix-forward.")
        return 4

    print("")
    ok("CI is GREEN for this exact commit. Safe to create and push release tags.")
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
    """Run release gates without changing the verified editable sources.

    Exit 0 means the release is consistent across LAB / 100号機 / mdx / hibino,
    the repo is release-ready, the retired non-Cubit PySide panel surface has
    not been reintroduced, AND LAB/100号機 still use the exact clean source
    verified by this command. Returning to the canonical development tree is a
    separate, explicit ``restore-editable`` operation after that tree catches
    up with the published release.
    """
    step("Definition-of-done check "
         "(preflight + editable tier + phase9 + retired standalone panel guard)")
    rc = cmd_preflight(args)
    if rc != 0:
        fail("preflight failed — repo state not release-ready.")
        return rc

    rc = _verify_local_release_source(_editable_repo_lab(), _release_head())
    if rc != 0:
        fail("active LAB editable source is not the exact clean release SHA.")
        return rc

    rc = _verify_head_release_tag()
    if rc != 0:
        fail("release HEAD is not anchored by its declared Radia version tag.")
        return rc

    drift = _verify_lab_editable()
    drift += _verify_100_editable()
    if drift > 0:
        fail(f"{drift} editable-tier check(s) drifted.  "
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
        rc = _verify_simulink_candidate_state(args.simulink_package)
        if rc != 0:
            fail("Simulink candidate did not satisfy the four-machine gate.")
            return rc

    print("")
    suffix = (" The supplied Simulink candidate also passed all four MATLAB "
              "machines." if getattr(args, "simulink_package", None) else "")
    ok("DEFINITION OF DONE met. Release is consistent across LAB / 100号機 / "
       "mdx / hibino, LAB/100号機 remain on the exact verified editable "
       "source, and the retired standalone PySide panel surface is absent. "
       "Run `release_quad restore-editable` explicitly after the canonical "
       "development tree catches up." + suffix)
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
_MOTOR_PYTEST = ("packages/radia-mcp/tests/"
                 "test_annular_motor_dual_lane_artifact.py")
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
    """Regenerate the motor MEX evidence artifact on HIBINO.

    Scripted form of the 2026-08-07 manual run: build the MEX on LAB,
    ship the snapshot closure to hibino over scp, run the generator in
    a SYNCHRONOUS ssh (Windows OpenSSH reaps detached children on
    session exit — Start-Process launches died twice before this was
    understood), fetch the artifact back, verify the SHA pins, and
    align the pytest test-count expectation.
    """
    step("evidence-motor: HIBINO MATLAB evidence regeneration")

    art, drifted = _motor_pin_status()
    if not drifted and not getattr(args, "force", False):
        ok(f"evidence already current (tests={art['test_count']}, "
           f"host={art['execution_environment']['hostname']}) — use --force to regenerate")
        return 0
    if getattr(args, "check", False):
        for rel in drifted:
            fail("pin drift: " + rel)
        return 4 if drifted else 0
    for rel in drifted:
        warn("pin drift: " + rel + " — regenerating")

    # 1) fresh MEX from current source
    # Preserve the caller's mapped-drive spelling on Windows. Path.resolve()
    # expands S:/W: to a UNC path, which cmd.exe cannot use as the working
    # directory of the CMake/Ninja post-build copy commands.
    build_script = Path.cwd() / "Build.ps1"
    if not build_script.is_file():
        build_script = REPO / "Build.ps1"
    p = run(["pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass",
             "-File", str(build_script), "-MatlabMexOnly"], check=False)
    if p.returncode != 0:
        fail("Build.ps1 -MatlabMexOnly failed")
        return 3

    # 2) snapshot zip
    zip_path = Path(r"C:\temp\radia_motor_evidence_quad.zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for f in _MOTOR_SNAPSHOT_FILES:
            z.write(REPO / f, f)
        for root in _MOTOR_SNAPSHOT_ROOTS:
            for dirpath, dirnames, filenames in os.walk(REPO / root):
                dirnames[:] = [d for d in dirnames
                               if d not in ("__pycache__", ".pytest_cache")]
                for fn in filenames:
                    full = Path(dirpath) / fn
                    z.write(full, full.relative_to(REPO).as_posix())
    info(f"snapshot: {zip_path} ({zip_path.stat().st_size/1e6:.1f} MB)")

    # 3) runner (ASCII only; synchronous MATLAB, refuses on busy host)
    gen_rel = r"validation_test\radia_mcp\generate_motor_angle_family_mex_artifact.m"
    runner = "\n".join([
        "$ErrorActionPreference = 'Stop'",
        "if (Get-Process -Name MATLAB -ErrorAction SilentlyContinue) {",
        "    Write-Output 'BUSY: MATLAB already running on hibino'; exit 2 }",
        f"$dest = '{_HIBINO_DEST}'",
        "if (Test-Path $dest) { Remove-Item $dest -Recurse -Force }",
        "Expand-Archive -Path 'C:\\temp\\radia_motor_evidence_quad.zip' -DestinationPath $dest -Force",
        "$env:RADIA_VALIDATION_HOST_ROLE = 'compute'",
        f"$gen = Join-Path $dest '{gen_rel}'",
        "$log = 'C:\\temp\\radia_motor_evidence_quad.log'",
        "Remove-Item $log -Force -ErrorAction SilentlyContinue",
        f"& '{MATLAB_EXE}' -wait -batch (\"run('\" + $gen + \"')\") -logfile $log",
        "Write-Output ('MATLAB_EXIT=' + $LASTEXITCODE)",
        "Get-Content $log -Tail 8 -ErrorAction SilentlyContinue",
        "exit $LASTEXITCODE",
    ]) + "\n"
    runner_path = Path(r"C:\temp\radia_motor_evidence_quad_runner.ps1")
    runner_path.write_text(runner, encoding="ascii")

    for src, dst in ((zip_path, "C:/temp/radia_motor_evidence_quad.zip"),
                     (runner_path, "C:/temp/radia_motor_evidence_quad_runner.ps1")):
        p = run(["scp", "-o", "ConnectTimeout=20", str(src),
                 f"{SSH_HIBINO}:{dst}"], check=False)
        if p.returncode != 0:
            fail(f"scp to hibino failed: {src}")
            return 3

    # 4) synchronous run (previous artifact: 79 tests in ~2 min + startup)
    p = run(["ssh", SSH_HIBINO, "pwsh -NoProfile -ExecutionPolicy Bypass "
             "-File C:\\temp\\radia_motor_evidence_quad_runner.ps1"],
            check=False, timeout=1800)
    if p.returncode != 0:
        fail(f"hibino generator failed (exit {p.returncode}) — see log above")
        return 3

    # 5) fetch artifact back + verify pins
    art_remote = (_HIBINO_DEST.replace("\\", "/") + "/" +
                  _MOTOR_ARTIFACT.replace("\\", "/"))
    p = run(["scp", "-o", "ConnectTimeout=20",
             f"{SSH_HIBINO}:{art_remote}", str(REPO / _MOTOR_ARTIFACT)],
            check=False)
    if p.returncode != 0:
        fail("scp of the regenerated artifact failed")
        return 3
    art, drifted = _motor_pin_status()
    if drifted:
        for rel in drifted:
            fail("pin STILL drifted after regeneration: " + rel)
        return 4
    ok(f"artifact regenerated: {art['status']}, "
       f"{art['passed_count']}/{art['test_count']} on "
       f"{art['execution_environment']['hostname']} ({art['matlab_release']})")

    # 6) align the pytest expectation with the measured suite size
    import re
    pytest_path = REPO / _MOTOR_PYTEST
    text = pytest_path.read_text(encoding="utf-8")
    new_text, n = re.subn(
        r'(assert native\["test_count"\] == native\["passed_count"\] == )\d+',
        rf"\g<1>{art['test_count']}", text)
    if n == 1 and new_text != text:
        pytest_path.write_text(new_text, encoding="utf-8", newline="\n")
        ok(f"pytest expectation aligned to {art['test_count']}")
    info("stage & commit:  git add " + _MOTOR_ARTIFACT + " " + _MOTOR_PYTEST)
    return 0


# ============================================================
# CLI
# ============================================================

def main():
    p = argparse.ArgumentParser(prog="release_quad",
                                 description="Enforce the release-quad flow.")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("preflight",
                    help="read-only state report (always safe)")
    sub.add_parser("phase0",
                    help="clean rebuild of the Cubit plugin")
    s8 = sub.add_parser("phase8",
                         help="deploy + verify + smoke on LAB / 100号機 / hibino")
    s8.add_argument("--target", default="lab,100",
                     help="comma list: lab, 100, hibino, all (default lab,100)")
    sub.add_parser("phase8e",
                    help="upgrade mdx from PyPI (after PyPI propagation)")
    sub.add_parser("phase9",
                    help="cross-machine consistency probe")
    ss = sub.add_parser(
        "simulink-candidate",
        help="verify one extracted Simulink package on four MATLAB machines")
    ss.add_argument("--package", required=True,
                    help="path to an IH preview or full Radia Simulink ZIP")
    ss.add_argument("--target", default="all",
                    help="comma list: lab, 100, mdx, hibino, all")
    optuna_candidate = sub.add_parser(
        "optuna-candidate",
        help="download one main-CI radia-optuna wheel and verify it on four MATLAB machines")
    optuna_candidate.add_argument(
        "--ci-run-id", required=True,
        help="successful main-push CI run containing radia-optuna-wheel")
    optuna_candidate.add_argument(
        "--target", default="all",
        help="comma list: lab, 100, mdx, hibino, all")
    optuna_done = sub.add_parser(
        "optuna-done",
        help="require the exact radia-optuna wheel to have passed all four machines")
    optuna_done.add_argument(
        "--wheel", required=True,
        help="retained wheel path emitted by optuna-candidate")
    sub.add_parser("all",
                    help="phase8 -> phase8e -> phase9 in one shot")
    sub.add_parser("verify-editable",
                    help="LAB/100号機 editable-install pointers check (read-only)")
    sub.add_parser(
        "restore-editable",
        help="stop MCP transports and restore LAB/100号機 to canonical editable sources")
    sub.add_parser("ci-verify",
                    help="Phase 5.5: SHA-bound CI-green gate (after push main, before tag)")
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
        help="also require a matching four-machine Simulink candidate pass")

    args = p.parse_args()
    handler = {
        "preflight":        cmd_preflight,
        "phase0":           cmd_phase0,
        "phase8":           cmd_phase8,
        "phase8e":          cmd_phase8e,
        "phase9":           cmd_phase9,
        "simulink-candidate": cmd_simulink_candidate,
        "optuna-candidate":  cmd_optuna_candidate,
        "optuna-done":       cmd_optuna_done,
        "all":              cmd_all,
        "verify-editable":  cmd_verify_editable,
        "restore-editable": cmd_restore_editable,
        "ci-verify":        cmd_ci_verify,
        "sync-main":        cmd_sync_main,
        "evidence-motor":   cmd_evidence_motor,
        "done":             cmd_done,
    }[args.cmd]
    raise SystemExit(handler(args))


if __name__ == "__main__":
    main()
