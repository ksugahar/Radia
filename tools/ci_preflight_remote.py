#!/usr/bin/env python3
"""Run Radia's fast pre-push contracts on 100 for an unpushed commit.

GitHub Actions can only inspect commits already on GitHub. This helper closes
that gap for ``pre-push``: it sends the candidate commit delta as a Git bundle
to 100, checks it out as an owned standalone export (no registered worktree), and runs the
deterministic fast gates there.

The build validation host is reached directly over SSH, without GitHub
runner registration. A host-local OS lock protects the shared preflight repo
and environment; isolated worktrees and environments preserve other sessions.
RADIA_PREFLIGHT_SSH_ALIAS may map 100 to an alternate SSH config alias.
"""

from __future__ import annotations

import argparse
import base64
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import uuid


ROOT = Path(__file__).resolve().parents[1]
ZERO = "0" * 40


def resolve_ssh_destination(host: str, mapping: str | None = None) -> str:
    """Resolve transport only; never change the selected validation host."""
    if host != "100":
        raise RuntimeError("Preflight requires the build target")
    if mapping is None:
        mapping = os.environ.get("RADIA_PREFLIGHT_SSH_ALIAS", "")
    if not mapping.strip():
        return host
    aliases = {}
    for entry in mapping.split(","):
        label, separator, alias = entry.strip().partition("=")
        label, alias = label.strip(), alias.strip()
        if (not separator or label != "100"
                or label in aliases
                or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", alias) is None):
            raise RuntimeError("Invalid RADIA_PREFLIGHT_SSH_ALIAS mapping")
        aliases[label] = alias
    if host not in aliases:
        raise RuntimeError(f"RADIA_PREFLIGHT_SSH_ALIAS has no mapping for selected host {host}")
    return aliases[host]


def run(command: list[str], *, cwd: Path | None = None,
        display: list[str] | None = None) -> None:
    """Run a command while preserving its output for the pre-push caller.

    ``display`` replaces the argv in failure diagnostics when it holds data
    that must not reach logs.
    """
    completed = subprocess.run(command, cwd=cwd, text=True)
    if completed.returncode:
        shown = ' '.join(command if display is None else display)
        raise RuntimeError(f"command failed ({completed.returncode}): {shown}")


def git_output(*args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=ROOT, text=True, encoding="utf-8"
    ).strip()


def resolve_candidate_base(base: str, head: str) -> str:
    """Resolve the comparison base, including a new remote branch push."""
    if git_output("rev-parse", "--verify", f"{head}^{{commit}}") != head:
        raise RuntimeError(f"candidate is not a commit: {head}")
    if base != ZERO:
        git_output("rev-parse", "--verify", f"{base}^{{commit}}")
        return base

    for main_ref in ("origin/main", "main"):
        try:
            git_output("rev-parse", "--verify", f"{main_ref}^{{commit}}")
            return git_output("merge-base", head, main_ref)
        except subprocess.CalledProcessError:
            continue
    raise RuntimeError("cannot resolve origin/main or main for a new branch push")


def create_bundle(base: str, head: str, bundle: Path, ref: str) -> str:
    """Create a bundle containing exactly the candidate side of the push."""
    effective_base = resolve_candidate_base(base, head)

    run(["git", "update-ref", ref, head], cwd=ROOT)
    try:
        run(
            ["git", "bundle", "create", str(bundle), f"{effective_base}..{ref}"],
            cwd=ROOT,
        )
    finally:
        run(["git", "update-ref", "-d", ref], cwd=ROOT)
    return effective_base


FORWARDED_ENV = ("RADIA_PUBLIC_NAME_RULES_FILE", "RADIA_OPTUNA_SOURCE")


def forwarded_env_lines(environ=os.environ) -> list[str]:
    """Pass private-policy and pinned-source locations to the remote preflight.

    Only locations are forwarded, never file contents; a missing path fails
    loudly on the remote side like a missing local setting does.
    """
    lines = []
    for name in FORWARDED_ENV:
        value = environ.get(name, "").strip()
        if value:
            quoted = value.replace("'", "''")
            lines.append(f"$env:{name} = '{quoted}'")
    return lines


def remote_command(script: str, host: str) -> None:
    encoded = base64.b64encode(script.encode("utf-16le")).decode("ascii")
    launcher = ["ssh", host, "pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass",
                "-EncodedCommand"]
    # The encoded script carries forwarded private locations; never echo it.
    run([*launcher, encoded], display=[*launcher, "<script>"])


# PowerShell helpers for the build-host candidate export.  Plain text (not an
# f-string) so the functions can also run on their own in an isolated test.
EXPORT_FUNCTIONS = r"""
function Assert-PlainParent($path) {
  # The owned parent (work or tmp) and the preflight root above it must be real
  # directories: a junction there would redirect creation and cleanup elsewhere.
  $dir = Split-Path -Parent ([IO.Path]::GetFullPath($path))
  foreach ($level in 1..2) {
    if (Test-Path -LiteralPath $dir) {
      $item = Get-Item -LiteralPath $dir -Force
      if ($item.Attributes -band [IO.FileAttributes]::ReparsePoint) {
        throw "refusing a reparse-point parent for owned scratch: $dir"
      }
    }
    $dir = Split-Path -Parent $dir
  }
}

function Export-Candidate($git, $repo, $export, $head) {
  # Owned standalone export of the exact candidate: objects come from the cached
  # repository through alternates; nothing is registered as a worktree there.
  Assert-PlainParent $export
  if (Test-Path -LiteralPath $export) { throw "candidate export path already exists: $export" }
  & $git init -q $export
  if ($LASTEXITCODE -ne 0) { throw "git init failed for the candidate export ($LASTEXITCODE)" }
  $script:exportOwned = $true  # created by this invocation; only then may cleanup remove it
  $objects = [IO.Path]::GetFullPath((Join-Path (Join-Path $repo '.git') 'objects'))
  $alternates = Join-Path (Join-Path (Join-Path (Join-Path $export '.git') 'objects') 'info') 'alternates'
  # Git reads alternates as LF-terminated UTF-8; a CR would become part of the path.
  [IO.File]::WriteAllText($alternates, $objects + "`n", (New-Object Text.UTF8Encoding($false)))
  & $git -c "safe.directory=$export" -C $export update-ref refs/radia-preflight/candidate $head
  if ($LASTEXITCODE -ne 0) { throw "candidate $head is not readable from the cached repository" }
  & $git -c "safe.directory=$export" -c core.autocrlf=false -C $export checkout -q --detach $head
  if ($LASTEXITCODE -ne 0) { throw "checkout of candidate $head failed ($LASTEXITCODE)" }
  $exported = (& $git -c "safe.directory=$export" -C $export rev-parse HEAD).Trim()
  if ($exported -ne $head) { throw "candidate export is at $exported, expected $head" }
}

function Remove-CandidateExport($root, $export) {
  if (-not $script:exportOwned) { return }  # never remove a path this run did not create
  if (-not (Test-Path -LiteralPath $export)) { return }
  Assert-PlainParent $export
  $item = Get-Item -LiteralPath $export -Force
  $owned = (Split-Path -Parent $item.FullName) -eq [IO.Path]::GetFullPath((Join-Path $root 'work'))
  if (-not $owned -or ($item.Attributes -band [IO.FileAttributes]::ReparsePoint)) {
    throw "refusing to remove candidate export outside the owned work root: $export"
  }
  $links = Get-ChildItem -LiteralPath $export -Recurse -Force -Attributes ReparsePoint -ErrorAction Stop
  if ($links) { throw "candidate export contains links, not removed: $export" }
  Remove-Item -LiteralPath $export -Recurse -Force -ErrorAction Stop
}

function New-PytestRoot($root, $path) {
  # Per-run pytest scratch: refuse an existing path, record ownership on creation.
  Assert-PlainParent $path
  if (Test-Path -LiteralPath $path) { throw "pytest temp root already exists: $path" }
  New-Item -ItemType Directory -Path $path -ErrorAction Stop | Out-Null
  $script:pytestRootOwned = $true
}

function Remove-PytestRoot($root, $path) {
  if (-not $script:pytestRootOwned) { return }  # never remove a path this run did not create
  if (-not (Test-Path -LiteralPath $path)) { return }
  Assert-PlainParent $path
  $item = Get-Item -LiteralPath $path -Force
  $owned = (Split-Path -Parent $item.FullName) -eq [IO.Path]::GetFullPath((Join-Path $root 'tmp'))
  if (-not $owned -or ($item.Attributes -band [IO.FileAttributes]::ReparsePoint)) {
    throw "refusing to remove pytest temp root outside the owned tmp root: $path"
  }
  # pytest leaves its own "current" links; unlink them without following.
  Get-ChildItem -LiteralPath $path -Recurse -Force -Attributes ReparsePoint -ErrorAction Stop |
    Sort-Object { $_.FullName.Length } -Descending |
    ForEach-Object { if ($_.PSIsContainer) { [IO.Directory]::Delete($_.FullName, $false) } else { [IO.File]::Delete($_.FullName) } }
  Remove-Item -LiteralPath $path -Recurse -Force -ErrorAction Stop
}
"""


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, help="remote SHA before the push")
    parser.add_argument("--head", required=True, help="local candidate SHA")
    parser.add_argument(
        "--host", default=os.environ.get("RADIA_PREFLIGHT_HOST", "100"),
        choices=("100",),
        help="Use the build validation host",
    )
    args = parser.parse_args(argv)
    try:
        ssh_destination = resolve_ssh_destination(args.host)
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"100 preflight unavai100le: {exc}", file=sys.stderr)
        return 1
    print(f"Preflight host: {args.host}")

    run_id = uuid.uuid4().hex
    remote_root = r"C:\temp\radia-preflight"
    remote_bundle = f"{remote_root}\\incoming\\{run_id}.bundle"
    remote_ref = f"refs/radia-preflight/{run_id}"
    remote_worktree = f"{remote_root}\\work\\{run_id}"

    try:
        with tempfile.TemporaryDirectory(prefix="radia-preflight-", dir=tempfile.gettempdir()) as temp:
            bundle = Path(temp) / "candidate.bundle"
            effective_base = create_bundle(args.base, args.head, bundle, remote_ref)
            remote_command(
                "\n".join(
                    [
                        "$ErrorActionPreference = 'Stop'",
                        f"New-Item -ItemType Directory -Force -Path '{remote_root}\\incoming', '{remote_root}\\work' | Out-Null",
                    ]
                ),
                ssh_destination,
            )
            run(["scp", str(bundle), f"{ssh_destination}:{remote_bundle.replace(chr(92), '/')}"])

        script = "\n".join(forwarded_env_lines()) + f"""
$ErrorActionPreference = 'Stop'
{EXPORT_FUNCTIONS}
$git = 'C:\\Program Files\\Git\\cmd\\git.exe'
$system_python = 'C:\\Program Files\\Python312\\python.exe'
$root = '{remote_root}'
$repo = Join-Path $root 'repo'
$bundle = '{remote_bundle}'
$worktree = '{remote_worktree}'
$ref = '{remote_ref}'
$head = '{args.head}'
$base = '{effective_base}'
$lock = $null
$script:exportOwned = $false
$script:pytestRootOwned = $false
$pytestRoot = Join-Path (Join-Path $root 'tmp') ([IO.Path]::GetFileName($worktree))
try {{
  # Fail fast on overlapping preflights; the OS releases the lock on process exit.
  $lock = [IO.File]::Open((Join-Path $root 'preflight.lock'),
    [IO.FileMode]::OpenOrCreate, [IO.FileAccess]::ReadWrite, [IO.FileShare]::None)
  if (Get-Process -Name Runner.Worker -ErrorAction SilentlyContinue) {{
    throw 'build runner is occupied; retry after its job finishes'
  }}
if (-not (Test-Path -LiteralPath $git)) {{ throw "Git is unavai100le at $git" }}
if (-not (Test-Path -LiteralPath $system_python)) {{ throw "Python is unavai100le at $system_python" }}
if (-not (Test-Path -LiteralPath (Join-Path $repo '.git'))) {{
  & $git clone --no-checkout https://github.com/ksugahar/Radia.git $repo
  if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
}} else {{
  & $git -C $repo fetch --prune origin
  if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
}}
& $git -C $repo fetch $bundle "${{ref}}:${{ref}}"
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
Export-Candidate $git $repo $worktree $head
  Push-Location $worktree
  try {{
    # ssh reports the account name in a different case than path resolution;
    # a session-owned temp root keeps pytest paths free of the account name.
    New-Item -ItemType Directory -Force -Path (Join-Path $root 'tmp') | Out-Null
    New-PytestRoot $root $pytestRoot
    $env:PYTEST_DEBUG_TEMPROOT = $pytestRoot
    $venv = Join-Path $root 'fast-venv'
    $python = Join-Path $venv 'Scripts\\python.exe'
    if (-not (Test-Path -LiteralPath $python)) {{
      & $system_python -m venv $venv
      if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
      & $python -m pip install --disable-pip-version-check --quiet numpy scipy pytest pyyaml setuptools wheel 'mcp>=1.0,<2'
      if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
    }}
    & $python -c "import mcp, numpy, pytest, scipy, setuptools, wheel, yaml"
    if ($LASTEXITCODE -ne 0) {{
      & $python -m pip install --disable-pip-version-check --quiet numpy scipy pytest pyyaml setuptools wheel 'mcp>=1.0,<2'
      if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
    }}
    & $python tools/audit_ci_no_system_install.py
    if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
    & $python tools/ci_preflight.py --since $base
    if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
    & $python tools/run_test_tier.py --profile fast-contracts --since $base
    if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
  }} finally {{
    Pop-Location
  }}
}} finally {{
  try {{
    if ($null -ne $lock -and (Test-Path -LiteralPath $git)) {{
      Remove-CandidateExport $root $worktree
      Remove-PytestRoot $root $pytestRoot
      if (Test-Path -LiteralPath (Join-Path $repo '.git')) {{ & $git -C $repo update-ref -d $ref }}
    }}
  }} finally {{
    try {{
      Remove-Item -LiteralPath $bundle -Force -ErrorAction SilentlyContinue
    }} finally {{
      if ($null -ne $lock) {{ $lock.Dispose() }}
    }}
  }}
}}
"""
        remote_command(script, ssh_destination)
    except (OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"100 preflight failed: {exc}", file=sys.stderr)
        return 1
    finally:
        # Also remove a partial upload when scp or the remote launcher fails.
        try:
            remote_command(
                f"if (Test-Path -LiteralPath '{remote_bundle}') {{ "
                f"Remove-Item -LiteralPath '{remote_bundle}' -Force -ErrorAction Stop }}\nexit 0",
                ssh_destination,
            )
        except (OSError, RuntimeError, subprocess.CalledProcessError) as exc:
            print(f"preflight upload cleanup needs attention: {exc}", file=sys.stderr)

    print(f"100 preflight passed for {args.head}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
