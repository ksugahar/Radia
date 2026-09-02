#!/usr/bin/env python3
"""Run Radia's fast pre-push contracts on mdx for an unpushed commit.

GitHub Actions can only inspect commits already on GitHub. This helper closes
that gap for ``pre-push``: it sends the candidate commit delta as a Git bundle
to mdx, checks it out in an isolated temporary worktree, and runs the
deterministic fast gates there.
"""

from __future__ import annotations

import argparse
import base64
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import uuid


ROOT = Path(__file__).resolve().parents[1]
ZERO = "0" * 40


def run(command: list[str], *, cwd: Path | None = None) -> None:
    """Run a command while preserving its output for the pre-push caller."""
    completed = subprocess.run(command, cwd=cwd, text=True)
    if completed.returncode:
        raise RuntimeError(f"command failed ({completed.returncode}): {' '.join(command)}")


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


def remote_command(script: str, host: str) -> None:
    encoded = base64.b64encode(script.encode("utf-16le")).decode("ascii")
    run(["ssh", host, "pwsh", "-NoProfile", "-ExecutionPolicy", "Bypass",
         "-EncodedCommand", encoded])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True, help="remote SHA before the push")
    parser.add_argument("--head", required=True, help="local candidate SHA")
    parser.add_argument(
        "--host", default=os.environ.get("RADIA_PREFLIGHT_HOST", "mdx"),
        help="SSH host for the dedicated preflight worker (default: mdx)",
    )
    args = parser.parse_args(argv)

    run_id = uuid.uuid4().hex
    remote_root = r"C:\temp\radia-preflight"
    remote_bundle = f"{remote_root}\\incoming\\{run_id}.bundle"
    remote_ref = f"refs/radia-preflight/{run_id}"
    remote_worktree = f"{remote_root}\\work\\{run_id}"

    try:
        with tempfile.TemporaryDirectory(prefix="radia-preflight-", dir=r"C:\temp") as temp:
            bundle = Path(temp) / "candidate.bundle"
            effective_base = create_bundle(args.base, args.head, bundle, remote_ref)
            remote_command(
                "\n".join(
                    [
                        "$ErrorActionPreference = 'Stop'",
                        f"New-Item -ItemType Directory -Force -Path '{remote_root}\\incoming', '{remote_root}\\work' | Out-Null",
                    ]
                ),
                args.host,
            )
            run(["scp", str(bundle), f"{args.host}:{remote_bundle.replace(chr(92), '/')}"])

        script = f"""
$ErrorActionPreference = 'Stop'
$git = 'C:\\actions-runner\\tools\\PortableGit\\bin\\git.exe'
$system_python = 'C:\\Program Files\\Python312\\python.exe'
$root = '{remote_root}'
$repo = Join-Path $root 'repo'
$bundle = '{remote_bundle}'
$worktree = '{remote_worktree}'
$ref = '{remote_ref}'
$head = '{args.head}'
$base = '{effective_base}'
if (-not (Test-Path -LiteralPath $git)) {{ throw "Git is unavailable at $git" }}
if (-not (Test-Path -LiteralPath $system_python)) {{ throw "Python is unavailable at $system_python" }}
if (-not (Test-Path -LiteralPath (Join-Path $repo '.git'))) {{
  & $git clone --no-checkout https://github.com/ksugahar/Radia.git $repo
  if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
}} else {{
  & $git -C $repo fetch --prune origin
  if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
}}
& $git -C $repo fetch $bundle "${{ref}}:${{ref}}"
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
& $git -C $repo worktree prune
& $git -C $repo worktree add --detach $worktree $head
if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
try {{
  Push-Location $worktree
  try {{
    $venv = Join-Path $root 'fast-venv'
    $python = Join-Path $venv 'Scripts\\python.exe'
    if (-not (Test-Path -LiteralPath $python)) {{
      & $system_python -m venv $venv
      if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
      & $python -m pip install --disable-pip-version-check --quiet pytest pyyaml 'mcp>=1.0,<2'
      if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
    }}
    & $python -c "import mcp, pytest, yaml"
    if ($LASTEXITCODE -ne 0) {{
      & $python -m pip install --disable-pip-version-check --quiet pytest pyyaml 'mcp>=1.0,<2'
      if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
    }}
    & $python tools/audit_ci_no_system_install.py
    if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
    & $python tools/ci_preflight.py --since $base
    if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
    & $python -m pytest tests/test_ci_execution_policy.py tests/test_docs_notebook_contract.py tests/test_release_workflow_ref_gate.py tests/test_wheel_package_policy.py tests/test_application_interface_manifest.py tests/axifem/test_docs_notebook_evidence.py -q
    if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
  }} finally {{
    Pop-Location
  }}
}} finally {{
  if (Test-Path -LiteralPath $worktree) {{ & $git -C $repo worktree remove --force $worktree }}
  Remove-Item -LiteralPath $bundle -Force -ErrorAction SilentlyContinue
  & $git -C $repo update-ref -d $ref
}}
"""
        remote_command(script, args.host)
    except (OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"mdx preflight failed: {exc}", file=sys.stderr)
        return 1

    print(f"mdx preflight passed for {args.head}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
