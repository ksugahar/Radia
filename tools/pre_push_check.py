#!/usr/bin/env python3
"""Check every branch candidate from Git pre-push stdin on the wheel test host."""
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).absolute().parents[1]
ZERO = "0" * 40


def run_local_candidate(base, head):
    """Local checks run on a detached candidate, never another session's WIP."""
    git = ["git", "-c", f"safe.directory={ROOT}"]
    try:
        if base == ZERO:
            base = subprocess.check_output(
                [*git, "merge-base", head, "origin/main"], cwd=ROOT,
                text=True, encoding="utf-8").strip()
        subprocess.run([*git, "cat-file", "-e", f"{base}^{{commit}}"],
                       cwd=ROOT, check=True, capture_output=True, timeout=30)
    except subprocess.SubprocessError as exc:
        print(f"pre-push: cannot resolve base {base}; fetch the remote base commit/history and retry: {exc}", file=sys.stderr)
        return 1
    parent = Path(tempfile.mkdtemp(prefix="radia-push-candidate-"))
    checkout = parent / "candidate"
    added = False
    result = 1
    try:
        subprocess.run([*git, "worktree", "add", "--detach", str(checkout), head],
                       cwd=ROOT, check=True, timeout=60)
        added = True
        result = subprocess.run(
            [sys.executable, "tools/ci_preflight.py", "--local-push-candidate", head,
             "--local-push-base", base], cwd=checkout, timeout=420).returncode
    except (ValueError, OSError, subprocess.SubprocessError) as exc:
        print(f"pre-push: local candidate check failed: {exc}", file=sys.stderr)
        result = 1
    finally:
        try:
            if added:
                if checkout.is_symlink() or checkout.resolve().parent != parent.resolve():
                    raise ValueError("owned candidate cleanup path escaped its temporary parent")
                cleanup = subprocess.run(
                    [*git, "worktree", "remove", "--force", str(checkout)], cwd=ROOT,
                    timeout=60)
                if cleanup.returncode:
                    raise ValueError(f"git worktree remove returned {cleanup.returncode}")
            if not checkout.exists():
                parent.rmdir()  # Empty exact parent only; no recursive fallback.
        except (ValueError, OSError, subprocess.SubprocessError) as exc:
            print(f"pre-push: owned candidate cleanup failed ({checkout}): {exc}", file=sys.stderr)
            result = 1
    return result


def candidates(lines):
    """Validate the entire input before running any remote preflight."""
    result = []
    for number, line in enumerate(lines, 1):
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(f"invalid pre-push record at line {number}")
        local_ref, head, remote_ref, base = fields
        if not all(re.fullmatch(r"[0-9a-f]{40}", sha) for sha in (head, base)):
            raise ValueError(f"invalid commit ID at line {number}")
        if not remote_ref.startswith("refs/heads/") or head == ZERO:
            continue  # Release tags have separate release gates; deletes have no candidate.
        if not local_ref:
            raise ValueError(f"missing local ref at line {number}")
        pair = (base, head)
        if pair not in result:
            result.append(pair)
    return result


def main(lines=None):
    try:
        updates = candidates(sys.stdin if lines is None else lines)
        helper = ROOT / "tools" / "ci_preflight_lab.py"
        if updates and not helper.is_file():
            raise ValueError("candidate preflight helper is missing")
        for base, head in updates:
            local_result = run_local_candidate(base, head)
            if local_result:
                return local_result
            print(f"pre-push: impact-scoped candidate preflight {head}", flush=True)
            completed = subprocess.run(
                [sys.executable, str(helper), "--base", base, "--head", head],
                cwd=ROOT,
            )
            if completed.returncode:
                return completed.returncode
        return 0
    except (ValueError, OSError, subprocess.SubprocessError) as exc:
        print(f"pre-push: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
