#!/usr/bin/env python3
"""Check every branch candidate from Git pre-push stdin locally, then on the build target.

The local stage exports the exact outgoing SHA into an owned standalone
repository (objects borrowed through alternates); it never registers a
worktree in this clone.
"""
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).absolute().parents[1]
ZERO = "0" * 40


# Overall budget for the local candidate checks (policy lint + compact contract
# tests).  The shared precheck measured 255.76 s for the policy tests alone;
# the budget leaves headroom for slower runs and still ends a hung check.
LOCAL_CANDIDATE_BUDGET_S = 2400


def child_env():
    """Environment for candidate processes: no inherited hook-local Git state."""
    return {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}


def _git(args, cwd, **kw):
    return subprocess.run(["git", *args], cwd=cwd, env=child_env(), check=True,
                          capture_output=True, text=True, encoding="utf-8",
                          timeout=kw.pop("timeout", 60), **kw)


def export_candidate(checkout, head, base):
    """Check out exact ``head`` in an owned standalone repository.

    The repository borrows objects from this clone through alternates, so the
    candidate and base history are readable without copying objects, and it
    is never registered as a worktree of this clone.
    """
    common = _git(["-c", f"safe.directory={ROOT}", "rev-parse", "--path-format=absolute",
                   "--git-common-dir"], ROOT).stdout.strip()
    objects = Path(common) / "objects"
    if not objects.is_dir():
        raise ValueError(f"object store not found: {objects}")
    _git(["init", "-q", str(checkout)], ROOT)
    (checkout / ".git" / "objects" / "info" / "alternates").write_text(
        str(objects) + "\n", encoding="utf-8", newline="\n")  # LF: Git rejects a CR here
    iso = ["-c", f"safe.directory={checkout}"]
    _git([*iso, "update-ref", "refs/radia-push/candidate", head], checkout)
    _git([*iso, "update-ref", "refs/radia-push/base", base], checkout)
    _git([*iso, "-c", "core.autocrlf=false", "checkout", "-q", "--detach", head], checkout,
         timeout=300)
    actual = _git([*iso, "rev-parse", "HEAD"], checkout).stdout.strip()
    if actual != head:
        raise ValueError(f"candidate export is at {actual}, expected {head}")


def _remove_owned(checkout, parent):
    """Delete the owned export only: exact path under its own temporary parent."""
    if (checkout.is_symlink() or os.path.isjunction(checkout)
            or checkout.resolve().parent != parent.resolve()):
        raise ValueError("owned candidate cleanup path escaped its temporary parent")
    for path in checkout.rglob("*"):
        if path.is_symlink() or (os.name == "nt" and os.path.isjunction(path)):
            raise ValueError(f"owned candidate contains a link, not removed: {path}")

    def writable(func, path, _exc):
        os.chmod(path, stat.S_IWRITE)
        func(path)

    shutil.rmtree(checkout, onexc=writable)


def run_local_candidate(base, head):
    """Local checks run on an exact exported candidate, never another session's WIP."""
    git = ["git", "-c", f"safe.directory={ROOT}"]
    try:
        if base == ZERO:
            base = subprocess.check_output(
                [*git, "merge-base", head, "origin/main"], cwd=ROOT, env=child_env(),
                text=True, encoding="utf-8").strip()
        subprocess.run([*git, "cat-file", "-e", f"{base}^{{commit}}"], cwd=ROOT,
                       env=child_env(), check=True, capture_output=True, timeout=30)
        subprocess.run([*git, "merge-base", "--is-ancestor", base, head], cwd=ROOT,
                       env=child_env(), check=True, capture_output=True, timeout=30)
    except subprocess.SubprocessError as exc:
        print(f"pre-push: cannot resolve base {base} as an ancestor of {head}; "
              f"fetch the remote base commit/history and retry: {exc}", file=sys.stderr)
        return 1
    parent = Path(tempfile.mkdtemp(prefix="radia-push-candidate-"))
    checkout = parent / "candidate"
    result = 1
    started = time.monotonic()
    try:
        export_candidate(checkout, head, base)
        result = subprocess.run(
            [sys.executable, "tools/ci_preflight.py", "--local-push-candidate", head,
             "--local-push-base", base], cwd=checkout, env=child_env(),
            timeout=LOCAL_CANDIDATE_BUDGET_S).returncode
    except subprocess.TimeoutExpired:
        print(f"pre-push: local candidate checks exceeded the {LOCAL_CANDIDATE_BUDGET_S} s "
              f"budget after {time.monotonic() - started:.0f} s", file=sys.stderr)
        result = 1
    except (ValueError, OSError, subprocess.SubprocessError) as exc:
        detail = getattr(exc, "stderr", "") or ""
        print(f"pre-push: local candidate check failed: {exc} {detail}".rstrip(), file=sys.stderr)
        result = 1
    finally:
        try:
            if checkout.exists():
                _remove_owned(checkout, parent)
            parent.rmdir()  # Empty exact parent only; no recursive fallback.
        except (ValueError, OSError) as exc:
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
        helper = ROOT / "tools" / "ci_preflight_remote.py"
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
