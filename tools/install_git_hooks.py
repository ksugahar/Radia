#!/usr/bin/env python
"""Install Radia's version-controlled git hooks into this clone's hooks dir.

The hooks live tracked under tools/git-hooks/ so the source is reviewable and
shared; this installer copies them into the real hooks directory (executable)
so they actually run.  Run ONCE per clone:

    python tools/install_git_hooks.py

Installs:
  pre-push  -- CI preflight gate (tools/ci_preflight.py) on pushes to main,
               plus the binary-upload-to-GitHub-Releases step.  Bypass the
               gate for an emergency push with:  CI_PREFLIGHT_SKIP=1 git push

Idempotent: re-run any time tools/git-hooks/ changes.  Worktree-safe (uses
git's common dir so hooks land in the shared .git, not a worktree stub).
"""
import os
import stat
import subprocess
import sys

# UNC-safe: do not Path.resolve() on the LAB S: drive (it canonicalises to a
# \\192.168.11.100\... form that some tools choke on).
_THIS = os.path.abspath(__file__)
REPO = os.path.dirname(os.path.dirname(_THIS))
SRC = os.path.join(REPO, "tools", "git-hooks")


def _hooks_dir():
    """The real hooks directory (handles worktrees via --git-common-dir)."""
    r = subprocess.run(["git", "-c", "safe.directory=*", "rev-parse",
                        "--git-common-dir"],
                       cwd=REPO, capture_output=True, text=True)
    if r.returncode == 0 and r.stdout.strip():
        gd = r.stdout.strip()
        if not os.path.isabs(gd):
            gd = os.path.join(REPO, gd)
        return os.path.join(gd, "hooks")
    return os.path.join(REPO, ".git", "hooks")


def main():
    if not os.path.isdir(SRC):
        print(f"ERROR: hook source dir not found: {SRC}", file=sys.stderr)
        return 1
    dst_dir = _hooks_dir()
    os.makedirs(dst_dir, exist_ok=True)

    installed = 0
    for name in sorted(os.listdir(SRC)):
        src = os.path.join(SRC, name)
        if not os.path.isfile(src):
            continue
        dst = os.path.join(dst_dir, name)
        # Read/write as bytes to preserve LF line endings (git bash hooks
        # must be LF even on Windows).
        with open(src, "rb") as f:
            data = f.read().replace(b"\r\n", b"\n")
        with open(dst, "wb") as f:
            f.write(data)
        os.chmod(dst, os.stat(dst).st_mode
                 | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        print(f"  installed {name} -> {dst}")
        installed += 1

    if installed == 0:
        print("no hooks found to install", file=sys.stderr)
        return 1
    print(f"{installed} hook(s) installed.")
    print("The pre-push gate runs tools/ci_preflight.py on pushes to main.")
    print("Bypass for an emergency push:  CI_PREFLIGHT_SKIP=1 git push")
    return 0


if __name__ == "__main__":
    sys.exit(main())
