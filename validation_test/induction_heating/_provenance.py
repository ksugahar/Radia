"""Provenance of an induction-heating validation run.

In a git checkout the commit is read (and must equal ``source_commit`` when
that is given).  A copy made with ``git archive`` has no .git and must name
its commit with ``--source-commit``; it is clean by construction.  The
native extension modules that ``import radia`` loads are not in git, so
their hashes are recorded, and radia must come from this tree.
"""
from __future__ import annotations

import datetime
import hashlib
import os
import platform
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))


def _sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def provenance(argv, source_commit, files):
    """``files``: the script and the modules it exercises (paths)."""
    def git(*args):
        return subprocess.run(["git", *args], cwd=ROOT, capture_output=True,
                              text=True, check=True).stdout.strip()

    hashes = {os.path.relpath(f, ROOT).replace("\\", "/"): _sha256(f)
              for f in files}
    pkg = os.path.join(ROOT, "src", "radia")
    native = {n: _sha256(os.path.join(pkg, n)) for n in sorted(os.listdir(pkg))
              if n.endswith((".pyd", ".dll", ".so"))}
    import ngsolve
    import radia
    import scipy
    if os.path.dirname(os.path.abspath(radia.__file__)) != os.path.abspath(pkg):
        raise SystemExit(f"radia is imported from {radia.__file__}, not from "
                         f"this tree's {pkg}")
    if os.path.exists(os.path.join(ROOT, ".git")):
        commit = git("rev-parse", "HEAD")
        dirty = bool(git("status", "--porcelain", "--", "src",
                         "validation_test/induction_heating"))
        if source_commit and source_commit != commit:
            raise SystemExit(f"--source-commit {source_commit} is not the "
                             f"checkout's HEAD {commit}")
        source = "git-checkout"
    else:
        if not source_commit:
            raise SystemExit("this copy has no .git; pass --source-commit "
                             "(the commit it was archived from)")
        commit, dirty, source = source_commit, False, "git-archive"
    return {"commit": commit, "dirty": dirty, "source": source,
            "sha256": hashes, "native_sha256": native,
            "command": [os.path.basename(sys.executable), *argv],
            "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(
                timespec="seconds"),
            "host": platform.node(), "python": platform.python_version(),
            "numpy": np.__version__, "scipy": scipy.__version__,
            "ngsolve": ngsolve.__version__, "cpu_count": os.cpu_count()}
