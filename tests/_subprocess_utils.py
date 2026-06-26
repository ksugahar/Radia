"""Shared subprocess + Cubit batch helpers for the Radia test suite.

Encapsulates the three lessons learned while building the Mesh
Evaluation p-convergence test (validation_test/panels/test_radia_export_menu.py)
so other tests do not re-discover them painfully.

Lesson 1 -- PySide6 (Qt6) ⊗ Cubit ``_cubit3.pyd`` (Qt6) DLL conflict.
    Both ship their own Qt6 binaries; loading both into the same
    Windows process resolves the wrong DLL chain and ``import cubit``
    raises ``ImportError: DLL load failed while importing _cubit3``.
    Mitigation: run Cubit batch in a clean child Python process
    (``run_cubit_python_script`` below).  Tests that ONLY use Cubit
    (no Qt) can ``import cubit`` directly.

Lesson 2 -- system Python 3.12 numpy/NGSolve ⊗ Cubit bundled
    Python 3.10 numpy.  ``cubit.init()`` prepends Cubit's bundled
    ``bin/python3/lib/site-packages`` to ``sys.path`` so any later
    ``import numpy`` resolves the wrong (3.10-ABI) module and raises
    ``Importing the numpy C-extensions failed``.  Mitigation:
    pre-import numpy / NGSolve / calc_mesh_eval / etc. BEFORE
    ``cubit.init()`` -- the modules become cached in ``sys.modules``
    and Cubit's later sys.path tweak no longer affects them.

Lesson 3 -- Cubit and NGSolve emit non-ASCII on stderr in Japanese
    Windows locale (cp932).  ``subprocess.run(..., text=True)`` uses
    the system default codec and crashes its reader thread with
    ``UnicodeDecodeError: 'cp932' codec can't decode byte 0x94``,
    losing the entire subprocess output INCLUDING the real failure
    cause.  Mitigation: read raw bytes and decode with
    ``errors="replace"``, plus set ``PYTHONIOENCODING=utf-8`` in the
    child env so its own print() emits UTF-8.

This module is intentionally dependency-free: stdlib only, no
PySide6 / numpy / NGSolve import at module top, so any test under
tests/ can import it without paying for unrelated DLL setup.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path


# ----------------------------------------------------------------------
# Lesson 3 -- UTF-8-safe subprocess.run for any Cubit / NGSolve invoker
# ----------------------------------------------------------------------

def run_subprocess_utf8(argv, timeout=300, extra_env=None,
                         cwd=None):
    """Drop-in replacement for ``subprocess.run(text=True)`` that is
    hardened against:
      * cp932 stderr from Cubit / NGSolve / Fortran MKL diagnostics
      * QT_QPA_PLATFORM / PYSIDE inheritance into the child process

    Returns (returncode, stdout_str, stderr_str).
    On ``subprocess.TimeoutExpired``, re-raises -- the caller decides
    whether to skip, fail, or retry.

    Args:
        argv: list of command-line tokens (same shape as
              subprocess.run's first arg).
        timeout: seconds (passed through to subprocess.run).
        extra_env: optional dict of additional env vars to set in the
              child (overrides inherited values).
        cwd: optional working directory for the child.
    """
    env = os.environ.copy()
    for k in list(env.keys()):
        if k.upper().startswith("QT_") or "PYSIDE" in k.upper():
            del env[k]
    env["PYTHONIOENCODING"] = "utf-8"
    if extra_env:
        env.update(extra_env)

    # text=False -> raw bytes; decode ourselves with errors="replace"
    proc = subprocess.run(argv, capture_output=True, timeout=timeout,
                          env=env, cwd=cwd)
    stdout = (proc.stdout or b"").decode("utf-8", errors="replace")
    stderr = (proc.stderr or b"").decode("utf-8", errors="replace")
    return proc.returncode, stdout, stderr


# ----------------------------------------------------------------------
# Lessons 1 + 2 + 3 -- run a script that needs Cubit batch in a CLEAN
# child Python interpreter
# ----------------------------------------------------------------------

# Wrapper template injected around the user script.  Pre-imports
# numpy + NGSolve + netgen BEFORE cubit.init() so Cubit's bundled
# Python 3.10 site-packages cannot shadow the system 3.12 ones
# (lesson 2).  The user script can then ``import cubit`` and call
# ``cubit.init(...)`` itself.
_CUBIT_BATCH_PRELUDE = """
# AUTO-INJECTED by tests/_subprocess_utils.py::run_cubit_python_script
# -- pre-binds system numpy / NGSolve so cubit.init() cannot poison
# their sys.path resolution.  Safe to skip individual imports that
# the user script does not need.
import os, sys

# Make src/radia + src/radia/panels importable for follow-on imports.
# Forward-slash so Python parses the same on Windows (no \\ escaping).
_REPO = "{repo}"
sys.path.insert(0, os.path.join(_REPO, "src", "radia"))
sys.path.insert(0, os.path.join(_REPO, "src", "radia", "panels"))

try:
    import numpy            # noqa: F401  -- bind system 3.12 numpy
except ImportError:
    pass
try:
    import netgen.meshing   # noqa: F401  -- bind NGSolve native deps
except ImportError:
    pass
try:
    from ngsolve import Mesh as _NGMesh  # noqa: F401
except ImportError:
    pass

# Help the user script locate Cubit's bin (idempotent).
try:
    from install_panels import find_cubit_bin as _find_cubit_bin
    _cb = _find_cubit_bin()
    if _cb:
        sys.path.insert(0, _cb)
        os.environ.setdefault(
            "CUBIT_PLUGIN_DIR",
            os.path.join(os.path.dirname(_cb), "bin", "plugins"))
except Exception:
    pass
# --------- END AUTO-INJECTED PRELUDE -------------------------------
"""


def run_cubit_python_script(script_src, args=(), timeout=300,
                              repo_root=None):
    """Run an arbitrary Python script in a CLEAN child process where:
      * PySide6 has NOT been loaded (avoids Qt6 DLL conflict with
        Cubit's _cubit3.pyd -- lesson 1);
      * system numpy / NGSolve / netgen are pre-bound BEFORE the
        script gets a chance to call cubit.init() (lesson 2);
      * stdout / stderr are decoded UTF-8 with errors="replace"
        (lesson 3).

    The script body is the caller's responsibility: it should do
    ``import cubit; cubit.init([...])`` and emit results on stdout
    (typically as a single JSON line wrapped in __BEGIN_JSON__ /
    __END_JSON__ markers -- see ``extract_json_from_output``).

    Args:
        script_src: string -- the Python source of the script body.
                    The prelude (numpy / NGSolve / find_cubit_bin) is
                    prepended automatically.
        args: extra argv passed AFTER the script path.
        timeout: seconds.
        repo_root: optional absolute path to the Radia repo root; used
                    to seed sys.path.  Defaults to the repo containing
                    THIS file.

    Returns: (returncode, stdout_str, stderr_str)
    """
    if repo_root is None:
        # tests/_subprocess_utils.py -> repo root is parent of tests/
        repo_root = str(Path(__file__).resolve().parent.parent)

    # Forward-slashes -> Python accepts them on Windows AND the source
    # we generate has no '\' escape characters that would mis-parse.
    prelude = _CUBIT_BATCH_PRELUDE.format(
        repo=repo_root.replace("\\", "/"))
    full_src = prelude + "\n" + script_src

    # Write to a temp file so the child has a real __file__ and
    # traceback line numbers point somewhere useful.
    fd, path = tempfile.mkstemp(prefix="cubit_batch_", suffix=".py")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(full_src)
        argv = [sys.executable, path, *map(str, args)]
        return run_subprocess_utf8(argv, timeout=timeout)
    finally:
        # Best-effort cleanup; ignore on Windows if the child still
        # has the file open.
        try:
            os.unlink(path)
        except OSError:
            pass


# ----------------------------------------------------------------------
# Convenience: extract a single JSON object from BEGIN/END markers
# ----------------------------------------------------------------------

def extract_json_from_output(stdout, begin_marker="__BEGIN_JSON__",
                              end_marker="__END_JSON__"):
    """Pull the JSON blob the child wrote between BEGIN / END markers.

    Returns the parsed dict, or raises ValueError with a diagnostic
    that includes the last bit of stdout so the test can pinpoint why
    the marker is missing.
    """
    import json
    if begin_marker not in stdout or end_marker not in stdout:
        tail = stdout[-1500:] if stdout else "(empty stdout)"
        raise ValueError(
            f"missing {begin_marker} / {end_marker} markers in "
            f"subprocess stdout.  Tail:\n{tail}")
    start = stdout.index(begin_marker) + len(begin_marker)
    end = stdout.index(end_marker)
    return json.loads(stdout[start:end].strip())
