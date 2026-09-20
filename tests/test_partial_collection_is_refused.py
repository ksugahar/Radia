"""A run that cannot collect the whole suite must say so, not look green.

`tests/conftest.py` skips a test file whose top-level imports name a module it
cannot import, which is right for a genuinely optional extra but wrong for
`ngsolve` and `netgen`: both are pinned in `pyproject.toml`, so a failed import
means the environment is broken.  On 2026-09-20 one full run collected 4784
items and the next collected 5313 -- the same command, 529 tests apart -- and
the short run still printed a green summary, because a transient import failure
had silently removed every NGSolve-dependent file.  A partial suite that is
indistinguishable from a complete one is worse than a red one, so it now has to
be requested by name.
"""

from __future__ import annotations

import os
import pathlib
import subprocess
import sys
import tempfile

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
TARGET = "tests/test_hdiv_environment_policy.py"


def _run_with_broken(module, *, allow_partial=False):
    """Run pytest with `module` shadowed by one that fails on import."""
    shadow = tempfile.mkdtemp()
    pathlib.Path(shadow, module + ".py").write_text(
        'raise OSError("simulated DLL load failure")\n', encoding="utf-8")
    env = dict(os.environ, PYTHONPATH=shadow)
    if allow_partial:
        env["RADIA_TESTS_ALLOW_PARTIAL"] = "1"
    else:
        env.pop("RADIA_TESTS_ALLOW_PARTIAL", None)
    return subprocess.run(
        [sys.executable, "-m", "pytest", TARGET, "-q"],
        cwd=str(ROOT), capture_output=True, text=True,
        encoding="utf-8", errors="replace", check=False, env=env)


@pytest.mark.parametrize("module", ["ngsolve", "netgen"])
def test_a_broken_pinned_dependency_refuses_the_run(module):
    out = _run_with_broken(module)
    assert out.returncode != 0, out.stdout[-2000:]
    text = out.stdout + out.stderr
    assert module in text
    assert "simulated DLL load failure" in text
    # The message has to carry the size of the loss, since the whole failure
    # mode is that the loss is invisible.
    assert "would be skipped" in text
    assert "RADIA_TESTS_ALLOW_PARTIAL" in text


def test_a_partial_run_is_available_but_has_to_be_named():
    out = _run_with_broken("ngsolve", allow_partial=True)
    assert out.returncode == 0, (out.stdout + out.stderr)[-2000:]
    assert "PARTIAL RUN" in out.stdout
    assert "ngsolve" in out.stdout


def test_the_header_states_how_many_files_were_left_out():
    out = subprocess.run(
        [sys.executable, "-m", "pytest", TARGET, "-q"],
        cwd=str(ROOT), capture_output=True, text=True,
        encoding="utf-8", errors="replace", check=False)
    assert out.returncode == 0, (out.stdout + out.stderr)[-2000:]
    # Either every file is collectable, or the count of skipped ones is printed;
    # what must never happen is an unreported reduction.
    from tests.conftest import collect_ignore
    if collect_ignore:
        assert "test files not collected: %d" % len(collect_ignore) in out.stdout
