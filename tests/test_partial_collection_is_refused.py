"""A run that quietly loses test files it was asked for must not look green.

`tests/conftest.py` skips a test file whose top-level imports name a module it
cannot import, which is right for a genuinely optional extra but wrong for
`ngsolve` and `netgen`: both are pinned in `pyproject.toml`, so a failed import
means the environment is broken.  On 2026-09-20 two consecutive full runs of
the same command collected 4784 and 5313 items -- 529 tests apart -- and the
short one still printed a green summary, because a transient import failure had
removed every NGSolve-dependent file.

The refusal is scoped to what the invocation asked for.  A tier that names its
files from the checked manifest cannot lose one this way, so the fast contract
lane on mdx keeps running without NGSolve installed; only a directory-wide
collection can shrink without saying so.
"""

from __future__ import annotations

import os
import pathlib
import subprocess
import sys
import tempfile

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
NAMED_FILE = "tests/test_hdiv_environment_policy.py"


def _pytest(args, broken=None, allow_partial=False):
    """Run pytest with `broken` shadowed by a module that fails on import."""
    env = dict(os.environ)
    if broken is not None:
        shadow = tempfile.mkdtemp()
        pathlib.Path(shadow, broken + ".py").write_text(
            'raise OSError("simulated DLL load failure")', encoding="utf-8")
        env["PYTHONPATH"] = shadow
    if allow_partial:
        env["RADIA_TESTS_ALLOW_PARTIAL"] = "1"
    else:
        env.pop("RADIA_TESTS_ALLOW_PARTIAL", None)
    return subprocess.run(
        [sys.executable, "-m", "pytest", *args],
        cwd=str(ROOT), capture_output=True, text=True,
        encoding="utf-8", errors="replace", check=False, env=env)


@pytest.mark.parametrize("module", ["ngsolve", "netgen"])
def test_a_broken_pinned_dependency_refuses_a_whole_directory_run(module):
    out = _pytest(["tests", "--collect-only", "-q"], broken=module)
    text = out.stdout + out.stderr
    assert out.returncode != 0, text[-2000:]
    assert module in text
    assert "simulated DLL load failure" in text
    # The message carries the size of the loss, because the whole failure mode
    # is that the loss is invisible.
    assert "would be skipped" in text
    assert "RADIA_TESTS_ALLOW_PARTIAL" in text


def test_a_named_file_still_runs_without_the_pinned_dependency():
    """The mdx fast-contracts lane: explicit files, no NGSolve, nothing lost."""
    out = _pytest([NAMED_FILE, "-q"], broken="ngsolve")
    assert out.returncode == 0, (out.stdout + out.stderr)[-2000:]


def test_a_partial_directory_run_is_available_but_has_to_be_named():
    """Naming it gets past the refusal; it does not silence anything else.

    A broken NGSolve still breaks the imports of the subdirectory tests that
    use it, and those stay loud -- the opt-out only says the missing files
    are expected, never that the run is sound.
    """
    out = _pytest(["tests", "--collect-only", "-q"], broken="ngsolve",
                  allow_partial=True)
    assert "PARTIAL RUN" in out.stdout, (out.stdout + out.stderr)[-2000:]
    assert "ngsolve" in out.stdout
    assert "would be skipped" not in out.stdout + out.stderr


def test_a_directory_run_states_how_much_of_the_directory_it_collected():
    out = _pytest(["tests/feec", "--collect-only", "-q"])
    assert out.returncode == 0, (out.stdout + out.stderr)[-2000:]
    line = next(l for l in out.stdout.splitlines() if l.startswith("test files:"))
    collected, on_disk = (int(n) for n in line.replace(":", " ").split()
                          if n.isdigit())
    assert collected == on_disk, line

def test_the_header_states_how_many_files_were_left_out():
    out = _pytest([NAMED_FILE, "-q"])
    assert out.returncode == 0, (out.stdout + out.stderr)[-2000:]
    from tests.conftest import collect_ignore
    if collect_ignore:
        assert "test files not collected: %d" % len(collect_ignore) in out.stdout
