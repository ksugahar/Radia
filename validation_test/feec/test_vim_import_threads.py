import subprocess
import sys
import textwrap


def test_import_radia_vim_does_not_override_taskmanager_threads():
    code = textwrap.dedent(
        """
        import ngsolve as ng

        ng.SetNumThreads(7)
        before = ng.ngsglobals.numthreads
        import radia.vim  # noqa: F401
        after = ng.ngsglobals.numthreads

        assert before == 7, before
        assert after == 7, (before, after)
        """
    )
    proc = subprocess.run(
        [sys.executable, "-c", code],
        text=True,
        capture_output=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
