"""The window itself, driven through its real WndProc.

Everything below radia.equation is ordinary functions with ordinary tests; the
WIN32 window on top of them -- key wiring, painting, the palette popups, the
mouse -- had none, and the editor's first crash in real use happened exactly
there: some fifty seconds of editing, an access violation, and no way to say
which operation did it.

`eqnedt64.exe --selftest` is the seam.  It injects window messages at the real
window (no keyboard, no foreground, so it runs beside a working user and on a
headless CI desktop), applies every published chord and every palette cell
from several caret states, then runs seeded random editing walks with a forced
repaint after every step.  Each step is journalled and flushed BEFORE it runs,
so a crash names its killer on the journal's last line, and the WER LocalDumps
registration for eqnedt64.exe catches the dump.

This test just runs it and reads the verdict.  A nonzero exit is either the
harness's failure count or, when a step crashed the process, the exception
code -- 0xC0000005 comes back as 3221225477.
"""

from __future__ import annotations

import pathlib
import subprocess

import pytest

_REPO = pathlib.Path(__file__).resolve().parents[2]


def _find_exe() -> pathlib.Path | None:
    """The build of THIS tree first: the post-build copy in src/radia.

    The installed radia package is second choice -- on a machine where the
    editable install points at another tree, testing that tree's exe would
    pass or fail for the wrong reasons.
    """
    local = _REPO / "src" / "radia" / "eqnedt64.exe"
    if local.exists():
        return local
    try:
        import radia
        packaged = pathlib.Path(radia.__file__).parent / "eqnedt64.exe"
        if packaged.exists():
            return packaged
    except ImportError:
        pass
    return None


_EXE = _find_exe()


@pytest.mark.skipif(_EXE is None, reason="eqnedt64.exe is not built")
def test_selftest_survives_every_operation(tmp_path):
    log = tmp_path / "selftest.log"
    proc = subprocess.run(
        [str(_EXE), "--selftest", "--log", str(log),
         "--walks", "2", "--steps", "800"],
        timeout=600,
        capture_output=True,
    )

    tail = ""
    if log.exists():
        lines = log.read_text(encoding="utf-8", errors="replace").splitlines()
        tail = "\n".join(lines[-15:])

    if proc.returncode != 0:
        hint = (
            "a crash dump, if any, is under C:\\temp\\wer_dumps\\eqnedt64; "
            "the journal's last line names the step that died"
        )
        pytest.fail(
            f"eqnedt64 --selftest exited {proc.returncode} "
            f"(0x{proc.returncode & 0xFFFFFFFF:08X}); {hint}\n"
            f"journal tail:\n{tail}"
        )

    assert log.exists(), "the selftest ran but wrote no journal"
    assert tail.splitlines()[-1].startswith("PASS"), tail
