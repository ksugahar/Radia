r"""What the EDITOR puts on the clipboard, byte for byte.

`tests/equation/test_paste_size.py` places the GVML package on the clipboard
itself and asks PowerPoint what size it pasted at.  That proves the package is
right.  It does not prove the editor puts that package on the clipboard, and
for a while the editor did not: its `put()` asked for one byte more than the
payload and wrote a NUL there.

For RTF and MathML that is correct -- they are read as C strings.  For the GVML
package, which is an OPC ZIP, it is fatal: PowerPoint refuses to open an
archive with a byte after the end, and the paste fails outright.  No error, no
picture, nothing on the slide.  Ctrl+C in the editor was broken for PowerPoint
while the Python API path, which writes exact-length buffers, worked -- and the
test above could not see the difference because it never went through the
editor.

So this drives the real gesture: open the window, press Ctrl+C, read the bytes
back.  Keys are POSTED with the input queues attached, which is how a modifier
registers in another process without touching the real keyboard.

It uses the clipboard, so it lives in validation_test rather than tests/.
"""

from __future__ import annotations

import ctypes
import ctypes.wintypes as wt
import pathlib
import subprocess
import sys
import time

import pytest

for _p in pathlib.Path(__file__).resolve().parents:
    if (_p / "src" / "radia").exists():
        _REPO = _p
        sys.path.insert(0, str(_p / "src"))
        break

equation = pytest.importorskip("radia._equation")

EXE = _REPO / "src" / "radia" / "eqnedt64.exe"
LATEX = r"\dfrac{a}{b}+\sqrt{c}"

u32 = ctypes.WinDLL("user32", use_last_error=True)
k32 = ctypes.WinDLL("kernel32", use_last_error=True)
u32.RegisterClipboardFormatW.argtypes = [ctypes.c_wchar_p]
u32.RegisterClipboardFormatW.restype = wt.UINT
u32.GetClipboardData.argtypes = [wt.UINT]
u32.GetClipboardData.restype = ctypes.c_void_p
u32.EnumWindows.argtypes = [ctypes.c_void_p, wt.LPARAM]
u32.GetClassNameW.argtypes = [ctypes.c_void_p, ctypes.c_wchar_p, ctypes.c_int]
u32.GetWindowThreadProcessId.argtypes = [ctypes.c_void_p,
                                         ctypes.POINTER(wt.DWORD)]
u32.PostMessageW.argtypes = [ctypes.c_void_p, wt.UINT, ctypes.c_void_p,
                             ctypes.c_void_p]
u32.AttachThreadInput.argtypes = [wt.DWORD, wt.DWORD, wt.BOOL]
k32.GlobalLock.argtypes = [ctypes.c_void_p]
k32.GlobalLock.restype = ctypes.c_void_p
k32.GlobalUnlock.argtypes = [ctypes.c_void_p]
k32.GlobalSize.argtypes = [ctypes.c_void_p]
k32.GlobalSize.restype = ctypes.c_size_t


def editor_window(pid):
    found = []
    PROC = ctypes.WINFUNCTYPE(wt.BOOL, ctypes.c_void_p, wt.LPARAM)

    def cb(h, _):
        got = wt.DWORD()
        u32.GetWindowThreadProcessId(h, ctypes.byref(got))
        if got.value != pid:
            return True
        name = ctypes.create_unicode_buffer(64)
        u32.GetClassNameW(h, name, 64)
        if name.value == "Eqnedt64Window":
            found.append(h)
            return False
        return True

    u32.EnumWindows(PROC(cb), 0)
    return found[0] if found else None


def press_ctrl_c(hwnd):
    """Ctrl+C, posted, with the modifier held across the post.

    SetKeyboardState writes THIS thread's input state; a window in another
    process reads its own, so the queues have to be attached first.  And the
    modifier must stay down until the other thread has drained the message --
    clearing it straight away races, and Ctrl+C arrives as a plain c.
    """
    got = wt.DWORD()
    theirs = u32.GetWindowThreadProcessId(hwnd, ctypes.byref(got))
    mine = k32.GetCurrentThreadId()
    u32.AttachThreadInput(mine, theirs, True)
    try:
        ks = (ctypes.c_ubyte * 256)()
        u32.GetKeyboardState(ks)
        ks[0x11] = ks[0xA2] = 0x80
        u32.SetKeyboardState(ks)
        u32.PostMessageW(hwnd, 0x0100, ctypes.c_void_p(0x43), None)
        time.sleep(1.5)
        u32.PostMessageW(hwnd, 0x0101, ctypes.c_void_p(0x43), None)
        time.sleep(0.8)
        ks[0x11] = ks[0xA2] = 0
        u32.SetKeyboardState(ks)
    finally:
        u32.AttachThreadInput(mine, theirs, False)


def clipboard_bytes(name):
    fmt = u32.RegisterClipboardFormatW(name)
    if not u32.OpenClipboard(None):
        return None
    try:
        h = u32.GetClipboardData(fmt)
        if not h:
            return None
        size = k32.GlobalSize(h)
        p = k32.GlobalLock(h)
        try:
            return ctypes.string_at(p, size)
        finally:
            k32.GlobalUnlock(h)
    finally:
        u32.CloseClipboard()


@pytest.mark.skipif(not EXE.exists(), reason="eqnedt64.exe is not built")
def test_ctrl_c_puts_the_gvml_package_on_the_clipboard_exactly():
    proc = subprocess.Popen([str(EXE), LATEX])
    try:
        hwnd = None
        for _ in range(20):
            time.sleep(0.5)
            hwnd = editor_window(proc.pid)
            if hwnd:
                break
        assert hwnd, "the editor window never appeared"

        press_ctrl_c(hwnd)
        got = clipboard_bytes("Art::GVML ClipFormat")
        assert got, "Ctrl+C put no GVML package on the clipboard"
    finally:
        proc.kill()
        # Reap it.  A Popen collected without wait() raises from __del__, and
        # pytest turns that into a failure of whatever test happens to be
        # running when the collector gets round to it.
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            pass

    want = equation.tex_to_gvml(LATEX, equation.PASTE_SIZE_PT, False)

    assert got[:4] == b"PK\x03\x04", (
        f"the package does not start as a ZIP: {got[:8]!r}")
    assert len(got) == len(want), (
        f"the editor wrote {len(got)} bytes where the package is {len(want)}; "
        "a byte after the end of an OPC archive makes PowerPoint refuse the "
        "paste outright")
    assert got == want, "the package on the clipboard is not the one we build"
